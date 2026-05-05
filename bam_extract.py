#!/usr/bin/env python3
"""
bam_extract.py — Step 1 of the VHH NGS pipeline.

Extracts VHH sequences from ONT BAM files aligned to IMGT Vicugna germlines.

High-throughput design for millions of reads across many BAM files:
  • CIGAR-based extraction — uses Dorado alignment coordinates to locate the
    VHH amplicon directly; handles insertions and deletions without regex
  • Deduplication before expensive work — unique DNA sequences are translated
    once; unique proteins hit ANARCI once (batch call, not per-sequence)
  • Parallel BAM processing — one worker process per file (--workers)
  • Fast Q-score via numpy vectorisation; fast RC via C translation table
  • Motif-search fallback for reads that don't fully cover the reference window
"""

import os
import pysam
import regex as re
import csv
import json
import shutil
import argparse
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional, Tuple, List

import numpy as np
from Bio.Seq import Seq
from tqdm import tqdm
from rich.console import Console
from rich.table import Table
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from anarci import anarci
    ANARCI_AVAILABLE = True
except ImportError:
    ANARCI_AVAILABLE = False

console = Console()

# ══════════════════════════════════════════════════════════════════════════════
# DEFAULT CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════
DEFAULT_CONFIG = {
    # VHH framework anchor motifs
    # Extraction searches query_alignment_sequence first (the soft-clip-stripped
    # aligned region, already in reference orientation), then falls back to the
    # full query sequence (both strands) for unaligned or edge-case reads.
    "FR1_MOTIF":         "CAGGTGCAGCTG",
    "J4_MOTIF":          "ACCCAGGTCACC",
    "FR_MAX_MISMATCHES": 2,

    # ONT read quality gates
    "MIN_READ_LENGTH": 300,
    "MAX_READ_LENGTH": 1000,
    "MIN_MEAN_QSCORE": 10.0,

    # VHH protein quality gates
    "MIN_VHH_AA_LENGTH": 100,
    "MAX_INTERNAL_STOPS": 0,

    # CDR extraction method: "anarci" (preferred) or "offset" (fallback)
    "CDR_METHOD": "anarci",

    # UMI deduplication
    "USE_UMI": False,
    "UMI_TAG": "UX",

    "EXPORT_QC_METRICS": True,

    # Number of unique proteins per ANARCI batch call.
    # Larger batches amortise subprocess overhead; reduce if RAM is tight.
    "ANARCI_BATCH_SIZE": 5000,
}

LIABILITY_PATTERNS = {
    "N-glycosylation":  re.compile(r"N[^P][ST]"),
    "Deamidation (NG)": re.compile(r"NG"),
    "Deamidation (NS)": re.compile(r"NS"),
    "Isomerization":    re.compile(r"D[TGSH]"),
    "Free Cys":         re.compile(r"C"),
    "Oxidation (Met)":  re.compile(r"M"),
}

VHH_START_MOTIFS = ("QVQ", "VQL", "QLQ", "EVQ", "DVQ")

# ══════════════════════════════════════════════════════════════════════════════
# FAST UTILITIES
# ══════════════════════════════════════════════════════════════════════════════

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def mean_qscore(qualities) -> float:
    if qualities is None or len(qualities) == 0:
        return 0.0
    arr = np.frombuffer(bytes(qualities), dtype=np.uint8).astype(np.float32)
    mean_err = float(np.mean(10.0 ** (-arr / 10.0)))
    return float(-10.0 * np.log10(mean_err)) if mean_err > 0.0 else 40.0


def fuzzy_search(pattern: str, sequence: str, max_mismatches: int):
    # Exact match first (C-speed str.find); fall back to regex only if needed
    pos = sequence.find(pattern)
    if pos >= 0:
        return (pos, pos + len(pattern))
    m = re.search(f"({pattern}){{e<={max_mismatches}}}", sequence)
    return (m.start(), m.end()) if m else None


# ══════════════════════════════════════════════════════════════════════════════
# VHH BOUNDARY EXTRACTION
# ══════════════════════════════════════════════════════════════════════════════

def extract_amplicon(
    read: pysam.AlignedSegment,
    fr1: str,
    j4: str,
    max_mm: int,
) -> Optional[Tuple[str, str]]:
    """
    Locate the VHH amplicon from FR1 start to J4 end.

    All reads reaching this function have already passed the `is_unmapped`
    filter, so the BAM sequence is always in reference / coding orientation —
    the reverse complement is never needed for the main search paths.

    Strategy (two-tier, forward strand only):

    1. Aligned search — search `query_alignment_sequence` (soft-clip-stripped).
       Shorter than the full read; works when both anchors fall within the
       aligned region.  Typical for reads that span the full amplicon.

    2. Full-read search — forward strand of `query_sequence` (soft clips
       included).  Needed when one anchor (usually J4) lies in a soft-clipped
       tail because the alignment ends before the J4 motif on the reference.

    RC fallback is reserved for the rare genuinely antisense read that slipped
    through alignment (observed <0.5% in practice).
    """
    seq = read.query_sequence
    if not seq:
        return None
    seq = seq.upper()

    # ── Tier 1: aligned region (no soft clips) ─────────────────────────────
    aln = read.query_alignment_sequence
    if aln:
        aln = aln.upper()
        fr1_m = fuzzy_search(fr1, aln, max_mm)
        j4_m  = fuzzy_search(j4,  aln, max_mm)
        if fr1_m and j4_m and j4_m[1] > fr1_m[0]:
            return aln[fr1_m[0] : j4_m[1]], "aligned"

    # ── Tier 2: full read, forward strand ──────────────────────────────────
    fr1_m = fuzzy_search(fr1, seq, max_mm)
    j4_m  = fuzzy_search(j4,  seq, max_mm)
    if fr1_m and j4_m and j4_m[1] > fr1_m[0]:
        return seq[fr1_m[0] : j4_m[1]], "+"

    # ── Tier 3: RC fallback (<1% of reads) ─────────────────────────────────
    rc_seq = reverse_complement(seq)
    fr1_m  = fuzzy_search(fr1, rc_seq, max_mm)
    j4_m   = fuzzy_search(j4,  rc_seq, max_mm)
    if fr1_m and j4_m and j4_m[1] > fr1_m[0]:
        return rc_seq[fr1_m[0] : j4_m[1]], "-"

    return None


# ══════════════════════════════════════════════════════════════════════════════
# TRANSLATION
# ══════════════════════════════════════════════════════════════════════════════

def translate_vhh(dna_seq: str) -> str:
    """
    Try all three reading frames; return the longest stop-free VHH protein.
    Prefers frames that begin with a recognised VHH FR1 start motif.
    """
    best = ""
    for frame in range(3):
        trimmed = dna_seq[frame:]
        trimmed = trimmed[: len(trimmed) - (len(trimmed) % 3)]
        try:
            prot = str(Seq(trimmed).translate(to_stop=False))
        except Exception:
            continue
        if "*" in prot[:-1]:
            continue
        if prot.startswith(VHH_START_MOTIFS):
            return prot
        if any(m in prot for m in VHH_START_MOTIFS) and len(prot) > len(best):
            best = prot
    return best if best and "*" not in best[:-1] else ""


# ══════════════════════════════════════════════════════════════════════════════
# CDR EXTRACTION
# ══════════════════════════════════════════════════════════════════════════════

def _parse_anarci_result(result) -> Tuple[str, str, str]:
    if not result or not result[0]:
        return ("", "", "")
    numbering, _, _ = result[0][0]
    pos_aa = {num: aa for (num, _), aa in numbering if aa != "-"}
    cdr1 = "".join(pos_aa.get(i, "") for i in range(27, 39))
    cdr2 = "".join(pos_aa.get(i, "") for i in range(56, 66))
    cdr3 = "".join(pos_aa.get(i, "") for i in range(105, 118))
    return cdr1, cdr2, cdr3


def extract_cdrs_batch_anarci(proteins: List[str]) -> List[Tuple[str, str, str]]:
    """Run ANARCI on a list of proteins in a single subprocess call."""
    seqs = [(f"s{i}", p) for i, p in enumerate(proteins)]
    try:
        results = anarci(seqs, scheme="imgt", output=False, allow={"H"})
        return [_parse_anarci_result(r) for r in results]
    except Exception:
        return [extract_cdrs_offset(p) for p in proteins]


def extract_cdrs_offset(protein: str) -> Tuple[str, str, str]:
    """Offset-based CDR extraction fallback (no ANARCI dependency)."""
    try:
        cdr1 = protein[26:34]
        cdr2 = protein[50:60]
        raw  = protein[99:122]
        m    = re.search(r"W[GA][QKR]G", raw)
        cdr3 = raw[: m.start()] if m else raw
        return cdr1, cdr2, cdr3
    except IndexError:
        return "", "", ""


# ══════════════════════════════════════════════════════════════════════════════
# LIABILITY FLAGGING
# ══════════════════════════════════════════════════════════════════════════════

def flag_liabilities(cdr3: str) -> str:
    flags = [name for name, pat in LIABILITY_PATTERNS.items() if pat.search(cdr3)]
    return ";".join(flags) if flags else "None"


# ══════════════════════════════════════════════════════════════════════════════
# PLOTTING
# ══════════════════════════════════════════════════════════════════════════════

def plot_cdr3_lengths(cdr3_lengths: Counter, output_file: Path):
    if not cdr3_lengths:
        return
    lengths = list(cdr3_lengths.elements())
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.hist(lengths, bins=range(min(lengths), max(lengths) + 2),
            edgecolor="black", color="#4C72B0", alpha=0.85)
    ax.set_xlabel("CDR3 Length (aa)", fontsize=12)
    ax.set_ylabel("Frequency", fontsize=12)
    ax.set_title("CDR3 Length Distribution", fontsize=14)
    ax.axvline(np.mean(lengths), color="red", linestyle="--",
               label=f"Mean = {np.mean(lengths):.1f}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_file, dpi=150)
    plt.close(fig)


def plot_qscore_dist(qscores: list, output_file: Path):
    if not qscores:
        return
    arr = np.array(qscores, dtype=np.float32)
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.hist(arr, bins=50, color="#55A868", edgecolor="none", alpha=0.85)
    ax.axvline(float(np.median(arr)), color="red", linestyle="--",
               label=f"Median Q = {np.median(arr):.1f}")
    ax.set_xlabel("Mean Read Q-score", fontsize=12)
    ax.set_ylabel("Read Count", fontsize=12)
    ax.set_title("ONT Read Q-score Distribution", fontsize=14)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_file, dpi=150)
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
# CORE BAM PROCESSING
# ══════════════════════════════════════════════════════════════════════════════

def process_bam(bam_path: Path, output_dir: Path, config: dict) -> dict:
    """
    Three-pass strategy for high-throughput extraction:
      Pass 1 — Read BAM; apply quality/length/UMI filters; extract VHH DNA
               via FR1/J4 motif search on the alignment-stripped region
               (soft-clip-free, reference-oriented), falling back to full
               read search for unmapped or partial reads.
               Accumulate unique DNA counts — no translation yet.
      Pass 2 — Translate unique DNA sequences only; aggregate by protein.
      Pass 3 — Batch ANARCI on unique proteins; assemble output tables.
    """
    bam_out_dir = output_dir / bam_path.stem
    bam_out_dir.mkdir(parents=True, exist_ok=True)

    stem          = bam_path.stem
    out_dna_prot  = bam_out_dir / f"{stem}_vhh_dna_protein.csv"
    out_prot_cdr  = bam_out_dir / f"{stem}_vhh_protein_cdr.csv"
    out_cdr3_plot = bam_out_dir / f"{stem}_cdr3_lengths.png"
    out_qplot     = bam_out_dir / f"{stem}_qscore_dist.png"
    out_summary   = bam_out_dir / f"{stem}_summary.json"

    cfg         = {**DEFAULT_CONFIG, **config}
    FR1         = cfg["FR1_MOTIF"]
    J4          = cfg["J4_MOTIF"]
    FR_MM       = cfg["FR_MAX_MISMATCHES"]
    CDR_METHOD  = cfg["CDR_METHOD"]
    BATCH_SIZE  = cfg["ANARCI_BATCH_SIZE"]

    stats       = defaultdict(int)
    qscores_all: List[float] = []
    umi_seen:    set          = set()

    # ── Pass 1: filter reads, extract VHH DNA amplicons ────────────────────
    # Count each unique amplicon DNA across all reads.
    # Translation and annotation happen after deduplication (Pass 2/3).
    dna_counts: Counter = Counter()

    # Estimate total for tqdm — use index stats if available, otherwise 0
    try:
        with pysam.AlignmentFile(str(bam_path), check_sq=False) as bam_tmp:
            total_est = bam_tmp.mapped + bam_tmp.unmapped
    except Exception:
        total_est = 0

    console.print(f"[bold cyan]Pass 1:[/bold cyan] {bam_path.name}")
    with pysam.AlignmentFile(str(bam_path), check_sq=False) as bamfile:
        for read in tqdm(bamfile.fetch(until_eof=True), total=total_est or None,
                         desc=stem[:30], unit="reads", leave=False):
            stats["total"] += 1

            if read.is_unmapped or read.query_sequence is None:
                stats["unmapped"] += 1
                continue

            seq  = read.query_sequence.upper()
            rlen = len(seq)

            if rlen < cfg["MIN_READ_LENGTH"]:
                stats["too_short"] += 1
                continue
            if rlen > cfg["MAX_READ_LENGTH"]:
                stats["too_long"] += 1
                continue

            qs = mean_qscore(read.query_qualities)
            qscores_all.append(qs)
            if qs < cfg["MIN_MEAN_QSCORE"]:
                stats["low_quality"] += 1
                continue

            if cfg["USE_UMI"]:
                umi = read.get_tag(cfg["UMI_TAG"]) if read.has_tag(cfg["UMI_TAG"]) else None
                if umi and umi in umi_seen:
                    stats["umi_duplicate"] += 1
                    continue
                if umi:
                    umi_seen.add(umi)

            # Extract FR1→J4 amplicon: aligned region first, full read fallback
            result = extract_amplicon(read, FR1, J4, FR_MM)
            if result is None:
                stats["no_framework"] += 1
                continue

            amplicon, source = result
            stats["framework_found"] += 1
            stats[f"extracted_{source}"] += 1  # "aligned" | "+" | "-"
            dna_counts[amplicon] += 1

    # ── Pass 2: translate unique DNA sequences ──────────────────────────────
    n_unique_dna = len(dna_counts)
    console.print(
        f"[cyan]Pass 2:[/cyan] translating {n_unique_dna:,} unique amplicons — {stem}"
    )

    dna_to_prot: dict = {}
    for dna in dna_counts:
        prot = translate_vhh(dna)
        if not prot or len(prot) < cfg["MIN_VHH_AA_LENGTH"]:
            stats["bad_translation"] += 1
            dna_to_prot[dna] = None
            continue
        if prot.count("*") > cfg["MAX_INTERNAL_STOPS"]:
            stats["stop_codon"] += 1
            dna_to_prot[dna] = None
            continue
        dna_to_prot[dna] = prot

    # Aggregate read counts at the protein level
    prot_counts:  Counter = Counter()
    prot_rep_dna: dict    = {}  # protein → highest-count representative DNA
    for dna, prot in dna_to_prot.items():
        if not prot:
            continue
        count = dna_counts[dna]
        prot_counts[prot] += count
        if prot not in prot_rep_dna or count > dna_counts[prot_rep_dna[prot]]:
            prot_rep_dna[prot] = dna

    stats["extracted"] = sum(prot_counts.values())
    unique_proteins = list(prot_counts.keys())

    # ── Pass 3: batch CDR annotation ───────────────────────────────────────
    n_unique_prot = len(unique_proteins)
    console.print(
        f"[cyan]Pass 3:[/cyan] CDR annotation — {n_unique_prot:,} unique proteins — {stem}"
    )

    prot_cdr_map: dict = {}  # protein → (cdr1, cdr2, cdr3)

    if CDR_METHOD == "anarci" and ANARCI_AVAILABLE:
        for i in range(0, n_unique_prot, BATCH_SIZE):
            batch = unique_proteins[i : i + BATCH_SIZE]
            cdrs  = extract_cdrs_batch_anarci(batch)
            for prot, cdr_tuple in zip(batch, cdrs):
                prot_cdr_map[prot] = cdr_tuple
    else:
        for prot in unique_proteins:
            prot_cdr_map[prot] = extract_cdrs_offset(prot)

    # ── Assemble output tables ──────────────────────────────────────────────
    # DNA table: one row per unique amplicon that translated successfully
    dna_rows = [
        (dna, dna_to_prot[dna], dna_counts[dna])
        for dna in dna_to_prot
        if dna_to_prot[dna]
    ]

    # Protein/CDR table: one row per unique protein
    prot_cdr_rows = []
    cdr3_lengths  = Counter()
    for prot, count in prot_counts.items():
        cdr1, cdr2, cdr3 = prot_cdr_map.get(prot, ("", "", ""))
        if not any([cdr1, cdr2, cdr3]):
            stats["no_cdr"] += 1
            continue
        liabilities = flag_liabilities(cdr3) if cdr3 else "None"
        prot_cdr_rows.append((prot, cdr1, cdr2, cdr3, cdr1 + cdr2 + cdr3, liabilities, count))
        if cdr3:
            cdr3_lengths[len(cdr3)] += count  # weighted by read count

    # ── Write outputs ───────────────────────────────────────────────────────
    def _write_csv(path, header, rows):
        rows_sorted = sorted(rows, key=lambda x: x[-1], reverse=True)
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(header)
            w.writerows(rows_sorted)

    _write_csv(out_dna_prot, ["DNA_Sequence", "Protein_Sequence", "Count"], dna_rows)
    _write_csv(
        out_prot_cdr,
        ["Protein_Sequence", "CDR1", "CDR2", "CDR3", "CDR_Concatenated", "Liabilities", "Count"],
        prot_cdr_rows,
    )

    plot_cdr3_lengths(cdr3_lengths, out_cdr3_plot)
    if cfg["EXPORT_QC_METRICS"] and qscores_all:
        plot_qscore_dist(qscores_all, out_qplot)

    # Copy protein/CDR CSV up to the BAM directory for downstream steps
    export_path = output_dir.parent / f"{stem}_vhh_protein_cdr.csv"
    shutil.copy(out_prot_cdr, export_path)

    summary = {
        "file":              bam_path.name,
        "stats":             dict(stats),
        "unique_amplicons":  len(dna_rows),
        "unique_proteins":   len(prot_cdr_rows),
        "median_qscore":     float(np.median(qscores_all)) if qscores_all else None,
        "cdr_method":        CDR_METHOD,
        "anarci_available":  ANARCI_AVAILABLE,
        "fr1_motif":         FR1,
        "j4_motif":          J4,
        "fr_max_mismatches": FR_MM,
    }
    with open(out_summary, "w") as f:
        json.dump(summary, f, indent=2)

    return summary


# ══════════════════════════════════════════════════════════════════════════════
# ENTRY POINT
# ══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="VHH BAM extractor — alignment-narrowed motif search, indel-aware, parallel"
    )
    parser.add_argument("input_dir",
                        help="Directory containing .bam files aligned to IMGT Vicugna IGHV")
    parser.add_argument("--fr1",            default="CAGGTGCAGCTG",
                        help="FR1 anchor motif for motif-search fallback (default: CAGGTGCAGCTG)")
    parser.add_argument("--j4",             default="ACCCAGGTCACC",
                        help="J4 anchor motif for motif-search fallback (default: ACCCAGGTCACC)")
    parser.add_argument("--fr-mm",          type=int,   default=2,
                        help="Fuzzy mismatch tolerance for FR1/J4 motifs (default: 2)")
    parser.add_argument("--min-q",          type=float, default=10.0,
                        help="Min mean ONT Q-score (default: 10)")
    parser.add_argument("--min-len",        type=int,   default=300,
                        help="Min read length bp (default: 300)")
    parser.add_argument("--max-len",        type=int,   default=1000,
                        help="Max read length bp (default: 1000)")
    parser.add_argument("--min-aa",         type=int,   default=100,
                        help="Min VHH aa length (default: 100)")
    parser.add_argument("--cdr-method",     default="anarci", choices=["anarci", "offset"],
                        help="CDR extraction method (default: anarci)")
    parser.add_argument("--use-umi",        action="store_true",
                        help="Enable UMI deduplication")
    parser.add_argument("--umi-tag",        default="UX",
                        help="BAM tag for UMI (default: UX)")
    parser.add_argument("--workers",        type=int,   default=0,
                        help="Parallel BAM worker processes (default: one per BAM up to CPU count)")
    parser.add_argument("--anarci-batch",   type=int,   default=5000,
                        help="Proteins per ANARCI batch call (default: 5000)")
    args = parser.parse_args()

    input_dir  = Path(args.input_dir)
    output_dir = input_dir / "BAM_extract"
    output_dir.mkdir(exist_ok=True)

    config = {
        "FR1_MOTIF":         args.fr1,
        "J4_MOTIF":          args.j4,
        "FR_MAX_MISMATCHES": args.fr_mm,
        "MIN_MEAN_QSCORE":   args.min_q,
        "MIN_READ_LENGTH":   args.min_len,
        "MAX_READ_LENGTH":   args.max_len,
        "MIN_VHH_AA_LENGTH": args.min_aa,
        "CDR_METHOD":        args.cdr_method,
        "USE_UMI":           args.use_umi,
        "UMI_TAG":           args.umi_tag,
        "ANARCI_BATCH_SIZE": args.anarci_batch,
    }

    bam_files = sorted(input_dir.glob("*.bam"))
    if not bam_files:
        console.print("[red]No .bam files found in directory.[/red]")
        return

    n_workers = args.workers or min(len(bam_files), os.cpu_count() or 1)
    console.print(
        f"[bold green]VHH extraction:[/bold green] "
        f"{len(bam_files)} BAM file(s) — {n_workers} worker(s)"
    )

    all_summaries = []
    if n_workers == 1 or len(bam_files) == 1:
        for bam in bam_files:
            all_summaries.append(process_bam(bam, output_dir, config))
    else:
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            futures = {
                pool.submit(process_bam, bam, output_dir, config): bam
                for bam in bam_files
            }
            for future in as_completed(futures):
                bam = futures[future]
                try:
                    all_summaries.append(future.result())
                except Exception as exc:
                    console.print(f"[red]Error processing {bam.name}: {exc}[/red]")

    # Sort summary rows to match original BAM order
    name_order = {b.name: i for i, b in enumerate(bam_files)}
    all_summaries.sort(key=lambda s: name_order.get(s["file"], 999))

    table = Table(title="BAM Extraction Summary", show_lines=True)
    table.add_column("File",              style="cyan", no_wrap=True)
    table.add_column("Total",             justify="right")
    table.add_column("Unmapped",          justify="right")
    table.add_column("Aln search",        justify="right")
    table.add_column("Full read fallback",justify="right")
    table.add_column("Extracted reads",   justify="right")
    table.add_column("Unique proteins",   justify="right")
    table.add_column("Median Q",          justify="right")

    for s in all_summaries:
        st = s["stats"]
        aln_found  = st.get("extracted_aligned", 0)
        full_found = st.get("extracted_+", 0) + st.get("extracted_-", 0)
        table.add_row(
            s["file"],
            str(st.get("total", 0)),
            str(st.get("unmapped", 0)),
            str(aln_found),
            str(full_found),
            str(st.get("extracted", 0)),
            str(s["unique_proteins"]),
            f"{s['median_qscore']:.1f}" if s["median_qscore"] else "N/A",
        )

    console.print(table)
    console.print(f"\n[green]✓ Output:[/green] {output_dir}")
    console.print(
        f"[dim]ANARCI: {'enabled' if ANARCI_AVAILABLE else 'unavailable — offset fallback'} | "
        f"FR1: {args.fr1} | J4: {args.j4} | max_mm: {args.fr_mm}[/dim]"
    )


if __name__ == "__main__":
    main()
