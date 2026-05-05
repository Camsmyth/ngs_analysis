#!/usr/bin/env python3
"""
bam_extract.py — Step 1 of the VHH NGS pipeline.

Extracts VHH sequences from ONT BAM files aligned to IMGT Vicugna germlines.
Uses alignment as a pre-filter (mapped reads only), then anchors extraction
to conserved VHH framework motifs (FR1 start, J4 end) rather than PCR primer
sequences — making it robust to variation in library prep primers.

  • Alignment pre-filter: skips ~5% unmapped reads immediately
  • FR1/J4 anchor extraction: strand-agnostic, primer-independent
  • ANARCI-based CDR annotation using IMGT numbering
  • ONT-aware quality filtering (Q-score, min/max read length)
  • UMI-based deduplication (if UMI tags present in BAM)
  • PTM / liability site flagging (Asn-glycosylation, free Cys, deamidation)
  • Rich progress display and structured JSON summary
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
from typing import Optional, Tuple

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
    # VHH framework anchor motifs (primer-independent)
    # FR1: conserved camelid VHH framework 1 start (covers QVQL/EVQL/DVQL)
    # J4:  conserved J-region end (WGQGTQVTVSS)
    "FR1_MOTIF": "CAGGTGCAGCTG",
    "J4_MOTIF":  "ACCCAGGTCACC",
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
}

LIABILITY_PATTERNS = {
    "N-glycosylation":  re.compile(r"N[^P][ST]"),
    "Deamidation (NG)": re.compile(r"NG"),
    "Deamidation (NS)": re.compile(r"NS"),
    "Isomerization":    re.compile(r"D[TGSH]"),
    "Free Cys":         re.compile(r"C"),
    "Oxidation (Met)":  re.compile(r"M"),
}

# ══════════════════════════════════════════════════════════════════════════════
# UTILITIES
# ══════════════════════════════════════════════════════════════════════════════

def mean_qscore(qualities) -> float:
    if qualities is None or len(qualities) == 0:
        return 0.0
    error_probs = [10 ** (-q / 10) for q in qualities]
    mean_err = sum(error_probs) / len(error_probs)
    return -10 * np.log10(mean_err) if mean_err > 0 else 40.0


def fuzzy_search(pattern: str, sequence: str, max_mismatches: int):
    m = re.search(f"({pattern}){{e<={max_mismatches}}}", sequence)
    return (m.start(), m.end()) if m else None


def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def find_vhh_boundaries(seq: str, fr1: str, j4: str, max_mm: int) -> Optional[Tuple[str, str]]:
    """
    Locate the VHH coding region using conserved framework anchors.

    Searches both the read and its reverse complement. Returns
    (amplicon_dna, strand) where strand is '+' or '-', or None if not found.
    The amplicon runs from the FR1 start codon through the end of the J4 motif.
    """
    for s, strand in [(seq, "+"), (reverse_complement(seq), "-")]:
        fr1_m = fuzzy_search(fr1, s, max_mm)
        j4_m  = fuzzy_search(j4,  s, max_mm)
        if fr1_m and j4_m and j4_m[1] > fr1_m[0]:
            return s[fr1_m[0]:j4_m[1]], strand
    return None


def detect_and_correct_frameshift(dna_seq: str) -> str:
    """Try all three frames; return the VHH protein with a recognised FR1 start."""
    VHH_START_MOTIFS = ("QVQ", "VQL", "QLQ", "EVQ", "DVQ")
    best_prot = ""
    for frame in range(3):
        trimmed = dna_seq[frame:]
        trimmed = trimmed[:len(trimmed) - (len(trimmed) % 3)]
        try:
            prot = str(Seq(trimmed).translate(to_stop=False))
        except Exception:
            continue
        if "*" in prot[:-1]:
            continue
        if prot.startswith(VHH_START_MOTIFS):
            return prot
        if any(m in prot for m in VHH_START_MOTIFS) and len(prot) > len(best_prot):
            best_prot = prot
    return best_prot if best_prot and "*" not in best_prot[:-1] else ""


# ══════════════════════════════════════════════════════════════════════════════
# CDR EXTRACTION
# ══════════════════════════════════════════════════════════════════════════════

def extract_cdrs_anarci(protein: str) -> Tuple[str, str, str]:
    if not ANARCI_AVAILABLE:
        return extract_cdrs_offset(protein)
    try:
        results = anarci(
            [("seq", protein)],
            scheme="imgt",
            output=False,
            allow={"H"},
        )
        if not results or not results[0]:
            return ("", "", "")
        numbering, _, _ = results[0][0][0]
        pos_aa = {num: aa for (num, _), aa in numbering if aa != "-"}
        cdr1 = "".join(pos_aa.get(i, "") for i in range(27, 39))
        cdr2 = "".join(pos_aa.get(i, "") for i in range(56, 66))
        cdr3 = "".join(pos_aa.get(i, "") for i in range(105, 118))
        return cdr1, cdr2, cdr3
    except Exception:
        return extract_cdrs_offset(protein)


def extract_cdrs_offset(protein: str) -> Tuple[str, str, str]:
    try:
        cdr1 = protein[26:34]
        cdr2 = protein[50:60]
        raw_cdr3 = protein[99:122]
        m = re.search(r"W[GA][QKR]G", raw_cdr3)
        cdr3 = raw_cdr3[:m.start()] if m else raw_cdr3
        return cdr1, cdr2, cdr3
    except IndexError:
        return "", "", ""


def extract_cdrs(protein: str, method: str = "anarci") -> Tuple[str, str, str]:
    if method == "anarci" and ANARCI_AVAILABLE:
        return extract_cdrs_anarci(protein)
    return extract_cdrs_offset(protein)


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
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.hist(qscores, bins=50, color="#55A868", edgecolor="none", alpha=0.85)
    ax.axvline(np.median(qscores), color="red", linestyle="--",
               label=f"Median Q = {np.median(qscores):.1f}")
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
    bam_out_dir = output_dir / bam_path.stem
    bam_out_dir.mkdir(parents=True, exist_ok=True)

    stem = bam_path.stem
    out_dna_prot  = bam_out_dir / f"{stem}_vhh_dna_protein.csv"
    out_prot_cdr  = bam_out_dir / f"{stem}_vhh_protein_cdr.csv"
    out_cdr3_plot = bam_out_dir / f"{stem}_cdr3_lengths.png"
    out_qplot     = bam_out_dir / f"{stem}_qscore_dist.png"
    out_summary   = bam_out_dir / f"{stem}_summary.json"

    cfg        = {**DEFAULT_CONFIG, **config}
    FR1        = cfg["FR1_MOTIF"]
    J4         = cfg["J4_MOTIF"]
    FR_MM      = cfg["FR_MAX_MISMATCHES"]
    CDR_METHOD = cfg["CDR_METHOD"]

    seq_counts      = Counter()
    prot_cdr_counts = Counter()
    cdr3_lengths    = Counter()
    qscores_all     = []
    umi_seen        = set()
    stats           = defaultdict(int)

    console.print(f"[bold cyan]Processing:[/bold cyan] {bam_path.name}")
    with pysam.AlignmentFile(str(bam_path), check_sq=False) as bam_tmp:
        total_est = sum(1 for _ in bam_tmp.fetch(until_eof=True))

    with pysam.AlignmentFile(str(bam_path), check_sq=False) as bamfile:
        for read in tqdm(bamfile.fetch(until_eof=True), total=total_est,
                         desc=stem[:30], unit="reads"):
            stats["total"] += 1

            # ── Alignment pre-filter ───────────────────────────────────────
            # Unmapped reads have no VHH germline hit — skip immediately.
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

            # ── ONT Q-score filter ─────────────────────────────────────────
            qs = mean_qscore(read.query_qualities)
            qscores_all.append(qs)
            if qs < cfg["MIN_MEAN_QSCORE"]:
                stats["low_quality"] += 1
                continue

            # ── UMI deduplication ──────────────────────────────────────────
            if cfg["USE_UMI"]:
                umi = read.get_tag(cfg["UMI_TAG"]) if read.has_tag(cfg["UMI_TAG"]) else None
                if umi and umi in umi_seen:
                    stats["umi_duplicate"] += 1
                    continue
                if umi:
                    umi_seen.add(umi)

            # ── FR1 → J4 anchor extraction (strand-agnostic) ──────────────
            # Operates on full query_sequence (including soft-clips) so the
            # VHH region is always reachable regardless of alignment boundaries.
            result = find_vhh_boundaries(seq, FR1, J4, FR_MM)
            if result is None:
                stats["no_framework"] += 1
                continue

            amplicon, strand = result
            stats["framework_found"] += 1

            # ── Translation / frameshift correction ───────────────────────
            prot = detect_and_correct_frameshift(amplicon)
            if not prot or len(prot) < cfg["MIN_VHH_AA_LENGTH"]:
                stats["bad_translation"] += 1
                continue
            if prot.count("*") > cfg["MAX_INTERNAL_STOPS"]:
                stats["stop_codon"] += 1
                continue

            # ── CDR extraction ─────────────────────────────────────────────
            cdr1, cdr2, cdr3 = extract_cdrs(prot, CDR_METHOD)
            if not any([cdr1, cdr2, cdr3]):
                stats["no_cdr"] += 1
                continue

            liabilities = flag_liabilities(cdr3) if cdr3 else "None"

            seq_counts[(amplicon, prot)] += 1
            prot_cdr_counts[(prot, cdr1, cdr2, cdr3,
                             cdr1 + cdr2 + cdr3, liabilities)] += 1
            if cdr3:
                cdr3_lengths[len(cdr3)] += 1
            stats["extracted"] += 1

    # ── Write outputs ──────────────────────────────────────────────────────────
    def _write_csv(path, header, rows):
        rows_sorted = sorted(rows, key=lambda x: x[-1], reverse=True)
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(header)
            w.writerows(rows_sorted)

    dna_rows = [(dna, prot, cnt) for (dna, prot), cnt in seq_counts.items()]
    cdr_rows = [
        (prot, c1, c2, c3, concat, liab, cnt)
        for (prot, c1, c2, c3, concat, liab), cnt in prot_cdr_counts.items()
    ]

    _write_csv(out_dna_prot, ["DNA_Sequence", "Protein_Sequence", "Count"], dna_rows)
    _write_csv(out_prot_cdr,
               ["Protein_Sequence", "CDR1", "CDR2", "CDR3",
                "CDR_Concatenated", "Liabilities", "Count"],
               cdr_rows)

    plot_cdr3_lengths(cdr3_lengths, out_cdr3_plot)
    if cfg["EXPORT_QC_METRICS"] and qscores_all:
        plot_qscore_dist(qscores_all, out_qplot)

    export_path = output_dir.parent / f"{stem}_vhh_protein_cdr.csv"
    shutil.copy(out_prot_cdr, export_path)

    summary = {
        "file": bam_path.name,
        "stats": dict(stats),
        "unique_vhh_sequences": len(dna_rows),
        "unique_proteins": len(cdr_rows),
        "median_qscore": float(np.median(qscores_all)) if qscores_all else None,
        "cdr_method": CDR_METHOD,
        "anarci_available": ANARCI_AVAILABLE,
        "fr1_motif": FR1,
        "j4_motif": J4,
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
        description="VHH BAM extractor — alignment-filtered, framework-anchored"
    )
    parser.add_argument("input_dir",    help="Directory containing .bam files aligned to IMGT VHH germlines")
    parser.add_argument("--fr1",        default="CAGGTGCAGCTG",  help="FR1 anchor motif (default: CAGGTGCAGCTG)")
    parser.add_argument("--j4",         default="ACCCAGGTCACC",  help="J4 anchor motif  (default: ACCCAGGTCACC)")
    parser.add_argument("--fr-mm",      type=int, default=2,     help="Max mismatches for framework motifs (default: 2)")
    parser.add_argument("--min-q",      type=float, default=10.0, help="Min mean ONT Q-score (default: 10)")
    parser.add_argument("--min-len",    type=int, default=300,    help="Min read length bp (default: 300)")
    parser.add_argument("--max-len",    type=int, default=1000,   help="Max read length bp (default: 1000)")
    parser.add_argument("--min-aa",     type=int, default=100,    help="Min VHH aa length (default: 100)")
    parser.add_argument("--cdr-method", default="anarci", choices=["anarci", "offset"],
                        help="CDR extraction method (default: anarci)")
    parser.add_argument("--use-umi",    action="store_true",      help="Enable UMI deduplication")
    parser.add_argument("--umi-tag",    default="UX",             help="BAM tag for UMI (default: UX)")
    args = parser.parse_args()

    input_dir  = Path(args.input_dir)
    output_dir = input_dir / "BAM_extract"
    output_dir.mkdir(exist_ok=True)

    config = {
        "FR1_MOTIF":        args.fr1,
        "J4_MOTIF":         args.j4,
        "FR_MAX_MISMATCHES": args.fr_mm,
        "MIN_MEAN_QSCORE":  args.min_q,
        "MIN_READ_LENGTH":  args.min_len,
        "MAX_READ_LENGTH":  args.max_len,
        "MIN_VHH_AA_LENGTH": args.min_aa,
        "CDR_METHOD":       args.cdr_method,
        "USE_UMI":          args.use_umi,
        "UMI_TAG":          args.umi_tag,
    }

    bam_files = list(input_dir.glob("*.bam"))
    if not bam_files:
        console.print("[red]No .bam files found in directory.[/red]")
        return

    all_summaries = []
    for bam_path in bam_files:
        all_summaries.append(process_bam(bam_path, output_dir, config))

    table = Table(title="BAM Extraction Summary", show_lines=True)
    table.add_column("File",             style="cyan", no_wrap=True)
    table.add_column("Total Reads",      justify="right")
    table.add_column("Unmapped",         justify="right")
    table.add_column("Framework Found",  justify="right")
    table.add_column("Extracted VHHs",   justify="right")
    table.add_column("Unique Proteins",  justify="right")
    table.add_column("Median Q",         justify="right")

    for s in all_summaries:
        st = s["stats"]
        table.add_row(
            s["file"],
            str(st.get("total", 0)),
            str(st.get("unmapped", 0)),
            str(st.get("framework_found", 0)),
            str(st.get("extracted", 0)),
            str(s["unique_proteins"]),
            f"{s['median_qscore']:.1f}" if s["median_qscore"] else "N/A",
        )

    console.print(table)
    console.print(f"\n[green]✓ Output written to:[/green] {output_dir}")
    console.print(
        f"[dim]ANARCI: {'enabled' if ANARCI_AVAILABLE else 'unavailable — using offset fallback'} | "
        f"FR1: {args.fr1} | J4: {args.j4} | max_mm: {args.fr_mm}[/dim]"
    )


if __name__ == "__main__":
    main()
