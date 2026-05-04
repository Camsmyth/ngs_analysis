#!/usr/bin/env python3
"""
bam_extract.py — Step 1 of the VHH NGS pipeline.

Redesigned for Oxford Nanopore (PromethION) BAM files from Dorado basecalling
with alpaca VHH genome alignment. Key improvements over original:

  • ONT-aware quality filtering (Q-score, min read length)
  • Both forward AND reverse-complement primer search (strand-agnostic)
  • ANARCI-based CDR annotation using IMGT numbering (replaces fixed-offset slicing)
  • UMI-based deduplication (if UMI tags present in BAM)
  • Per-read quality metrics exported for QC
  • PTM / liability site flagging (Asn-glycosylation, free Cys, deamidation)
  • Rich progress display and structured JSON summary
  • Configurable via CLI flags or config dict
"""

import pysam
import regex as re
import csv
import json
import os
import shutil
import time
import logging
import argparse
from collections import defaultdict, Counter
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

# ── Try to import ANARCI for numbered CDR extraction ──────────────────────────
try:
    from anarci import anarci
    ANARCI_AVAILABLE = True
except ImportError:
    ANARCI_AVAILABLE = False

console = Console()

# ══════════════════════════════════════════════════════════════════════════════
# DEFAULT CONFIGURATION  (override via CLI or pass config dict to process_bam)
# ══════════════════════════════════════════════════════════════════════════════
DEFAULT_CONFIG = {
    # Primer sequences (alpaca VHH amplicon)
    "VHH_F": "CAGGTACAGCTGCA",
    "VHH_R": "CGGTGTCTAGCACT",

    # Fuzzy matching
    "MAX_MISMATCHES": 1,

    # ONT read quality gates
    "MIN_READ_LENGTH": 300,        # bp — discard very short reads
    "MAX_READ_LENGTH": 1000,       # bp — discard chimeric/concatenated reads
    "MIN_MEAN_QSCORE": 10.0,       # Phred Q-score threshold (Q10 = 90% acc)

    # VHH protein quality gates
    "MIN_VHH_AA_LENGTH": 100,      # minimum translated VHH length
    "MAX_INTERNAL_STOPS": 0,       # tolerated internal stop codons

    # CDR extraction method: "anarci" (preferred) or "offset" (fallback)
    "CDR_METHOD": "anarci",

    # UMI deduplication
    "USE_UMI": False,              # set True if UMI tags present (umi:Z: in BAM)
    "UMI_TAG": "UX",               # BAM tag for UMI (Dorado: UX; custom: MI)

    # Output
    "PROGRESS_UPDATE_INTERVAL": 3,
    "EXPORT_QC_METRICS": True,
}

# ── Liability motifs to flag ───────────────────────────────────────────────────
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
    """Convert Phred quality array to mean Q-score (error-probability space)."""
    if qualities is None or len(qualities) == 0:
        return 0.0
    error_probs = [10 ** (-q / 10) for q in qualities]
    mean_err = sum(error_probs) / len(error_probs)
    return -10 * np.log10(mean_err) if mean_err > 0 else 40.0


def fuzzy_search(pattern: str, sequence: str, max_mismatches: int = 1):
    """Fuzzy regex search; returns (start, end) or None."""
    m = re.search(f"({pattern}){{e<={max_mismatches}}}", sequence)
    return (m.start(), m.end()) if m else None


def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def find_primers_both_strands(seq: str, fwd: str, rev: str, max_mm: int):
    """
    Search for primer pair on both orientations.
    Returns (start, end, orientation) or None.
    Orientation: '+' means fwd primer leads; '-' means found on RC strand.
    """
    # Forward orientation
    fwd_m = fuzzy_search(fwd, seq, max_mm)
    rev_m = fuzzy_search(rev, seq, max_mm)
    if fwd_m and rev_m and rev_m[1] > fwd_m[0]:
        return fwd_m[0], rev_m[1], "+"

    # Try reverse-complement orientation (ONT reads are strand-random)
    rc_seq = reverse_complement(seq)
    fwd_m2 = fuzzy_search(fwd, rc_seq, max_mm)
    rev_m2 = fuzzy_search(rev, rc_seq, max_mm)
    if fwd_m2 and rev_m2 and rev_m2[1] > fwd_m2[0]:
        return fwd_m2[0], rev_m2[1], "-"

    return None


def mask_indels(sequence: str, read) -> Tuple[str, bool]:
    """Mask insertions and mark deletions as N using CIGAR string."""
    if read.cigartuples is None:
        return sequence, False
    masked_seq, seq_pos, was_masked = [], 0, False
    for op, length in read.cigartuples:
        if op == 0:   # M
            masked_seq.extend(sequence[seq_pos:seq_pos + length])
            seq_pos += length
        elif op == 1: # I
            seq_pos += length
            was_masked = True
        elif op == 2: # D
            masked_seq.extend(["N"] * length)
            was_masked = True
        elif op in (4, 5):  # S/H clips
            seq_pos += length
    return "".join(masked_seq), was_masked


def detect_and_correct_frameshift(dna_seq: str) -> str:
    """Try all three frames; return best VHH protein (starts with known VHH motif or longest valid)."""
    VHH_START_MOTIFS = ("QVQ", "VQL", "QLQ", "EVQ", "DVQ")  # expanded motif list
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
# CDR EXTRACTION — two backends
# ══════════════════════════════════════════════════════════════════════════════

def extract_cdrs_anarci(protein: str) -> Tuple[str, str, str]:
    """
    Use ANARCI with IMGT numbering to extract CDR1 (27–38), CDR2 (56–65),
    CDR3 (105–117) for VHH sequences. Returns ("","","") on failure.
    """
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

        # Build position→AA dict
        pos_aa = {num: aa for (num, _), aa in numbering if aa != "-"}

        # IMGT CDR definitions for VHH (nanobody heavy chain)
        cdr1 = "".join(pos_aa.get(i, "") for i in range(27, 39))
        cdr2 = "".join(pos_aa.get(i, "") for i in range(56, 66))
        cdr3 = "".join(pos_aa.get(i, "") for i in range(105, 118))
        return cdr1, cdr2, cdr3
    except Exception:
        return extract_cdrs_offset(protein)


def extract_cdrs_offset(protein: str) -> Tuple[str, str, str]:
    """
    Fallback fixed-offset CDR extraction.
    Improved CDR3 trimming: searches for conserved WGxG motif.
    """
    try:
        cdr1 = protein[26:34]
        cdr2 = protein[50:60]
        raw_cdr3 = protein[99:122]
        # Trim CDR3 at WGxG (conserved VHH framework 4 boundary)
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
# LIABILITY / PTM FLAGGING
# ══════════════════════════════════════════════════════════════════════════════

def flag_liabilities(cdr3: str) -> str:
    """Return comma-separated liability flags found in CDR3."""
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
    ax.hist(lengths, bins=range(min(lengths), max(lengths) + 2), edgecolor="black",
            color="#4C72B0", alpha=0.85)
    ax.set_xlabel("CDR3 Length (aa)", fontsize=12)
    ax.set_ylabel("Frequency", fontsize=12)
    ax.set_title("CDR3 Length Distribution", fontsize=14)
    ax.axvline(np.mean(lengths), color="red", linestyle="--", label=f"Mean = {np.mean(lengths):.1f}")
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
    """
    Process a single BAM file and write per-sample output CSVs.
    Returns a summary dict.
    """
    bam_out_dir = output_dir / bam_path.stem
    bam_out_dir.mkdir(parents=True, exist_ok=True)

    stem = bam_path.stem
    out_dna_prot  = bam_out_dir / f"{stem}_vhh_dna_protein.csv"
    out_prot_cdr  = bam_out_dir / f"{stem}_vhh_protein_cdr.csv"
    out_cdr3_plot = bam_out_dir / f"{stem}_cdr3_lengths.png"
    out_qplot     = bam_out_dir / f"{stem}_qscore_dist.png"
    out_summary   = bam_out_dir / f"{stem}_summary.json"

    cfg = {**DEFAULT_CONFIG, **config}
    VHH_F, VHH_R = cfg["VHH_F"], cfg["VHH_R"]
    MAX_MM = cfg["MAX_MISMATCHES"]
    CDR_METHOD = cfg["CDR_METHOD"]

    seq_counts      = Counter()
    prot_cdr_counts = Counter()
    cdr3_lengths    = Counter()
    qscores_all     = []
    umi_seen        = set()

    stats = defaultdict(int)

    # Count total reads for progress bar
    console.print(f"[bold cyan]Processing:[/bold cyan] {bam_path.name}")
    with pysam.AlignmentFile(str(bam_path), check_sq=False) as bam_tmp:
        total_est = sum(1 for _ in bam_tmp.fetch(until_eof=True))

    with pysam.AlignmentFile(str(bam_path), check_sq=False) as bamfile:
        for read in tqdm(bamfile.fetch(until_eof=True), total=total_est,
                         desc=stem[:30], unit="reads"):
            stats["total"] += 1

            # ── Basic filters ──────────────────────────────────────────────
            if read.is_unmapped or read.query_sequence is None:
                stats["unmapped"] += 1
                continue

            seq = read.query_sequence.upper()
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

            # ── Primer search (strand-agnostic) ────────────────────────────
            hit = find_primers_both_strands(seq, VHH_F, VHH_R, MAX_MM)
            if hit is None:
                stats["no_primer"] += 1
                continue

            start, end, strand = hit
            working_seq = seq if strand == "+" else reverse_complement(seq)
            # Re-search on oriented sequence for exact boundaries
            fwd_m = fuzzy_search(VHH_F, working_seq, MAX_MM)
            rev_m = fuzzy_search(VHH_R, working_seq, MAX_MM)
            if not (fwd_m and rev_m and rev_m[1] > fwd_m[0]):
                stats["no_primer"] += 1
                continue

            amplicon = working_seq[fwd_m[0]:rev_m[1]]
            stats["primer_found"] += 1

            # ── Indel masking ──────────────────────────────────────────────
            masked_seq, was_masked = mask_indels(amplicon, read) if strand == "+" \
                else (amplicon, False)  # can't use CIGAR on RC

            # ── Translation / frameshift correction ───────────────────────
            prot = detect_and_correct_frameshift(masked_seq)
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
            tag = "Masked" if was_masked else "Clean"

            seq_counts[(masked_seq, prot, tag)] += 1
            prot_cdr_counts[(prot, cdr1, cdr2, cdr3,
                             cdr1 + cdr2 + cdr3, liabilities, tag)] += 1
            if cdr3:
                cdr3_lengths[len(cdr3)] += 1
            stats["extracted"] += 1

    # ── Write outputs ──────────────────────────────────────────────────────────
    def _write_csv(path, header, rows):
        rows_sorted = sorted(rows, key=lambda x: x[-1], reverse=True)
        with open(path, "w", newline="") as f:
            csv.writer(f).writerow(header)
            csv.writer(f).writerows(rows_sorted)

    dna_rows = [
        (dna, prot, tag, cnt)
        for (dna, prot, tag), cnt in seq_counts.items()
    ]
    cdr_rows = [
        (prot, c1, c2, c3, concat, liab, tag, cnt)
        for (prot, c1, c2, c3, concat, liab, tag), cnt in prot_cdr_counts.items()
    ]

    _write_csv(out_dna_prot, ["DNA_Sequence", "Protein_Sequence", "Status", "Count"], dna_rows)
    _write_csv(out_prot_cdr,
               ["Protein_Sequence", "CDR1", "CDR2", "CDR3",
                "CDR_Concatenated", "Liabilities", "Status", "Count"],
               cdr_rows)

    plot_cdr3_lengths(cdr3_lengths, out_cdr3_plot)
    if cfg["EXPORT_QC_METRICS"] and qscores_all:
        plot_qscore_dist(qscores_all, out_qplot)

    # ── Copy protein_cdr CSV to parent dir (for clustering step) ─────────────
    export_path = output_dir.parent / f"{stem}_vhh_protein_cdr.csv"
    shutil.copy(out_prot_cdr, export_path)

    # ── JSON summary ───────────────────────────────────────────────────────────
    summary = {
        "file": bam_path.name,
        "stats": dict(stats),
        "unique_vhh_sequences": len(dna_rows),
        "unique_proteins": len(cdr_rows),
        "median_qscore": float(np.median(qscores_all)) if qscores_all else None,
        "cdr_method": CDR_METHOD,
        "anarci_available": ANARCI_AVAILABLE,
    }
    with open(out_summary, "w") as f:
        json.dump(summary, f, indent=2)

    return summary


# ══════════════════════════════════════════════════════════════════════════════
# ENTRY POINT
# ══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="VHH BAM extractor — ONT-optimised with ANARCI CDR numbering"
    )
    parser.add_argument("input_dir", help="Directory containing .bam files")
    parser.add_argument("--max-mm",       type=int,   default=1,       help="Max primer mismatches (default: 1)")
    parser.add_argument("--min-q",        type=float, default=10.0,    help="Min mean ONT Q-score (default: 10)")
    parser.add_argument("--min-len",      type=int,   default=300,     help="Min read length bp (default: 300)")
    parser.add_argument("--max-len",      type=int,   default=1000,    help="Max read length bp (default: 1000)")
    parser.add_argument("--min-aa",       type=int,   default=100,     help="Min VHH aa length (default: 100)")
    parser.add_argument("--cdr-method",   default="anarci",
                        choices=["anarci", "offset"],                  help="CDR extraction method")
    parser.add_argument("--use-umi",      action="store_true",         help="Enable UMI deduplication")
    parser.add_argument("--umi-tag",      default="UX",                help="BAM tag for UMI (default: UX)")
    parser.add_argument("--vhh-f",        default="CAGGTACAGCTGCA",    help="Forward primer sequence")
    parser.add_argument("--vhh-r",        default="CGGTGTCTAGCACT",    help="Reverse primer sequence")
    args = parser.parse_args()

    input_dir  = Path(args.input_dir)
    output_dir = input_dir / "BAM_extract"
    output_dir.mkdir(exist_ok=True)

    config = {
        "MAX_MISMATCHES":   args.max_mm,
        "MIN_MEAN_QSCORE":  args.min_q,
        "MIN_READ_LENGTH":  args.min_len,
        "MAX_READ_LENGTH":  args.max_len,
        "MIN_VHH_AA_LENGTH": args.min_aa,
        "CDR_METHOD":       args.cdr_method,
        "USE_UMI":          args.use_umi,
        "UMI_TAG":          args.umi_tag,
        "VHH_F":            args.vhh_f,
        "VHH_R":            args.vhh_r,
    }

    all_summaries = []
    bam_files = list(input_dir.glob("*.bam"))
    if not bam_files:
        console.print("[red]No .bam files found in directory.[/red]")
        return

    for bam_path in bam_files:
        summary = process_bam(bam_path, output_dir, config)
        all_summaries.append(summary)

    # ── Print summary table ────────────────────────────────────────────────────
    table = Table(title="BAM Extraction Summary", show_lines=True)
    table.add_column("File", style="cyan", no_wrap=True)
    table.add_column("Total Reads", justify="right")
    table.add_column("Primer Found", justify="right")
    table.add_column("Extracted VHHs", justify="right")
    table.add_column("Unique Proteins", justify="right")
    table.add_column("Median Q", justify="right")

    for s in all_summaries:
        st = s["stats"]
        table.add_row(
            s["file"],
            str(st.get("total", 0)),
            str(st.get("primer_found", 0)),
            str(st.get("extracted", 0)),
            str(s["unique_proteins"]),
            f"{s['median_qscore']:.1f}" if s["median_qscore"] else "N/A",
        )
    console.print(table)
    console.print(f"\n[green]✓ Output written to:[/green] {output_dir}")
    console.print(f"[dim]ANARCI CDR numbering: {'enabled' if ANARCI_AVAILABLE else 'unavailable — using offset fallback'}[/dim]")


if __name__ == "__main__":
    main()
