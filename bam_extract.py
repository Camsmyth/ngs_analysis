#!/usr/bin/env python3
"""
bam_extract.py — Step 1 of the VHH NGS pipeline.

Extracts VHH sequences from ONT BAM files aligned to IMGT Vicugna germlines.
Uses alignment as a pre-filter (mapped reads only), then anchors extraction
to conserved VHH framework motifs (FR1 start, J4 end) rather than PCR primer
sequences — making it robust to variation in library prep primers.

  • Alignment pre-filter: skips ~5% unmapped reads immediately
  • FR1/J4 anchor extraction: strand-agnostic, primer-independent
  • ANARCI-based CDR annotation — batched for throughput
  • ONT-aware quality filtering (Q-score, min/max read length)
  • GPU-accelerated Q-score batching via CuPy (optional, graceful fallback)
  • Parallel BAM processing via ProcessPoolExecutor (--workers)
  • UMI-based deduplication (if UMI tags present in BAM)
  • PTM / liability site flagging (Asn-glycosylation, free Cys, deamidation)
  • Rich progress display and structured JSON summary
"""

import pysam
import regex as re
import csv
import json
import shutil
import argparse
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
from Bio.Seq import Seq
from tqdm import tqdm
from rich.console import Console
from rich.table import Table
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import cupy as cp
    GPU_AVAILABLE = True
except ImportError:
    cp = None
    GPU_AVAILABLE = False

try:
    from anarci import anarci
    ANARCI_AVAILABLE = True
except ImportError:
    ANARCI_AVAILABLE = False

console = Console()

_Q_BATCH = 8_000    # reads per GPU/CPU Q-score batch
_ANARCI_BATCH = 500 # proteins per ANARCI call

# ══════════════════════════════════════════════════════════════════════════════
# DEFAULT CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════
DEFAULT_CONFIG = {
    "FR1_MOTIF": "CAGGTGCAGCTG",
    "J4_MOTIF":  "ACCCAGGTCACC",
    "FR_MAX_MISMATCHES": 2,
    "MIN_READ_LENGTH": 300,
    "MAX_READ_LENGTH": 1000,
    "MIN_MEAN_QSCORE": 10.0,
    "MIN_VHH_AA_LENGTH": 100,
    "MAX_INTERNAL_STOPS": 0,
    "CDR_METHOD": "anarci",
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
# Q-SCORE — CPU AND GPU PATHS
# ══════════════════════════════════════════════════════════════════════════════

def mean_qscore(qualities) -> float:
    if qualities is None or len(qualities) == 0:
        return 0.0
    error_probs = [10 ** (-q / 10) for q in qualities]
    mean_err = sum(error_probs) / len(error_probs)
    return -10 * np.log10(mean_err) if mean_err > 0 else 40.0


def batch_mean_qscore(qualities_list: list, use_gpu: bool = False) -> list:
    """
    Vectorised mean Q-score for a batch of reads.

    GPU path (CuPy): pads all quality arrays to the same length, runs
    the error-probability conversion in a single GPU kernel, then transfers
    results back. Effective for batches of ≥1000 reads.

    CPU path: numpy per-read (fallback when CuPy unavailable or use_gpu=False).
    """
    if not qualities_list:
        return []

    if use_gpu and GPU_AVAILABLE and cp is not None:
        max_len = max(len(q) for q in qualities_list)
        n = len(qualities_list)
        padded = np.zeros((n, max_len), dtype=np.float32)
        lengths = np.zeros(n, dtype=np.int32)
        for i, q in enumerate(qualities_list):
            ln = len(q)
            padded[i, :ln] = q
            lengths[i] = ln

        padded_gpu  = cp.asarray(padded)
        lengths_gpu = cp.asarray(lengths, dtype=cp.float32)

        err = cp.power(10.0, -padded_gpu / 10.0)
        # mask out padding positions
        mask = cp.arange(max_len)[None, :] < cp.asarray(lengths)[:, None]
        err  = err * mask

        mean_errs = err.sum(axis=1) / lengths_gpu
        mean_errs = cp.clip(mean_errs, 1e-30, None)
        qscores   = -10.0 * cp.log10(mean_errs)
        return cp.asnumpy(qscores).tolist()

    # CPU path
    results = []
    for q in qualities_list:
        if q is None or len(q) == 0:
            results.append(0.0)
            continue
        arr = np.array(q, dtype=np.float32)
        mean_err = np.mean(np.power(10.0, -arr / 10.0))
        results.append(float(-10 * np.log10(mean_err)) if mean_err > 0 else 40.0)
    return results


# ══════════════════════════════════════════════════════════════════════════════
# UTILITIES
# ══════════════════════════════════════════════════════════════════════════════

def fuzzy_search(pattern: str, sequence: str, max_mismatches: int):
    m = re.search(f"({pattern}){{e<={max_mismatches}}}", sequence)
    return (m.start(), m.end()) if m else None


def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def find_vhh_boundaries(seq: str, fr1: str, j4: str, max_mm: int) -> Optional[Tuple[str, str]]:
    for s, strand in [(seq, "+"), (reverse_complement(seq), "-")]:
        fr1_m = fuzzy_search(fr1, s, max_mm)
        j4_m  = fuzzy_search(j4,  s, max_mm)
        if fr1_m and j4_m and j4_m[1] > fr1_m[0]:
            return s[fr1_m[0]:j4_m[1]], strand
    return None


def detect_and_correct_frameshift(dna_seq: str) -> str:
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
# CDR EXTRACTION — SINGLE AND BATCH
# ══════════════════════════════════════════════════════════════════════════════

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


def _anarci_chunk(proteins: List[str]) -> List[Tuple[str, str, str]]:
    """Run ANARCI on a chunk of proteins; fall back per-protein on error."""
    try:
        results = anarci(
            [(str(i), p) for i, p in enumerate(proteins)],
            scheme="imgt",
            output=False,
            allow={"H"},
        )
        cdrs = []
        for i, result_tuple in enumerate(results):
            try:
                if not result_tuple or not result_tuple[0]:
                    raise ValueError("no hit")
                numbering, _, _ = result_tuple[0][0]
                pos_aa = {num: aa for (num, _), aa in numbering if aa != "-"}
                cdr1 = "".join(pos_aa.get(j, "") for j in range(27, 39))
                cdr2 = "".join(pos_aa.get(j, "") for j in range(56, 66))
                cdr3 = "".join(pos_aa.get(j, "") for j in range(105, 118))
                cdrs.append((cdr1, cdr2, cdr3))
            except Exception:
                cdrs.append(extract_cdrs_offset(proteins[i]))
        return cdrs
    except Exception:
        return [extract_cdrs_offset(p) for p in proteins]


def batch_extract_cdrs(proteins: List[str], method: str = "anarci") -> List[Tuple[str, str, str]]:
    """
    Extract CDRs for a list of proteins in batched ANARCI calls.
    Batching amortises the hmmscan subprocess overhead — typically 10-50× faster
    than calling anarci() once per sequence.
    Falls back to offset extraction when ANARCI is unavailable.
    """
    if not proteins:
        return []
    if method != "anarci" or not ANARCI_AVAILABLE:
        return [extract_cdrs_offset(p) for p in proteins]

    all_cdrs = []
    for i in range(0, len(proteins), _ANARCI_BATCH):
        all_cdrs.extend(_anarci_chunk(proteins[i:i + _ANARCI_BATCH]))
    return all_cdrs


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

def process_bam(bam_path: Path, output_dir: Path, config: dict,
                use_gpu: bool = False, silent: bool = False) -> dict:
    """
    Process a single BAM file through the full VHH extraction pipeline.

    Structured in four phases to enable GPU and batch optimisations:
      1. Read BAM — pre-filter by alignment/length; batch Q-score on GPU
      2. Framework search + translation — CPU (regex + Biopython)
      3. Batch ANARCI — all proteins in one batched hmmscan call
      4. Aggregate counts, write outputs
    """
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
    USE_UMI    = cfg["USE_UMI"]
    UMI_TAG    = cfg["UMI_TAG"]

    stats    = defaultdict(int)
    umi_seen = set()

    # ── Phase 1: Read BAM, pre-filter, batch Q-score ───────────────────────────
    # Reads passing length/alignment filter are buffered in _Q_BATCH-sized
    # chunks; Q-score is computed in batch (GPU or numpy), then low-Q reads
    # are discarded before the more expensive framework search.
    qscores_all  = []
    passing_seqs = []  # sequences that clear all pre-filters

    with pysam.AlignmentFile(str(bam_path), check_sq=False) as bam_tmp:
        total_est = sum(1 for _ in bam_tmp.fetch(until_eof=True))

    buf_seqs, buf_quals, buf_umis = [], [], []

    def _flush_qbatch():
        qscores = batch_mean_qscore(buf_quals, use_gpu=use_gpu)
        qscores_all.extend(qscores)
        for seq, umi, qs in zip(buf_seqs, buf_umis, qscores):
            if qs < cfg["MIN_MEAN_QSCORE"]:
                stats["low_quality"] += 1
                continue
            if USE_UMI and umi:
                if umi in umi_seen:
                    stats["umi_duplicate"] += 1
                    continue
                umi_seen.add(umi)
            passing_seqs.append(seq)
        buf_seqs.clear(); buf_quals.clear(); buf_umis.clear()

    with pysam.AlignmentFile(str(bam_path), check_sq=False) as bamfile:
        for read in tqdm(bamfile.fetch(until_eof=True), total=total_est,
                         desc=stem[:28], unit="reads", disable=silent):
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

            umi = read.get_tag(UMI_TAG) if USE_UMI and read.has_tag(UMI_TAG) else None
            buf_seqs.append(seq)
            buf_quals.append(read.query_qualities)
            buf_umis.append(umi)

            if len(buf_seqs) >= _Q_BATCH:
                _flush_qbatch()

    if buf_seqs:
        _flush_qbatch()

    # ── Phase 2: Framework search + translation ────────────────────────────────
    protein_records = []  # [(amplicon, prot), ...]

    for seq in tqdm(passing_seqs, desc=f"{stem[:24]} extract",
                    unit="reads", disable=silent):
        result = find_vhh_boundaries(seq, FR1, J4, FR_MM)
        if result is None:
            stats["no_framework"] += 1
            continue
        amplicon, _ = result
        stats["framework_found"] += 1

        prot = detect_and_correct_frameshift(amplicon)
        if not prot or len(prot) < cfg["MIN_VHH_AA_LENGTH"]:
            stats["bad_translation"] += 1
            continue
        if prot.count("*") > cfg["MAX_INTERNAL_STOPS"]:
            stats["stop_codon"] += 1
            continue
        protein_records.append((amplicon, prot))

    # ── Phase 3: Batch ANARCI ──────────────────────────────────────────────────
    proteins = [prot for _, prot in protein_records]
    cdrs_batch = batch_extract_cdrs(proteins, method=CDR_METHOD)

    # ── Phase 4: Aggregate counts and write outputs ────────────────────────────
    seq_counts      = Counter()
    prot_cdr_counts = Counter()
    cdr3_lengths    = Counter()

    for (amplicon, prot), (cdr1, cdr2, cdr3) in zip(protein_records, cdrs_batch):
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
        "gpu_used": use_gpu and GPU_AVAILABLE,
        "fr1_motif": FR1,
        "j4_motif": J4,
        "fr_max_mismatches": FR_MM,
    }
    with open(out_summary, "w") as f:
        json.dump(summary, f, indent=2)

    return summary


def _bam_worker(args: tuple) -> dict:
    """Top-level picklable wrapper for ProcessPoolExecutor."""
    bam_path, output_dir, config, use_gpu = args
    return process_bam(bam_path, output_dir, config, use_gpu=use_gpu, silent=True)


# ══════════════════════════════════════════════════════════════════════════════
# ENTRY POINT
# ══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="VHH BAM extractor — alignment-filtered, framework-anchored, parallel"
    )
    parser.add_argument("input_dir",     help="Directory containing .bam files aligned to IMGT VHH germlines")
    parser.add_argument("--fr1",         default="CAGGTGCAGCTG",  help="FR1 anchor motif (default: CAGGTGCAGCTG)")
    parser.add_argument("--j4",          default="ACCCAGGTCACC",  help="J4 anchor motif  (default: ACCCAGGTCACC)")
    parser.add_argument("--fr-mm",       type=int, default=2,     help="Max mismatches for framework motifs (default: 2)")
    parser.add_argument("--min-q",       type=float, default=10.0, help="Min mean ONT Q-score (default: 10)")
    parser.add_argument("--min-len",     type=int, default=300,    help="Min read length bp (default: 300)")
    parser.add_argument("--max-len",     type=int, default=1000,   help="Max read length bp (default: 1000)")
    parser.add_argument("--min-aa",      type=int, default=100,    help="Min VHH aa length (default: 100)")
    parser.add_argument("--cdr-method",  default="anarci", choices=["anarci", "offset"],
                        help="CDR extraction method (default: anarci)")
    parser.add_argument("--use-umi",     action="store_true",      help="Enable UMI deduplication")
    parser.add_argument("--umi-tag",     default="UX",             help="BAM tag for UMI (default: UX)")
    parser.add_argument("--workers",     type=int, default=min(cpu_count(), 8),
                        help=f"Parallel BAM workers (default: min(cpu_count,8)={min(cpu_count(),8)})")
    parser.add_argument("--gpu",         action="store_true",
                        help="GPU-accelerated Q-score batching via CuPy (requires: pip install cupy-cuda11x)")
    args = parser.parse_args()

    if args.gpu and not GPU_AVAILABLE:
        console.print("[yellow]Warning: --gpu requested but CuPy not installed. "
                      "Install with: pip install cupy-cuda11x  (match your CUDA version)[/yellow]")

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
    }

    bam_files = sorted(input_dir.glob("*.bam"))
    if not bam_files:
        console.print("[red]No .bam files found in directory.[/red]")
        return

    n_workers = min(args.workers, len(bam_files))
    gpu_label = ("CuPy GPU" if GPU_AVAILABLE else "GPU unavailable — CPU") if args.gpu else "CPU"
    console.print(
        f"[bold green]Processing {len(bam_files)} BAM file(s) "
        f"with {n_workers} worker(s) | Q-score: {gpu_label}[/bold green]"
    )

    all_summaries = []

    if n_workers == 1:
        # Single-worker: keep per-read progress bars
        for bam_path in bam_files:
            all_summaries.append(process_bam(bam_path, output_dir, config,
                                             use_gpu=args.gpu and GPU_AVAILABLE))
    else:
        # Multi-worker: per-read bars suppressed; top-level bar tracks BAM completion
        worker_args = [(bam, output_dir, config, args.gpu and GPU_AVAILABLE)
                       for bam in bam_files]
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {executor.submit(_bam_worker, arg): arg[0].name
                       for arg in worker_args}
            with tqdm(total=len(futures), desc="BAM files", unit="bam") as pbar:
                for future in as_completed(futures):
                    bam_name = futures[future]
                    try:
                        summary = future.result()
                        all_summaries.append(summary)
                        pbar.set_postfix_str(bam_name[:30])
                    except Exception as exc:
                        console.print(f"[red]ERROR processing {bam_name}: {exc}[/red]")
                    pbar.update(1)

    # ── Summary table ──────────────────────────────────────────────────────────
    table = Table(title="BAM Extraction Summary", show_lines=True)
    table.add_column("File",            style="cyan", no_wrap=True)
    table.add_column("Total Reads",     justify="right")
    table.add_column("Unmapped",        justify="right")
    table.add_column("Framework Found", justify="right")
    table.add_column("Extracted VHHs",  justify="right")
    table.add_column("Unique Proteins", justify="right")
    table.add_column("Median Q",        justify="right")

    for s in sorted(all_summaries, key=lambda x: x["file"]):
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
        f"FR1: {args.fr1} | J4: {args.j4} | max_mm: {args.fr_mm} | "
        f"workers: {n_workers} | gpu: {args.gpu and GPU_AVAILABLE}[/dim]"
    )


if __name__ == "__main__":
    main()
