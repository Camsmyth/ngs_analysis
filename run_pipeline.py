#!/usr/bin/env python3
"""
run_pipeline.py — Orchestrator for the full VHH NGS enrichment pipeline.

Usage examples:
  # Full pipeline: extract → cluster → enrich
  python run_pipeline.py \\
    --bam-dir /data/round2/ \\
    --r1-consensus /data/round1/R1_cluster_consensus.csv \\
    --min-q 12 \\
    --cdr-method anarci \\
    --min-r2-count 10

  # Extract + cluster only (no enrichment)
  python run_pipeline.py \\
    --bam-dir /data/round2/ \\
    --no-enrich

  # Enrichment only (pre-computed consensus files)
  python run_pipeline.py \\
    --enrich-only \\
    --r1-consensus round1_cluster_consensus.csv \\
    --r2-consensus round2_cluster_consensus.csv
"""

import argparse
import sys
from pathlib import Path
from rich.console import Console
from rich.rule import Rule

console = Console()


def run_extraction(args):
    console.print(Rule("[bold cyan]Step 1 — BAM Extraction[/bold cyan]"))
    from bam_extract import process_bam, DEFAULT_CONFIG

    input_dir  = Path(args.bam_dir)
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
    }

    bam_files = sorted(input_dir.glob("*.bam"))
    if not bam_files:
        console.print(f"[red]No .bam files found in {input_dir}[/red]")
        sys.exit(1)

    use_gpu = getattr(args, "gpu", False)
    n_workers = min(getattr(args, "workers", 1), len(bam_files))

    if n_workers == 1:
        for bam in bam_files:
            process_bam(bam, output_dir, config, use_gpu=use_gpu)
    else:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        from bam_extract import _bam_worker
        from tqdm import tqdm
        worker_args = [(bam, output_dir, config, use_gpu) for bam in bam_files]
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {executor.submit(_bam_worker, arg): arg[0].name for arg in worker_args}
            with tqdm(total=len(futures), desc="BAM files", unit="bam") as pbar:
                for future in as_completed(futures):
                    try:
                        future.result()
                    except Exception as exc:
                        console.print(f"[red]ERROR: {exc}[/red]")
                    pbar.update(1)

    protein_cdr_files = []
    for bam in bam_files:
        pf = input_dir / f"{bam.stem}_vhh_protein_cdr.csv"
        if pf.exists():
            protein_cdr_files.append(pf)
    return protein_cdr_files


def run_clustering(protein_cdr_files, args):
    console.print(Rule("[bold cyan]Step 2 — Levenshtein Clustering[/bold cyan]"))
    from cluster_levenshtein import run_levenshtein_clustering

    consensus_files = []
    for csv_path in protein_cdr_files:
        output_prefix = str(csv_path.parent / csv_path.stem)
        consensus_csv = run_levenshtein_clustering(
            str(csv_path),
            output_prefix=output_prefix,
            threshold=args.threshold,
            seq_column="CDR3",
            min_cluster_count=args.min_cluster_count,
            use_vgene=args.use_vgene,
        )
        if consensus_csv and Path(consensus_csv).exists():
            consensus_files.append(Path(consensus_csv))

    return consensus_files


def run_enrichment(r1_file, r2_file, args):
    console.print(Rule("[bold cyan]Step 3 — Enrichment Analysis[/bold cyan]"))
    from cluster_enrichment import calculate_enrichment

    output = args.enrich_output or str(Path(r2_file).parent)
    calculate_enrichment(
        str(r1_file), str(r2_file),
        output_dir=output,
        threshold=args.match_threshold,
        use_fuzzy=not args.no_fuzzy,
        log2_cutoff=args.log2_cutoff,
        fdr_cutoff=args.fdr_cutoff,
        min_r2_count=args.min_r2_count,
        entropy_flag=args.entropy_flag,
    )


def main():
    p = argparse.ArgumentParser(description="VHH NGS Pipeline Orchestrator")

    # ── Mode ──────────────────────────────────────────────────────────────────
    p.add_argument("--enrich-only", action="store_true",
                   help="Skip extraction and clustering; run enrichment only")
    p.add_argument("--no-enrich",   action="store_true",
                   help="Skip enrichment step")

    # ── BAM extraction ────────────────────────────────────────────────────────
    p.add_argument("--bam-dir",  help="Directory with .bam files")
    p.add_argument("--fr1",      default="CAGGTGCAGCTG",
                   help="FR1 anchor motif (default: CAGGTGCAGCTG)")
    p.add_argument("--j4",       default="ACCCAGGTCACC",
                   help="J4 anchor motif (default: ACCCAGGTCACC)")
    p.add_argument("--fr-mm",    type=int, default=2,
                   help="Max mismatches for FR1/J4 anchor search (default: 2)")
    p.add_argument("--min-q",    type=float, default=10.0, help="Min ONT Q-score")
    p.add_argument("--min-len",  type=int,   default=300)
    p.add_argument("--max-len",  type=int,   default=1000)
    p.add_argument("--min-aa",   type=int,   default=100)
    p.add_argument("--cdr-method", default="anarci", choices=["anarci", "offset"])
    p.add_argument("--use-umi",  action="store_true")
    p.add_argument("--workers",  type=int, default=1,
                   help="Parallel BAM workers for Step 1 (default: 1)")
    p.add_argument("--gpu",      action="store_true",
                   help="GPU Q-score batching via CuPy in Step 1")

    # ── Clustering ────────────────────────────────────────────────────────────
    p.add_argument("--threshold",        type=float, default=0.85,
                   help="Levenshtein similarity threshold (default: 0.85)")
    p.add_argument("--min-cluster-count", type=int, default=5)
    p.add_argument("--use-vgene",        action="store_true",
                   help="Stratify clustering by V_gene column if present")

    # ── Enrichment ────────────────────────────────────────────────────────────
    p.add_argument("--r1-consensus",  help="Round 1 cluster consensus .csv (or .xlsx)")
    p.add_argument("--r2-consensus",  help="Round 2 cluster consensus .csv (enrich-only mode)")
    p.add_argument("--enrich-output", default=".",
                   help="Output directory for enrichment CSV and plots (default: R2 consensus directory)")
    p.add_argument("--match-threshold", type=float, default=0.85,
                   help="Levenshtein similarity threshold for fuzzy CDR3 matching (default: 0.85)")
    p.add_argument("--no-fuzzy",  action="store_true")
    p.add_argument("--log2-cutoff",   type=float, default=1.0)
    p.add_argument("--fdr-cutoff",    type=float, default=0.05)
    p.add_argument("--min-r2-count",  type=int,   default=0,
                   help="Exclude R2 clusters with fewer than N reads from enrichment (default: 0)")
    p.add_argument("--entropy-flag",  type=float, default=1.5,
                   help="Shannon entropy threshold for heterogeneous cluster flag (default: 1.5)")

    args = p.parse_args()

    console.print(Rule("[bold green]VHH NGS Pipeline[/bold green]"))

    if args.enrich_only:
        if not args.r1_consensus or not args.r2_consensus:
            console.print("[red]--enrich-only requires --r1-consensus and --r2-consensus[/red]")
            sys.exit(1)
        run_enrichment(args.r1_consensus, args.r2_consensus, args)
        return

    if not args.bam_dir:
        console.print("[red]--bam-dir is required unless using --enrich-only[/red]")
        sys.exit(1)

    # Step 1
    protein_cdr_files = run_extraction(args)

    # Step 2
    consensus_files = run_clustering(protein_cdr_files, args)

    # Step 3
    if not args.no_enrich:
        if not args.r1_consensus:
            console.print("[yellow]Skipping enrichment: --r1-consensus not provided.[/yellow]")
        elif not consensus_files:
            console.print("[yellow]Skipping enrichment: no R2 consensus file produced.[/yellow]")
        else:
            r2_file = consensus_files[0]
            run_enrichment(args.r1_consensus, r2_file, args)

    console.print(Rule("[bold green]Pipeline complete[/bold green]"))


if __name__ == "__main__":
    main()
