#!/usr/bin/env python3
"""
run_pipeline.py — Orchestrator for the full VHH NGS enrichment pipeline.

Usage examples:
  # Full pipeline: extract → cluster → enrich
  python run_pipeline.py \\
    --bam-dir /data/round2/ \\
    --r1-consensus /data/round1/R1_cluster_consensus.xlsx \\
    --cluster-method hdbscan \\
    --min-q 12 \\
    --cdr-method anarci

  # Extract + cluster only (no enrichment)
  python run_pipeline.py \\
    --bam-dir /data/round2/ \\
    --cluster-method blosum \\
    --blosum-threshold 30

  # Enrichment only (pre-computed consensus files)
  python run_pipeline.py \\
    --enrich-only \\
    --r1-consensus round1_consensus.xlsx \\
    --r2-consensus round2_consensus.xlsx
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
        "MAX_MISMATCHES":    args.max_mm,
        "MIN_MEAN_QSCORE":   args.min_q,
        "MIN_READ_LENGTH":   args.min_len,
        "MAX_READ_LENGTH":   args.max_len,
        "MIN_VHH_AA_LENGTH": args.min_aa,
        "CDR_METHOD":        args.cdr_method,
        "USE_UMI":           args.use_umi,
        "VHH_F":             args.vhh_f,
        "VHH_R":             args.vhh_r,
    }

    bam_files = list(input_dir.glob("*.bam"))
    if not bam_files:
        console.print(f"[red]No .bam files found in {input_dir}[/red]")
        sys.exit(1)

    protein_cdr_files = []
    for bam in bam_files:
        summary = process_bam(bam, output_dir, config)
        pf = input_dir / f"{bam.stem}_vhh_protein_cdr.csv"
        if pf.exists():
            protein_cdr_files.append(pf)
    return protein_cdr_files


def run_clustering(protein_cdr_files, args):
    console.print(Rule(f"[bold cyan]Step 2 — Clustering ({args.cluster_method.upper()})[/bold cyan]"))

    consensus_files = []
    if args.cluster_method == "blosum":
        from cluster_blosum import run_blosum_clustering
        for csv_path in protein_cdr_files:
            run_blosum_clustering(
                str(csv_path),
                threshold=args.blosum_threshold,
                seq_column="CDR3",
                min_cluster_count=args.min_cluster_count,
                n_jobs=-1,
            )
            xlsx = str(csv_path).replace("_vhh_protein_cdr.csv", "_cluster_consensus.xlsx")
            if Path(xlsx).exists():
                consensus_files.append(Path(xlsx))

    elif args.cluster_method == "hdbscan":
        from hdbscan_cluster import process_csv
        cfg = {
            "min_cluster_counts": args.min_cluster_count,
            "min_cluster_size":   args.hdbscan_min_size,
            "min_samples":        args.hdbscan_min_samples,
            "use_umap":           not args.no_umap,
        }
        for csv_path in protein_cdr_files:
            process_csv(str(csv_path), cfg)
            # Find most recent output dir
            parent    = csv_path.parent
            stem      = csv_path.stem
            out_dirs  = sorted(parent.glob(f"{stem}_clustering_output_*"))
            if out_dirs:
                xlsx = out_dirs[-1] / f"{stem}_cluster_consensus.xlsx"
                if xlsx.exists():
                    consensus_files.append(xlsx)

    return consensus_files


def run_enrichment(r1_file, r2_file, args):
    console.print(Rule("[bold cyan]Step 3 — Enrichment Analysis[/bold cyan]"))
    from cluster_enrichment import calculate_enrichment

    output = args.enrich_output or str(
        Path(r2_file).parent / "VHH_enrichment.xlsx"
    )
    calculate_enrichment(
        str(r1_file), str(r2_file),
        output_file=output,
        threshold=args.match_threshold,
        use_fuzzy=not args.no_fuzzy,
        log2_cutoff=args.log2_cutoff,
        fdr_cutoff=args.fdr_cutoff,
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
    p.add_argument("--max-mm",   type=int,   default=1)
    p.add_argument("--min-q",    type=float, default=10.0, help="Min ONT Q-score")
    p.add_argument("--min-len",  type=int,   default=300)
    p.add_argument("--max-len",  type=int,   default=1000)
    p.add_argument("--min-aa",   type=int,   default=100)
    p.add_argument("--cdr-method", default="anarci", choices=["anarci", "offset"])
    p.add_argument("--use-umi",  action="store_true")
    p.add_argument("--vhh-f",    default="CAGGTACAGCTGCA")
    p.add_argument("--vhh-r",    default="CGGTGTCTAGCACT")

    # ── Clustering ────────────────────────────────────────────────────────────
    p.add_argument("--cluster-method", default="hdbscan",
                   choices=["blosum", "hdbscan"],
                   help="Clustering algorithm (default: hdbscan)")
    p.add_argument("--min-cluster-count", type=int, default=10)
    p.add_argument("--blosum-threshold",  type=float, default=30.0,
                   help="BLOSUM distance threshold (default: 30)")
    p.add_argument("--hdbscan-min-size",    type=int, default=2)
    p.add_argument("--hdbscan-min-samples", type=int, default=2)
    p.add_argument("--no-umap", action="store_true")

    # ── Enrichment ────────────────────────────────────────────────────────────
    p.add_argument("--r1-consensus",  help="Round 1 cluster consensus .xlsx")
    p.add_argument("--r2-consensus",  help="Round 2 cluster consensus .xlsx (enrich-only mode)")
    p.add_argument("--enrich-output", help="Output enrichment Excel path")
    p.add_argument("--match-threshold", type=float, default=0.90)
    p.add_argument("--no-fuzzy",  action="store_true")
    p.add_argument("--log2-cutoff", type=float, default=1.0)
    p.add_argument("--fdr-cutoff",  type=float, default=0.05)

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
