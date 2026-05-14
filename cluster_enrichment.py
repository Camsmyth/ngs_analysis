#!/usr/bin/env python3
"""
cluster_enrichment.py — Step 3 of the VHH NGS pipeline.

Compare clustering outputs from two phage display selection rounds and
compute enrichment scores.

  • CPM-normalised log2 fold change with pseudocount
  • Two-tailed binomial test (theoretically correct for sequencing count data
    where only total reads per round are fixed, not per-cluster totals)
  • Laplace-smoothed expected proportion handles novel clusters (Count_R1=0)
  • Benjamini-Hochberg FDR correction
  • 1:1 R1→R2 matching enforced; duplicate fuzzy matches treated as novel
  • Shannon entropy quality flag for heterogeneous clusters
  • Volcano plot with bubble size proportional to R2 cluster depth
  • Accepts both CSV and Excel consensus inputs (auto-detects by extension)
  • Outputs CSV + plots to a user-specified directory
"""

import sys
import os
import argparse
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests
from rapidfuzz.distance import Levenshtein
from rapidfuzz import process as fuzz_process
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from rich.console import Console
from rich.table import Table

console = Console()
warnings.filterwarnings("ignore", category=FutureWarning)

# ══════════════════════════════════════════════════════════════════════════════
# CDR3 MATCHING
# ══════════════════════════════════════════════════════════════════════════════

def build_exact_lookup(df: pd.DataFrame, cdr_col: str = "CDR3") -> dict:
    """Map CDR3 sequence → row index for O(1) exact lookup."""
    return {str(seq).upper().strip(): i
            for i, seq in df[cdr_col].items() if pd.notna(seq)}


def match_cdr3(
    query:     str,
    lookup:    dict,
    ref_seqs:  list,
    threshold: float = 0.85,
    use_fuzzy: bool  = True,
) -> Optional[str]:
    """
    Try exact match first; fall back to normalised Levenshtein similarity via
    rapidfuzz (SIMD-accelerated). extractOne applies score_cutoff for early exit,
    making this O(n) with a small constant vs O(n * L^2) for BLOSUM alignment.
    """
    q = str(query).upper().strip()
    if q in lookup:
        return q

    if not use_fuzzy:
        return None

    result = fuzz_process.extractOne(
        q, ref_seqs,
        scorer=Levenshtein.normalized_similarity,
        score_cutoff=threshold,
    )
    return result[0] if result else None


# ══════════════════════════════════════════════════════════════════════════════
# STATISTICS
# ══════════════════════════════════════════════════════════════════════════════

def compute_enrichment(df: pd.DataFrame, min_r2_count: int = 0) -> pd.DataFrame:
    """
    Add enrichment metrics using the binomial test.

    Normalisation: % frequency = cluster_count / sum(all cluster counts).
    A fixed pseudocount of 1e-3 is added solely to prevent log(0); it has
    negligible effect on any cluster with real counts. Library totals are
    computed from raw counts before min_r2_count filtering so the denominator
    is stable across thresholds.
    """
    pseudocount = 1e-3

    # Stable library totals — computed before filtering
    total_r1_int = int(df["Count_R1"].sum())
    total_r2_int = int(df["Count_R2"].sum())

    if min_r2_count > 0:
        df = df[df["Count_R2"] >= min_r2_count].copy().reset_index(drop=True)

    df["Freq_R1"] = (df["Count_R1"] + pseudocount) / total_r1_int * 100
    df["Freq_R2"] = (df["Count_R2"] + pseudocount) / total_r2_int * 100
    df["Log2_Enrichment"] = np.log2(df["Freq_R2"] / df["Freq_R1"])

    # Two-tailed binomial test: detects both enrichment and depletion
    pvals = []
    for _, row in df.iterrows():
        c2 = int(row["Count_R2"])
        c1 = int(row["Count_R1"])
        p_expected = (c1 + pseudocount) / (total_r1_int + pseudocount)
        result = binomtest(c2, total_r2_int, p=p_expected, alternative="two-sided")
        pvals.append(result.pvalue)

    df["Pvalue"] = pvals
    _, fdrs, _, _ = multipletests(pvals, alpha=0.05, method="fdr_bh")
    df["FDR"]         = fdrs
    df["Neg_log10_FDR"] = -np.log10(df["FDR"].clip(lower=1e-300))

    return df.sort_values("Log2_Enrichment", ascending=False).reset_index(drop=True)


# ══════════════════════════════════════════════════════════════════════════════
# PLOTTING
# ══════════════════════════════════════════════════════════════════════════════

def plot_volcano(df: pd.DataFrame, output_path: Path,
                 log2_threshold: float = 1.0, fdr_threshold: float = 0.05):
    """Volcano plot: log2 enrichment vs -log10 FDR. Dot size = sqrt(R2 count)."""
    matched  = df["Match_Type"].isin(["exact", "fuzzy"])
    enriched = (df["Log2_Enrichment"] >= log2_threshold) & (df["FDR"] < fdr_threshold)
    depleted = (df["Log2_Enrichment"] <= -log2_threshold) & (df["FDR"] < fdr_threshold)

    enriched_matched   = enriched &  matched
    enriched_unmatched = enriched & ~matched
    neutral            = ~enriched & ~depleted

    def _sizes(mask):
        return np.sqrt(df.loc[mask, "Count_R2"].clip(lower=1)).clip(5, 60)

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.scatter(df.loc[neutral,            "Log2_Enrichment"], df.loc[neutral,            "Neg_log10_FDR"],
               c="grey",    alpha=0.4, s=_sizes(neutral),            label="Not significant")
    ax.scatter(df.loc[depleted,           "Log2_Enrichment"], df.loc[depleted,           "Neg_log10_FDR"],
               c="#1f77b4", alpha=0.8, s=_sizes(depleted),           label=f"Depleted (n={depleted.sum()})")
    ax.scatter(df.loc[enriched_unmatched, "Log2_Enrichment"], df.loc[enriched_unmatched, "Neg_log10_FDR"],
               c="#d62728", alpha=0.8, s=_sizes(enriched_unmatched), label=f"Enriched – novel (n={enriched_unmatched.sum()})")
    ax.scatter(df.loc[enriched_matched,   "Log2_Enrichment"], df.loc[enriched_matched,   "Neg_log10_FDR"],
               c="#2ca02c", alpha=0.9, s=_sizes(enriched_matched),   label=f"Enriched – matched (n={enriched_matched.sum()})")

    ax.axvline(log2_threshold,  color="red",   linestyle="--", linewidth=0.8, alpha=0.7)
    ax.axvline(-log2_threshold, color="blue",  linestyle="--", linewidth=0.8, alpha=0.7)
    ax.axhline(-np.log10(fdr_threshold), color="black", linestyle="--",
               linewidth=0.8, alpha=0.7, label=f"FDR={fdr_threshold}")

    # Label all enriched+matched, then top depleted
    for sub_df, ha in [(df[enriched_matched], "left"),
                       (df[depleted].nsmallest(5, "Log2_Enrichment"), "right")]:
        for _, row in sub_df.iterrows():
            ax.annotate(str(row.get("CDR3", ""))[:12],
                        (row["Log2_Enrichment"], row["Neg_log10_FDR"]),
                        fontsize=6, ha=ha, va="bottom", alpha=0.9)

    ax.set_xlabel("log₂(Enrichment)", fontsize=12)
    ax.set_ylabel("-log₁₀(FDR)", fontsize=12)
    ax.set_title("VHH Cluster Enrichment Volcano", fontsize=14)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    console.print(f"[green]✓ Volcano plot:[/green] {output_path}")


def plot_rank_enrichment(df: pd.DataFrame, output_path: Path):
    """Rank-ordered enrichment plot."""
    vals = df["Log2_Enrichment"].sort_values(ascending=False).values
    fig, ax = plt.subplots(figsize=(10, 5))
    colors = ["#d62728" if v > 0 else "#1f77b4" for v in vals]
    ax.bar(range(len(vals)), vals, color=colors, width=1.0, edgecolor="none")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Cluster rank", fontsize=12)
    ax.set_ylabel("log₂(Enrichment)", fontsize=12)
    ax.set_title("Rank-ordered Cluster Enrichment", fontsize=14)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    console.print(f"[green]✓ Rank enrichment plot:[/green] {output_path}")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def _read_consensus(filepath: str) -> pd.DataFrame:
    """Read a cluster consensus file — accepts .csv or .xlsx."""
    ext = Path(filepath).suffix.lower()
    if ext == ".csv":
        return pd.read_csv(filepath)
    elif ext in (".xlsx", ".xls"):
        return pd.read_excel(filepath)
    else:
        console.print(f"[red]ERROR: Unsupported file extension '{ext}'. Use .csv or .xlsx.[/red]")
        sys.exit(1)


def calculate_enrichment(
    file_r1:      str,
    file_r2:      str,
    output_dir:   str   = ".",
    threshold:    float = 0.85,
    use_fuzzy:    bool  = True,
    log2_cutoff:  float = 0.7,
    fdr_cutoff:   float = 0.05,
    min_r2_count: int   = 100,
    entropy_flag: float = 1.5,
    novel_cutoff: int   = 25,
) -> pd.DataFrame:

    df_r1 = _read_consensus(file_r1)
    df_r2 = _read_consensus(file_r2)

    for df in [df_r1, df_r2]:
        df.columns = df.columns.str.strip()

    for name, df in [("R1", df_r1), ("R2", df_r2)]:
        missing = [c for c in ["CDR3", "Cluster_Count"] if c not in df.columns]
        if missing:
            console.print(f"[red]ERROR: {name} file missing columns: {missing}[/red]")
            sys.exit(1)

    console.print(f"  R1 clusters: {len(df_r1):,}  |  R2 clusters: {len(df_r2):,}")

    lookup_r1 = build_exact_lookup(df_r1, "CDR3")
    seqs_r1   = list(lookup_r1.keys())

    # ── Match each R2 cluster to R1 (1:1 enforced) ───────────────────────────
    # Greedy first-match: once an R1 CDR3 is claimed, subsequent R2 clusters
    # that would fuzzy-match the same R1 sequence are treated as novel.
    # Prevents shared-R1-count inflation across multiple R2 clusters.
    claimed_r1: set = set()
    rows = []

    for _, row_r2 in df_r2.iterrows():
        cdr3_r2 = str(row_r2["CDR3"]).upper().strip()
        matched = match_cdr3(cdr3_r2, lookup_r1, seqs_r1,
                             threshold, use_fuzzy)

        if matched is not None and matched in claimed_r1:
            matched = None   # R1 already consumed by a better-ranked R2 cluster

        if matched is not None:
            claimed_r1.add(matched)
            r1_row   = df_r1.iloc[lookup_r1[matched]]
            count_r1 = int(r1_row["Cluster_Count"])
            meta_r1  = {f"{c}_R1": r1_row[c] for c in df_r1.columns
                        if c not in ["CDR3", "Cluster_Count"]}
        else:
            count_r1 = 0
            meta_r1  = {}

        combined = {
            "CDR3":       cdr3_r2,
            "Count_R1":   count_r1,
            "Count_R2":   int(row_r2["Cluster_Count"]),
            "Match_Type": "exact" if cdr3_r2 in lookup_r1 else ("fuzzy" if matched else "novel"),
            **{f"{c}_R2": row_r2[c] for c in df_r2.columns
               if c not in ["CDR3", "Cluster_Count"]},
            **meta_r1,
        }
        rows.append(combined)

    merged = pd.DataFrame(rows)
    merged = compute_enrichment(merged, min_r2_count=min_r2_count)

    # ── Shannon entropy quality flag ──────────────────────────────────────────
    if "Shannon_Entropy_R2" in merged.columns:
        merged["Quality_Flag"] = np.where(
            merged["Shannon_Entropy_R2"] > entropy_flag, "heterogeneous",
            np.where(merged["Count_R2"] < 20, "low_depth", "")
        )

    # ── Outputs ───────────────────────────────────────────────────────────────
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = "VHH_enrichment"

    csv_path = out_dir / f"{stem}.csv"
    merged.to_csv(csv_path, index=False)
    console.print(f"[green]✓ Enrichment CSV:[/green] {csv_path}")

    known_enriched = merged[
        (merged["Log2_Enrichment"] >= log2_cutoff) &
        (merged["FDR"] < fdr_cutoff) &
        (merged["Match_Type"].isin(["exact", "fuzzy"]))
    ]
    known_path = out_dir / f"{stem}_known_enriched.csv"
    known_enriched.to_csv(known_path, index=False)
    console.print(f"[green]✓ Known enriched CSV:[/green] {known_path} ({len(known_enriched)} clusters)")

    novel_enriched = merged[
        (merged["Log2_Enrichment"] >= log2_cutoff) &
        (merged["FDR"] < fdr_cutoff) &
        (merged["Match_Type"] == "novel") &
        (merged["Count_R2"] >= novel_cutoff)
    ]
    novel_path = out_dir / f"{stem}_novel_enriched.csv"
    novel_enriched.to_csv(novel_path, index=False)
    console.print(f"[green]✓ Novel enriched CSV:[/green] {novel_path} ({len(novel_enriched)} clusters)")

    plot_volcano(merged, out_dir / f"{stem}_volcano.png",
                 log2_threshold=log2_cutoff, fdr_threshold=fdr_cutoff)
    plot_rank_enrichment(merged, out_dir / f"{stem}_rank_enrichment.png")

    # ── Summary table ─────────────────────────────────────────────────────────
    sig_enriched = ((merged["Log2_Enrichment"] >= log2_cutoff) & (merged["FDR"] < fdr_cutoff)).sum()
    sig_depleted = ((merged["Log2_Enrichment"] <= -log2_cutoff) & (merged["FDR"] < fdr_cutoff)).sum()
    novel        = (merged["Match_Type"] == "novel").sum()

    summary_rows = [
        ("R2 clusters total",                                       str(len(merged))),
        ("Exact matches",                                           str((merged["Match_Type"] == "exact").sum())),
        ("Fuzzy matches",                                           str((merged["Match_Type"] == "fuzzy").sum())),
        ("Novel (absent in R1)",                                    str(novel)),
        (f"Enriched (log2≥{log2_cutoff}, FDR<{fdr_cutoff})",          str(sig_enriched)),
        (f"Enriched + matched (known_enriched.csv)",                    str(len(known_enriched))),
        (f"Enriched + novel ≥{novel_cutoff} reads (novel_enriched.csv)", str(len(novel_enriched))),
        (f"Depleted (log2≤-{log2_cutoff}, FDR<{fdr_cutoff})",         str(sig_depleted)),
    ]

    table = Table(title="Enrichment Summary", show_lines=True)
    table.add_column("Metric",  style="cyan")
    table.add_column("Value",   justify="right")
    for metric, value in summary_rows:
        table.add_row(metric, value)
    console.print(table)

    # ── Summary log ───────────────────────────────────────────────────────────
    import datetime
    log_path = out_dir / f"{stem}_summary.txt"
    col_w = max(len(m) for m, _ in summary_rows) + 2
    with open(log_path, "w") as fh:
        fh.write("Enrichment Summary\n")
        fh.write("=" * (col_w + 12) + "\n")
        fh.write(f"Generated : {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        fh.write(f"R1 input  : {file_r1}\n")
        fh.write(f"R2 input  : {file_r2}\n")
        fh.write(f"min_r2_count : {min_r2_count}\n")
        fh.write(f"log2_cutoff  : {log2_cutoff}\n")
        fh.write(f"fdr_cutoff   : {fdr_cutoff}\n")
        fh.write(f"threshold    : {threshold}\n")
        fh.write("-" * (col_w + 12) + "\n")
        for metric, value in summary_rows:
            fh.write(f"{metric:<{col_w}}{value}\n")
    console.print(f"[green]✓ Summary log:[/green] {log_path}")

    console.print(f"\n[green]✓ Enrichment written to:[/green] {out_dir}/")

    return merged


# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="VHH phage display enrichment analysis (R1 vs R2 cluster comparison)"
    )
    parser.add_argument("file_r1",    help="Round 1 cluster consensus file (.csv or .xlsx)")
    parser.add_argument("file_r2",    help="Round 2 cluster consensus file (.csv or .xlsx)")
    parser.add_argument("--output", default=".",
                        help="Output directory for enrichment CSV and plots (default: current directory)")
    parser.add_argument("--threshold", type=float, default=0.85,
                        help="Normalised Levenshtein similarity threshold for fuzzy CDR3 matching (default: 0.85)")
    parser.add_argument("--no-fuzzy", action="store_true",
                        help="Disable fuzzy matching; exact CDR3 only")
    parser.add_argument("--log2-cutoff", type=float, default=0.7,
                        help="Log2 enrichment cutoff for volcano (default: 0.7)")
    parser.add_argument("--fdr-cutoff",  type=float, default=0.05,
                        help="FDR cutoff for significance (default: 0.05)")
    parser.add_argument("--min-r2-count", type=int, default=0,
                        help="Exclude R2 clusters with fewer than N total reads (default: 0 = off)")
    parser.add_argument("--entropy-flag", type=float, default=1.5,
                        help="Shannon entropy above which clusters are flagged as heterogeneous "
                             "(default: 1.5 bits; requires Shannon_Entropy column in R2 input)")
    parser.add_argument("--novel-cutoff", type=int, default=25,
                        help="Min R2 Count for a novel (unmatched) enriched cluster to appear in "
                             "*_novel_enriched.csv (default: 25)")
    args = parser.parse_args()

    if not os.path.exists(args.file_r1) or not os.path.exists(args.file_r2):
        console.print("[red]ERROR: Input file(s) not found.[/red]")
        sys.exit(1)

    calculate_enrichment(
        args.file_r1, args.file_r2,
        output_dir=args.output,
        threshold=args.threshold,
        use_fuzzy=not args.no_fuzzy,
        log2_cutoff=args.log2_cutoff,
        fdr_cutoff=args.fdr_cutoff,
        min_r2_count=args.min_r2_count,
        entropy_flag=args.entropy_flag,
        novel_cutoff=args.novel_cutoff,
    )


if __name__ == "__main__":
    main()
