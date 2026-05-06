# ngs_analysis

VHH nanobody NGS pipeline for Oxford Nanopore (PromethION) phage display enrichment analysis.

Processes Dorado-basecalled BAM files (pre-aligned to IMGT Vicugna germlines) through CDR extraction, Levenshtein-based clonotyping, and round-over-round enrichment scoring.

---

## Pipeline overview

```
BAM files (Dorado/PromethION)
pre-aligned to IMGT Vicugna IGHV germlines
          │
          ▼
  bam_extract.py              Step 1 — alignment filter, FR1/J4 anchor extraction,
                              ANARCI CDR annotation, quality filtering
                              → *_vhh_protein_cdr.csv  (unique proteins + read Counts)
          │
          ▼
  cluster_levenshtein.py      Step 2 — graph-based clonotyping via normalised
                              Levenshtein distance; Count-weighted consensus,
                              biophysical metrics, liability flags
                              → *_clonotypes.csv
                              → *_cluster_consensus.csv
          │
          ▼
  cluster_enrichment.py       Step 3 — log2 % frequency enrichment (any two rounds),
                              two-tailed binomial test + BH FDR correction
                              → {sample}_{comparison}.csv + volcano/rank plots + summary log
```

Or run all three steps in sequence:
```
  run_pipeline.py             Orchestrates Steps 1–3 from a single CLI call
```

---

## Quick start

```bash
git clone https://github.com/<you>/ngs_analysis.git
cd ngs_analysis
bash setup.sh
source ngs/bin/activate

# Full pipeline — Round 2 BAMs vs pre-computed Round 1 consensus
python run_pipeline.py \
  --bam-dir /data/R2/ \
  --r1-consensus /data/R1/R1_cluster_consensus.csv \
  --min-q 12 \
  --threshold 0.85 \
  --min-r2-count 10
```

---

## Requirements

- Python ≥ 3.10
- [HMMER](http://hmmer.org/) — enables full ANARCI IMGT CDR numbering
  ```bash
  sudo apt install hmmer        # Ubuntu/Debian
  brew install hmmer            # macOS
  conda install -c bioconda hmmer
  ```
- BAM files must be aligned to the IMGT Vicugna IGHV germline FASTA **before** running `bam_extract.py`:
  ```bash
  minimap2 -ax map-ont imgt_vicugna_ighv.fa sample.bam > sample_aligned.bam
  samtools sort -o sample_aligned_sorted.bam sample_aligned.bam
  samtools index sample_aligned_sorted.bam
  ```

All Python dependencies installed by `setup.sh`.

---

## Scripts

### `bam_extract.py` — Step 1

Extracts VHH sequences from BAM files aligned to IMGT Vicugna germlines.
Designed for millions of reads across 10+ BAM files.

**Extraction approach (three-tier):**
1. **Aligned region search** — searches `query_alignment_sequence` (the soft-clip-stripped aligned region). Already in reference/coding orientation for both strands; no reverse complement needed. Fast for reads where both FR1 and J4 fall within the aligned span.
2. **Full read forward search** — searches the complete `query_sequence` (including soft clips) on the forward strand. Needed when J4 is soft-clipped beyond the alignment end, which is common when the reference has flanking sequence before FR1.
3. **RC fallback** — reverse complement search for the rare genuinely antisense read (<1% in practice).

All motif searches try an exact `str.find` first; fuzzy regex only fires on mismatches, keeping high-Q-score reads fast.

**Deduplication before annotation:**
Unique DNA amplicons are counted across all reads in pass 1. Translation and ANARCI CDR annotation run only once per unique sequence — not once per read — dramatically reducing work on high-depth libraries. ANARCI is called in batches (`--anarci-batch`) to amortise hmmer subprocess overhead.

**Parallelism:**
Multiple BAM files are processed in parallel worker processes (`--workers`). Default is one worker per BAM up to the available CPU count.

**Output:** One row per unique protein sequence with a `Count` column (total ONT reads). This Count propagates through all downstream steps and drives consensus weighting.

| Flag | Default | Description |
|------|---------|-------------|
| `--min-q` | 10.0 | Minimum mean ONT Q-score |
| `--min-len` | 300 | Minimum read length (bp) |
| `--max-len` | 1000 | Maximum read length (bp) |
| `--min-aa` | 100 | Minimum VHH protein length (aa) |
| `--cdr-method` | anarci | CDR extraction: `anarci` or `offset` |
| `--fr1` | CAGGTGCAGCTG | FR1 anchor motif |
| `--j4` | ACCCAGGTCACC | J4 anchor motif |
| `--fr-mm` | 2 | Fuzzy mismatch tolerance |
| `--use-umi` | off | Enable UMI deduplication (Dorado `UX` tag) |
| `--umi-tag` | UX | BAM tag containing the UMI |
| `--workers` | #BAMs | Parallel BAM worker processes |
| `--anarci-batch` | 5000 | Proteins per ANARCI batch call |

```bash
python bam_extract.py /path/to/bams/ --min-q 12 --workers 10
```

The summary table shows **Aln search** vs **Full read fallback** counts per file. A high full-read fallback rate is normal when the IGHV reference contains flanking sequence before FR1 (common with IMGT germlines).

---

### `cluster_levenshtein.py` — Step 2

Graph-based clonotyping using normalised Levenshtein (edit) distance.

**Approach:** Sequences are nodes; edges connect pairs with similarity ≥ threshold. Connected components define clonotypes. The graph operates on **unique CDR3 sequences**; all aggregation (counts, consensus, entropy) is **weighted by the `Count` column** from Step 1, ensuring that a sequence observed 5,000 times contributes proportionally more than one observed once. Uses [rapidfuzz](https://github.com/maxbachmann/RapidFuzz) SIMD acceleration; switches to length-binned BK-tree for > 5,000 unique sequences.

**Why Levenshtein over BLOSUM for VHH?** VHH CDR3s vary substantially in length between clones (germline D-segment usage, exonuclease trimming, somatic hypermutation indels). Levenshtein normalises by length naturally; BLOSUM gap penalties are arbitrary and length-blind.

**Outputs:**
- `*_clonotypes.csv` — every input sequence with its cluster assignment and `Cluster_Count` (total reads in that clonotype)
- `*_cluster_consensus.csv` — one row per clonotype; `Cluster_Count` = total reads (sum of member Counts); consensus sequences and biophysical metrics are Count-weighted

| Flag | Default | Description |
|------|---------|-------------|
| `--input` | — | Input CSV from Step 1 (`*_vhh_protein_cdr.csv`) |
| `--output` | input stem | Output file prefix/stem |
| `--threshold` | 0.80 | Similarity threshold 0–1 (higher = stricter) |
| `--column` | CDR3 | Sequence column to cluster |
| `--min-count` | 5 | Min total reads to report a clonotype |
| `--use-vgene` | off | Stratify by `V_gene` column if present |

```bash
python cluster_levenshtein.py \
  --input sample_vhh_protein_cdr.csv \
  --output results/sample \
  --threshold 0.85
```

**Threshold guide:**

| Threshold | Edits per 14 aa CDR3 | Use case |
|-----------|----------------------|----------|
| 0.95 | ~1 | Exact-clone tracking across runs |
| 0.85 | ~2 | Standard enrichment analysis (recommended) |
| 0.80 | ~3 | Captures affinity-matured variants |
| 0.70 | ~4 | Groups distant CDR3 families |

**Consensus output columns:**

| Column | Description |
|--------|-------------|
| `Cluster` | Clonotype ID (1 = most abundant) |
| `CDR3` | Count-weighted majority-vote consensus CDR3 |
| `CDR1`, `CDR2` | From the highest-Count member sequence |
| `CDR_Concatenated` | CDR1+CDR2+CDR3 from highest-Count member |
| `Representative_CDR3` | Highest-Count CDR3 in cluster |
| `Cluster_Count` | **Total reads** across all cluster members (sum of Count) |
| `Unique_Sequences` | Number of distinct CDR3 strings in cluster |
| `Mean_CDR3_Length` | Count-weighted mean CDR3 length |
| `CDR1/2/3_Length` | Lengths of consensus CDR sequences |
| `Shannon_Entropy` | Count-weighted diversity within clonotype |
| `pI` | Isoelectric point of CDR_Concatenated |
| `GRAVY` | Hydrophobicity (grand average of hydropathicity) |
| `Charge_pH7` | Net charge at pH 7.0 |
| `Aromaticity` | Aromatic residue fraction |
| `MW_kDa` | Molecular weight (kDa) |
| `Liabilities` | Union of PTM/liability flags across all members |

---

### `cluster_enrichment.py` — Step 3

Compares any two selection rounds (R0 vs R1, R0 vs R2, R1 vs R2) using cluster consensus files. Computes log2 % frequency enrichment, two-tailed binomial test p-values, and Benjamini-Hochberg FDR correction. Outputs a named CSV, summary log, and volcano and rank-enrichment plots.

**Normalisation:** Enrichment is calculated as log2(Freq_R2 / Freq_R1), where `Freq = cluster_count / sum(all cluster counts) × 100`. A fixed pseudocount of 1e-3 is added to each cluster count solely to prevent log(0) for clusters absent in the earlier round — it has negligible effect on any cluster with real counts. Library totals are computed from raw counts before any `--min-r2-count` filtering so the denominator is stable across thresholds.

**Statistical approach:** Two-tailed binomial test — the correct model for sequencing proportion data where only total reads per round are experimentally fixed. Laplace smoothing `(c1+0.5)/(total_r1+0.5)` gives a well-behaved expected proportion for novel clusters. BH FDR correction applied across all tested clusters.

**CDR3 matching:** Exact identity first, falling back to normalised Levenshtein similarity (rapidfuzz, SIMD-accelerated). 1:1 matching enforced — if multiple late-round clusters fuzzy-match the same early-round cluster, the first (highest-ranked) match wins and the rest are treated as novel, preventing shared-count inflation.

**Volcano plot:** Four categories, colour-coded:

| Colour | Category |
|--------|----------|
| Green | Enriched AND matched (exact or fuzzy) to early round — labelled with CDR3 |
| Red | Enriched but novel (absent in early round) |
| Blue | Depleted (top 5 labelled) |
| Grey | Not significant |

Accepts `.csv` or `.xlsx` consensus input.

```bash
python cluster_enrichment.py \
  R1_cluster_consensus.csv \
  R2_cluster_consensus.csv \
  --output results/ \
  --min-r2-count 50 \
  --log2-cutoff 1.0 \
  --fdr-cutoff 0.05
```

**Multi-comparison shell loop** (R0 vs R2, R0 vs R1, R1 vs R2 per sample, with output files named by sample and comparison):

```bash
for DIR in vhh_extract_05052026/*/; do
    SAMPLE=$(basename "$DIR")
    [[ "$SAMPLE" == "BAM_extract" ]] && continue

    R0=$(ls "$DIR"Alpaca*_cluster_consensus.csv 2>/dev/null | head -1)
    R1=$(ls "$DIR"*-1_vhh_protein_cdr_cluster_consensus.csv 2>/dev/null | head -1)
    R2=$(ls "$DIR"*-2_vhh_protein_cdr_cluster_consensus.csv 2>/dev/null | head -1)

    [[ -z "$R2" ]] && continue

    run_enrichment() {
        local early="$1" late="$2" label="$3" outdir="$4"
        python3 cluster_enrichment.py "$early" "$late" \
            --output "$outdir" \
            --min-r2-count 50
        for f in "$outdir"/VHH_enrichment*; do
            mv "$f" "${f/VHH_enrichment/${SAMPLE}_${label}}"
        done
    }

    [[ -n "$R0" ]]             && run_enrichment "$R0" "$R2" "R0vsR2" "${DIR}enrichment_R0vsR2"
    [[ -n "$R0" && -n "$R1" ]] && run_enrichment "$R0" "$R1" "R0vsR1" "${DIR}enrichment_R0vsR1"
    [[ -n "$R1" ]]             && run_enrichment "$R1" "$R2" "R1vsR2" "${DIR}enrichment_R1vsR2"
done
```

Output files per comparison (e.g. `h19/enrichment_R0vsR2/`):

| File | Description |
|------|-------------|
| `{sample}_{comparison}.csv` | Enrichment table: Freq_R1, Freq_R2, Log2_Enrichment, p-value, FDR |
| `{sample}_{comparison}_volcano.png` | Volcano plot with matched/novel/depleted colour coding |
| `{sample}_{comparison}_rank_enrichment.png` | Rank-ordered enrichment bar chart |
| `{sample}_{comparison}_summary.txt` | Plain-text summary statistics + run parameters |

| Flag | Default | Description |
|------|---------|-------------|
| `file_r1` | — | Early round consensus file (.csv or .xlsx) |
| `file_r2` | — | Late round consensus file (.csv or .xlsx) |
| `--output` | `.` | Output directory |
| `--threshold` | 0.85 | Normalised Levenshtein similarity for fuzzy CDR3 matching |
| `--no-fuzzy` | off | Exact CDR3 matching only |
| `--log2-cutoff` | 1.0 | log2 enrichment cutoff for volcano significance threshold |
| `--fdr-cutoff` | 0.05 | FDR significance threshold |
| `--min-r2-count` | 0 | Exclude late-round clusters with fewer than N reads |
| `--entropy-flag` | 1.5 | Shannon entropy above which clusters are flagged as `heterogeneous` |

**Quality_Flag column** (present when `Shannon_Entropy` is in the late-round input):

| Value | Meaning |
|-------|---------|
| `heterogeneous` | Shannon entropy > `--entropy-flag`; cluster contains many divergent CDR3s |
| `low_depth` | Late-round cluster has < 20 total reads; fold change unreliable |
| _(empty)_ | No quality concern |

---

### `run_pipeline.py` — Orchestrator

Chains Steps 1–3 from a single command.

```bash
# Full pipeline
python run_pipeline.py \
  --bam-dir /data/R2/ \
  --r1-consensus /data/R1/R1_cluster_consensus.csv \
  --min-q 12 \
  --threshold 0.85

# Enrichment only (consensus CSVs already computed)
python run_pipeline.py \
  --enrich-only \
  --r1-consensus R1_cluster_consensus.csv \
  --r2-consensus R2_cluster_consensus.csv
```

---

## Output files

| File | Step | Description |
|------|------|-------------|
| `*_vhh_protein_cdr.csv` | 1 | Unique proteins: CDR1/2/3, liability flags, read Count |
| `*_vhh_dna_protein.csv` | 1 | DNA + protein sequences with Count |
| `*_cdr3_lengths.png` | 1 | CDR3 length distribution |
| `*_qscore_dist.png` | 1 | ONT read Q-score distribution |
| `*_summary.json` | 1 | Per-BAM read fate statistics |
| `*_clonotypes.csv` | 2 | Per-sequence clonotype assignments + Cluster_Count |
| `*_cluster_consensus.csv` | 2 | One row per clonotype; Count-weighted consensus + biophysics |
| `{sample}_{comparison}.csv` | 3 | Enrichment table: Freq_R1, Freq_R2, Log2_Enrichment, p-value, FDR |
| `{sample}_{comparison}_summary.txt` | 3 | Plain-text summary statistics + run parameters |
| `*_volcano.png` | 3 | Volcano plot (log2 enrichment vs −log10 FDR) |
| `*_rank_enrichment.png` | 3 | Rank-ordered enrichment bar chart |

---

## Directory structure

```
ngs_analysis/
├── ngs/                        ← Python virtual environment (not committed)
├── bam_extract.py              ← Step 1: ONT BAM → unique VHH proteins + Counts
├── cluster_levenshtein.py      ← Step 2: clonotyping + biophysical annotation
├── cluster_enrichment.py       ← Step 3: R1 vs R2 enrichment statistics
├── run_pipeline.py             ← Orchestrator
├── setup.sh                    ← Environment installer
├── requirements.txt
├── .gitignore
└── README.md
```

---

## Upstream notes

- Align BAM files to the IMGT Vicugna IGHV germline FASTA before Step 1 (`minimap2 -ax map-ont`)
- Use HAC or SUP Dorado basecalling model (`dna_r10.4.1_e8.2_400bps_hac@v4.3` minimum)
- Use `--no-trim` in Dorado if VHH primers are part of the amplicon
- Add dual UMIs in library prep for accurate PCR deduplication — enable with `--use-umi`
- Aim for ≥ 10,000 reads per barcode post-filtering

---

## License

MIT
