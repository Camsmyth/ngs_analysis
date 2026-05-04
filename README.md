# ngs_analysis

VHH nanobody NGS pipeline for Oxford Nanopore (PromethION) phage display enrichment analysis.

Processes Dorado-basecalled BAM files (aligned to IMGT Vicugna germlines) through CDR extraction, sequence clustering, and round-over-round enrichment scoring.

---

## Pipeline overview

```
BAM files (Dorado/PromethION)
aligned to IMGT Vicugna germlines
        │
        ▼
  bam_extract.py            ← alignment filter, FR1/J4 anchor extraction, CDR annotation
        │
        ▼
  cluster_blosum.py         ← BLOSUM62 hierarchical clustering        ┐
        or                                                             ├─ choose one
  hdbscan_cluster.py        ← physicochemical HDBSCAN clustering      ┤
        or                                                             │
  cluster_levenshtein.py    ← graph-based clonotyping (edit distance) ┘
        │
        ▼
  cluster_enrichment.py     ← log2 CPM enrichment + FDR statistics
```

Or run everything via the orchestrator:
```
  run_pipeline.py           ← chains all steps from a single CLI call
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
  --r1-consensus /data/R1/R1_cluster_consensus.xlsx \
  --cluster-method hdbscan \
  --min-q 12 \
  --cdr-method anarci
```

---

## Requirements

- Python ≥ 3.10
- [HMMER](http://hmmer.org/) (optional but recommended — enables ANARCI IMGT numbering)
  ```bash
  brew install hmmer          # macOS
  sudo apt install hmmer      # Ubuntu/Debian
  ```

All Python dependencies are installed automatically by `setup.sh`.

---

## Scripts

### `bam_extract.py`
Extracts VHH sequences from BAM files aligned to IMGT Vicugna germlines.

Uses the alignment as a pre-filter (mapped reads only), then locates the VHH
coding region using conserved framework motifs — FR1 start (`CAGGTGCAGCTG`) and
J4 end (`ACCCAGGTCACC`) — rather than PCR primer sequences. This makes extraction
robust to variation in library prep primers and recovers sequences regardless of
which strand the aligner placed the read on.

| Flag | Default | Description |
|------|---------|-------------|
| `--min-q` | 10.0 | Minimum mean ONT Q-score |
| `--min-len` | 300 | Minimum read length (bp) |
| `--max-len` | 1000 | Maximum read length (bp) |
| `--cdr-method` | anarci | CDR extraction: `anarci` or `offset` |
| `--use-umi` | off | Enable UMI deduplication (Dorado `UX` tag) |
| `--fr1` | CAGGTGCAGCTG | FR1 anchor motif |
| `--j4` | ACCCAGGTCACC | J4 anchor motif |
| `--fr-mm` | 2 | Fuzzy mismatch tolerance for framework motifs |

```bash
python bam_extract.py /path/to/bams --min-q 12 --cdr-method anarci
```

> **Note:** BAM files must be aligned to the IMGT Vicugna VHH germline database
> before running this script (e.g. with `minimap2 -ax map-ont`). Unaligned BAMs
> will yield ~5% recovery at best.

---

### `cluster_blosum.py`
BLOSUM62-distance hierarchical clustering of CDR3 sequences. Best suited to
grouping sequences by amino-acid physicochemical similarity without regard to
CDR3 length variation.

| Flag | Default | Description |
|------|---------|-------------|
| `--input` | — | Input CSV (`*_vhh_protein_cdr.csv`) |
| `--threshold` | — | Distance cutoff 0–100 (lower = stricter, **required**) |
| `--column` | CDR3 | Sequence column to cluster |
| `--min-count` | 5 | Min total reads to report a cluster |
| `--linkage` | average | Linkage method: `single`, `complete`, `average`, `ward` |
| `--jobs` | -1 | CPU cores for distance matrix (-1 = all) |

```bash
python cluster_blosum.py \
  --input sample_vhh_protein_cdr.csv \
  --threshold 30 \
  --min-count 5 \
  --jobs -1
```

**Linkage method guide:**

| Method | Behaviour | Best for |
|--------|-----------|----------|
| `average` (UPGMA) | Cluster distance = mean of all pairwise distances | General use — balanced, outlier-resistant |
| `complete` | Cluster distance = max pairwise (furthest neighbour) | Tight homogeneous clusters; many small groups |
| `single` | Cluster distance = min pairwise (nearest neighbour) | Rarely useful — prone to chaining artefacts |
| `ward` | Minimise within-cluster variance | Avoid with BLOSUM distances (not Euclidean) |

---

### `cluster_levenshtein.py`
Graph-based clonotyping using normalised Levenshtein (edit) distance. The
preferred method for VHH CDR3 clustering when length variation between clones
is high — somatic hypermutation and junction exonuclease activity can introduce
insertions/deletions that BLOSUM gap penalties handle arbitrarily, whereas
Levenshtein naturally scales with sequence length.

Two CDR3s are linked if their normalised Levenshtein similarity ≥ `--threshold`
(default 0.80, i.e. ≤ 20% edits). Connected components of the resulting graph
define clonotypes. Uses [rapidfuzz](https://github.com/maxbachmann/RapidFuzz)
SIMD acceleration; switches to a BK-tree approximation for > 5,000 unique sequences.

| Flag | Default | Description |
|------|---------|-------------|
| `--input` | — | Input CSV (`*_vhh_protein_cdr.csv`) |
| `--threshold` | 0.80 | Similarity threshold 0–1 (higher = stricter) |
| `--column` | CDR3 | Sequence column to cluster |
| `--min-count` | 5 | Min total reads to report a clonotype |
| `--use-vgene` | off | Stratify by `V_gene` column if present |

```bash
python cluster_levenshtein.py \
  --input sample_vhh_protein_cdr.csv \
  --threshold 0.85 \
  --min-count 5
```

**Threshold guide:**

| Threshold | Edits per 14 aa CDR3 | Use case |
|-----------|----------------------|----------|
| 0.95 | ~1 | Exact-clone tracking across runs |
| 0.85 | ~2 | Standard enrichment analysis (recommended) |
| 0.80 | ~3 | Captures affinity-matured variants |
| 0.70 | ~4 | Groups distant CDR3 families |

**BLOSUM vs Levenshtein — when to choose which:**

- **Levenshtein** — when CDR3 length varies across your library (typical for
  naive or early-selection datasets); when you want B-cell repertoire-style
  clonotyping; when speed matters on large datasets.
- **BLOSUM** — when sequences are length-normalised or CDR3 lengths are
  homogeneous; when amino-acid physicochemical similarity should drive grouping
  (e.g. grouping HCDR3s that bind the same epitope via different sequences).

---

### `hdbscan_cluster.py`
Physicochemical embedding (CDR1+2+3) + UMAP + HDBSCAN density clustering.

```bash
python hdbscan_cluster.py /path/to/csvs \
  --min-cluster-counts 10 \
  --embed-cdrs CDR1 CDR2 CDR3
```

---

### `cluster_enrichment.py`
Compare two rounds of selection; computes log2(CPM) enrichment with Fisher's
exact test and BH FDR correction. Outputs an annotated Excel file, volcano plot,
and rank enrichment chart. Compatible with consensus files from all three
clustering methods.

```bash
python cluster_enrichment.py R1_consensus.xlsx R2_consensus.xlsx \
  --output VHH_enrichment.xlsx \
  --log2-cutoff 1.0 \
  --fdr-cutoff 0.05
```

---

### `run_pipeline.py`
Orchestrates all steps.

```bash
# Full pipeline
python run_pipeline.py \
  --bam-dir /data/R2/ \
  --r1-consensus /data/R1/R1_consensus.xlsx \
  --cluster-method hdbscan \
  --min-q 12

# Enrichment only (consensus files already computed)
python run_pipeline.py \
  --enrich-only \
  --r1-consensus R1.xlsx \
  --r2-consensus R2.xlsx \
  --log2-cutoff 1.5 \
  --fdr-cutoff 0.05
```

---

## Output files

| File | Script | Description |
|------|--------|-------------|
| `*_vhh_dna_protein.csv` | bam_extract | DNA + protein sequences with read counts |
| `*_vhh_protein_cdr.csv` | bam_extract | CDR1/2/3, liability flags, read counts |
| `*_cdr3_lengths.png` | bam_extract | CDR3 length distribution histogram |
| `*_qscore_dist.png` | bam_extract | ONT read Q-score distribution |
| `*_summary.json` | bam_extract | Per-BAM read fate statistics |
| `*_cluster_consensus.xlsx` | blosum / hdbscan | Cluster consensus with physicochemical properties |
| `*_lev_clonotypes.csv` | levenshtein | Per-sequence clonotype assignments |
| `*_lev_cluster_consensus.csv` | levenshtein | One row per clonotype with consensus CDR3 |
| `*_lev_cluster_consensus.xlsx` | levenshtein | Formatted Excel with count heatmap |
| `VHH_enrichment.xlsx` | enrichment | Enrichment table: log2, CPM, p-value, FDR |
| `*_volcano.png` | enrichment | Volcano plot (log2 enrichment vs −log10 FDR) |
| `*_rank_enrichment.png` | enrichment | Rank-ordered enrichment bar chart |

---

## Directory structure

```
ngs_analysis/
├── ngs/                        ← Python virtual environment (not committed)
├── bam_extract.py
├── cluster_blosum.py
├── cluster_levenshtein.py
├── hdbscan_cluster.py
├── cluster_enrichment.py
├── run_pipeline.py
├── setup.sh
├── requirements.txt
├── .gitignore
└── README.md
```

---

## Upstream notes (Dorado / pre-pipeline)

- Align BAM files to the IMGT Vicugna IGHV germline FASTA before running `bam_extract.py`; `minimap2 -ax map-ont` works well for ONT reads
- Use `--no-trim` in Dorado if VHH primers are part of the amplicon
- Use HAC or SUP basecalling model (`dna_r10.4.1_e8.2_400bps_hac@v4.3` minimum)
- Add dual UMIs in library prep for accurate PCR deduplication — enable with `--use-umi`
- Aim for ≥ 10,000 reads per barcode post-filtering

---

## License

MIT
