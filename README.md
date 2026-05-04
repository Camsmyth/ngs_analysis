# ngs_analysis

VHH nanobody NGS pipeline for Oxford Nanopore (PromethION) phage display enrichment analysis.

Processes Dorado-basecalled BAM files through CDR extraction, sequence clustering, and round-over-round enrichment scoring.

---

## Pipeline overview

```
BAM files (Dorado/PromethION)
        │
        ▼
  bam_extract.py          ← primer search, quality filter, CDR annotation
        │
        ▼
  cluster_blosum.py       ← BLOSUM62 hierarchical clustering   ┐
        or                                                      ├─ choose one
  hdbscan_cluster.py      ← physicochemical HDBSCAN clustering ┘
        │
        ▼
  cluster_enrichment.py   ← log2 CPM enrichment + FDR statistics
```

Or run everything via the orchestrator:
```
  run_pipeline.py         ← chains all steps from a single CLI call
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
Extracts VHH sequences from BAM files. Key options:

| Flag | Default | Description |
|------|---------|-------------|
| `--min-q` | 10.0 | Minimum mean ONT Q-score |
| `--min-len` | 300 | Minimum read length (bp) |
| `--max-len` | 1000 | Maximum read length (bp) |
| `--cdr-method` | anarci | CDR extraction: `anarci` or `offset` |
| `--use-umi` | off | Enable UMI deduplication (Dorado `UX` tag) |
| `--max-mm` | 1 | Fuzzy primer mismatch tolerance |
| `--vhh-f` | CAGGTACAGCTGCA | Forward primer |
| `--vhh-r` | CGGTGTCTAGCACT | Reverse primer |

```bash
python bam_extract.py /path/to/bams --min-q 12 --cdr-method anarci
```

### `cluster_blosum.py`
BLOSUM62 hierarchical clustering of CDR3 sequences.

```bash
python cluster_blosum.py \
  --input sample_vhh_protein_cdr.csv \
  --threshold 30 \
  --min-count 5 \
  --jobs -1
```

### `hdbscan_cluster.py`
Physicochemical embedding (CDR1+2+3) + UMAP + HDBSCAN density clustering.

```bash
python hdbscan_cluster.py /path/to/csvs \
  --min-cluster-counts 10 \
  --embed-cdrs CDR1 CDR2 CDR3
```

### `cluster_enrichment.py`
Compare two rounds of selection; computes log2(CPM) enrichment with Fisher's exact test and BH FDR correction. Outputs an annotated Excel file, volcano plot, and rank enrichment chart.

```bash
python cluster_enrichment.py R1_consensus.xlsx R2_consensus.xlsx \
  --output VHH_enrichment.xlsx \
  --log2-cutoff 1.0 \
  --fdr-cutoff 0.05
```

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

| File | Description |
|------|-------------|
| `*_vhh_dna_protein.csv` | DNA + protein sequences with read counts |
| `*_vhh_protein_cdr.csv` | CDR1/2/3, liability flags, read counts |
| `*_cdr3_lengths.png` | CDR3 length distribution histogram |
| `*_qscore_dist.png` | ONT read Q-score distribution |
| `*_cluster_consensus.xlsx` | Cluster consensus sequences with physicochemical properties |
| `VHH_enrichment.xlsx` | Enrichment table with log2, CPM, p-value, FDR |
| `*_volcano.png` | Volcano plot (log2 enrichment vs −log10 FDR) |
| `*_rank_enrichment.png` | Rank-ordered enrichment bar chart |
| `*_summary.json` | Per-BAM read fate statistics |

---

## Directory structure

```
ngs_analysis/
├── ngs/                        ← Python virtual environment (not committed)
├── bam_extract.py
├── cluster_blosum.py
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

- Use `--no-trim` in Dorado if VHH primers are part of the amplicon
- Use HAC or SUP basecalling model (`dna_r10.4.1_e8.2_400bps_hac@v4.3` minimum)
- Add dual UMIs in library prep for accurate PCR deduplication — enable with `--use-umi`
- Aim for ≥10,000 reads per barcode post-filtering

---

## License

MIT
