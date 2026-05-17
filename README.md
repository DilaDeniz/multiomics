# BioMultiOmics

[![CI](https://github.com/diladeniz/multiomics/actions/workflows/ci.yml/badge.svg)](https://github.com/diladeniz/multiomics/actions/workflows/ci.yml)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)

A parallel multi-omics analysis engine written in Rust. Ingests genomics (VCF), transcriptomics (TSV/BAM), epigenomics (BED), ATAC-seq (narrowPeak), and copy-number variation data simultaneously and produces an integrated HTML report, MultiQC-compatible JSON, and a live terminal dashboard.

---

## Features

**Genomics**
- VCF parsing with Ti/Tv classification (A↔G and C↔T transitions; all other SNP pairs are transversions)
- Per-chromosome variant density tracking with HyperLogLog cardinality estimation
- Allele-frequency histogram (20 bins)
- High-impact variant extraction (QUAL > 30) and gene annotation

**Copy-Number Variation**
- Parses `SVTYPE`, `CN`, `CNA`, `TCN` INFO fields from VCF
- Classifies segments: HomozygousDeletion / HeterozygousDeletion / Diploid / LowAmplification / HighAmplification
- Weighted-mean ploidy estimation and fraction-genome-altered (FGA) calculation

**Transcriptomics**
- Expression matrix TSV parsing (genes × samples)
- BAM input with noodles-bam
- DESeq2 size-factor normalization for raw integer count matrices (median-of-ratios, Anders & Huber 2010)
- Differential expression via Welch t-test on log₂(normalized + 0.5) with Benjamini-Hochberg FDR correction
- `--raw-counts` flag switches from TPM mode to DESeq2-normalized counts

**Epigenomics**
- BED methylation parsing (ENCODE bisulfite and 4-column formats)
- CpG island detection using the Gardiner-Garden criterion: GC% > 50%, length ≥ 200 bp, CpO/E ≥ 0.6 (O(n) two-pointer algorithm)
- Long-read methylation from Nanopore/PacBio BAM via MM/ML tags with CIGAR-based query-to-reference mapping (`feature = longread`)

**ATAC-seq**
- ENCODE narrowPeak (BED6+4) format parsing
- Per-chromosome peak density, signal distribution, and open-chromatin bp quantification

**Integration**
- Pearson and Spearman cross-modality correlation matrix (genomics × transcriptomics × epigenomics)
- PCA projection to 2D via linfa-reduction
- Pathway enrichment with Fisher's exact test (hypergeometric upper-tail p-value) and BH FDR correction across 75+ built-in KEGG pathways
- GSEA pre-ranked enrichment (Subramanian 2005 classic KS statistic) with 1000-permutation empirical p-values, NES normalization, and leading-edge extraction
- Custom GMT pathway file support (`--gmt`)
- Rule-based biological insight engine (Ti/Tv anomalies, global hypomethylation, correlated modalities, pathway hits)

**Output**
- Single-file HTML report with Chart.js visualizations and an inline Circos-style SVG genomic overview (no external dependencies at runtime)
- MultiQC-compatible JSON (`multiqc_bioomics.json`)
- Live terminal dashboard (ratatui TUI) with per-modality progress bars and streaming insight panel
- Python bindings via PyO3 / maturin (`pybioomics`)

---

## Installation

**From source**

```bash
git clone https://github.com/diladeniz/multiomics
cd multiomics
cargo build --release --bin bioomics
# Binary at target/release/bioomics
```

**Docker**

```bash
docker pull ghcr.io/diladeniz/bioomics:latest
docker run --rm -v $PWD/data:/data ghcr.io/diladeniz/bioomics \
  --genomics /data/variants.vcf \
  --transcriptomics /data/expression.tsv \
  --epigenomics /data/methylation.bed \
  --output /data/out
```

**Singularity**

```bash
singularity pull bioomics.sif docker://ghcr.io/diladeniz/bioomics:latest
singularity run bioomics.sif \
  --genomics variants.vcf \
  --transcriptomics expression.tsv \
  --epigenomics methylation.bed
```

**Python bindings**

```bash
pip install maturin
cd pybioomics
maturin develop --release
python -c "import pybioomics; print(pybioomics.__version__)"
```

---

## Usage

```
bioomics [OPTIONS] --genomics <FILE> --transcriptomics <FILE> --epigenomics <FILE>

Options:
  --genomics <FILE>         VCF for variant analysis
  --transcriptomics <FILE>  Expression TSV or raw count matrix
  --epigenomics <FILE>      BED methylation file
  --atac <FILE>             ENCODE narrowPeak for ATAC-seq
  --cnv <FILE>              VCF with CN INFO field for copy-number analysis
  --fastq <FILE>            FASTQ for sequence-level QC
  --gmt <FILE>              Custom GMT pathway file
  --raw-counts              Apply DESeq2 normalization to count matrix
  --output <DIR>            Output directory [default: ./bioomics_out]
  --threads <N>             Worker threads [default: all cores]
  --no-ml                   Skip PCA and correlation
  --json                    JSON output only, no TUI or HTML
```

**Minimal run**

```bash
bioomics \
  --genomics variants.vcf \
  --transcriptomics expression.tsv \
  --epigenomics methylation.bed
```

**With ATAC-seq, CNV, and custom pathways**

```bash
bioomics \
  --genomics somatic.vcf \
  --transcriptomics counts.tsv --raw-counts \
  --epigenomics methyl.bed \
  --atac peaks.narrowPeak \
  --cnv cnv_segments.vcf \
  --gmt MSigDB_Hallmarks.gmt \
  --threads 32 \
  --output results/
```

**JSON only (no TUI, no HTML — suitable for CI pipelines)**

```bash
bioomics --genomics v.vcf --transcriptomics t.tsv --epigenomics e.bed \
  --json --output pipeline_out/
```

---

## Presets

Named threshold bundles for common use cases. Use `--preset` to load one, optionally overriding specific fields with `--config`.

| Preset | Use case | Key changes |
|--------|----------|-------------|
| `cancer` | Somatic tumor/normal | QUAL≥20, strict Ti/Tv, hypomethylation alerts |
| `plant` | Plant/agricultural genomics | Relaxed Ti/Tv, lower expressed TPM |
| `rnaseq` | Bulk RNA-seq DE | Strict padj=0.01, top 200 genes |
| `wgbs` | Whole-genome bisulfite | Tight CpG island criteria |
| `atac` | ATAC-seq focus | Signal threshold=5, top 500 peaks |
| `clinical` | Clinical/translational | Conservative thresholds throughout |

```bash
bioomics --preset cancer --genomics tumor.vcf ...
bioomics --list-presets
bioomics --preset cancer --config my_overrides.toml --genomics ...
```

---

## Input Formats

### VCF (`--genomics` / `--cnv`)

Standard VCF 4.x. For CNV analysis the INFO field should contain at least one of:
- `SVTYPE=CNV|DEL|DUP|GAIN|LOSS`
- `CN=<integer>` (copy number)
- `CNA=<float>` or `TCN=<float>` (tumour copy number)

### Expression TSV (`--transcriptomics`)

Tab-delimited, first row is a header with sample names, first column is the gene identifier.

```
gene_id    sample_A    sample_B    sample_C
BRCA1      12.4        0.3         45.1
TP53       8.9         22.1        9.0
```

With `--raw-counts`, values are treated as raw integer counts and DESeq2 size-factor normalization is applied before differential expression.

### Methylation BED (`--epigenomics`)

4-column BED (chr, start, end, methylation%) or ENCODE bisulfite BED. The fifth column, if present, is interpreted as coverage depth.

### narrowPeak (`--atac`)

ENCODE BED6+4 format: chr, start, end, name, score, strand, signalValue, pValue, qValue, peak.

### GMT (`--gmt`)

Standard Gene Matrix Transposed format: one pathway per line, tab-delimited. First field is the pathway name, second is description, remaining fields are gene symbols.

---

## Output

After a successful run `--output` contains:

| File | Description |
|------|-------------|
| `report.html` | Self-contained HTML with all charts and the Circos overview |
| `multiqc_bioomics.json` | MultiQC-compatible JSON for pipeline integration |

The HTML report includes:
- Summary cards (variants, expression, methylation)
- Per-chromosome variant density bar chart
- Allele-frequency distribution histogram
- Volcano plot of differential expression (log₂FC vs −log₁₀ padj)
- Top-20 expressed genes
- Per-chromosome methylation levels
- Circos-style SVG genomic overview (three concentric rings: variant density, methylation, chromosome backbone)
- 3×3 cross-modality correlation heatmap
- PCA scatter plot
- Pathway enrichment table with Fisher's p-values and BH-adjusted q-values
- Insight list colour-coded by severity (INFO / WARN / CRIT)

---

## Architecture

```
biomics_core            — BatchAccum trait, parallel_fold, statistics (Welch, BH, Fisher)
genomics_core           — VCF parser, Ti/Tv accumulator, CNV extractor
transcriptomics_core    — TSV/BAM parser, DESeq2 normalizer, Welch t-test DE
epigenomics_core        — BED parser, CpG island detector, long-read MM/ML parser
atacseq_core            — narrowPeak parser, open-chromatin accumulator
integration_layer       — Pearson/Spearman correlation, PCA, KEGG enrichment, GSEA, GMT
cli                     — clap CLI, ratatui TUI, HTML/JSON output, Circos SVG
pybioomics              — PyO3 Python bindings (maturin)
```

### Parallelism

All three primary modalities (genomics, transcriptomics, epigenomics) are analysed concurrently via `std::thread::scope`. Within each modality, records are split into chunks of 64 000 and distributed across the rayon thread pool. Each worker builds a local accumulator, which is merged at the end via a lock-free reduce — no Mutex touches the hot path.

GSEA permutations are also parallelized with rayon, using a deterministic per-permutation xorshift64 seed derived from the permutation index.

### Memory

Files are memory-mapped with `memmap2` and `madvise(Sequential)`. Parsing is zero-copy: `ByteLines` and `TabFields` borrow directly from the mapped slice. Numeric fields use `fast-float`. `AHashMap` / `AHashSet` with AES-NI hardware hashing replace `std::HashMap` throughout.

---

## Statistics

**Differential expression**

Welch's t-test (unequal variance) on log₂(normalized + 0.5) expression values. Samples are split into two groups by index. Multiple testing correction uses the Benjamini-Hochberg step-up procedure. Genes with `padj < 0.05` are reported as significant.

With `--raw-counts`, DESeq2 size factors are estimated first:
1. Compute the geometric mean of each gene across samples in log-space.
2. Divide each sample's counts by the gene geometric means to get per-gene ratios.
3. Take the median ratio per sample as the size factor.
4. Normalize by dividing each count by the sample's size factor.

**Pathway enrichment**

Fisher's exact test (hypergeometric upper-tail p-value) using the Lanczos approximation of the log-gamma function. The query set is the union of high-impact variant genes and significantly DE genes (padj < 0.05). Results are sorted by Benjamini-Hochberg adjusted p-value.

**GSEA pre-ranked**

Classic KS enrichment score from Subramanian et al. 2005. Hit increment = √((N−Nh)/Nh), miss decrement = √(Nh/(N−Nh)). Empirical p-values use 1000 permutations computed in parallel via rayon. NES is ES / mean(|ES_null|). Leading-edge genes are extracted from the subset contributing to the ES peak.

**CpG island detection**

Gardiner-Garden & Frommer (1987) criterion: GC content > 50%, length ≥ 200 bp, observed/expected CpG ratio ≥ 0.6 where CpO/E = (n_CpG × window_length) / (n_C × n_G). Detection uses an O(n) two-pointer sliding window over sorted per-chromosome sites.

### Numerical Validation

The DESeq2 normalization and GSEA implementations are validated against
pre-computed reference values in `tests/`. Size factors match R DESeq2
to within 5% relative tolerance; log₂FC directions are verified against
manually confirmed reference datasets. GSEA enrichment scores satisfy
the theoretical bounds from Subramanian et al. 2005.

Run validation tests:

```bash
cargo test -p transcriptomics_core deseq2_validation -- --nocapture
cargo test -p integration_layer gsea_validation -- --nocapture
```

---

## Python API

```python
import pybioomics

# Analyse a VCF file
result = pybioomics.analyze_vcf("variants.vcf")
print(result.total_variants, result.titv_ratio)

# Pathway enrichment on a gene list
hits = pybioomics.enrich_pathways(["BRCA1", "TP53", "KRAS", "PIK3CA"])
for h in hits[:3]:
    print(h.pathway_name, h.padj)

# Full pipeline
out = pybioomics.run_full_pipeline(
    genomics="variants.vcf",
    transcriptomics="expression.tsv",
    epigenomics="methylation.bed",
)
```

---

## Bioconda

A recipe is provided in `conda-recipe/`. To build locally:

```bash
conda install -c conda-forge conda-build
conda build conda-recipe/
```

---

## CI

GitHub Actions runs on every push and pull request:

| Job | Platform | Description |
|-----|----------|-------------|
| `test` | ubuntu-latest, macos-latest | `cargo test --workspace --exclude pybioomics` |
| `clippy` | ubuntu-latest | `cargo clippy -- -D warnings` |
| `release` | ubuntu-latest (musl), macos-latest (arm64) | Cross-compiled static binaries attached to releases |
| `docker` | ubuntu-latest | Build and push to `ghcr.io` on version tags |

---

## License

Apache-2.0. See [LICENSE](LICENSE).
