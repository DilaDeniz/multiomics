# multiomics

[![License: Apache-2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/Rust-1.78%2B-orange.svg)](https://www.rust-lang.org)

A fast, parallel multi-omics analysis CLI written in Rust. Feed it genomics (VCF), transcriptomics (TSV), epigenomics (BED), proteomics (mzML), single-cell RNA (10x MEX), ATAC-seq, CNV, and more — get back an integrated HTML report, a MultiQC-compatible JSON file, and a live terminal dashboard.

```
┌─ multiomics ─────────────────────────────────────────────────────────────────┐
│  Phase: Genomics Analysis                               Elapsed: 00:01:23    │
├─────────────────────────────────┬────────────────────────────────────────────┤
│  GENOMICS      [████████░░] 82% │  LIVE INSIGHTS                             │
│  58,432 rec/s  ETA: 00:00:12   │  [INFO]  Ti/Tv = 2.14 (normal range)       │
│                                 │  [WARN]  12 high-impact KRAS variants      │
│  TRANSCRIPTOMICS [░░░░░░░░]  0% │  [INFO]  8,942 expressed genes (TPM≥1)    │
│  EPIGENOMICS     [░░░░░░░░]  0% │  [CRIT]  Global methylation: 38.2%        │
│  INTEGRATION     [░░░░░░░░]  0% │  [INFO]  PI3K-Akt pathway enriched        │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

## Features

| Module | Input | What it computes |
|---|---|---|
| **Genomics** | VCF | Ti/Tv ratio, variant density per chromosome, allele-frequency histogram, high-impact variant genes |
| **Transcriptomics** | TSV (genes × samples) | Expressed gene count, top-100 by TPM, log₂ fold-change DE (when ≥ 2 samples) |
| **Epigenomics** | BED (bisulfite) | Global methylation %, CpG island detection, hyper/hypomethylated regions |
| **Proteomics** | mzML + FASTA | Database search, PSM/peptide/protein counts at user-defined FDR, phosphoproteomics |
| **Single-cell RNA** | 10x MEX directory | QC, log-normalisation, HVG selection, GPU-accelerated UMAP |
| **ATAC-seq** | narrowPeak BED6+4 | Peak count, median signal, per-chromosome accessibility |
| **CNV** | VCF with CN field | Segment count, fraction genome altered, estimated ploidy |
| **Gene quantification** | BAM + GTF | Read-level counting, assignment rate, top expressed genes |
| **FASTQ QC** | FASTQ | Read count, GC%, Q30% |
| **Integration** | all of the above | PCA, cross-modality Pearson correlation, KEGG pathway enrichment, rule-based insights |

---

## Installation

### Requirements

- Rust 1.78 or newer — install via [rustup.rs](https://rustup.rs)
- A C linker (usually already present; on Windows install [Build Tools for Visual Studio](https://visualstudio.microsoft.com/visual-cpp-build-tools/))

### Build (CPU only — works everywhere)

```bash
git clone https://github.com/diladeniz/multiomics.git
cd multiomics
cargo build --release --bin multiomics
# binary at: target/release/multiomics
```

### Build with GPU acceleration

GPU support uses [wgpu](https://wgpu.rs) and runs on any modern GPU via Vulkan (Linux/Windows), Metal (macOS), or DX12 (Windows). No CUDA required.

It accelerates the k-nearest-neighbour step of UMAP for single-cell datasets with 5 000+ cells.

```bash
cargo build --release --bin multiomics --features gpu
```

### Build everything (ATAC + CNV + longread + GPU)

```bash
cargo build --release --bin multiomics --features full,gpu
```

### Install to PATH

```bash
cargo install --path cli                    # CPU only
cargo install --path cli --features gpu     # with GPU
```

---

## Quick start

```bash
# Genomics + transcriptomics + epigenomics — HTML report + JSON
multiomics \
  --genomics sample.vcf \
  --transcriptomics expression.tsv \
  --epigenomics methylation.bed \
  --output ./results

# JSON only (no TUI, no HTML) — good for automated pipelines
multiomics \
  --genomics sample.vcf \
  --transcriptomics expression.tsv \
  --epigenomics methylation.bed \
  --json \
  --output ./results

# Single-cell RNA with GPU UMAP
multiomics \
  --scrna /path/to/10x_mex_dir \
  --output ./sc_results

# Proteomics database search at 1% FDR
multiomics \
  --proteomics run1.mzML run2.mzML \
  --fasta human_proteome.fasta \
  --proteomics-fdr 0.01 \
  --output ./prot_results

# Full multi-omics run
multiomics \
  --genomics tumor.vcf \
  --transcriptomics rna.tsv \
  --epigenomics wgbs.bed \
  --atac peaks.narrowPeak \
  --cnv cnv.vcf \
  --scrna mex/ \
  --proteomics ms1.mzML ms2.mzML \
  --fasta proteome.fasta \
  --output ./full_results
```

---

## Input formats

| Flag | Format |
|---|---|
| `--genomics` | VCF 4.x (gzipped or plain). QUAL, INFO/AF, and gene name from INFO/ANN are used when present. |
| `--transcriptomics` | Tab-separated matrix: first row = sample names, first column = gene ID, values = TPM (or raw counts with `--raw-counts`). |
| `--epigenomics` | BED with columns: chrom, start, end, name, score, strand, methylation% (ENCODE bisulfite). A simpler 4-column BED (chrom, start, end, methylation%) is also accepted. |
| `--atac` | ENCODE narrowPeak (BED6+4). Requires `--features atac`. |
| `--cnv` | VCF with `CN=<int>` in INFO. Requires `--features cnv`. |
| `--scrna` | 10x Genomics MEX directory: `matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`. |
| `--proteomics` | One or more mzML files. |
| `--proteomics-dir` | Directory of `*.mzML` files — scanned automatically (alternative to listing files individually). |
| `--fasta` | Protein FASTA for proteomics search. |
| `--bam` | BAM for gene quantification (requires `--gtf`). |
| `--gtf` | GTF or GFF3 annotation. |
| `--fastq` | FASTQ for read-level QC. |
| `--gmt` | Gene-set file in GMT format (name, description, genes…). |

---

## Output files

All files are written to `--output` (default: `./multiomics_out`).

| File | Description |
|---|---|
| `report.html` | Self-contained HTML report with Chart.js charts — open in any browser, no server needed. |
| `multiqc_bioomics.json` | MultiQC-compatible JSON with all summary statistics. |

---

## All CLI flags

### Primary inputs

| Flag | Description |
|---|---|
| `--genomics FILE` | VCF for variant analysis |
| `--transcriptomics FILE` | Expression matrix TSV |
| `--epigenomics FILE` | Bisulfite BED for methylation |
| `--atac FILE` | ATAC-seq narrowPeak |
| `--cnv FILE` | VCF with CN field |
| `--fastq FILE` | FASTQ for QC |
| `--scrna DIR` | 10x MEX directory for single-cell analysis |
| `--bam FILE` | BAM for gene quantification |
| `--gtf FILE` | GTF/GFF3 annotation |

### Proteomics

| Flag | Default | Description |
|---|---|---|
| `--proteomics FILE…` | — | One or more mzML files |
| `--proteomics-dir DIR` | — | Directory of mzML files |
| `--fasta FILE` | — | Protein database FASTA |
| `--proteomics-fdr N` | `0.01` | FDR threshold (0.01 = 1%) |
| `--phospho-max-sites N` | `0` | Max phospho sites per peptide; 0 disables phosphoproteomics |

### Somatic variant calling

| Flag | Default | Description |
|---|---|---|
| `--tumor-bam FILE` | — | Tumor BAM (requires `--normal-bam`) |
| `--normal-bam FILE` | — | Matched normal BAM |
| `--somatic-min-lod N` | `6.3` | Minimum log-odds score for PASS calls |

### Comparison mode (tumor vs. normal / treatment vs. control)

| Flag | Description |
|---|---|
| `--compare-genomics FILE` | Control VCF (enables comparison mode) |
| `--compare-transcriptomics FILE` | Control expression TSV |
| `--compare-epigenomics FILE` | Control methylation BED |
| `--compare-atac FILE` | Control ATAC narrowPeak |

### Single-cell options

| Flag | Default | Description |
|---|---|---|
| `--umap-neighbors N` | `15` | Number of UMAP neighbours |
| `--no-umap` | off | Skip UMAP embedding (clustering still runs) |

### Output and performance

| Flag | Default | Description |
|---|---|---|
| `--output DIR` | `./multiomics_out` | Output directory (created if absent) |
| `--threads N` | all cores | Worker thread count |
| `--json` | off | JSON output only — no TUI, no HTML |
| `--no-ml` | off | Skip PCA and cross-modality correlation |
| `--no-gpu` | off | Disable GPU acceleration; use CPU for UMAP |
| `--raw-counts` | off | Treat `--transcriptomics` as raw counts and apply DESeq2 normalisation |

### Skip flags — disable individual modules

These let you pass input files (e.g. for path validation or future use) without actually running that analysis module.

| Flag | What it skips |
|---|---|
| `--skip-genomics` | Variant analysis |
| `--skip-transcriptomics` | Expression analysis |
| `--skip-epigenomics` | Methylation analysis |
| `--skip-proteomics` | Proteomics database search |
| `--skip-scrna` | Single-cell analysis |

### Presets and configuration

| Flag | Description |
|---|---|
| `--preset NAME` | Load a predefined threshold set |
| `--list-presets` | Print all available presets and exit |
| `--config FILE` | TOML configuration file |
| `--dump-config` | Print the default configuration as TOML and exit |

---

## Skip flags in practice

```bash
# Provide VCF but skip genomics (run only transcriptomics + epigenomics)
multiomics \
  --genomics mutations.vcf \
  --transcriptomics rna.tsv \
  --epigenomics wgbs.bed \
  --skip-genomics

# Full pipeline but skip the proteomics database search
multiomics \
  --genomics tumor.vcf \
  --transcriptomics rna.tsv \
  --epigenomics wgbs.bed \
  --proteomics mass_spec.mzML \
  --fasta proteome.fasta \
  --skip-proteomics

# Single-cell without UMAP (faster QC + normalisation only)
multiomics --scrna mex/ --no-umap

# Force CPU UMAP even when GPU is available
multiomics --scrna mex/ --no-gpu

# Skip ML integration layer (no PCA, no correlation matrix)
multiomics \
  --genomics sample.vcf \
  --transcriptomics rna.tsv \
  --epigenomics wgbs.bed \
  --no-ml
```

---

## Presets

Presets adjust analysis thresholds for common study types without editing a config file.

| Preset | Use case |
|---|---|
| `cancer` | Somatic mutation analysis; lower Ti/Tv warning threshold |
| `plant` | Plant genome; adjusted GC and methylation norms |
| `rna-seq` | Expression-focused; DE-centric insights |
| `wgbs` | Whole-genome bisulfite; strict CpG island parameters |
| `atac` | Chromatin accessibility; ATAC signal thresholds |
| `clinical` | Clinical WGS; conservative QUAL and FDR cutoffs |

```bash
multiomics --genomics tumor.vcf --preset cancer --output ./results
multiomics --list-presets
```

---

## Configuration file

The full set of thresholds and performance knobs can be saved and edited as TOML:

```bash
multiomics --dump-config > my_config.toml
# edit my_config.toml to taste
multiomics --config my_config.toml --genomics sample.vcf ...
```

---

## GPU acceleration details

When compiled with `--features gpu`, UMAP uses a WebGPU compute shader for the k-NN step.

**Algorithm — GPU k-selection**: one GPU thread per cell. Each thread keeps a private max-heap of size k (max 64), scans every other cell, and writes only n×k indices and n×k distances. No n×n distance matrix is ever allocated or transferred.

- Readback size: n × k × 8 bytes — about **12 MB for 100 000 cells at k=15**
- Falls back to CPU automatically when n < 5 000 or no GPU adapter is found

**Platforms**: Vulkan (Linux, Windows), Metal (macOS / Apple Silicon), DX12 (Windows). No CUDA dependency.

**Benchmark on RTX 4050 (6 GB VRAM)**:

| Cells | CPU | GPU (k-selection) |
|---|---|---|
| 2 000 | ~730 ms | CPU fallback |
| 5 000 | ~1.8 s | ~1.6 s |
| 10 000 | ~7 s | ~5 s |
| 20 000 | ~30 s | ~18 s |
| 50 000 | ~200 s | ~90 s |

---

## Feature flags at build time

| Cargo feature | What it adds |
|---|---|
| _(none / default)_ | Genomics, transcriptomics, epigenomics, ATAC-seq, CNV |
| `gpu` | GPU-accelerated UMAP (wgpu / Vulkan / Metal / DX12) |
| `longread` | Long-read (PacBio / Nanopore) epigenomics parsing |
| `full` | ATAC + CNV + longread |
| `full,gpu` | Everything |

---

## Development

```bash
# Verify all crates compile
cargo check --workspace

# Lint — zero warnings policy
cargo clippy --workspace -- -D warnings

# Tests
cargo test --workspace

# UMAP benchmarks (CPU vs GPU)
cd scrna_core
cargo bench                    # CPU only
cargo bench --features gpu     # CPU + GPU comparison
```

---

## License

Apache-2.0 — see [LICENSE](LICENSE).
