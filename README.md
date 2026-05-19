# multiomics

[![License: Apache-2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/Rust-1.80%2B-orange.svg)](https://www.rust-lang.org)
[![Version](https://img.shields.io/badge/version-0.1.0-green.svg)](Cargo.toml)

A parallel multi-omics analysis engine written in Rust. Ingests genomics (VCF/BAM), transcriptomics (TSV/BAM/GTF), epigenomics (BED), ATAC-seq (narrowPeak), proteomics (mzML), single-cell (10x MEX), spatial transcriptomics (Visium), and copy-number variation data simultaneously and produces an integrated HTML report, MultiQC-compatible JSON, and a live terminal dashboard.

2.6× faster than bcftools on VCF parsing. No Python or R runtime required.

---

## Features

**Genomics**
- VCF parsing with Ti/Tv classification and HyperLogLog cardinality estimation (2.6× faster than bcftools)
- Bayesian SNP genotyping — Li 2011 pileup model with Phred LazyLock table and Hardy-Weinberg prior
- GATK-grade local reassembly — De Bruijn k-mer graph (multi-k: 10/15/20/25), Smith-Waterman affine-gap alignment
- Pair-HMM forward algorithm — M/I/D states in log-space for accurate haplotype likelihoods
- Somatic tumor/normal calling — Mutect2-style LOD scoring, strand-bias filter, COSMIC mutation spectrum
- Short-read alignment — BWA-MEM-inspired seed/chain/extend, banded Smith-Waterman, FASTA reference index
- Splice-aware alignment — STAR-style two-pass with exon junction index for RNA-seq reads

**Copy-Number Variation**
- VCF-based CNV: parses `SVTYPE`, `CN`, `CNA`, `TCN`, `LOG2` INFO fields (PURPLE/Sequenza/GATK-SV compatible)
- Coverage-based CNV: CNVkit-style windowed BAM depth, GC correction, change-point segmentation
- Weighted-mean ploidy estimation and fraction-genome-altered (FGA)

**Transcriptomics**
- BAM + GTF gene quantification — featureCounts-style sweep-line overlap, strandedness-aware
- Full DESeq2 NB-GLM — MLE dispersions, parametric trend, MAP shrinkage, IRLS, Wald test, Cook's distance, independent filtering, apeglm LFC shrinkage (Love et al. 2014)
- Multi-factor design — general n×p design matrix, Cholesky IRLS, contrast-based Wald test, LRT

**Epigenomics**
- BED methylation parsing with CpG island detection (Gardiner-Garden criterion, O(n) two-pointer)
- Long-read methylation from Nanopore/PacBio BAM via MM/ML tags (`--features longread`)

**ATAC-seq**
- ENCODE narrowPeak parsing, per-chromosome peak density, open-chromatin quantification
- De novo peak calling — MACS2-style Poisson model, local background lambda, FRiP scoring

**Single-cell RNA-seq**
- 10x Genomics MEX format I/O (matrix.mtx/.gz, barcodes, features)
- MAD-based QC filtering, scran-inspired pooling normalization (Lun 2016)
- Seurat v3 HVG selection, random-projection KNN graph, Leiden clustering (Traag 2019)
- Wilcoxon rank-sum cluster marker detection
- UMAP dimensionality reduction (McInnes 2018) — fuzzy simplicial set + SGD
- Scrublet doublet detection (Wolock 2019), diffusion pseudotime (Haghverdi 2016)
- Harmony batch correction (Korsunsky 2019)
- RNA velocity — spliced/unspliced ratio model (La Manno 2018)

**Spatial Transcriptomics**
- Visium / Slide-seq spot-level analysis
- Spatial autocorrelation (Moran's I), spatially variable gene detection

**Multi-modal Single-cell**
- CITE-seq: antibody-derived tag (ADT) + RNA joint analysis
- Weighted nearest neighbor (WNN) combined embedding (Hao et al. 2021)
- Cell-cell communication: ligand-receptor scoring, communication probability matrix

**Proteomics**
- mzML parsing — streaming quick-xml, base64 + zlib binary array decoding (32-bit and 64-bit float)
- In-silico tryptic digest — K/R|not-P cleavage, 0–2 missed cleavages, length 6–50
- Database search — hyperscore (Craig & Beavis 2004, X!Tandem), precursor 10 ppm, fragment 20 ppm
- 1-Da mass-bin peptide index — O(1) candidate lookup per spectrum (same principle as Sage)
- Target-decoy competition FDR — reversed-sequence decoys, monotone q-value (Elias & Gygi 2007)
- Label-free quantification — MS1 XIC extraction (±5 ppm, ±30 s), trapezoidal peak integration
- Protein inference — PSM → peptide → protein group rollup with protein-level FDR

**Integration**
- Pearson cross-modality correlation, PCA via linfa-reduction
- fgsea multilevel GSEA (Korotkevich 2021) with importance sampling for p < 1e-4
- Pathway enrichment: 75+ built-in KEGG pathways + custom GMT file support
- Rule-based biological insight engine

**Output**
- Single-file HTML report with native SVG plots (no external CDN dependencies)
- MultiQC-compatible JSON (`multiqc_bioomics.json`)
- Live terminal dashboard (ratatui TUI) with per-modality progress bars
- Python bindings via PyO3 / maturin (`pybioomics`)

---

## Installation

**From source**

```bash
git clone https://github.com/diladeniz/multiomics
cd multiomics
cargo build --release --bin multiomics
# Binary at target/release/multiomics
```

**Docker**

```bash
docker pull ghcr.io/diladeniz/multiomics:latest
docker run --rm -v $PWD/data:/data ghcr.io/diladeniz/multiomics \
  --genomics /data/variants.vcf \
  --transcriptomics /data/expression.tsv \
  --epigenomics /data/methylation.bed \
  --output /data/out
```

**Singularity**

```bash
singularity pull multiomics.sif docker://ghcr.io/diladeniz/multiomics:latest
singularity run multiomics.sif \
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
multiomics [OPTIONS]

Primary inputs:
  --genomics <FILE>         VCF for variant analysis
  --transcriptomics <FILE>  Expression TSV or raw count matrix
  --epigenomics <FILE>      BED methylation file
  --atac <FILE>             ENCODE narrowPeak for ATAC-seq
  --cnv <FILE>              VCF with CN INFO field for copy-number analysis
  --bam <FILE>              BAM for gene quantification (requires --gtf)
  --gtf <FILE>              GTF/GFF3 annotation for gene quantification
  --scrna <DIR>             10x MEX directory for single-cell analysis
  --fastq <FILE>            FASTQ for sequence-level QC
  --reference <FILE>        Reference FASTA for read alignment
  --proteomics <FILE>       mzML for proteomics database search (requires --fasta)
  --fasta <FILE>            Protein database FASTA for proteomics search
  --proteomics-fdr <F>      FDR threshold for proteomics reporting [default: 0.01]

Somatic variant calling:
  --tumor-bam <FILE>        Tumor BAM (requires --normal-bam)
  --normal-bam <FILE>       Matched normal BAM
  --somatic-min-lod <F>     Minimum tumor LOD score [default: 6.3]

Enrichment:
  --gmt <FILE>              Custom GMT pathway file

Options:
  --raw-counts              Apply DESeq2 normalization to count matrix
  --preset <NAME>           Threshold preset: cancer, plant, rnaseq, wgbs, atac, clinical
  --output <DIR>            Output directory [default: ./multiomics_out]
  --threads <N>             Worker threads [default: all cores]
  --no-ml                   Skip PCA and correlation
  --json                    JSON output only, no TUI or HTML
```

**Minimal run**

```bash
multiomics \
  --genomics variants.vcf \
  --transcriptomics expression.tsv \
  --epigenomics methylation.bed
```

**Tumor/normal somatic calling**

```bash
multiomics \
  --tumor-bam tumor.bam \
  --normal-bam normal.bam \
  --genomics somatic.vcf \
  --somatic-min-lod 6.3 \
  --preset cancer \
  --output results/
```

**Single-cell with UMAP**

```bash
multiomics \
  --scrna /path/to/10x_mex/ \
  --umap-neighbors 15 \
  --output sc_results/
```

**Gene quantification from BAM**

```bash
multiomics \
  --bam aligned.bam \
  --gtf annotation.gtf \
  --raw-counts \
  --output quant_results/
```

**JSON only (suitable for pipelines)**

```bash
multiomics --genomics v.vcf --transcriptomics t.tsv --epigenomics e.bed \
  --json --output pipeline_out/
```

---

## Presets

| Preset | Use case | Key changes |
|--------|----------|-------------|
| `cancer` | Somatic tumor/normal | QUAL≥20, strict Ti/Tv, hypomethylation alerts |
| `plant` | Plant/agricultural genomics | Relaxed Ti/Tv, lower expressed TPM |
| `rnaseq` | Bulk RNA-seq DE | Strict padj=0.01, top 200 genes |
| `wgbs` | Whole-genome bisulfite | Tight CpG island criteria |
| `atac` | ATAC-seq focus | Signal threshold=5, top 500 peaks |
| `clinical` | Clinical/translational | Conservative thresholds throughout |

---

## Input Formats

### VCF (`--genomics` / `--cnv`)

Standard VCF 4.x. For CNV analysis the INFO field should contain at least one of:
- `SVTYPE=CNV|DEL|DUP|GAIN|LOSS`
- `CN=<integer>`, `CNA=<float>`, or `TCN=<float>`

### Expression TSV (`--transcriptomics`)

Tab-delimited, first row is a header with sample names, first column is gene identifier.

```
gene_id    sample_A    sample_B    sample_C
BRCA1      12.4        0.3         45.1
TP53       8.9         22.1        9.0
```

### Methylation BED (`--epigenomics`)

4-column BED (chr, start, end, methylation%) or ENCODE bisulfite BED. Fifth column = coverage depth.

### narrowPeak (`--atac`)

ENCODE BED6+4 format: chr, start, end, name, score, strand, signalValue, pValue, qValue, peak.

### GMT (`--gmt`)

Standard Gene Matrix Transposed: pathway name, description, gene symbols — one pathway per line.

---

## Output

| File | Description |
|------|-------------|
| `report.html` | Self-contained HTML with native SVG charts and Circos overview |
| `multiqc_bioomics.json` | MultiQC-compatible JSON for pipeline integration |

The HTML report includes: summary cards, variant density bar chart, allele-frequency histogram, volcano plot, top-20 expressed genes, methylation levels, Circos SVG genomic overview, 3×3 correlation heatmap, PCA scatter, pathway enrichment table, insight list.

---

## Architecture

```
biomics_core            — BatchAccum trait, parallel_fold, statistics (Welch, BH, Fisher, HLL)
genomics_core           — VCF parser, Ti/Tv, CNV (VCF+coverage), Bayesian SNP, GATK assembly,
                          pair-HMM, somatic calling, short-read aligner, splice-aware aligner
transcriptomics_core    — TSV/BAM parser, DESeq2 NB-GLM, multi-factor design, BAM+GTF quant
epigenomics_core        — BED parser, CpG island detector, long-read MM/ML parser
atacseq_core            — narrowPeak parser, de novo peak calling (MACS2-style Poisson)
proteomics_core         — mzML parser, tryptic digest, hyperscore search, target-decoy FDR, XIC quant
scrna_core              — 10x MEX I/O, QC, scran normalization, HVG, Leiden, Wilcoxon DE,
                          UMAP, doublets, pseudotime, Harmony, RNA velocity,
                          spatial transcriptomics, CITE-seq WNN, cell-cell communication
integration_layer       — Pearson/Spearman, PCA, fgsea multilevel GSEA, GMT, pathway enrichment
cli                     — clap CLI, ratatui TUI, native SVG plots, HTML/JSON output, Circos SVG
pybioomics              — PyO3 Python bindings (maturin)
```

### Parallelism

All modalities are analysed concurrently via `std::thread::scope`. Within each modality, records are split into chunks of 64 000 distributed across the rayon thread pool. Each worker builds a local accumulator merged at the end via lock-free reduce — no Mutex on the hot path.

### Memory

Files are memory-mapped with `memmap2` + `madvise(Sequential)`. Parsing is zero-copy via `ByteLines`/`TabFields` borrowing directly from the mapped slice. `AHashMap` with AES-NI hardware hashing replaces `std::HashMap` throughout.

---

## Statistics

**Differential expression** — Full DESeq2 NB-GLM (Love et al. 2014): median-of-ratios size factors, gene-wise MLE dispersions, parametric trend α=a₀+a₁/μ, empirical Bayes MAP shrinkage, IRLS, Wald z-test, Cook's distance outlier flagging, Bourgon independent filtering, apeglm LFC shrinkage.

**GSEA** — fgsea multilevel Monte Carlo (Korotkevich 2021): adaptive doubling until SE < 10% of p-estimate, null distribution cache keyed by pathway size. Importance sampling (Owen & Zhou 2000) kicks in automatically for p < 1e-4 for accurate tail estimation.

**Somatic calling** — Mutect2 LOD model: LOD_tumor = Σ log(AF/ε) for alt reads, LOD_normal = Σ log(0.5/ε) for germline-het test. Strand-bias filter via chi-squared approximation of Fisher's exact test.

**CpG island detection** — Gardiner-Garden & Frommer (1987): GC% > 50%, length ≥ 200 bp, CpO/E ≥ 0.6. O(n) two-pointer sliding window.

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

## License

Apache-2.0. See [LICENSE](LICENSE).
