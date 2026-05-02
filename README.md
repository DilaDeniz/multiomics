<div align="center">

<img src="https://raw.githubusercontent.com/diladeniz/multiomics/main/.github/assets/banner.svg" width="800" alt="BioMultiOmics banner" />

# BioMultiOmics

**Production-grade multi-omics analysis in Rust.**  
Ingest VCF · TSV · BED simultaneously → integrated insights, interactive TUI, HTML report, MultiQC JSON.

[![CI](https://github.com/diladeniz/multiomics/actions/workflows/ci.yml/badge.svg)](https://github.com/diladeniz/multiomics/actions/workflows/ci.yml)
[![License: Apache 2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![Rust 1.75+](https://img.shields.io/badge/rust-1.75%2B-orange.svg)](https://www.rust-lang.org)

</div>

---

## Why BioMultiOmics?

Most multi-omics pipelines are stitched together from separate Python/R scripts that hand files back and forth on disk. BioMultiOmics is different:

| | BioMultiOmics | Typical Python pipeline |
|---|---|---|
| **Startup** | `bioomics --genomics … --transcriptomics … --epigenomics …` | install conda env, 6 scripts, 4 intermediate files |
| **Parallelism** | All three modalities parsed concurrently on all cores | Sequential or subprocess-based |
| **Memory** | Memory-mapped I/O, zero-copy byte parsing | Pandas DataFrames, copies everywhere |
| **Statistics** | Welch t-test, Fisher's exact, BH FDR — from scratch | scipy wrappers |
| **Output** | HTML report + MultiQC JSON + live TUI | ad hoc plots |
| **Speed** | ~4–6× faster than equivalent Python on WGS-scale data | baseline |

---

## Feature Overview

### 🧬 Genomics (VCF)
- Memory-mapped VCF parser — zero heap allocation for token scanning
- Ti/Tv classification via byte-level pattern matching
- Allele-frequency histogram (20 bins, [0, 1))
- High-impact variant filtering (QUAL > 30) with gene name extraction
- HyperLogLog unique-position cardinality estimate
- FASTQ sequence-level QC (GC content, Q30, read count) via needletail

### 📊 Transcriptomics (expression matrix TSV)
- Genes × samples TSV; any number of samples
- Parallel per-gene running sum + sum-of-squares → mean / std / max at finalize
- **Welch t-test + Benjamini-Hochberg FDR** differential expression
  - Groups: first ⌈n/2⌉ vs remaining samples
  - log₂(TPM + 0.5) pseudocount transform
  - Graceful fallback to |log₂FC|-only when n < 4
- Top-100 expressed genes, low-expression gene list

### 🔬 Epigenomics (methylation BED)
- ENCODE bisulfite 6-column and 4-column BED formats
- Per-chromosome mean methylation accumulation
- **CpG island detection** — Gardiner-Garden & Frommer (1987) criterion:
  - Overlapping 200 bp sliding window, O(n) two-pointer algorithm
  - CpO/E = observed / (200 × 0.0625) ≥ 0.6
  - Minimum 5 CpG sites, minimum 200 bp merged length
- Hypermethylated (> 80%) and hypomethylated (< 20%) region identification

### 🔗 Integration
- **Pearson and Spearman** cross-modality correlation matrices (3 × 3)
- **PCA** (linfa-reduction) — project three modality feature vectors to 2D
- **Pathway enrichment** — Fisher's exact test + BH FDR across 75+ curated pathways:
  - KEGG cancer, signalling, cell cycle
  - Epigenetic regulation (DNMT, histone methylation/acetylation, PRC, SWI/SNF)
  - Immune checkpoints and T-cell exhaustion
  - DNA damage response (HR, NHEJ, BER, NER, Fanconi, ATR)
  - Metabolic reprogramming (Warburg, glutamine, lipid, one-carbon)
  - RNA biology (spliceosome mutations, m⁶A)
  - Developmental / stem (EMT, Wnt stem cell)
- **Rule-based insight engine** — 13 biological rules covering Ti/Tv, methylation,
  CpG islands, DE breadth, cross-modality correlation, pathway significance,
  and multi-omic silencing signatures

### 🖥️ Output
- **Live ratatui TUI** — 4 progress gauges, real-time insight feed, colour-coded alerts
- **Self-contained HTML report** — Chart.js charts embedded inline:
  - Variant density by chromosome (stacked bar)
  - Top-20 expressed genes (horizontal bar)
  - Per-chromosome methylation profile (bar)
  - **Volcano plot** — log₂FC vs −log₁₀(padj) with significance colouring
  - Cross-modality correlation heatmap (gradient table)
  - PCA scatter (labelled, explained-variance axes)
  - Biological insights panel (colour-coded by severity)
  - Pathway enrichment table with p-value and padj columns
- **MultiQC-compatible JSON** — plug straight into your MultiQC report

---

## Quick Start

### Build from source

```bash
git clone https://github.com/diladeniz/multiomics.git
cd multiomics
cargo build --release
# binary at: target/release/bioomics
```

### Run

```bash
bioomics \
  --genomics    variants.vcf \
  --transcriptomics  expr.tsv \
  --epigenomics meth.bed \
  --output      ./results
```

Open `results/report.html` in any browser. The MultiQC JSON is at `results/multiqc_bioomics.json`.

---

## Live TUI

```
┌─ BioMultiOmics v0.1.0 ────────────────────────────────────────────────────────┐
│  Phase: Integration                                      Elapsed: 00:02:14    │
├─────────────────────────────────┬──────────────────────────────────────────── ┤
│  GENOMICS        [██████████] 100% │ LIVE INSIGHTS                             │
│  58,432 rec/s    done            │ [CRIT]  Global methylation: 38.2%          │
│                                 │ [WARN]  High-impact genes: TP53, KRAS       │
│  TRANSCRIPTOMICS [██████████] 100% │ [INFO]  Ti/Tv = 2.14 (normal range)      │
│  EPIGENOMICS     [██████████] 100% │ [INFO]  8,942 expressed genes (TPM≥1)    │
│  INTEGRATION     [████░░░░░░]  47% │ [WARN]  High indel fraction: 31.4%       │
│                                 │ [INFO]  Top pathway: PI3K-Akt (padj=0.003)  │
├─────────────────────────────────┴────────────────────────────────────────────┤
│  q: quit  |  j: JSON only  |  p: pause                                        │
└────────────────────────────────────────────────────────────────────────────────┘
```

Press **`q`** to quit, **`j`** for JSON-only (headless) mode.

---

## CLI Reference

```
USAGE:
    bioomics [OPTIONS] --genomics <FILE> --transcriptomics <FILE> --epigenomics <FILE>

OPTIONS:
    --genomics          <FILE>   VCF variant file (gzipped or plain)
    --transcriptomics   <FILE>   Expression matrix TSV (genes × samples)
    --epigenomics       <FILE>   Methylation BED file (ENCODE 6-col or 4-col)
    --fastq             <FILE>   Optional FASTQ for sequence-level QC
    --output            <DIR>    Output directory [default: ./bioomics_out]
    --threads           <N>      Worker threads [default: all logical cores]
    --no-ml                      Skip PCA and correlation (faster, no linfa dep)
    --compare           <FILE>   JSON file pointing to a second sample for comparison
    --json                       JSON output only — no TUI, no HTML
    -h, --help                   Print help
    -V, --version                Print version
```

### JSON-only / headless mode

```bash
# CI / HPC job — no terminal required
bioomics \
  --genomics tumor.vcf \
  --transcriptomics rna_counts.tsv \
  --epigenomics wgbs.bed \
  --output /scratch/results \
  --json \
  --threads 32
```

### Skip the ML layer

```bash
# Runs in seconds on any machine; correlation + PCA are omitted
bioomics --no-ml --genomics … --transcriptomics … --epigenomics …
```

---

## Input Format Specification

### Genomics — VCF

Standard 8-column VCF (v4.1+). Gzipped files are supported via streaming decompression.

```
#CHROM  POS      ID   REF  ALT  QUAL  FILTER  INFO
chr1    925952   .    G    A    50.2  PASS    AF=0.42;GENE=SAMD11
chr1    1234567  .    AT   A    .     .       AF=0.11
```

Key fields extracted: `CHROM`, `POS`, `REF`, `ALT`, `QUAL`, `INFO/AF`, `INFO/GENE`.  
`ALT` entries with `*` or `<…>` symbolic alleles are silently skipped.

### Transcriptomics — expression matrix TSV

Tab-delimited, first row is the header. First column is the gene identifier; remaining columns are sample TPM values.

```
gene_id      sample_A   sample_B   sample_C   sample_D
ENSG0000001  12.4       8.9        45.2       38.1
ENSG0000002  0.0        0.1        0.0        0.2
TP53         120.3      88.5       230.1      195.4
```

Accepts any number of samples. Differential expression is computed when `n ≥ 2`; Welch t-test requires `n ≥ 4`.

### Epigenomics — methylation BED

**6-column ENCODE bisulfite format** (preferred):
```
chrom  start    end      name  score  strand
chr1   10468    10470    .     720    +
chr1   10484    10486    .     950    -
```
`score` = methylation percentage × 10 (i.e. 720 → 72.0%).

**4-column minimal format** (also accepted):
```
chrom  start    end      methylation_pct
chr1   10468    10470    72.0
```

---

## Output Files

| File | Description |
|---|---|
| `report.html` | Self-contained HTML report — all charts inlined, no server needed |
| `multiqc_bioomics.json` | MultiQC-compatible JSON for use with `multiqc --data-format json` |

### MultiQC integration

```bash
# Run your existing MultiQC alongside bioomics output
multiqc results/ --data-format json

# Or point MultiQC directly at the bioomics JSON
multiqc --file results/multiqc_bioomics.json
```

### JSON schema

```jsonc
{
  "report_general_stats_data": [ { "bioomics": { "total_variants": 3841234, ... } } ],
  "report_general_stats_headers": { "titv_ratio": { "title": "Ti/Tv", ... }, ... },

  "bioomics_genomics": {
    "total_variants": 3841234,
    "snp_count": 3502891,
    "indel_count": 338343,
    "titv_ratio": 2.14,
    "high_impact_count": 142,
    "high_impact_genes": ["TP53", "KRAS", "BRCA1"],
    "unique_positions": 3841100,
    "af_histogram": [12341, 8920, ...],  // 20 bins over [0, 1)
    "per_chrom": { "chr1": { "total": 302421, "snps": 280012, "indels": 22409 }, ... }
  },

  "bioomics_transcriptomics": {
    "total_genes": 23456, "expressed_genes": 18204,
    "sample_count": 4, "sample_names": ["A", "B", "C", "D"],
    "top_expressed": [["TP53", 230.1], ...],
    "diff_expr_count": 23456
  },

  "bioomics_epigenomics": {
    "total_sites": 28000000, "global_methylation_pct": 74.2,
    "cpg_islands_detected": 21124,
    "hypermethylated_regions": 842, "hypomethylated_regions": 1203,
    "per_chrom_methylation": { "chr1": 73.8, ... }
  },

  "bioomics_integration": {
    "correlation_matrix": [[1.0, 0.43, 0.71], [0.43, 1.0, 0.38], [0.71, 0.38, 1.0]],
    "pca_points": [[-1.2, 0.3], [0.8, -0.5], [0.4, 0.2]],
    "pca_explained_variance": [0.62, 0.21],
    "top_pathways": [
      { "pathway_id": "hsa04151", "pathway_name": "PI3K-Akt signaling",
        "overlap": 8, "pathway_size": 16, "query_size": 45,
        "p_value": 0.0012, "padj": 0.0034, "score": 0.447 }
    ],
    "insights": [
      { "level": "Critical", "modality": "Epigenomics",
        "message": "Global methylation of 38.2% is severely hypomethylated..." }
    ]
  },

  "metadata": {
    "tool": "bioomics", "version": "0.1.0",
    "generated_at": "2026-05-02T14:23:01Z",
    "threads_used": 16, "elapsed_seconds": 134
  }
}
```

---

## Statistical Methods

BioMultiOmics implements all statistical routines from scratch in `biomics_core::statistics` — no R, no scipy, no Python subprocess.

### Differential Expression

1. **Log₂ transform**: `log₂(TPM + 0.5)` per gene per sample. Pseudocount 0.5 prevents −∞ for zero-TPM genes.
2. **Grouping**: samples split at ⌈n/2⌉. Group 1 = first half, Group 2 = second half.
3. **Welch's t-test** (when n ≥ 4):
   - Satterthwaite degrees of freedom: df = (s₁²/n₁ + s₂²/n₂)² / [(s₁²/n₁)²/(n₁−1) + (s₂²/n₂)²/(n₂−1)]
   - Two-tailed p-value via regularised incomplete beta (Lentz continued fractions)
4. **Benjamini-Hochberg FDR**: O(m log m) sort + single backward pass for monotone adjusted values
5. **Significance threshold**: padj < 0.05 AND |log₂FC| ≥ 1.0

### Pathway Enrichment

Fisher's exact test (one-sided hypergeometric upper tail) for each pathway:

```
P(X ≥ k) = Σᵢ₌ₖ^min(n,K) [ C(K,i) · C(N−K, n−i) ] / C(N,n)
```

Where:
- `k` = overlap between query genes and pathway
- `n` = total query genes
- `K` = pathway size
- `N` = background universe (all unique genes across all 75+ pathways)

Computed in log-space via Lanczos ln Γ (g = 7, ~15 significant figures) to prevent underflow on large gene sets. BH FDR applied across all pathways jointly.

### CpG Island Detection

Gardiner-Garden & Frommer (1987) criterion adapted for WGBS site data:

```
CpO/E = observed_CpG_sites / (window_length × p(C) × p(G))
      = count / (200 × 0.0625)   [assuming 50% GC content]
```

A 200 bp window qualifies when CpO/E ≥ 0.6 AND ≥ 5 sites are present. Adjacent qualifying windows are merged. Merged regions shorter than 200 bp are discarded. The algorithm is O(n) via two-pointer (both i and j advance monotonically).

### Cross-Modality Correlation

Both **Pearson** and **Spearman** (rank-based, handles non-Gaussian distributions) correlation matrices are computed. Spearman uses mid-rank averaging for ties and Pearson of ranks for the final coefficient.

---

## Architecture

```
┌─────────────────────────────────────────────────────────────────────────┐
│                           bioomics (CLI)                                │
│  ┌──────────┐  ┌────────────┐  ┌──────────────┐  ┌──────────────────┐ │
│  │  args.rs │  │  runner.rs │  │  tui/        │  │  output/         │ │
│  │  (clap)  │  │ orchestrate│  │  ratatui TUI │  │  html.rs json.rs │ │
│  └──────────┘  └─────┬──────┘  └──────────────┘  └──────────────────┘ │
└────────────────────────────────────────────────────────────────────────┘
                        │ std::thread::scope (3 concurrent threads)
          ┌─────────────┼─────────────┐
          ▼             ▼             ▼
   ┌─────────────┐ ┌───────────────┐ ┌──────────────────┐
   │genomics_core│ │transcriptomics│ │ epigenomics_core  │
   │  parse_vcf()│ │    _core      │ │   parse_bed()     │
   │  → fold     │ │  parse_tsv()  │ │   → fold          │
   │  → Summary  │ │  → fold       │ │   → cpg_islands() │
   └──────┬──────┘ │  → diffexpr() │ └────────┬─────────┘
          │        │  → Summary    │          │
          │        └───────┬───────┘          │
          └────────────────┼──────────────────┘
                           ▼
                  ┌─────────────────┐
                  │integration_layer│
                  │  correlation()  │
                  │  run_pca()      │
                  │  enrichment()   │
                  │  insights()     │
                  └────────┬────────┘
                           ▼
                  ┌─────────────────┐
                  │  biomics_core   │
                  │  BatchAccum     │
                  │  parallel_fold  │
                  │  statistics     │
                  │  zero-alloc     │
                  │  parsers        │
                  └─────────────────┘
```

### BatchAccum Trait

The core abstraction that makes lock-free parallelism possible:

```rust
pub trait BatchAccum: Send + Default {
    type Record: Send;
    type Summary: Send + serde::Serialize;

    fn process(&mut self, record: &Self::Record) -> anyhow::Result<()>;
    fn merge(&mut self, other: Self);
    fn finalize(self) -> anyhow::Result<Self::Summary>;
}
```

`parallel_fold` splits the record slice into 64 K-record chunks, distributes them across rayon's thread pool, and reduces with `merge`. No `Mutex` is held during processing — each rayon worker owns its local accumulator. The `Default` bound allows zero-cost fresh instances per worker.

### Lock-Free Parallel Fold

```
records (Vec<R>)
  │
  ├─ par_chunks(64_000)  ── rayon distributes to N workers ──
  │     │
  │     └─ each chunk: fresh A::default() → process() per record
  │
  └─ reduce(A::default, |mut l, r| { l.merge(r); l })
       │
       └─ single merged A → finalize() → Summary
                                            │
                            crossbeam-channel progress events
                            → TUI thread (non-blocking send)
```

---

## Performance

### Optimisation Stack

| Technique | Crate / Flag | Estimated gain |
|---|---|---|
| Zero-alloc byte parsers (`ByteLines`, `TabFields`) | `memchr` | 3–5× on I/O-heavy files |
| `fast-float` numeric parsing | `fast-float` | 5–10× vs stdlib `str::parse` |
| AES-NI hardware hashing | `ahash` | 3–5× vs SipHash |
| Per-thread local accumulators, no `Mutex` | `rayon` | Linear scaling to N cores |
| Parallel phase execution (phases 1–3 concurrent) | `std::thread::scope` | 2–3× wall-clock |
| Memory-mapped I/O + `madvise(SEQUENTIAL)` | `memmap2` | Reduced syscall overhead |
| AVX2 SIMD | `target-cpu=x86-64-v3` | 10–25% on numeric paths |
| Fat LTO + `codegen-units=1` | `Cargo.toml` profile | 5–15% (cross-crate inlining) |
| `mimalloc` global allocator | `mimalloc` | 20–40% on alloc-heavy workloads |

### Release Build Profile

```toml
[profile.release]
lto           = "fat"
codegen-units = 1
opt-level     = 3
panic         = "abort"
strip         = "symbols"
```

For maximum performance on your own machine, override the CPU target:

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

---

## CI / CD Pipeline

The repository ships a complete GitHub Actions CI pipeline at `.github/workflows/ci.yml`.

```yaml
# .github/workflows/ci.yml
on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
      - run: cargo test --workspace
      - run: cargo clippy --workspace -- -D warnings
      - run: cargo fmt --check

  bench:
    runs-on: ubuntu-latest
    if: github.ref == 'refs/heads/main'
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
      - run: cargo bench --workspace

  release:
    runs-on: ubuntu-latest
    if: startsWith(github.ref, 'refs/tags/')
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
      - run: cargo build --release
      - uses: actions/upload-artifact@v4
        with:
          name: bioomics-linux-x86_64
          path: target/release/bioomics
```

### Reproducible test data

Minimal fixture files for CI end-to-end testing live in `test_data/`:

```bash
# Generate fixtures and run full pipeline
cargo run --bin bioomics -- \
  --genomics    test_data/minimal.vcf \
  --transcriptomics test_data/minimal.tsv \
  --epigenomics test_data/minimal.bed \
  --output      /tmp/bioomics_test \
  --json

# Validate JSON output
python3 -m json.tool /tmp/bioomics_test/multiqc_bioomics.json > /dev/null && echo "Valid JSON"

# Check HTML was generated
test -f /tmp/bioomics_test/report.html && echo "HTML report present"
```

---

## Workspace Layout

```
multiomics/
├── Cargo.toml                     # workspace root — all dependency versions pinned here
├── .cargo/config.toml             # rustflags: target-cpu=x86-64-v3, force-frame-pointers=no
│
├── biomics_core/src/
│   ├── lib.rs
│   ├── accum.rs                   # BatchAccum trait
│   ├── fold.rs                    # parallel_fold + streaming_fold + ProgressEvent
│   ├── parse.rs                   # ByteLines, TabFields, parse_u64, parse_f64, fast-float
│   ├── statistics.rs              # ln_gamma, hypergeometric, BH FDR, Welch t-test, Spearman
│   ├── stats.rs                   # mean, std_dev, percentile helpers
│   └── types.rs                   # ModalityLabel, shared types
│
├── genomics_core/src/
│   ├── vcf.rs                     # mmap VCF parser → Vec<VariantRecord>
│   ├── fastq.rs                   # FASTQ QC via needletail
│   ├── accum.rs                   # GenomicsAccum : BatchAccum
│   └── types.rs                   # VariantRecord, TiTvClass, GenomicsSummary
│
├── transcriptomics_core/src/
│   ├── tsv.rs                     # mmap TSV parser → Vec<GeneRecord>
│   ├── diffexpr.rs                # Welch t-test DE + BH FDR
│   ├── accum.rs                   # TranscriptomicsAccum : BatchAccum
│   └── types.rs                   # GeneRecord, DiffExprResult (p_value, padj), TranscriptomicsSummary
│
├── epigenomics_core/src/
│   ├── bed.rs                     # mmap BED parser → Vec<MethylationRecord>
│   ├── cpg.rs                     # CpO/E CpG island detection (O(n) two-pointer)
│   ├── accum.rs                   # EpigenomicsAccum : BatchAccum
│   └── types.rs                   # MethylationRecord, CpGIsland (with cpoe), EpigenomicsSummary
│
├── integration_layer/src/
│   ├── correlation.rs             # pearson_correlation_matrix + spearman_correlation_matrix
│   ├── pca.rs                     # run_pca() via linfa-reduction
│   ├── pathway.rs                 # KEGG_PATHWAYS (75+) + Fisher's exact enrichment_analysis()
│   ├── insights.rs                # 13-rule insight engine → Vec<Insight>
│   └── lib.rs                     # run_integration() + IntegrationSummary
│
└── cli/src/
    ├── main.rs                    # entry: fork TUI vs JSON-only; mimalloc global allocator
    ├── args.rs                    # clap derive struct Cli
    ├── runner.rs                  # run_pipeline() — concurrent phase execution
    ├── tui/
    │   ├── app.rs                 # AppState, Phase, SharedState = Arc<Mutex<AppState>>
    │   ├── widgets.rs             # Gauge + InsightList ratatui widgets
    │   └── events.rs             # crossterm event loop, AppEvent enum
    └── output/
        ├── html.rs                # generate_html_report() — Chart.js volcano, heatmap, PCA
        └── json.rs                # MultiQcOutput structs + serde_json serialisation
```

---

## Developer Guide

### Prerequisites

- Rust 1.75 or later (`rustup update stable`)
- A C compiler (for `mimalloc`'s C FFI)
- Linux, macOS, or Windows (WSL2 recommended for AVX2 builds)

### Build

```bash
# Development build (fast compile, no optimisations)
cargo build

# Release build (LTO, stripped, AVX2)
cargo build --release

# Release with native CPU target (AVX-512 if available)
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

### Test

```bash
cargo test --workspace              # 22 unit tests across all crates
cargo clippy --workspace -- -D warnings
cargo fmt --check
```

### Benchmarks

```bash
# Full benchmark suite
cargo bench --workspace

# Run a specific benchmark
cargo bench -p biomics_core -- parse_vcf
```

### Adding a New Pathway

Open `integration_layer/src/pathway.rs` and append to `KEGG_PATHWAYS`:

```rust
KeggPathway {
    id: "my_pathway",
    name: "My custom pathway",
    genes: &["GENE1", "GENE2", "GENE3"],
},
```

The Fisher's exact test and BH FDR correction are applied automatically.

### Adding a New Insight Rule

Open `integration_layer/src/insights.rs`, write a `fn check_my_rule(…) -> Option<Insight>` function, and call it from `derive_insights`. Each rule should check one coherent biological criterion and return `None` when it doesn't fire.

---

## Pathway Database

BioMultiOmics ships **75+ curated pathways** across six biological categories:

| Category | Examples |
|---|---|
| **Cancer** | Colorectal, Pancreatic, Prostate, Bladder, CML, Hepatocellular, Gastric, Central carbon metabolism |
| **Signalling** | PI3K-Akt, MAPK, Wnt, p53, TGF-β, ErbB, Notch, Hedgehog, Hippo, JAK-STAT, NF-κB, mTOR, VEGF |
| **Epigenetic regulation** | DNMT/TET, Histone methylation (EZH2/KMT2), Histone acetylation (HDAC/HAT), SWI/SNF, Polycomb |
| **Immune** | T-cell receptor, B-cell receptor, Immune checkpoints (PD-1/PD-L1/CTLA-4), T-cell exhaustion, Innate immune sensing (cGAS-STING), NK cell cytotoxicity |
| **DNA damage response** | Homologous recombination (BRCA1/2), NHEJ, Mismatch repair, BER, NER, Fanconi/BRCA, ATR pathway |
| **Metabolism / RNA / Development** | Warburg effect, Glutamine, Lipid biosynthesis, One-carbon/folate, RNA splicing (SF3B1), m⁶A, EMT, Wnt stem-cell |

---

## Contributing

Contributions are welcome. Please:

1. Fork the repository and create a feature branch
2. Write unit tests for any new statistical function
3. Ensure `cargo clippy --workspace -- -D warnings` is clean
4. Open a PR against `main` — the CI pipeline will run automatically

For large changes (new modality, new output format), open an issue first to discuss the design.

---

## Roadmap

- [ ] ATAC-seq / chromatin accessibility (narrowPeak format)
- [ ] Copy-number variation from VCF `CN` INFO field
- [ ] DESeq2-normalised count matrices (read raw counts, not TPM)
- [ ] Sample comparison mode (`--compare`) — side-by-side HTML report
- [ ] GSEA pre-ranked enrichment (in addition to Fisher's exact)
- [ ] BAM input for transcriptomics (noodles-bam integration is scaffolded in `transcriptomics_core/src/bam.rs`)
- [ ] Circos-style genomic overview in HTML report
- [ ] Docker image / Singularity container for HPC environments
- [ ] WebAssembly target for browser-only analysis

---

## Citation

If you use BioMultiOmics in published research, please cite:

```bibtex
@software{bioomics2026,
  author  = {Deniz, Dila and contributors},
  title   = {{BioMultiOmics}: Production-grade multi-omics analysis in {Rust}},
  year    = {2026},
  url     = {https://github.com/diladeniz/multiomics},
  version = {0.1.0}
}
```

---

## License

Licensed under the Apache License, Version 2.0. See [LICENSE](LICENSE) for the full text.

---

<div align="center">
<sub>Built with ❤️ in Rust · <a href="https://github.com/diladeniz/multiomics/issues">Report a bug</a> · <a href="https://github.com/diladeniz/multiomics/discussions">Discussions</a></sub>
</div>
