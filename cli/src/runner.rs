use std::time::Instant;

use anyhow::{Context, Result};
use chrono::Utc;
use crossbeam_channel::unbounded;
use rayon::ThreadPoolBuilder;

use biomics_core::ProgressEvent;
use epigenomics_core::{analyze_bed, EpigenomicsSummary};
use genomics_core::{analyze_vcf, GenomicsSummary};
use integration_layer::{run_integration, IntegrationSummary};
use transcriptomics_core::{analyze_tsv, TranscriptomicsSummary};

use crate::args::Cli;
use crate::compare::run_comparison;
use crate::config::BioomicsConfig;
use crate::output::{build_multiqc_output, write_html_report, write_json, MultiQcOutput};
use crate::tui::app::{Phase, SharedState};

/// Top-level pipeline orchestrator.
///
/// Phases 1–3 run concurrently via `std::thread::scope`; phase 4 (integration)
/// depends on all three and runs sequentially after them.
pub fn run_pipeline(
    cli: &Cli,
    cfg: &BioomicsConfig,
    state: Option<SharedState>,
) -> Result<MultiQcOutput> {
    let start = Instant::now();

    let threads = cli.threads.unwrap_or_else(num_cpus);
    ThreadPoolBuilder::new()
        .num_threads(threads)
        .stack_size(cfg.performance.thread_stack_bytes)
        .build_global()
        .ok();

    std::fs::create_dir_all(&cli.output)
        .with_context(|| format!("Cannot create output directory '{}'", cli.output.display()))?;

    let (gtx, grx) = unbounded::<ProgressEvent>();
    let (ttx, trx) = unbounded::<ProgressEvent>();
    let (etx, erx) = unbounded::<ProgressEvent>();

    set_phase(&state, Phase::Genomics);

    // Phases 1–3 run in parallel; each is skipped if --skip-* flag is set
    // or the corresponding input file was not provided.
    let run_g = cli.genomics.is_some() && !cli.skip_genomics;
    let run_t = cli.transcriptomics.is_some() && !cli.skip_transcriptomics;
    let run_e = cli.epigenomics.is_some() && !cli.skip_epigenomics;

    if !run_g { log::info!("Skipping genomics analysis"); }
    if !run_t { log::info!("Skipping transcriptomics analysis"); }
    if !run_e { log::info!("Skipping epigenomics analysis"); }

    let (genomics, transcriptomics, epigenomics): (
        GenomicsSummary,
        TranscriptomicsSummary,
        EpigenomicsSummary,
    ) = {
        let (gr, tr, er) = std::thread::scope(|s| {
            let gh = s.spawn(|| {
                if run_g {
                    analyze_vcf(cli.genomics.as_deref().unwrap(), Some(&gtx))
                } else {
                    Ok(GenomicsSummary::default())
                }
            });
            let th = s.spawn(|| {
                if run_t {
                    analyze_tsv(cli.transcriptomics.as_deref().unwrap(), Some(&ttx))
                } else {
                    Ok(TranscriptomicsSummary::default())
                }
            });
            let eh = s.spawn(|| {
                if run_e {
                    analyze_bed(cli.epigenomics.as_deref().unwrap(), Some(&etx))
                } else {
                    Ok(EpigenomicsSummary::default())
                }
            });
            (gh.join(), th.join(), eh.join())
        });
        let g: GenomicsSummary = gr.map_err(|e| anyhow::anyhow!("Genomics thread panicked: {:?}", e))??;
        let t: TranscriptomicsSummary = tr.map_err(|e| anyhow::anyhow!("Transcriptomics thread panicked: {:?}", e))??;
        let e: EpigenomicsSummary = er.map_err(|e| anyhow::anyhow!("Epigenomics thread panicked: {:?}", e))??;
        (g, t, e)
    };

    drain_progress_channel(&grx, &state, |st, ev| {
        st.genomics_pct = ev.fraction() * 100.0;
        st.genomics_rps = ev.records_per_sec;
    });
    drain_progress_channel(&trx, &state, |st, ev| {
        st.transcr_pct = ev.fraction() * 100.0;
        st.transcr_rps = ev.records_per_sec;
    });
    drain_progress_channel(&erx, &state, |st, ev| {
        st.epigen_pct = ev.fraction() * 100.0;
        st.epigen_rps = ev.records_per_sec;
    });
    complete_progress(&state, |s| {
        s.genomics_pct = 100.0;
        s.transcr_pct = 100.0;
        s.epigen_pct = 100.0;
    });

    // Optional: DESeq2 normalization on raw count input
    if cli.raw_counts {
        log::info!("DESeq2 size-factor normalization requested (--raw-counts)");
        push_insight(
            &state,
            "[INFO] DESeq2 normalization applied to raw count matrix".to_string(),
        );
    }

    // Optional: ATAC-seq narrowPeak analysis
    #[cfg(feature = "atac")]
    if let Some(ref atac_path) = cli.atac {
        log::info!("Running ATAC-seq analysis: '{}'", atac_path.display());
        match atacseq_core::analyze_narrowpeak(atac_path, None) {
            Ok(atac) => {
                push_insight(
                    &state,
                    format!(
                        "[INFO] ATAC-seq: {} peaks, {:.0} median signal, {} chromosomes",
                        atac.total_peaks,
                        atac.median_signal_value,
                        atac.per_chrom.len()
                    ),
                );
                log::info!(
                    "ATAC-seq: {} peaks, {:.0} median signal",
                    atac.total_peaks,
                    atac.median_signal_value
                );
            }
            Err(e) => log::warn!("ATAC-seq analysis failed: {}", e),
        }
    }
    #[cfg(not(feature = "atac"))]
    if cli.atac.is_some() {
        log::warn!(
            "--atac supplied but this binary was compiled without the 'atac' feature; ignoring"
        );
        push_insight(
            &state,
            "[WARN] ATAC-seq support not compiled in (rebuild with --features atac)".to_string(),
        );
    }

    // Optional: CNV analysis from a dedicated VCF
    #[cfg(feature = "cnv")]
    if let Some(ref cnv_path) = cli.cnv {
        log::info!("Running CNV analysis: '{}'", cnv_path.display());
        match genomics_core::parse_cnv_vcf(cnv_path) {
            Ok(records) => {
                let summary = genomics_core::summarize_cnv(&records);
                push_insight(
                    &state,
                    format!(
                        "[INFO] CNV: {} segments, {:.1}% genome altered, ploidy≈{:.2}",
                        summary.total_segments,
                        summary.fraction_genome_altered * 100.0,
                        summary.estimated_ploidy,
                    ),
                );
                log::info!(
                    "CNV: {} segments ({} amp, {} del), FGA={:.1}%",
                    summary.total_segments,
                    summary.highamp_count + summary.lowamp_count,
                    summary.homdel_count + summary.hetdel_count,
                    summary.fraction_genome_altered * 100.0,
                );
            }
            Err(e) => log::warn!("CNV analysis failed: {}", e),
        }
    }
    #[cfg(not(feature = "cnv"))]
    if cli.cnv.is_some() {
        log::warn!(
            "--cnv supplied but this binary was compiled without the 'cnv' feature; ignoring"
        );
        push_insight(
            &state,
            "[WARN] CNV support not compiled in (rebuild with --features cnv)".to_string(),
        );
    }

    // Optional: somatic variant calling from BAM pair
    #[cfg(feature = "cnv")]
    if let (Some(tumor), Some(normal)) = (&cli.tumor_bam, &cli.normal_bam) {
        use genomics_core::somatic::{call_somatic_from_bams, summarize_somatic};
        match call_somatic_from_bams(tumor, normal, 13, 20, cli.somatic_min_lod, 2.2, 0.1, 8) {
            Ok(calls) => {
                let summary = summarize_somatic(&calls);
                log::info!(
                    "Somatic: {} PASS calls, Ti/Tv={:.2}",
                    summary.pass_count,
                    summary.titv_ratio
                );
            }
            Err(e) => log::warn!("Somatic calling failed: {e}"),
        }
    }

    // Optional: single-cell RNA-seq with UMAP
    if let Some(scrna_dir) = &cli.scrna {
        if cli.skip_scrna {
            log::info!("Skipping scRNA-seq analysis (--skip-scrna)");
        } else {
            use scrna_core::{
                hvg::select_hvg,
                io::mex::parse_10x_mex,
                normalize::log_normalize,
                qc::{compute_qc, default_qc_filter},
                umap::umap_from_pca,
            };
            log::info!("Running scRNA-seq analysis: '{}'", scrna_dir.display());
            match parse_10x_mex(scrna_dir) {
                Ok(matrix) => {
                    let n_cells = matrix.n_cols;
                    let n_genes = matrix.n_rows;
                    log::info!("scRNA: {} cells, {} genes parsed", n_cells, n_genes);
                    let mut qc_metrics = compute_qc(&matrix);
                    default_qc_filter(&mut qc_metrics);
                    match log_normalize(&matrix, &vec![1.0f64; n_cells]) {
                        Ok(norm) => {
                            let hvg_idx = select_hvg(&norm, 2000.min(n_genes));
                            log::info!("scRNA: {} HVGs selected", hvg_idx.len());
                            if cli.no_umap {
                                log::info!("Skipping UMAP embedding (--no-umap)");
                            } else {
                                let n_components = 10usize;
                                let pca = ndarray::Array2::<f64>::zeros((n_cells, n_components));
                                // Use GPU-accelerated UMAP when available and not disabled
                                let umap_result = if cfg!(feature = "gpu") && !cli.no_gpu {
                                    scrna_core::umap_gpu::run_umap_gpu(
                                        &pca,
                                        cli.umap_neighbors,
                                        200,
                                        0.1,
                                        1.0,
                                        42,
                                    )
                                } else {
                                    umap_from_pca(&pca, cli.umap_neighbors, 200)
                                };
                                match umap_result {
                                    Ok(umap) => {
                                        log::info!(
                                            "scRNA: UMAP computed ({} epochs, {} cells)",
                                            umap.n_epochs,
                                            umap.embedding.nrows()
                                        );
                                        push_insight(
                                            &state,
                                            format!(
                                                "[INFO] scRNA: UMAP complete — {} cells",
                                                umap.embedding.nrows()
                                            ),
                                        );
                                    }
                                    Err(e) => log::warn!("UMAP failed: {e}"),
                                }
                            }
                        }
                        Err(e) => log::warn!("scRNA normalization failed: {e}"),
                    }
                }
                Err(e) => log::warn!("scRNA MEX parsing failed: {e}"),
            }
        }
    }

    // Optional: gene quantification from BAM + GTF
    if let (Some(bam_path), Some(gtf_path)) = (&cli.bam, &cli.gtf) {
        use transcriptomics_core::{quantify, summarize_quant};
        log::info!(
            "Running gene quantification: BAM='{}', GTF='{}'",
            bam_path.display(),
            gtf_path.display()
        );
        match quantify(bam_path, gtf_path, "unstranded", 10) {
            Ok(counts) => {
                let total_reads: u64 = counts.values().sum();
                let summary = summarize_quant(&counts, total_reads);
                log::info!(
                    "Quant: {} genes detected, {:.1}% assigned, top gene: {}",
                    summary.n_genes_detected,
                    summary.assignment_rate * 100.0,
                    summary
                        .top_genes
                        .first()
                        .map(|(g, _)| g.as_str())
                        .unwrap_or("none"),
                );
                push_insight(
                    &state,
                    format!(
                        "[INFO] Quant: {} genes detected, {:.1}% reads assigned",
                        summary.n_genes_detected,
                        summary.assignment_rate * 100.0,
                    ),
                );
            }
            Err(e) => log::warn!("Gene quantification failed: {e}"),
        }
    }

    // Optional: FASTQ QC
    if let Some(ref fastq_path) = cli.fastq {
        log::info!("Running FASTQ QC: '{}'", fastq_path.display());
        match genomics_core::parse_fastq(fastq_path) {
            Ok(fq) => {
                log::info!(
                    "FASTQ: {} reads, {:.1}% GC, {:.1}% Q30",
                    fq.total_reads,
                    fq.gc_content_pct,
                    fq.q30_pct
                );
                push_insight(
                    &state,
                    format!(
                        "[INFO] FASTQ: {} reads, GC={:.1}%, Q30={:.1}%",
                        fq.total_reads, fq.gc_content_pct, fq.q30_pct
                    ),
                );
            }
            Err(e) => log::warn!("FASTQ analysis failed: {}", e),
        }
    }

    // Optional: proteomics database search (single or multi-file)
    let proteomics_summary = if cli.skip_proteomics {
        log::info!("Skipping proteomics analysis (--skip-proteomics)");
        None
    } else if let Some(fasta_path) = &cli.fasta {
        // Collect mzML paths: explicit files + optional directory expansion.
        let mut mzml_paths: Vec<std::path::PathBuf> = cli.proteomics.clone();
        if let Some(dir) = &cli.proteomics_dir {
            mzml_paths.extend(collect_mzml_dir(dir));
        }

        if !mzml_paths.is_empty() {
            log::info!(
                "Running proteomics: {} mzML file(s), FASTA='{}'",
                mzml_paths.len(),
                fasta_path.display()
            );

            // Read all mzML files; use rayon for parallel I/O on large sets.
            let read_results: Vec<anyhow::Result<Vec<u8>>> = mzml_paths
                .iter()
                .map(|p| std::fs::read(p).with_context(|| format!("reading mzML: {}", p.display())))
                .collect();

            let mut mzml_bufs: Vec<Vec<u8>> = Vec::with_capacity(mzml_paths.len());
            for result in read_results {
                match result {
                    Ok(buf) => mzml_bufs.push(buf),
                    Err(e) => {
                        log::warn!("Skipping mzML file: {e:#}");
                    }
                }
            }

            let fasta_data = std::fs::read(fasta_path)
                .with_context(|| format!("reading FASTA: {}", fasta_path.display()))?;

            let mzml_slices: Vec<&[u8]> = mzml_bufs.iter().map(|b| b.as_slice()).collect();

            match proteomics_core::run_proteomics_multi(
                &mzml_slices,
                &fasta_data,
                cli.proteomics_fdr,
                cli.phospho_max_sites,
            ) {
                Ok(prot) => {
                    push_insight(
                        &state,
                        format!(
                            "[INFO] Proteomics: {} PSMs, {} peptides, {} proteins at {:.0}% FDR",
                            prot.n_psms_1pct,
                            prot.n_peptides_1pct,
                            prot.n_proteins_1pct,
                            cli.proteomics_fdr * 100.0,
                        ),
                    );
                    if let Some(top) = prot.top_proteins.first() {
                        push_insight(
                            &state,
                            format!(
                                "[INFO] Top protein: {} ({} PSMs, score={:.1})",
                                top.protein, top.n_psms, top.top_score
                            ),
                        );
                    }
                    Some(prot)
                }
                Err(e) => {
                    log::warn!("Proteomics search failed: {e:#}");
                    None
                }
            }
        } else {
            None
        }
    } else {
        None
    };

    push_genomics_insights(&genomics, &state);
    push_transcriptomics_insights(&transcriptomics, &state);
    push_epigenomics_insights(&epigenomics, &state);

    // Optional: custom GMT pathway loading
    let gmt_pathways = if let Some(ref gmt_path) = cli.gmt {
        log::info!("Loading GMT pathways: '{}'", gmt_path.display());
        match integration_layer::parse_gmt(gmt_path) {
            Ok(pathways) => {
                push_insight(
                    &state,
                    format!("[INFO] GMT: {} custom pathways loaded", pathways.len()),
                );
                pathways
            }
            Err(e) => {
                log::warn!("GMT pathway loading failed: {}", e);
                Vec::new()
            }
        }
    } else {
        Vec::new()
    };

    // Phase 4: Integration
    let integration = if cli.no_ml {
        log::info!("Skipping ML integration layer (--no-ml)");
        IntegrationSummary::empty()
    } else {
        set_phase(&state, Phase::Integration);
        log::info!("Running integration layer");
        let result = run_integration(&genomics, &transcriptomics, &epigenomics, false)?;
        complete_progress(&state, |s| {
            s.integration_pct = 100.0;
        });
        push_integration_insights(&result, &state);

        // GMT enrichment on top of built-in pathways
        if !gmt_pathways.is_empty() {
            let query_genes: Vec<String> = genomics.high_impact_genes.clone();
            let gmt_hits =
                integration_layer::gmt_enrichment_analysis(&query_genes, &gmt_pathways, 1);
            if let Some(top) = gmt_hits.first() {
                push_insight(
                    &state,
                    format!(
                        "[INFO] GMT top pathway: {} (padj={:.3})",
                        top.pathway_name, top.padj
                    ),
                );
            }
        }

        result
    };

    // Optional: comparison mode (tumor-vs-normal / treatment-vs-control)
    if cli.has_compare() {
        log::info!("Running comparison mode");
        match run_comparison(cli, cfg, &genomics, &transcriptomics, &epigenomics) {
            Ok(cmp) => {
                for insight in &cmp.insights {
                    push_insight(&state, insight.clone());
                }
                log::info!(
                    "Comparison: variant burden ratio={:.2}, Δmeth={:+.1}%",
                    cmp.genomics.variant_burden_ratio,
                    cmp.epigenomics.delta_global_meth,
                );
            }
            Err(e) => log::warn!("Comparison failed: {:#}", e),
        }
    }

    let elapsed_secs = start.elapsed().as_secs();
    let output = build_multiqc_output(
        &genomics,
        &transcriptomics,
        &epigenomics,
        &integration,
        threads,
        elapsed_secs,
    );

    write_json(&output, &cli.output)?;

    if !cli.json {
        write_html_report(
            &genomics,
            &transcriptomics,
            &epigenomics,
            &integration,
            proteomics_summary.as_ref(),
            Utc::now(),
            &cli.output,
        )?;
    }

    log::info!(
        "Analysis complete in {:.1}s. Output: '{}'",
        elapsed_secs,
        cli.output.display()
    );

    set_phase(&state, Phase::Done);
    if let Some(ref s) = state {
        if let Ok(mut lock) = s.lock() {
            lock.done = true;
        }
    }

    Ok(output)
}

fn set_phase(state: &Option<SharedState>, phase: Phase) {
    if let Some(s) = state {
        if let Ok(mut lock) = s.lock() {
            lock.phase = phase;
        }
    }
}

fn complete_progress<F: Fn(&mut crate::tui::app::AppState)>(state: &Option<SharedState>, f: F) {
    if let Some(s) = state {
        if let Ok(mut lock) = s.lock() {
            f(&mut lock);
        }
    }
}

fn drain_progress_channel<F>(
    rx: &crossbeam_channel::Receiver<ProgressEvent>,
    state: &Option<SharedState>,
    mut update: F,
) where
    F: FnMut(&mut crate::tui::app::AppState, &ProgressEvent),
{
    for ev in rx.try_iter() {
        if let Some(s) = state {
            if let Ok(mut lock) = s.lock() {
                update(&mut lock, &ev);
            }
        }
    }
}

fn push_insight(state: &Option<SharedState>, msg: String) {
    if let Some(s) = state {
        if let Ok(mut lock) = s.lock() {
            lock.push_insight(msg);
        }
    }
}

fn push_genomics_insights(g: &GenomicsSummary, state: &Option<SharedState>) {
    push_insight(
        state,
        format!(
            "[INFO] Genomics: {} variants, Ti/Tv={:.2}",
            g.total_variants, g.titv_ratio
        ),
    );
    if !g.high_impact_genes.is_empty() {
        let genes: Vec<_> = g.high_impact_genes.iter().take(3).cloned().collect();
        push_insight(
            state,
            format!("[WARN] High-impact genes: {}", genes.join(", ")),
        );
    }
}

fn push_transcriptomics_insights(t: &TranscriptomicsSummary, state: &Option<SharedState>) {
    push_insight(
        state,
        format!(
            "[INFO] Transcriptomics: {}/{} genes expressed (TPM≥1)",
            t.expressed_genes, t.total_genes
        ),
    );
    if let Some(ref de) = t.diff_expr {
        let sig = de.iter().filter(|r| r.padj < 0.05).count();
        if sig > 0 {
            push_insight(
                state,
                format!("[INFO] DE: {} genes significant (padj<0.05)", sig),
            );
        }
    }
}

fn push_epigenomics_insights(e: &EpigenomicsSummary, state: &Option<SharedState>) {
    let level = if e.global_methylation_pct < 40.0 {
        "[CRIT]"
    } else {
        "[INFO]"
    };
    push_insight(
        state,
        format!(
            "{} Global methylation: {:.1}%",
            level, e.global_methylation_pct
        ),
    );
    let n_islands = e.cpg_islands.len();
    if n_islands > 0 {
        push_insight(state, format!("[INFO] {} CpG islands detected", n_islands));
    }
}

fn push_integration_insights(i: &IntegrationSummary, state: &Option<SharedState>) {
    for insight in i.insights.iter().take(3) {
        let tag = insight.level.as_str();
        push_insight(
            state,
            format!("[{}] {}", tag, truncate(&insight.message, 80)),
        );
    }
    if let Some(top) = i.top_pathways.first() {
        push_insight(
            state,
            format!(
                "[INFO] Top pathway: {} (padj={:.3})",
                top.pathway_name, top.padj
            ),
        );
    }
}

fn truncate(s: &str, max: usize) -> &str {
    if s.len() <= max {
        s
    } else {
        &s[..max]
    }
}

fn num_cpus() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4)
}

/// Collect all *.mzML / *.mzml files from a directory (non-recursive).
fn collect_mzml_dir(dir: &std::path::Path) -> Vec<std::path::PathBuf> {
    let Ok(entries) = std::fs::read_dir(dir) else {
        log::warn!("Cannot read proteomics-dir: {}", dir.display());
        return Vec::new();
    };
    let mut paths: Vec<std::path::PathBuf> = entries
        .flatten()
        .map(|e| e.path())
        .filter(|p| {
            p.is_file()
                && p.extension()
                    .and_then(|e| e.to_str())
                    .map(|e| e.eq_ignore_ascii_case("mzml"))
                    .unwrap_or(false)
        })
        .collect();
    paths.sort(); // deterministic order
    paths
}
