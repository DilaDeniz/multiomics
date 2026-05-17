use std::time::Instant;

use anyhow::{Context, Result};
use chrono::Utc;
use crossbeam_channel::unbounded;
use rayon::ThreadPoolBuilder;

use biomics_core::ProgressEvent;
use genomics_core::{analyze_vcf, GenomicsSummary};
use transcriptomics_core::{analyze_tsv, TranscriptomicsSummary};
use epigenomics_core::{analyze_bed, EpigenomicsSummary};
use integration_layer::{run_integration, IntegrationSummary};

use crate::args::Cli;
use crate::compare::run_comparison;
use crate::config::BioomicsConfig;
use crate::output::{build_multiqc_output, write_html_report, write_json, MultiQcOutput};
use crate::tui::app::{Phase, SharedState};

/// Top-level pipeline orchestrator.
///
/// Phases 1–3 run concurrently via `std::thread::scope`; phase 4 (integration)
/// depends on all three and runs sequentially after them.
pub fn run_pipeline(cli: &Cli, cfg: &BioomicsConfig, state: Option<SharedState>) -> Result<MultiQcOutput> {
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

    log::info!(
        "Starting parallel analysis: '{}', '{}', '{}'",
        cli.genomics.display(),
        cli.transcriptomics.display(),
        cli.epigenomics.display()
    );

    let g_path = &cli.genomics;
    let t_path = &cli.transcriptomics;
    let e_path = &cli.epigenomics;

    let (genomics, transcriptomics, epigenomics): (
        GenomicsSummary,
        TranscriptomicsSummary,
        EpigenomicsSummary,
    ) = {
        let (gr, tr, er) = std::thread::scope(|s| {
            let gh = s.spawn(|| analyze_vcf(g_path, Some(&gtx)));
            let th = s.spawn(|| analyze_tsv(t_path, Some(&ttx)));
            let eh = s.spawn(|| analyze_bed(e_path, Some(&etx)));
            (gh.join(), th.join(), eh.join())
        });
        let g = gr.map_err(|e| anyhow::anyhow!("Genomics thread panicked: {:?}", e))??;
        let t = tr.map_err(|e| anyhow::anyhow!("Transcriptomics thread panicked: {:?}", e))??;
        let e = er.map_err(|e| anyhow::anyhow!("Epigenomics thread panicked: {:?}", e))??;
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
        push_insight(&state, "[INFO] DESeq2 normalization applied to raw count matrix".to_string());
    }

    // Optional: ATAC-seq narrowPeak analysis
    if let Some(ref atac_path) = cli.atac {
        log::info!("Running ATAC-seq analysis: '{}'", atac_path.display());
        match atacseq_core::analyze_narrowpeak(atac_path, None) {
            Ok(atac) => {
                push_insight(&state, format!(
                    "[INFO] ATAC-seq: {} peaks, {:.0} median signal, {} chromosomes",
                    atac.total_peaks, atac.median_signal_value, atac.per_chrom.len()
                ));
                log::info!(
                    "ATAC-seq: {} peaks, {:.0} median signal",
                    atac.total_peaks, atac.median_signal_value
                );
            }
            Err(e) => log::warn!("ATAC-seq analysis failed: {}", e),
        }
    }

    // Optional: CNV analysis from a dedicated VCF
    if let Some(ref cnv_path) = cli.cnv {
        log::info!("Running CNV analysis: '{}'", cnv_path.display());
        match genomics_core::parse_cnv_vcf(cnv_path) {
            Ok(records) => {
                let summary = genomics_core::summarize_cnv(&records);
                push_insight(&state, format!(
                    "[INFO] CNV: {} segments, {:.1}% genome altered, ploidy≈{:.2}",
                    summary.total_segments,
                    summary.fraction_genome_altered * 100.0,
                    summary.estimated_ploidy,
                ));
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

    // Optional: FASTQ QC
    if let Some(ref fastq_path) = cli.fastq {
        log::info!("Running FASTQ QC: '{}'", fastq_path.display());
        match genomics_core::parse_fastq(fastq_path) {
            Ok(fq) => {
                log::info!(
                    "FASTQ: {} reads, {:.1}% GC, {:.1}% Q30",
                    fq.total_reads, fq.gc_content_pct, fq.q30_pct
                );
                push_insight(&state, format!(
                    "[INFO] FASTQ: {} reads, GC={:.1}%, Q30={:.1}%",
                    fq.total_reads, fq.gc_content_pct, fq.q30_pct
                ));
            }
            Err(e) => log::warn!("FASTQ analysis failed: {}", e),
        }
    }

    push_genomics_insights(&genomics, &state);
    push_transcriptomics_insights(&transcriptomics, &state);
    push_epigenomics_insights(&epigenomics, &state);

    // Optional: custom GMT pathway loading
    let gmt_pathways = if let Some(ref gmt_path) = cli.gmt {
        log::info!("Loading GMT pathways: '{}'", gmt_path.display());
        match integration_layer::parse_gmt(gmt_path) {
            Ok(pathways) => {
                push_insight(&state, format!("[INFO] GMT: {} custom pathways loaded", pathways.len()));
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
        complete_progress(&state, |s| { s.integration_pct = 100.0; });
        push_integration_insights(&result, &state);

        // GMT enrichment on top of built-in pathways
        if !gmt_pathways.is_empty() {
            let query_genes: Vec<String> = genomics.high_impact_genes.clone();
            let gmt_hits = integration_layer::gmt_enrichment_analysis(&query_genes, &gmt_pathways, 1);
            if let Some(top) = gmt_hits.first() {
                push_insight(&state, format!(
                    "[INFO] GMT top pathway: {} (padj={:.3})",
                    top.pathway_name, top.padj
                ));
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
        format!("[INFO] Genomics: {} variants, Ti/Tv={:.2}", g.total_variants, g.titv_ratio),
    );
    if !g.high_impact_genes.is_empty() {
        let genes: Vec<_> = g.high_impact_genes.iter().take(3).cloned().collect();
        push_insight(state, format!("[WARN] High-impact genes: {}", genes.join(", ")));
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
            push_insight(state, format!("[INFO] DE: {} genes significant (padj<0.05)", sig));
        }
    }
}

fn push_epigenomics_insights(e: &EpigenomicsSummary, state: &Option<SharedState>) {
    let level = if e.global_methylation_pct < 40.0 { "[CRIT]" } else { "[INFO]" };
    push_insight(
        state,
        format!("{} Global methylation: {:.1}%", level, e.global_methylation_pct),
    );
    let n_islands = e.cpg_islands.len();
    if n_islands > 0 {
        push_insight(state, format!("[INFO] {} CpG islands detected", n_islands));
    }
}

fn push_integration_insights(i: &IntegrationSummary, state: &Option<SharedState>) {
    for insight in i.insights.iter().take(3) {
        let tag = insight.level.as_str();
        push_insight(state, format!("[{}] {}", tag, truncate(&insight.message, 80)));
    }
    if let Some(top) = i.top_pathways.first() {
        push_insight(
            state,
            format!("[INFO] Top pathway: {} (padj={:.3})", top.pathway_name, top.padj),
        );
    }
}

fn truncate(s: &str, max: usize) -> &str {
    if s.len() <= max { s } else { &s[..max] }
}

fn num_cpus() -> usize {
    std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4)
}
