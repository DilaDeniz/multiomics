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
use crate::output::{build_multiqc_output, write_html_report, write_json, MultiQcOutput};
use crate::tui::app::{Phase, SharedState};

/// Top-level pipeline orchestrator.
///
/// ## Parallelism strategy
/// - **Phases 1–3 (Genomics / Transcriptomics / Epigenomics)** run concurrently
///   on three OS threads via `std::thread::scope`. Each spawned thread submits
///   work to the shared rayon pool, so all available CPU cores stay saturated
///   while all three files are being parsed and folded simultaneously.
/// - **Phase 4 (Integration)** runs after phases 1-3 complete — it depends on
///   all three summaries.
pub fn run_pipeline(cli: &Cli, state: Option<SharedState>) -> Result<MultiQcOutput> {
    let start = Instant::now();

    let threads = cli.threads.unwrap_or_else(num_cpus);
    ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .ok(); // ignore error if pool already initialised

    std::fs::create_dir_all(&cli.output)
        .with_context(|| format!("Cannot create output directory '{}'", cli.output.display()))?;

    // Build progress channels before spawning threads (channels are Send).
    let (gtx, grx) = unbounded::<ProgressEvent>();
    let (ttx, trx) = unbounded::<ProgressEvent>();
    let (etx, erx) = unbounded::<ProgressEvent>();

    // Mark all three modality phases as active at once
    set_phase(&state, Phase::Genomics);

    // ── Phases 1–3: parallel I/O + fold ─────────────────────────────────────
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

    // Drain accumulated progress events to TUI state
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

    // Optional FASTQ analysis (genomics augmentation)
    if let Some(ref fastq_path) = cli.fastq {
        log::info!("Running FASTQ QC: '{}'", fastq_path.display());
        match genomics_core::parse_fastq(fastq_path) {
            Ok(fq) => {
                log::info!(
                    "FASTQ: {} reads, {:.1}% GC, {:.1}% Q30",
                    fq.total_reads, fq.gc_content_pct, fq.q30_pct
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

    push_genomics_insights(&genomics, &state);
    push_transcriptomics_insights(&transcriptomics, &state);
    push_epigenomics_insights(&epigenomics, &state);

    // ── Phase 4: Integration ─────────────────────────────────────────────────
    let integration = if cli.no_ml {
        log::info!("Skipping ML integration layer (--no-ml)");
        IntegrationSummary::empty()
    } else {
        set_phase(&state, Phase::Integration);
        log::info!("Running integration layer");
        let result = run_integration(&genomics, &transcriptomics, &epigenomics, false)?;
        complete_progress(&state, |s| { s.integration_pct = 100.0; });
        push_integration_insights(&result, &state);
        result
    };

    // ── Output ───────────────────────────────────────────────────────────────
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

// ── Progress / state helpers ─────────────────────────────────────────────────

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
}

fn push_epigenomics_insights(e: &EpigenomicsSummary, state: &Option<SharedState>) {
    let level = if e.global_methylation_pct < 40.0 { "[CRIT]" } else { "[INFO]" };
    push_insight(
        state,
        format!("{} Global methylation: {:.1}%", level, e.global_methylation_pct),
    );
}

fn push_integration_insights(i: &IntegrationSummary, state: &Option<SharedState>) {
    for insight in i.insights.iter().take(3) {
        let tag = insight.level.as_str();
        push_insight(state, format!("[{}] {}", tag, truncate(&insight.message, 80)));
    }
    if let Some(top) = i.top_pathways.first() {
        push_insight(
            state,
            format!("[INFO] Top pathway: {} (score={:.3})", top.pathway_name, top.score),
        );
    }
}

fn truncate(s: &str, max: usize) -> &str {
    if s.len() <= max { s } else { &s[..max] }
}

fn num_cpus() -> usize {
    std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4)
}
