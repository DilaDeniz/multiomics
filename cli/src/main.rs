mod args;
mod compare;
mod config;
mod context_detect;
mod output;
mod runner;
mod tui;

use anyhow::Result;
use clap::Parser;

use args::Cli;
use config::{dump_default_config, load_config_onto, preset_config, BioomicsConfig};
use tui::new_shared_state;

// Replace the default system allocator with mimalloc.
// Benchmarks on allocation-heavy workloads show 20-40% wall-clock improvement.
#[global_allocator]
static ALLOC: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// Print which optional modalities are compiled into this binary.
fn print_compiled_features() {
    println!("Compiled-in optional modalities:");
    #[cfg(feature = "atac")]
    println!("  [+] atac     — ATAC-seq narrowPeak analysis");
    #[cfg(not(feature = "atac"))]
    println!("  [-] atac     — not compiled (rebuild with --features atac)");

    #[cfg(feature = "cnv")]
    println!("  [+] cnv      — Copy-number variation from VCF");
    #[cfg(not(feature = "cnv"))]
    println!("  [-] cnv      — not compiled (rebuild with --features cnv)");

    #[cfg(feature = "longread")]
    println!("  [+] longread — Long-read (PacBio/Nanopore) epigenomics");
    #[cfg(not(feature = "longread"))]
    println!("  [-] longread — not compiled (rebuild with --features longread)");
}

/// Log warnings for optional flags that need a feature not compiled in.
fn warn_disabled_features(_cli: &Cli) {
    #[cfg(not(feature = "atac"))]
    if _cli.atac.is_some() {
        log::warn!(
            "--atac supplied but multiomics was compiled without the 'atac' feature (rebuild with --features atac)"
        );
    }

    #[cfg(not(feature = "cnv"))]
    if _cli.cnv.is_some() {
        log::warn!(
            "--cnv supplied but multiomics was compiled without the 'cnv' feature (rebuild with --features cnv)"
        );
    }
}

/// Check that at least one analysis input was supplied and that every supplied
/// path exists. Each modality is optional; the runner skips absent ones.
fn validate_inputs(cli: &Cli) -> Result<()> {
    let mut paths: Vec<(&str, &std::path::Path)> = [
        ("--genomics", cli.genomics.as_deref()),
        ("--transcriptomics", cli.transcriptomics.as_deref()),
        ("--epigenomics", cli.epigenomics.as_deref()),
        ("--atac", cli.atac.as_deref()),
        ("--cnv", cli.cnv.as_deref()),
        ("--fastq", cli.fastq.as_deref()),
        ("--scrna", cli.scrna.as_deref()),
        ("--bam", cli.bam.as_deref()),
        ("--tumor-bam", cli.tumor_bam.as_deref()),
        ("--proteomics-dir", cli.proteomics_dir.as_deref()),
    ]
    .into_iter()
    .filter_map(|(label, p)| p.map(|p| (label, p)))
    .collect();
    paths.extend(cli.proteomics.iter().map(|p| ("--proteomics", p.as_path())));

    if paths.is_empty() {
        anyhow::bail!(
            "no input supplied: pass at least one of --genomics, --transcriptomics, \
             --epigenomics, --atac, --cnv, --fastq, --scrna, --proteomics/--proteomics-dir, \
             --bam (with --gtf) or --tumor-bam (with --normal-bam)"
        );
    }

    let companions = [
        ("--gtf", cli.gtf.as_deref()),
        ("--normal-bam", cli.normal_bam.as_deref()),
        ("--fasta", cli.fasta.as_deref()),
        ("--reference", cli.reference.as_deref()),
        ("--gmt", cli.gmt.as_deref()),
        ("--config", cli.config.as_deref()),
        ("--compare-genomics", cli.compare_genomics.as_deref()),
        (
            "--compare-transcriptomics",
            cli.compare_transcriptomics.as_deref(),
        ),
        ("--compare-epigenomics", cli.compare_epigenomics.as_deref()),
        ("--compare-atac", cli.compare_atac.as_deref()),
    ];
    for (label, p) in paths.into_iter().chain(
        companions
            .into_iter()
            .filter_map(|(l, p)| p.map(|p| (l, p))),
    ) {
        if !p.exists() {
            anyhow::bail!("Input for {} not found: '{}'", label, p.display());
        }
    }

    let has_proteomics = !cli.proteomics.is_empty() || cli.proteomics_dir.is_some();
    if has_proteomics && cli.fasta.is_none() {
        anyhow::bail!("--proteomics/--proteomics-dir requires --fasta");
    }
    if cli.bam.is_some() != cli.gtf.is_some() {
        anyhow::bail!("--bam and --gtf must be supplied together");
    }
    if cli.tumor_bam.is_some() != cli.normal_bam.is_some() {
        anyhow::bail!("--tumor-bam and --normal-bam must be supplied together");
    }
    if let Some(mb) = cli.tmb_genome_mb {
        if !(mb.is_finite() && mb > 0.0) {
            anyhow::bail!("--tmb-genome-mb must be a positive number of megabases (got {mb})");
        }
    }
    Ok(())
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"))
        .format_timestamp_millis()
        .init();

    let cli = Cli::parse();

    // --dump-config: print defaults and exit (no input files required)
    if cli.dump_config {
        print!("{}", dump_default_config());
        return Ok(());
    }

    // --list-presets: show available presets and exit
    if cli.list_presets {
        println!("Available presets:\n");
        println!(
            "  {:<10} Somatic mutation analysis (tumor/normal, strict thresholds)",
            "cancer"
        );
        println!(
            "  {:<10} Plant/agricultural genomics (relaxed Ti/Tv, plant-specific)",
            "plant"
        );
        println!(
            "  {:<10} Bulk RNA-seq DE focus (strict padj, high gene count)",
            "rna-seq"
        );
        println!(
            "  {:<10} Whole-genome bisulfite sequencing (tight CpG island criteria)",
            "wgbs"
        );
        println!(
            "  {:<10} ATAC-seq chromatin accessibility (peak signal tuned)",
            "atac"
        );
        println!(
            "  {:<10} Clinical/translational (conservative, fewer false positives)",
            "clinical"
        );
        println!();
        println!("Use: multiomics --preset cancer --genomics ...");
        println!("Combine with --config to override specific fields:");
        println!("  multiomics --preset cancer --config my_tweaks.toml --genomics ...");
        println!();
        print_compiled_features();
        return Ok(());
    }

    // Load configuration: start with preset (if any), then overlay the config
    // file so that explicit file settings always win over preset defaults.
    let cfg = {
        let base = if let Some(preset) = cli.preset {
            preset_config(preset)
        } else {
            BioomicsConfig::default()
        };
        if let Some(ref path) = cli.config {
            load_config_onto(path, base)?
        } else {
            base
        }
    };

    validate_inputs(&cli)?;

    // Warn at startup about any disabled-feature flags that were passed
    warn_disabled_features(&cli);

    // Feature #3: Automatic sample context detection.
    // Scan the first N records of each input file to infer species, assay,
    // and suggest a preset — without reading the whole file.
    let sample_context = context_detect::detect_sample_context(
        cli.genomics.as_deref(),
        cli.transcriptomics.as_deref(),
        cli.epigenomics.as_deref(),
        cli.atac.as_deref(),
    );

    // Log auto-detected context
    if let Some(ref suggestion) = sample_context.suggested_preset {
        let assay_str = sample_context
            .genomics_assay
            .as_ref()
            .map(|a| a.as_str())
            .unwrap_or("");
        log::info!(
            "Auto-detected context: species={}, assay={} — suggested preset: {} (confidence={:.0}%)",
            sample_context.species.as_str(),
            assay_str,
            suggestion.preset,
            suggestion.confidence * 100.0
        );
    }
    for w in &sample_context.concordance.warnings {
        log::warn!("{}", w);
    }

    if cli.json {
        runner::run_pipeline(&cli, &cfg, None, Some(&sample_context))?;
        eprintln!("Done. Output written to '{}'", cli.output.display());
    } else {
        let state = new_shared_state();
        let state_pipeline = state.clone();
        let cli_clone = cli.clone();
        let cfg_clone = cfg.clone();
        let ctx_clone = sample_context.clone();

        let pipeline_handle = std::thread::spawn(move || {
            match runner::run_pipeline(
                &cli_clone,
                &cfg_clone,
                Some(state_pipeline.clone()),
                Some(&ctx_clone),
            ) {
                Ok(_) => {}
                Err(e) => {
                    if let Ok(mut lock) = state_pipeline.lock() {
                        lock.error = Some(format!("{:#}", e));
                        lock.done = true;
                    }
                    log::error!("Pipeline error: {:#}", e);
                }
            }
        });

        tui::run_tui(state)?;

        if let Err(e) = pipeline_handle.join() {
            log::error!("Pipeline thread panicked: {:?}", e);
        }

        eprintln!("Done. Output written to '{}'", cli.output.display());
    }

    Ok(())
}
