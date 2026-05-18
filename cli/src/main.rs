mod args;
mod compare;
mod config;
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
            "rnaseq"
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

    // Validate required inputs (they are Option only so --dump-config can skip them)
    for (label, path) in [
        ("--genomics", cli.genomics.as_deref()),
        ("--transcriptomics", cli.transcriptomics.as_deref()),
        ("--epigenomics", cli.epigenomics.as_deref()),
    ] {
        match path {
            None => anyhow::bail!("{} is required", label),
            Some(p) if !p.exists() => {
                anyhow::bail!("Input file for {} not found: '{}'", label, p.display());
            }
            _ => {}
        }
    }

    // Warn at startup about any disabled-feature flags that were passed
    warn_disabled_features(&cli);

    if cli.json {
        runner::run_pipeline(&cli, &cfg, None)?;
        eprintln!("Done. Output written to '{}'", cli.output.display());
    } else {
        let state = new_shared_state();
        let state_pipeline = state.clone();
        let cli_clone = cli.clone();
        let cfg_clone = cfg.clone();

        let pipeline_handle = std::thread::spawn(move || {
            match runner::run_pipeline(&cli_clone, &cfg_clone, Some(state_pipeline.clone())) {
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
