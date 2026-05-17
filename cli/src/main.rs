mod args;
mod compare;
mod config;
mod output;
mod runner;
mod tui;

use anyhow::Result;
use clap::Parser;

use args::Cli;
use config::{dump_default_config, load_config, BioomicsConfig};
use tui::new_shared_state;

// Replace the default system allocator with mimalloc.
// Benchmarks on allocation-heavy workloads show 20-40% wall-clock improvement.
#[global_allocator]
static ALLOC: mimalloc::MiMalloc = mimalloc::MiMalloc;

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

    // Load configuration (config file merged with defaults; CLI flags take priority)
    let cfg = if let Some(ref path) = cli.config {
        load_config(path)?
    } else {
        BioomicsConfig::default()
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
