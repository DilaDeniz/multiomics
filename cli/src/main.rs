mod args;
mod output;
mod runner;
mod tui;

use anyhow::Result;
use clap::Parser;

use args::Cli;
use tui::new_shared_state;

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"))
        .format_timestamp_millis()
        .init();

    let cli = Cli::parse();

    // Validate required input files exist before doing any work
    for (label, path) in [
        ("--genomics", &cli.genomics),
        ("--transcriptomics", &cli.transcriptomics),
        ("--epigenomics", &cli.epigenomics),
    ] {
        if !path.exists() {
            anyhow::bail!(
                "Input file for {} not found: '{}'",
                label,
                path.display()
            );
        }
    }

    if cli.json {
        // JSON-only mode: run pipeline directly on main thread, no TUI
        runner::run_pipeline(&cli, None)?;
        eprintln!("Done. Output written to '{}'", cli.output.display());
    } else {
        // TUI mode: pipeline on a background thread, TUI renders on main thread
        let state = new_shared_state();
        let state_pipeline = state.clone();
        let cli_clone = cli.clone();

        let pipeline_handle = std::thread::spawn(move || {
            match runner::run_pipeline(&cli_clone, Some(state_pipeline.clone())) {
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

        // Block main thread on TUI render loop
        tui::run_tui(state)?;

        // Wait for pipeline to finish
        if let Err(e) = pipeline_handle.join() {
            log::error!("Pipeline thread panicked: {:?}", e);
        }

        eprintln!("Done. Output written to '{}'", cli.output.display());
    }

    Ok(())
}
