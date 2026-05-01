pub mod app;
pub mod events;
pub mod widgets;

pub use app::{new_shared_state, AppState, SharedState};
pub use events::{poll_event, AppEvent};
pub use widgets::render;

use anyhow::Result;
use crossterm::{
    execute,
    terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
};
use ratatui::{backend::CrosstermBackend, Terminal};
use std::{
    io,
    time::{Duration, Instant},
};

/// Run the TUI render loop, blocking until the pipeline sets `state.done` or the user quits.
///
/// The loop renders at ~10 fps (100 ms tick). State is read under the Mutex at each tick.
pub fn run_tui(state: SharedState) -> Result<()> {
    enable_raw_mode()?;
    let mut stdout = io::stdout();
    execute!(stdout, EnterAlternateScreen)?;

    let backend = CrosstermBackend::new(stdout);
    let mut terminal = Terminal::new(backend)?;

    let tick = Duration::from_millis(100);
    let start = Instant::now();

    loop {
        // Read state snapshot under lock
        let snapshot = {
            let mut s = state.lock().expect("state mutex poisoned");
            s.elapsed_secs = start.elapsed().as_secs();
            AppState {
                phase: s.phase.clone(),
                genomics_pct: s.genomics_pct,
                genomics_rps: s.genomics_rps,
                transcr_pct: s.transcr_pct,
                transcr_rps: s.transcr_rps,
                epigen_pct: s.epigen_pct,
                epigen_rps: s.epigen_rps,
                integration_pct: s.integration_pct,
                elapsed_secs: s.elapsed_secs,
                eta_secs: s.eta_secs,
                insights_live: s.insights_live.clone(),
                error: s.error.clone(),
                done: s.done,
            }
        };

        terminal.draw(|frame| render(frame, &snapshot))?;

        match poll_event(tick)? {
            Some(AppEvent::Quit) => break,
            _ => {}
        }

        if snapshot.done || snapshot.error.is_some() {
            // Give user a moment to read the final state
            std::thread::sleep(Duration::from_millis(500));
            break;
        }
    }

    disable_raw_mode()?;
    execute!(terminal.backend_mut(), LeaveAlternateScreen)?;
    terminal.show_cursor()?;

    Ok(())
}
