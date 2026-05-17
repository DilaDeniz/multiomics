use crossterm::event::{self, Event, KeyCode, KeyEvent, KeyModifiers};
use std::time::Duration;

/// Events produced by the crossterm event loop.
#[derive(Debug, Clone)]
pub enum AppEvent {
    /// User pressed 'q' or Ctrl-C — quit requested.
    Quit,
    /// User pressed 'p' — toggle pause (future use).
    Pause,
    /// Resize or tick — redraw the TUI.
    Tick,
}

/// Poll for the next crossterm event, waiting at most `timeout`.
///
/// Returns `None` when the poll times out with no event.
pub fn poll_event(timeout: Duration) -> anyhow::Result<Option<AppEvent>> {
    if !event::poll(timeout)? {
        return Ok(None);
    }

    match event::read()? {
        Event::Key(KeyEvent {
            code, modifiers, ..
        }) => match (code, modifiers) {
            (KeyCode::Char('q'), _) | (KeyCode::Char('c'), KeyModifiers::CONTROL) => {
                Ok(Some(AppEvent::Quit))
            }
            (KeyCode::Char('p'), _) => Ok(Some(AppEvent::Pause)),
            _ => Ok(Some(AppEvent::Tick)),
        },
        Event::Resize(_, _) => Ok(Some(AppEvent::Tick)),
        _ => Ok(Some(AppEvent::Tick)),
    }
}
