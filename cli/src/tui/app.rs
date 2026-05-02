use std::sync::{Arc, Mutex};

/// Current phase of the analysis pipeline.
#[derive(Debug, Clone, PartialEq, Default)]
#[allow(dead_code)]
pub enum Phase {
    #[default]
    Idle,
    Genomics,
    Transcriptomics,
    Epigenomics,
    Integration,
    Done,
    Error(String),
}

impl Phase {
    pub fn label(&self) -> &str {
        match self {
            Phase::Idle => "Idle",
            Phase::Genomics => "Genomics Analysis",
            Phase::Transcriptomics => "Transcriptomics Analysis",
            Phase::Epigenomics => "Epigenomics Analysis",
            Phase::Integration => "Integration Layer",
            Phase::Done => "Complete",
            Phase::Error(_) => "Error",
        }
    }
}

/// Shared mutable state written by the pipeline thread and read by the TUI renderer.
///
/// The `Mutex` hold time is always sub-microsecond (copy of primitives + push to insight vec),
/// so there is no meaningful contention.
#[derive(Debug, Default)]
pub struct AppState {
    pub phase: Phase,

    pub genomics_pct: f64,
    pub genomics_rps: f64,

    pub transcr_pct: f64,
    pub transcr_rps: f64,

    pub epigen_pct: f64,
    pub epigen_rps: f64,

    pub integration_pct: f64,

    pub elapsed_secs: u64,
    pub eta_secs: Option<u64>,

    /// Last 8 insight strings for the live feed, newest first.
    pub insights_live: Vec<String>,

    pub error: Option<String>,
    pub done: bool,
}

impl AppState {
    pub fn push_insight(&mut self, msg: String) {
        self.insights_live.insert(0, msg);
        self.insights_live.truncate(8);
    }
}

/// Thread-safe handle to the shared application state.
pub type SharedState = Arc<Mutex<AppState>>;

/// Create a new shared state handle.
pub fn new_shared_state() -> SharedState {
    Arc::new(Mutex::new(AppState::default()))
}
