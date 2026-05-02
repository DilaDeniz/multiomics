use serde::{Deserialize, Serialize};

/// Labels for the four analysis modalities, used in progress events and log messages.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ModalityLabel {
    Genomics,
    Transcriptomics,
    Epigenomics,
    Integration,
}

impl ModalityLabel {
    pub fn as_str(self) -> &'static str {
        match self {
            ModalityLabel::Genomics => "genomics",
            ModalityLabel::Transcriptomics => "transcriptomics",
            ModalityLabel::Epigenomics => "epigenomics",
            ModalityLabel::Integration => "integration",
        }
    }
}

impl std::fmt::Display for ModalityLabel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}
