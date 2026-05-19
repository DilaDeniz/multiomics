use serde::{Deserialize, Serialize};

pub const PROTON: f64 = 1.007_276_466_621;
pub const WATER: f64 = 18.010_564_684;

/// Monoisotopic residue mass indexed by ASCII byte (`b'A'` = 65).
/// Returns 0.0 for unknown residues.
#[inline(always)]
pub fn aa_mass(aa: u8) -> f64 {
    // Indexed as (aa - b'A') for A-Z range.
    const T: [f64; 26] = [
        71.037_113_8,  // A
        0.0,           // B (not standard)
        103.009_184_5, // C
        115.026_943_1, // D
        129.042_593_1, // E
        147.068_413_9, // F
        57.021_463_7,  // G
        137.058_911_9, // H
        113.084_064_0, // I
        0.0,           // J
        128.094_963_0, // K
        113.084_064_0, // L
        131.040_484_6, // M
        114.042_927_5, // N
        0.0,           // O
        97.052_763_9,  // P
        128.058_577_5, // Q
        156.101_111_0, // R
        87.032_028_4,  // S
        101.047_678_5, // T
        0.0,           // U
        99.068_413_9,  // V
        186.079_312_9, // W
        0.0,           // X
        163.063_328_5, // Y
        0.0,           // Z
    ];
    if aa.is_ascii_uppercase() {
        T[(aa - b'A') as usize]
    } else {
        0.0
    }
}

/// Raw mass-spectrometry spectrum (MS1 or MS2).
#[derive(Debug, Clone, Default)]
pub struct Spectrum {
    pub scan: u32,
    pub ms_level: u8,
    /// Retention time in seconds.
    pub rt: f32,
    /// Precursor m/z (0 for MS1).
    pub precursor_mz: f32,
    /// Precursor charge state (0 if unknown or MS1).
    pub precursor_z: u8,
    pub mz: Vec<f32>,
    pub intensity: Vec<f32>,
}

impl Spectrum {
    /// Monoisotopic neutral precursor mass.
    #[inline]
    pub fn precursor_mass(&self) -> f64 {
        (self.precursor_mz as f64 - PROTON) * self.precursor_z.max(1) as f64
    }

    /// Normalize intensities to [0, 1] by the base peak.
    pub fn normalize(&mut self) {
        let max = self.intensity.iter().cloned().fold(0.0f32, f32::max);
        if max > 0.0 {
            for v in self.intensity.iter_mut() {
                *v /= max;
            }
        }
    }

    /// Remove peaks below `threshold` fraction of the base peak, then re-sort.
    pub fn filter_noise(&mut self, threshold: f32) {
        let max = self.intensity.iter().cloned().fold(0.0f32, f32::max);
        let cutoff = max * threshold;
        let mut pairs: Vec<(f32, f32)> = self
            .mz
            .iter()
            .copied()
            .zip(self.intensity.iter().copied())
            .filter(|&(_, int)| int >= cutoff)
            .collect();
        pairs.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        self.mz = pairs.iter().map(|p| p.0).collect();
        self.intensity = pairs.iter().map(|p| p.1).collect();
    }
}

/// Tryptic peptide with pre-computed monoisotopic mass.
#[derive(Debug, Clone)]
pub struct Peptide {
    pub sequence: String,
    /// Monoisotopic neutral mass (sum of residues + H₂O).
    pub mass: f64,
    pub protein_idx: u32,
    pub is_decoy: bool,
    pub missed_cleavages: u8,
}

/// Peptide-spectrum match.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Psm {
    pub scan: u32,
    pub rt: f32,
    pub peptide: String,
    pub protein: String,
    pub charge: u8,
    pub hyperscore: f64,
    /// Score gap to second-best candidate (higher = more confident).
    pub delta_score: f64,
    pub mass_error_ppm: f64,
    pub n_matched_b: u32,
    pub n_matched_y: u32,
    pub q_value: f64,
    pub is_decoy: bool,
}

/// Protein-level result after protein inference and FDR control.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProteinGroup {
    pub protein: String,
    pub n_psms: u32,
    pub n_unique_peptides: u32,
    pub top_score: f64,
    pub q_value: f64,
    pub is_decoy: bool,
}

/// Top-level summary emitted to the HTML report and JSON output.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProteomicsSummary {
    pub n_spectra_total: u32,
    pub n_ms2: u32,
    /// PSMs passing 1 % FDR.
    pub n_psms_1pct: u32,
    /// Unique peptide sequences passing 1 % FDR.
    pub n_peptides_1pct: u32,
    /// Protein groups passing 1 % FDR.
    pub n_proteins_1pct: u32,
    pub median_hyperscore: f64,
    pub score_histogram: Vec<u64>,
    pub top_proteins: Vec<ProteinGroup>,
    /// Capped at 5 000 PSMs for the HTML report.
    pub psms: Vec<Psm>,
}
