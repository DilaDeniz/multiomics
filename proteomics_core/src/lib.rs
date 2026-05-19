//! Proteomics analysis pipeline.
//!
//! Implements mzML parsing, in-silico tryptic digest, database search
//! (hyperscore, target-decoy FDR), and label-free XIC quantification.
//!
//! # Speed vs competitors
//! - mzML parsing: streaming quick-xml, no DOM allocation.
//! - Candidate lookup: 1-Da mass-bin index — O(1) per spectrum.
//! - Scoring: binary-search fragment matching, no hash tables on the hot path.
//! - Parallelism: rayon over spectra; each worker is stateless.
//!
//! # References
//! - Craig R & Beavis RC (2004) TANDEM: matching proteins with tandem mass spectra.
//!   Bioinformatics 20:1466-1467. (hyperscore)
//! - Elias JE & Gygi SP (2007) Nature Methods 4:207-214. (target-decoy FDR)

pub mod fasta;
pub mod fdr;
pub mod index;
pub mod mzml;
pub mod quant;
pub mod score;
pub mod search;
pub mod types;

pub use fasta::{digest, parse_fasta, peptide_mass};
pub use fdr::{assign_qvalues, filter_psms};
pub use index::PeptideIndex;
pub use mzml::parse_mzml;
pub use quant::{quantify_psms, PeptideQuant};
pub use score::{b_ions, hyperscore, y_ions, FRAG_TOL_PPM};
pub use search::{infer_proteins, search_spectra, PRECURSOR_TOL_PPM};
pub use types::{aa_mass, Peptide, ProteinGroup, ProteomicsSummary, Psm, Spectrum, PROTON, WATER};

use anyhow::Result;

/// Run the complete proteomics pipeline.
///
/// - `mzml_data`: raw mzML file bytes.
/// - `fasta_data`: protein database FASTA bytes (target only; decoys generated internally).
/// - `fdr_threshold`: FDR cutoff for reporting (e.g. 0.01 for 1 %).
pub fn run_proteomics(
    mzml_data: &[u8],
    fasta_data: &[u8],
    fdr_threshold: f64,
) -> Result<ProteomicsSummary> {
    // 1. Parse spectra.
    log::info!("Parsing mzML...");
    let mut spectra = parse_mzml(mzml_data)?;
    let n_spectra_total = spectra.len() as u32;
    let n_ms2 = spectra.iter().filter(|s| s.ms_level == 2).count() as u32;
    log::info!("{n_spectra_total} spectra ({n_ms2} MS2)");

    // 2. Pre-process MS2 spectra: filter noise, normalize.
    for spec in spectra.iter_mut().filter(|s| s.ms_level == 2) {
        spec.filter_noise(0.01);
        spec.normalize();
    }

    // 3. Build peptide index.
    log::info!("Digesting protein database...");
    let proteins = parse_fasta(fasta_data);
    let protein_names: Vec<String> = proteins.iter().map(|p| first_word(&p.header)).collect();
    let peptides = digest(&proteins, 2, 6, 50);
    log::info!("{} peptides (target+decoy)", peptides.len());
    let index = PeptideIndex::build(peptides);

    // 4. Database search.
    log::info!("Searching {n_ms2} MS2 spectra...");
    let all_psms = search_spectra(
        &spectra,
        &index,
        &protein_names,
        FRAG_TOL_PPM,
        PRECURSOR_TOL_PPM,
    );

    // 5. Filter to FDR threshold.
    let passing = filter_psms(&all_psms, fdr_threshold);
    let n_psms_1pct = passing.len() as u32;
    let mut peptide_set: std::collections::HashSet<&str> = std::collections::HashSet::new();
    for p in &passing {
        peptide_set.insert(p.peptide.as_str());
    }
    let n_peptides_1pct = peptide_set.len() as u32;

    // 6. Protein inference.
    let protein_groups = infer_proteins(&passing);
    let n_proteins_1pct = protein_groups
        .iter()
        .filter(|p| !p.is_decoy && p.q_value <= fdr_threshold)
        .count() as u32;

    // 7. Compute summary statistics.
    let median_hyperscore = median_f64(&passing.iter().map(|p| p.hyperscore).collect::<Vec<_>>());
    let score_histogram = build_histogram(&passing, 20);
    let top_proteins: Vec<ProteinGroup> = protein_groups
        .into_iter()
        .filter(|p| !p.is_decoy && p.q_value <= fdr_threshold)
        .take(20)
        .collect();

    // Cap PSMs for the HTML report.
    let psms_report: Vec<Psm> = passing.into_iter().take(5000).collect();

    log::info!(
        "Proteomics: {n_psms_1pct} PSMs, {n_peptides_1pct} peptides, {n_proteins_1pct} proteins at {:.0}% FDR",
        fdr_threshold * 100.0
    );

    Ok(ProteomicsSummary {
        n_spectra_total,
        n_ms2,
        n_psms_1pct,
        n_peptides_1pct,
        n_proteins_1pct,
        median_hyperscore,
        score_histogram,
        top_proteins,
        psms: psms_report,
    })
}

fn first_word(s: &str) -> String {
    s.split_whitespace().next().unwrap_or(s).to_string()
}

fn median_f64(v: &[f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    let mut s = v.to_vec();
    s.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    let n = s.len();
    if n % 2 == 1 {
        s[n / 2]
    } else {
        (s[n / 2 - 1] + s[n / 2]) / 2.0
    }
}

fn build_histogram(psms: &[Psm], n_bins: usize) -> Vec<u64> {
    if psms.is_empty() {
        return vec![0; n_bins];
    }
    let max_score = psms.iter().map(|p| p.hyperscore).fold(0.0f64, f64::max);
    let bin_width = (max_score / n_bins as f64).max(1.0);
    let mut hist = vec![0u64; n_bins];
    for p in psms {
        let bin = ((p.hyperscore / bin_width) as usize).min(n_bins - 1);
        hist[bin] += 1;
    }
    hist
}
