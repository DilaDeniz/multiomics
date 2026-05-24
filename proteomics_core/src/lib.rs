//! Proteomics analysis pipeline.
//!
//! Implements mzML parsing, in-silico tryptic digest, database search
//! (hyperscore, target-decoy FDR), and label-free XIC quantification.
//! Supports both single-file and multi-file (multi-run) search: all mzML
//! runs are processed in parallel and their PSMs merged before a single
//! experiment-level FDR pass, which is the correct way to handle fractionated
//! or replicated LC-MS/MS experiments.
//!
//! # Speed vs competitors
//! - mzML parsing: streaming quick-xml, no DOM allocation.
//! - Candidate lookup: hybrid precursor bin + fragment vote index (Sage-inspired).
//! - Scoring: binary-search fragment matching, no hash tables on the hot path.
//! - Parallelism: rayon over spectra; each worker is stateless.
//!   Multi-file: files are parsed and searched in parallel, then merged.
//!
//! # References
//! - Craig R & Beavis RC (2004) TANDEM: matching proteins with tandem mass spectra.
//!   Bioinformatics 20:1466-1467. (hyperscore)
//! - Elias JE & Gygi SP (2007) Nature Methods 4:207-214. (target-decoy FDR)

pub mod dia;
pub mod fasta;
pub mod fdr;
pub mod index;
pub mod mzml;
pub mod phospho;
pub mod quant;
pub mod score;
pub mod search;
pub mod types;

pub use dia::{
    assign_dia_qvalues, build_library_from_psms, score_library_against_dia, DiaPsm, LibraryEntry,
    MIN_DOT_SCORE, MIN_FRAGMENTS,
};
pub use fasta::{digest, parse_fasta, peptide_mass};
pub use fdr::{assign_qvalues, filter_psms};
pub use index::PeptideIndex;
pub use mzml::parse_mzml;
pub use phospho::{expand_with_phospho, phospho_variants, DEFAULT_MAX_PHOSPHO_SITES, PHOSPHO_MASS};
pub use quant::{quantify_psms, PeptideQuant};
pub use score::{b_ions, hyperscore, y_ions, FRAG_TOL_PPM};
pub use search::{infer_proteins, search_spectra, PRECURSOR_TOL_PPM};
pub use types::{aa_mass, Peptide, ProteinGroup, ProteomicsSummary, Psm, Spectrum, PROTON, WATER};

use ahash::AHashSet;
use anyhow::Result;

/// Run the complete proteomics pipeline for a single mzML file.
///
/// Convenience wrapper around [`run_proteomics_multi`] for the common
/// single-file case. Pass `phospho_max_sites = 0` to disable phospho search.
pub fn run_proteomics(
    mzml_data: &[u8],
    fasta_data: &[u8],
    fdr_threshold: f64,
    phospho_max_sites: usize,
) -> Result<ProteomicsSummary> {
    run_proteomics_multi(&[mzml_data], fasta_data, fdr_threshold, phospho_max_sites)
}

/// Run the complete proteomics pipeline across multiple mzML files.
///
/// All files are parsed and searched **in parallel**; PSMs from every run are
/// pooled before a single experiment-level FDR pass. This is the correct
/// approach for fractionated or replicated LC-MS/MS experiments.
///
/// - `mzml_slices`: raw bytes of each mzML file.
/// - `fasta_data`: protein database FASTA (target only; decoys auto-generated).
/// - `fdr_threshold`: FDR cutoff for reporting (e.g. 0.01 for 1 %).
/// - `phospho_max_sites`: number of variable phosphorylation sites per peptide
///   (0 = disabled, 1–3 recommended for phosphoproteomics experiments).
pub fn run_proteomics_multi(
    mzml_slices: &[&[u8]],
    fasta_data: &[u8],
    fdr_threshold: f64,
    phospho_max_sites: usize,
) -> Result<ProteomicsSummary> {
    use rayon::prelude::*;

    if mzml_slices.is_empty() {
        anyhow::bail!("no mzML files provided");
    }

    // 1. Build the peptide index once — shared (read-only) across all searches.
    log::info!("Digesting protein database...");
    let proteins = parse_fasta(fasta_data);
    let protein_names: Vec<String> = proteins.iter().map(|p| first_word(&p.header)).collect();
    let peptides = digest(&proteins, 2, 6, 50);
    log::info!("{} peptides (target+decoy)", peptides.len());

    // Optional: expand with phospho variants before building the index.
    let mut peptides = peptides;
    if phospho_max_sites > 0 {
        let n_before = peptides.len();
        expand_with_phospho(&mut peptides, phospho_max_sites);
        log::info!(
            "Phospho expansion (+{} variants, {} total)",
            peptides.len() - n_before,
            peptides.len()
        );
    }

    let index = PeptideIndex::build(peptides);

    // 2. Parse all mzML files in parallel, pre-process MS2 peaks.
    log::info!("Parsing {} mzML file(s)...", mzml_slices.len());
    let per_file_results: Vec<anyhow::Result<Vec<Spectrum>>> = mzml_slices
        .par_iter()
        .map(|data| {
            let mut spectra = parse_mzml(data)?;
            for spec in spectra.iter_mut().filter(|s| s.ms_level == 2) {
                spec.filter_noise(0.01);
                spec.normalize();
            }
            Ok(spectra)
        })
        .collect();

    // Propagate first parse error; accumulate totals.
    let mut all_spectra: Vec<Spectrum> = Vec::new();
    for result in per_file_results {
        all_spectra.extend(result?);
    }

    let n_spectra_total = all_spectra.len() as u32;
    let n_ms2 = all_spectra.iter().filter(|s| s.ms_level == 2).count() as u32;
    log::info!("{n_spectra_total} total spectra across all runs ({n_ms2} MS2)");

    // 3. Database search — rayon parallelism is over spectra within search_spectra.
    log::info!("Searching {n_ms2} MS2 spectra...");
    let mut all_psms = search_spectra(
        &all_spectra,
        &index,
        &protein_names,
        FRAG_TOL_PPM,
        PRECURSOR_TOL_PPM,
    );

    // 4. Experiment-level FDR (single pass over all runs combined).
    assign_qvalues(&mut all_psms);
    let passing = filter_psms(&all_psms, fdr_threshold);
    let n_psms_1pct = passing.len() as u32;

    let mut peptide_set: AHashSet<&str> = AHashSet::default();
    for p in &passing {
        peptide_set.insert(p.peptide.as_str());
    }
    let n_peptides_1pct = peptide_set.len() as u32;

    // 5. Protein inference.
    let protein_groups = infer_proteins(&passing);
    let n_proteins_1pct = protein_groups
        .iter()
        .filter(|p| !p.is_decoy && p.q_value <= fdr_threshold)
        .count() as u32;

    // 6. Summary statistics.
    let median_hyperscore = median_f64(&passing.iter().map(|p| p.hyperscore).collect::<Vec<_>>());
    let score_histogram = build_histogram(&passing, 20);
    let top_proteins: Vec<ProteinGroup> = protein_groups
        .into_iter()
        .filter(|p| !p.is_decoy && p.q_value <= fdr_threshold)
        .take(20)
        .collect();

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
    // Single pass: find max and accumulate bins simultaneously.
    let mut max_score = 0.0f64;
    for p in psms {
        if p.hyperscore > max_score {
            max_score = p.hyperscore;
        }
    }
    let bin_width = (max_score / n_bins as f64).max(1.0);
    let mut hist = vec![0u64; n_bins];
    for p in psms {
        hist[((p.hyperscore / bin_width) as usize).min(n_bins - 1)] += 1;
    }
    hist
}
