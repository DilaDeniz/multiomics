pub mod accum;
pub mod cancer;
pub mod cnv;
pub mod cnv_coverage;
pub mod fastq;
pub mod types;
pub mod vcf;

pub use accum::GenomicsAccum;
pub use cancer::{
    compute_hrd_score, compute_hrd_score_with_reference, detect_kataegis, detect_loh,
    estimate_tumor_purity, HrdScore, KataegisLocus, LohChromosome, TumorPurityResult,
};
pub use cnv::{parse_cnv_vcf, summarize_cnv, CnvRecord, CnvSummary};
pub use cnv_coverage::{
    compute_depth_bins, segment_bins, summarize_coverage_cnv, CnvSegment, CoverageCnvSummary,
    DepthBin,
};
pub use fastq::{parse_fastq, FastqSummary};
pub use types::{ChromDensity, GenomicsSummary, TiTvClass, VariantRecord};
pub use vcf::parse_vcf;

pub mod assembly;
pub mod genotyper;
pub mod pairhmm;
pub mod pileup;
pub use assembly::{assemble_haplotypes, find_active_regions, ActiveRead, ActiveRegion};
pub use genotyper::{
    call_variants, call_variants_assembled, call_variants_from_bam, Genotype, GenotypeCall,
};
pub use pairhmm::{log_sum_exp, pair_hmm_log_prob};
pub use pileup::{build_pileup, PileupBase, PileupColumn};

pub mod aligner;
pub use aligner::{Alignment, ReferenceIndex};

pub mod somatic;
pub use somatic::{
    call_somatic_from_bams, call_somatic_variants, summarize_somatic, SomaticCall, SomaticSummary,
};

pub mod splice_aligner;
pub use splice_aligner::{SpliceAlignment, SpliceIndex, SpliceJunction};

use anyhow::Result;
use crossbeam_channel::Sender;
use std::path::Path;

use biomics_core::{parallel_fold, ProgressEvent};

/// Run the full genomics analysis pipeline on a VCF file.
pub fn analyze_vcf(
    path: &Path,
    progress_tx: Option<&Sender<ProgressEvent>>,
) -> Result<GenomicsSummary> {
    let records = parse_vcf(path)?;
    parallel_fold::<GenomicsAccum>(&records, "genomics", progress_tx)
}
