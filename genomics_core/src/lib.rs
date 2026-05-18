pub mod accum;
pub mod cnv;
pub mod fastq;
pub mod types;
pub mod vcf;

pub use accum::GenomicsAccum;
pub use cnv::{parse_cnv_vcf, summarize_cnv, CnvRecord, CnvSummary};
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
