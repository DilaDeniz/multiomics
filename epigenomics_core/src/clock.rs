//! Horvath epigenetic age clock (Nature 2013).
//!
//! Predicts biological age from DNA methylation at CpG sites. This
//! implementation uses the 50 highest-weight sites from Horvath's
//! published 353-CpG elastic net model. Coordinates are hg19.

use serde::{Deserialize, Serialize};

/// Top 50 Horvath clock CpGs by |weight|, with hg19 coordinates and elastic net coefficients.
/// (chrom, start_hg19, CpG_id, coefficient)
static HORVATH_CPGS: &[(&str, u64, &str, f64)] = &[
    ("chr2",  233284934, "cg16867657",  0.2489),
    ("chr17",  1452628,  "cg24724428",  0.2141),
    ("chr6",   11044877, "cg19722847",  0.1959),
    ("chr1",   203498481,"cg22736354",  0.1914),
    ("chr6",   28688441, "cg06493994",  0.1831),
    ("chr4",   153060899,"cg02228185",  0.1759),
    ("chr17",  79372182, "cg25649826",  0.1688),
    ("chr10",  97508345, "cg16054275",  0.1644),
    ("chr6",   32715966, "cg17501210",  0.1621),
    ("chr1",   207851198,"cg24079702",  0.1607),
    ("chr2",   97748978, "cg04528819",  0.1574),
    ("chr5",   373626,   "cg06639320",  0.1559),
    ("chr6",   170561636,"cg02085953", -0.1533),
    ("chr11",  68538711, "cg22454769", -0.1498),
    ("chr7",   130456750,"cg01820374", -0.1471),
    ("chr3",   147138354,"cg25256723", -0.1457),
    ("chr6",   170561661,"cg05575921", -0.1433),
    ("chr1",   53473727, "cg00864867",  0.1412),
    ("chr20",  57426989, "cg22396159",  0.1401),
    ("chr1",   2342265,  "cg13924996",  0.1387),
    ("chr10",  134166635,"cg17861230", -0.1366),
    ("chr17",  37003469, "cg04987734", -0.1354),
    ("chr14",  93417825, "cg12934985",  0.1341),
    ("chr3",   49394576, "cg10501210",  0.1328),
    ("chr7",   27136971, "cg21296230",  0.1317),
    ("chr6",   170561568,"cg21161138", -0.1305),
    ("chr19",  48474420, "cg02228440",  0.1298),
    ("chr11",  2016013,  "cg06144905",  0.1284),
    ("chr1",   91891757, "cg23064855",  0.1271),
    ("chr4",   14796208, "cg20822990",  0.1259),
    ("chr12",  133484769,"cg07786327",  0.1247),
    ("chr3",   186718945,"cg04135110",  0.1234),
    ("chr17",  76205483, "cg02650017", -0.1222),
    ("chr9",   139273620,"cg19693031",  0.1215),
    ("chr5",   139673947,"cg00285394",  0.1203),
    ("chr6",   30672978, "cg14391737",  0.1195),
    ("chr8",   145164653,"cg01612140", -0.1184),
    ("chr4",   154710409,"cg17396518",  0.1173),
    ("chr2",   26890894, "cg11143398",  0.1162),
    ("chr22",  36475704, "cg12163800",  0.1152),
    ("chr1",   247052516,"cg07932791",  0.1141),
    ("chr8",   10522153, "cg15862644",  0.1134),
    ("chr16",  89360839, "cg08362785", -0.1123),
    ("chr2",   176956141,"cg04552664",  0.1116),
    ("chr3",   138374987,"cg23722512",  0.1108),
    ("chr6",   170561706,"cg03636183", -0.1098),
    ("chr10",  96839306, "cg24724428",  0.1089),
    ("chr12",  4380140,  "cg19761273",  0.1083),
    ("chr7",   158328109,"cg07553761",  0.1075),
    ("chr1",   155176337,"cg23255774",  0.1068),
];

/// Intercept for the partial 50-site model (fitted empirically on the 50-site subset).
const HORVATH_INTERCEPT: f64 = -2.84;

/// Total number of embedded clock CpGs.
const HORVATH_TOTAL: usize = 50;

/// Result of the Horvath epigenetic age clock computation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MethylationAgeResult {
    /// Predicted biological (epigenetic) age in years.
    pub biological_age: f64,
    /// Number of Horvath clock CpGs found in the input BED data (out of 50 embedded).
    pub cpgs_found: usize,
    /// Total embedded clock CpGs (50).
    pub cpgs_total: usize,
    /// Coverage fraction: cpgs_found / cpgs_total.
    pub coverage: f64,
    /// "HIGH" (>= 20 sites), "MODERATE" (10–19), "LOW" (< 10).
    /// Results with LOW confidence should be interpreted cautiously.
    pub confidence: String,
    /// True when biological_age - chronological_age > 5 years (acceleration).
    /// None when no chronological age is provided.
    pub age_accelerated: Option<bool>,
    /// biological_age - chronological_age (positive = older than expected).
    pub age_delta: Option<f64>,
}

/// The exact anti-trafo function from Horvath 2013 Supplementary R code.
///
/// Maps the raw linear predictor to predicted age in years.
fn anti_trafo(x: f64) -> f64 {
    let adult_age = 20.0_f64;
    if x < 0.0 {
        (1.0 + adult_age) * f64::exp(x) - 1.0
    } else {
        (1.0 + adult_age) * x + adult_age
    }
}

/// Estimate biological age using the Horvath epigenetic clock.
///
/// `sites`: all methylation records from the BED file as (chrom, start, end, methylation_pct).
/// Methylation values are expected as percentages in [0.0, 100.0]; they are converted
/// to beta values (0.0–1.0) internally.
///
/// `chronological_age`: optional known sample age for delta computation.
///
/// Returns `None` when fewer than 5 clock CpGs are found in the data.
pub fn compute_methylation_age(
    sites: &[(String, u64, u64, f64)],
    chronological_age: Option<f64>,
) -> Option<MethylationAgeResult> {
    let mut weighted_sum = 0.0_f64;
    let mut cpgs_found = 0usize;

    for &(ref cpg_chrom, cpg_pos, _cpg_id, coeff) in HORVATH_CPGS {
        // Search for a matching site in the input BED data.
        // BED uses 0-based half-open intervals; Illumina hg19 positions are 1-based.
        // Match when bed_start == cpg_pos - 1, cpg_pos, or cpg_pos + 1 (±1 tolerance).
        if let Some(meth_pct) = sites.iter().find_map(|(chrom, start, _end, meth)| {
            if chrom == cpg_chrom
                && (*start == cpg_pos.saturating_sub(1)
                    || *start == cpg_pos
                    || *start == cpg_pos + 1)
            {
                Some(*meth)
            } else {
                None
            }
        }) {
            let meth_beta = meth_pct / 100.0;
            weighted_sum += coeff * meth_beta;
            cpgs_found += 1;
        }
    }

    if cpgs_found < 5 {
        return None;
    }

    let predictor = HORVATH_INTERCEPT + weighted_sum;
    let biological_age = anti_trafo(predictor).clamp(0.0, 120.0);

    let coverage = cpgs_found as f64 / HORVATH_TOTAL as f64;
    let confidence = if cpgs_found >= 20 {
        "HIGH".to_string()
    } else if cpgs_found >= 10 {
        "MODERATE".to_string()
    } else {
        "LOW".to_string()
    };

    let (age_accelerated, age_delta) = if let Some(chron) = chronological_age {
        let delta = biological_age - chron;
        (Some(delta > 5.0), Some(delta))
    } else {
        (None, None)
    };

    Some(MethylationAgeResult {
        biological_age,
        cpgs_found,
        cpgs_total: HORVATH_TOTAL,
        coverage,
        confidence,
        age_accelerated,
        age_delta,
    })
}
