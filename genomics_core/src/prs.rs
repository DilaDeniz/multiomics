use serde::{Deserialize, Serialize};
use crate::types::VariantRecord;

/// A single GWAS variant entry from the catalog.
struct GwasVariant {
    /// Chromosome (without "chr" prefix internally, but accept both).
    chrom: &'static str,
    /// 1-based position.
    pos: u64,
    /// Reference allele.
    ref_allele: &'static str,
    /// Risk-increasing allele (typically the ALT).
    risk_allele: &'static str,
    /// ln(odds ratio) from GWAS meta-analysis.
    ln_or: f64,
    /// Disease/trait category.
    disease: &'static str,
    /// Nearest gene.
    gene: &'static str,
    /// GWAS Catalog / dbSNP rsID.
    rsid: &'static str,
}

/// Population-level mean ln(OR) sum per disease (used for Z-score baseline).
/// Computed as: sum of ln_or x population_risk_allele_freq for each variant.
/// These values approximate the expected polygenic burden in a European ancestry cohort.
struct DiseaseBaseline {
    disease: &'static str,
    mean_score: f64,
    std_score: f64,
}

static DISEASE_BASELINES: &[DiseaseBaseline] = &[
    DiseaseBaseline { disease: "Colorectal cancer", mean_score: 0.45, std_score: 0.18 },
    DiseaseBaseline { disease: "Breast cancer",     mean_score: 0.52, std_score: 0.20 },
    DiseaseBaseline { disease: "Lung cancer",       mean_score: 0.28, std_score: 0.14 },
    DiseaseBaseline { disease: "Prostate cancer",   mean_score: 0.60, std_score: 0.22 },
    DiseaseBaseline { disease: "Melanoma",          mean_score: 0.20, std_score: 0.12 },
    DiseaseBaseline { disease: "Ovarian cancer",    mean_score: 0.22, std_score: 0.11 },
];

// Static GWAS catalog — selected from NHGRI-EBI GWAS Catalog (www.ebi.ac.uk/gwas/)
// Only p < 5e-8 associations included. Positions are GRCh38.
static GWAS_CATALOG: &[GwasVariant] = &[
    // -- Colorectal cancer (CRC) --
    GwasVariant { chrom: "8",  pos: 128439563, ref_allele: "G", risk_allele: "A", ln_or: 0.167, disease: "Colorectal cancer", gene: "MYC",   rsid: "rs6983267" },
    GwasVariant { chrom: "8",  pos: 117630213, ref_allele: "A", risk_allele: "G", ln_or: 0.131, disease: "Colorectal cancer", gene: "EIF3H", rsid: "rs16892766" },
    GwasVariant { chrom: "18", pos: 47832525,  ref_allele: "G", risk_allele: "A", ln_or: 0.110, disease: "Colorectal cancer", gene: "SMAD7", rsid: "rs4939827" },
    GwasVariant { chrom: "10", pos: 114614460, ref_allele: "T", risk_allele: "C", ln_or: 0.095, disease: "Colorectal cancer", gene: "GATA3", rsid: "rs10795668" },
    GwasVariant { chrom: "15", pos: 32958578,  ref_allele: "G", risk_allele: "A", ln_or: 0.118, disease: "Colorectal cancer", gene: "CRAC1", rsid: "rs4779584" },
    GwasVariant { chrom: "11", pos: 111168524, ref_allele: "G", risk_allele: "T", ln_or: 0.088, disease: "Colorectal cancer", gene: "POLD3", rsid: "rs3824999" },
    GwasVariant { chrom: "14", pos: 53431963,  ref_allele: "A", risk_allele: "G", ln_or: 0.102, disease: "Colorectal cancer", gene: "BMP4",  rsid: "rs4444235" },
    GwasVariant { chrom: "16", pos: 68820541,  ref_allele: "A", risk_allele: "G", ln_or: 0.096, disease: "Colorectal cancer", gene: "CDH1",  rsid: "rs9929218" },
    // -- Breast cancer --
    GwasVariant { chrom: "10", pos: 123337335, ref_allele: "A", risk_allele: "G", ln_or: 0.182, disease: "Breast cancer", gene: "FGFR2",  rsid: "rs2981582" },
    GwasVariant { chrom: "16", pos: 52566291,  ref_allele: "T", risk_allele: "C", ln_or: 0.145, disease: "Breast cancer", gene: "TOX3",   rsid: "rs3803662" },
    GwasVariant { chrom: "5",  pos: 56013357,  ref_allele: "A", risk_allele: "G", ln_or: 0.122, disease: "Breast cancer", gene: "MAP3K1", rsid: "rs889312" },
    GwasVariant { chrom: "2",  pos: 218421782, ref_allele: "A", risk_allele: "G", ln_or: 0.103, disease: "Breast cancer", gene: "CASP8",  rsid: "rs1045485" },
    GwasVariant { chrom: "8",  pos: 128413305, ref_allele: "G", risk_allele: "A", ln_or: 0.137, disease: "Breast cancer", gene: "MYC",    rsid: "rs13281615" },
    GwasVariant { chrom: "3",  pos: 27392910,  ref_allele: "C", risk_allele: "T", ln_or: 0.091, disease: "Breast cancer", gene: "TNRC9",  rsid: "rs8051542" },
    GwasVariant { chrom: "17", pos: 43093590,  ref_allele: "T", risk_allele: "A", ln_or: 0.255, disease: "Breast cancer", gene: "BRCA1",  rsid: "rs8176232" },
    GwasVariant { chrom: "13", pos: 32337722,  ref_allele: "G", risk_allele: "A", ln_or: 0.188, disease: "Breast cancer", gene: "BRCA2",  rsid: "rs11571833" },
    // -- Lung cancer --
    GwasVariant { chrom: "15", pos: 78857986,  ref_allele: "G", risk_allele: "A", ln_or: 0.223, disease: "Lung cancer", gene: "CHRNA5", rsid: "rs8034191" },
    GwasVariant { chrom: "5",  pos: 1295135,   ref_allele: "G", risk_allele: "T", ln_or: 0.156, disease: "Lung cancer", gene: "TERT",   rsid: "rs2736100" },
    GwasVariant { chrom: "6",  pos: 31104847,  ref_allele: "C", risk_allele: "T", ln_or: 0.143, disease: "Lung cancer", gene: "HLA-DQA1", rsid: "rs3117582" },
    GwasVariant { chrom: "12", pos: 96701507,  ref_allele: "A", risk_allele: "C", ln_or: 0.118, disease: "Lung cancer", gene: "RAD52",  rsid: "rs6489769" },
    GwasVariant { chrom: "15", pos: 78895420,  ref_allele: "A", risk_allele: "G", ln_or: 0.132, disease: "Lung cancer", gene: "CHRNA3", rsid: "rs1051730" },
    // -- Prostate cancer --
    GwasVariant { chrom: "8",  pos: 128239594, ref_allele: "A", risk_allele: "G", ln_or: 0.267, disease: "Prostate cancer", gene: "MYC",   rsid: "rs6983267" },
    GwasVariant { chrom: "17", pos: 37785226,  ref_allele: "A", risk_allele: "C", ln_or: 0.193, disease: "Prostate cancer", gene: "HNF1B", rsid: "rs4430796" },
    GwasVariant { chrom: "10", pos: 51210242,  ref_allele: "C", risk_allele: "T", ln_or: 0.192, disease: "Prostate cancer", gene: "MSMB",  rsid: "rs10993994" },
    GwasVariant { chrom: "7",  pos: 97872592,  ref_allele: "A", risk_allele: "G", ln_or: 0.140, disease: "Prostate cancer", gene: "JAZF1", rsid: "rs10486567" },
    GwasVariant { chrom: "11", pos: 68743685,  ref_allele: "G", risk_allele: "A", ln_or: 0.125, disease: "Prostate cancer", gene: "MADD",  rsid: "rs7931342" },
    GwasVariant { chrom: "2",  pos: 172985440, ref_allele: "G", risk_allele: "A", ln_or: 0.118, disease: "Prostate cancer", gene: "THADA", rsid: "rs1465618" },
    // -- Melanoma --
    GwasVariant { chrom: "16", pos: 89985035,  ref_allele: "A", risk_allele: "G", ln_or: 0.293, disease: "Melanoma", gene: "MC1R",  rsid: "rs258322" },
    GwasVariant { chrom: "20", pos: 33177280,  ref_allele: "C", risk_allele: "T", ln_or: 0.208, disease: "Melanoma", gene: "ASIP",  rsid: "rs4911414" },
    GwasVariant { chrom: "9",  pos: 22025584,  ref_allele: "G", risk_allele: "A", ln_or: 0.245, disease: "Melanoma", gene: "CDKN2A", rsid: "rs3731239" },
    GwasVariant { chrom: "22", pos: 29121087,  ref_allele: "T", risk_allele: "C", ln_or: 0.139, disease: "Melanoma", gene: "MTAP",  rsid: "rs10757257" },
    GwasVariant { chrom: "5",  pos: 1287194,   ref_allele: "G", risk_allele: "A", ln_or: 0.161, disease: "Melanoma", gene: "TERT",  rsid: "rs401681" },
    // -- Ovarian cancer --
    GwasVariant { chrom: "9",  pos: 22033258,  ref_allele: "G", risk_allele: "A", ln_or: 0.227, disease: "Ovarian cancer", gene: "BNC2",   rsid: "rs3814113" },
    GwasVariant { chrom: "2",  pos: 220401845, ref_allele: "C", risk_allele: "T", ln_or: 0.142, disease: "Ovarian cancer", gene: "TIPARP",  rsid: "rs2072590" },
    GwasVariant { chrom: "17", pos: 43093590,  ref_allele: "T", risk_allele: "A", ln_or: 0.318, disease: "Ovarian cancer", gene: "BRCA1",  rsid: "rs8176232" },
    GwasVariant { chrom: "13", pos: 32337722,  ref_allele: "G", risk_allele: "A", ln_or: 0.271, disease: "Ovarian cancer", gene: "BRCA2",  rsid: "rs11571833" },
    GwasVariant { chrom: "8",  pos: 75769584,  ref_allele: "A", risk_allele: "G", ln_or: 0.118, disease: "Ovarian cancer", gene: "GPR39",  rsid: "rs10088218" },
];

/// PRS result for a single disease.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrsResult {
    /// Disease category.
    pub disease: String,
    /// Raw polygenic risk score (sum of ln_OR for matched risk alleles).
    pub raw_score: f64,
    /// Z-score relative to population baseline (positive = above average risk).
    pub z_score: f64,
    /// "HIGH" (z>1.5), "ABOVE_AVERAGE" (z 0.5-1.5), "AVERAGE" (z -0.5 to 0.5), "BELOW_AVERAGE" (z<-0.5).
    pub risk_class: String,
    /// Number of GWAS catalog variants matched in VCF.
    pub n_variants_matched: usize,
    /// Total GWAS catalog variants for this disease.
    pub n_variants_total: usize,
    /// Matched risk variant details: (rsid, gene, ln_or).
    pub matched_variants: Vec<(String, String, f64)>,
}

/// Compute polygenic risk scores across all disease categories.
pub fn compute_prs(variants: &[VariantRecord]) -> Vec<PrsResult> {
    // Build a fast lookup: (chrom_norm, pos) -> VariantRecord refs
    use std::collections::HashMap;
    let mut variant_map: HashMap<(String, u64), Vec<&VariantRecord>> = HashMap::new();
    for v in variants {
        let chrom_norm = v.chrom.trim_start_matches("chr").to_ascii_uppercase();
        variant_map.entry((chrom_norm, v.pos)).or_default().push(v);
    }

    // Group GWAS catalog by disease
    let mut disease_scores: HashMap<&str, (f64, usize, Vec<(String, String, f64)>)> = HashMap::new();
    let mut disease_totals: HashMap<&str, usize> = HashMap::new();

    for gv in GWAS_CATALOG {
        *disease_totals.entry(gv.disease).or_insert(0) += 1;

        let chrom_norm = gv.chrom.to_ascii_uppercase();
        if let Some(vcf_variants) = variant_map.get(&(chrom_norm, gv.pos)) {
            for vcf_v in vcf_variants {
                if vcf_v.ref_allele.eq_ignore_ascii_case(gv.ref_allele)
                    && vcf_v.alt_allele.eq_ignore_ascii_case(gv.risk_allele)
                {
                    let entry = disease_scores.entry(gv.disease).or_insert((0.0, 0, Vec::new()));
                    entry.0 += gv.ln_or;
                    entry.1 += 1;
                    entry.2.push((gv.rsid.to_string(), gv.gene.to_string(), gv.ln_or));
                }
            }
        }
    }

    // Build results for all diseases, including those with 0 matches
    let mut results: Vec<PrsResult> = Vec::new();
    let diseases: Vec<&str> = {
        let mut d: Vec<&str> = GWAS_CATALOG.iter().map(|g| g.disease).collect();
        d.sort_unstable();
        d.dedup();
        d
    };

    for disease in diseases {
        let n_total = *disease_totals.get(disease).unwrap_or(&0);
        let (raw_score, n_matched, matched_variants) = disease_scores
            .get(disease)
            .cloned()
            .unwrap_or((0.0, 0, Vec::new()));

        let (z_score, risk_class) = DISEASE_BASELINES
            .iter()
            .find(|b| b.disease == disease)
            .map(|b| {
                let z = if b.std_score > 0.0 {
                    (raw_score - b.mean_score) / b.std_score
                } else {
                    0.0
                };
                let cls = if z > 1.5 {
                    "HIGH".to_string()
                } else if z > 0.5 {
                    "ABOVE_AVERAGE".to_string()
                } else if z >= -0.5 {
                    "AVERAGE".to_string()
                } else {
                    "BELOW_AVERAGE".to_string()
                };
                (z, cls)
            })
            .unwrap_or((0.0, "AVERAGE".to_string()));

        results.push(PrsResult {
            disease: disease.to_string(),
            raw_score,
            z_score,
            risk_class,
            n_variants_matched: n_matched,
            n_variants_total: n_total,
            matched_variants,
        });
    }

    // Sort by z_score descending
    results.sort_by(|a, b| b.z_score.partial_cmp(&a.z_score).unwrap_or(std::cmp::Ordering::Equal));
    results
}
