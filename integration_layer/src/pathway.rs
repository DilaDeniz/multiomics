use serde::{Deserialize, Serialize};

/// A KEGG pathway with its associated gene set for enrichment testing.
pub struct KeggPathway {
    pub id: &'static str,
    pub name: &'static str,
    pub genes: &'static [&'static str],
}

/// Top-50 KEGG pathways covering cancer, signaling, metabolism, and immune function.
pub static KEGG_PATHWAYS: &[KeggPathway] = &[
    KeggPathway {
        id: "hsa05210",
        name: "Colorectal cancer",
        genes: &["APC", "CTNNB1", "TP53", "KRAS", "SMAD4", "PIK3CA", "BRAF", "MLH1", "MSH2", "MSH6", "AXIN1", "AXIN2", "TCF7L2"],
    },
    KeggPathway {
        id: "hsa05212",
        name: "Pancreatic cancer",
        genes: &["KRAS", "TP53", "CDKN2A", "SMAD4", "BRCA2", "PALB2", "ATM", "BRAF", "MKK4", "RB1CC1"],
    },
    KeggPathway {
        id: "hsa05215",
        name: "Prostate cancer",
        genes: &["AR", "PTEN", "TP53", "RB1", "NKX3-1", "CDKN1B", "PIK3CA", "AKT1", "BRCA2", "CDK1"],
    },
    KeggPathway {
        id: "hsa05219",
        name: "Bladder cancer",
        genes: &["TP53", "RB1", "FGFR3", "HRAS", "CDKN2A", "PIK3CA", "EGFR", "VEGFA", "MMP1", "MMP2"],
    },
    KeggPathway {
        id: "hsa05220",
        name: "Chronic myeloid leukemia",
        genes: &["BCR", "ABL1", "GRB2", "SOS1", "KRAS", "RAF1", "MAP2K1", "MAPK1", "MYC", "STAT5A", "STAT5B"],
    },
    KeggPathway {
        id: "hsa04110",
        name: "Cell cycle",
        genes: &["CDK1", "CDK2", "CDK4", "CDK6", "CCND1", "CCNE1", "CCNA2", "CCNB1", "RB1", "E2F1", "TP53", "CDKN1A", "CDKN2A", "ATM", "CHEK1", "CHEK2"],
    },
    KeggPathway {
        id: "hsa04151",
        name: "PI3K-Akt signaling",
        genes: &["PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PTEN", "AKT1", "AKT2", "AKT3", "MTOR", "TSC1", "TSC2", "GSK3B", "FOXO1", "MDM2", "TP53", "CCND1"],
    },
    KeggPathway {
        id: "hsa04010",
        name: "MAPK signaling",
        genes: &["KRAS", "BRAF", "RAF1", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "MAPK8", "MAPK9", "MAPK14", "JUN", "FOS", "MYC", "HRAS", "NRAS"],
    },
    KeggPathway {
        id: "hsa04310",
        name: "Wnt signaling",
        genes: &["WNT1", "WNT2", "WNT3A", "FZD1", "LRP5", "LRP6", "AXIN1", "APC", "CTNNB1", "TCF7L2", "MYC", "CCND1", "DVL1", "GSK3B"],
    },
    KeggPathway {
        id: "hsa04115",
        name: "p53 signaling",
        genes: &["TP53", "MDM2", "CDKN1A", "BBC3", "BAX", "CASP3", "CASP9", "APAF1", "CHEK1", "CHEK2", "ATM", "ATR", "GADD45A", "PUMA"],
    },
    KeggPathway {
        id: "hsa04210",
        name: "Apoptosis",
        genes: &["BCL2", "BCL2L1", "BAX", "BAK1", "CASP3", "CASP8", "CASP9", "FADD", "FAS", "TNFRSF10A", "TNFRSF10B", "APAF1", "CYCS", "BID"],
    },
    KeggPathway {
        id: "hsa04350",
        name: "TGF-beta signaling",
        genes: &["TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4", "SMAD7", "ID1", "ID2", "ID3", "MYC"],
    },
    KeggPathway {
        id: "hsa04012",
        name: "ErbB signaling",
        genes: &["EGFR", "ERBB2", "ERBB3", "ERBB4", "KRAS", "PIK3CA", "AKT1", "MAP2K1", "MAPK1", "SRC", "STAT3", "CCND1", "MYC"],
    },
    KeggPathway {
        id: "hsa04020",
        name: "Calcium signaling",
        genes: &["CALM1", "CALM2", "CAMK2A", "CAMK4", "PRKCA", "PRKCB", "PLCB1", "PLCG1", "ATP2A1", "ATP2A2", "RYR1", "RYR2", "IP3R1"],
    },
    KeggPathway {
        id: "hsa04630",
        name: "JAK-STAT signaling",
        genes: &["JAK1", "JAK2", "JAK3", "TYK2", "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6", "IL6", "IL10", "IFNG", "PIK3CA"],
    },
    KeggPathway {
        id: "hsa04620",
        name: "Toll-like receptor signaling",
        genes: &["TLR1", "TLR2", "TLR3", "TLR4", "TLR7", "TLR9", "MYD88", "TRIF", "IRF3", "NF-KB1", "RELA", "TNF", "IL6", "IL12A", "IFNB1"],
    },
    KeggPathway {
        id: "hsa04660",
        name: "T cell receptor signaling",
        genes: &["CD3D", "CD3E", "CD3G", "ZAP70", "LCK", "FYN", "LAT", "PLCG1", "PIK3CA", "AKT1", "NFATC1", "JUN", "FOS", "RELA", "IL2"],
    },
    KeggPathway {
        id: "hsa04662",
        name: "B cell receptor signaling",
        genes: &["CD19", "CD79A", "CD79B", "SYK", "BTK", "BLNK", "PLCG2", "PIK3CA", "AKT1", "NFKB1", "RELA", "JUN", "FOS"],
    },
    KeggPathway {
        id: "hsa04668",
        name: "TNF signaling",
        genes: &["TNF", "TNFRSF1A", "TNFRSF1B", "TRADD", "TRAF2", "RIPK1", "MAP3K7", "NFKB1", "RELA", "JUN", "FOS", "CASP8", "FADD", "BID"],
    },
    KeggPathway {
        id: "hsa04370",
        name: "VEGF signaling",
        genes: &["VEGFA", "VEGFB", "VEGFC", "VEGFD", "FLT1", "KDR", "FLT4", "PLCG1", "PIK3CA", "AKT1", "SRC", "FAK1", "MAPK1", "NOS3"],
    },
    KeggPathway {
        id: "hsa04520",
        name: "Adherens junction",
        genes: &["CDH1", "CDH2", "CTNNB1", "CSNK1A1", "APC", "AXIN1", "PTPN11", "FYN", "SRC", "EGFR", "FGFR1", "MET", "SNAI1"],
    },
    KeggPathway {
        id: "hsa04540",
        name: "Gap junction",
        genes: &["GJA1", "GJA4", "GJA5", "GJB1", "GJB2", "PRKG1", "PRKG2", "PRKCA", "MAPK1", "TUBA1A", "TUBB"],
    },
    KeggPathway {
        id: "hsa04810",
        name: "Regulation of actin cytoskeleton",
        genes: &["ACTB", "ACTA1", "ARPC2", "ARPC3", "CDC42", "RAC1", "RHOA", "PFN1", "VCL", "TLN1", "ITGA1", "ITGB1", "PTK2", "SRC"],
    },
    KeggPathway {
        id: "hsa04920",
        name: "Adipocytokine signaling",
        genes: &["ADIPOQ", "ADIPOR1", "ADIPOR2", "PPARA", "PPARG", "LEPR", "LEP", "SOCS3", "IRS1", "PIK3CA", "AKT1", "AMPK", "PRKAA1"],
    },
    KeggPathway {
        id: "hsa04910",
        name: "Insulin signaling",
        genes: &["INS", "INSR", "IRS1", "IRS2", "PIK3CA", "AKT1", "FOXO1", "GSK3B", "PHKG1", "GYS1", "PYGL", "SLC2A4", "KRAS", "RAF1"],
    },
    KeggPathway {
        id: "hsa00190",
        name: "Oxidative phosphorylation",
        genes: &["MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND5", "MT-CO1", "MT-CO2", "MT-CO3", "MT-ATP6", "MT-ATP8", "NDUFS1", "SDHA", "UQCRC1", "COX4I1", "ATP5A1"],
    },
    KeggPathway {
        id: "hsa00010",
        name: "Glycolysis / Gluconeogenesis",
        genes: &["HK1", "HK2", "GPI", "PFKM", "ALDOA", "GAPDH", "PGK1", "PGAM1", "ENO1", "PKM", "LDHA", "G6PC", "PCK1", "FBP1"],
    },
    KeggPathway {
        id: "hsa00020",
        name: "Citrate cycle (TCA cycle)",
        genes: &["CS", "ACO2", "IDH1", "IDH2", "OGDH", "SUCLA2", "SUCLG1", "SDHA", "FH", "MDH2", "PC", "PDHB", "PDHA1"],
    },
    KeggPathway {
        id: "hsa00071",
        name: "Fatty acid degradation",
        genes: &["ACSL1", "ACSL4", "CPT1A", "CPT2", "HADHA", "HADHB", "ACADM", "ACADL", "ACADVL", "ECHS1", "HADH", "ACAA2"],
    },
    KeggPathway {
        id: "hsa00240",
        name: "Pyrimidine metabolism",
        genes: &["UMPS", "UPP1", "UPP2", "TYMP", "TK1", "TK2", "DCTD", "CDA", "DTYMK", "NME1", "NME2", "RRM1", "RRM2"],
    },
    KeggPathway {
        id: "hsa04062",
        name: "Chemokine signaling",
        genes: &["CCL2", "CCL5", "CXCL8", "CXCL12", "CCR2", "CCR5", "CXCR3", "CXCR4", "GNB1", "GNG2", "PIK3CA", "AKT1", "KRAS", "RAF1"],
    },
    KeggPathway {
        id: "hsa04150",
        name: "mTOR signaling",
        genes: &["MTOR", "RPTOR", "RICTOR", "TSC1", "TSC2", "RHEB", "AKT1", "PTEN", "PIK3CA", "RPS6KB1", "EIF4EBP1", "DEPTOR", "MLST8"],
    },
    KeggPathway {
        id: "hsa04330",
        name: "Notch signaling",
        genes: &["NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "DLL1", "DLL3", "DLL4", "JAG1", "JAG2", "RBPJ", "HES1", "HES5", "MYC", "CCND1"],
    },
    KeggPathway {
        id: "hsa04340",
        name: "Hedgehog signaling",
        genes: &["SHH", "IHH", "DHH", "PTCH1", "PTCH2", "SMO", "GLI1", "GLI2", "GLI3", "SUFU", "HHIP", "CDON", "BOC"],
    },
    KeggPathway {
        id: "hsa04390",
        name: "Hippo signaling",
        genes: &["MST1", "MST2", "LATS1", "LATS2", "MOB1A", "SAV1", "YAP1", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "NF2", "CTGF", "CYR61"],
    },
    KeggPathway {
        id: "hsa04064",
        name: "NF-kB signaling",
        genes: &["NFKB1", "NFKB2", "RELA", "RELB", "REL", "IKBKA", "IKBKB", "IKBKG", "TNFRSF1A", "IL1B", "TNF", "LTA", "TRAF2", "TRAF3", "TRAF5", "TRAF6"],
    },
    KeggPathway {
        id: "hsa04217",
        name: "Necroptosis",
        genes: &["RIPK1", "RIPK3", "MLKL", "CASP8", "FADD", "TNFRSF1A", "TNF", "TLR3", "TLR4", "DAI", "TRIF", "PGAM5"],
    },
    KeggPathway {
        id: "hsa04141",
        name: "Protein processing in endoplasmic reticulum",
        genes: &["HSP90B1", "DNAJB11", "CALR", "CANX", "PDIA3", "ERP29", "ATF6", "IRE1", "PERK", "XBP1", "DDIT3", "BIP", "GRP78"],
    },
    KeggPathway {
        id: "hsa03030",
        name: "DNA replication",
        genes: &["POLA1", "POLA2", "POLB", "POLD1", "POLE", "RFC1", "RFC2", "RFC3", "RFC4", "RFC5", "PCNA", "MCMC2", "CDC6", "CDT1", "ORC1"],
    },
    KeggPathway {
        id: "hsa03430",
        name: "Mismatch repair",
        genes: &["MLH1", "MLH3", "MSH2", "MSH3", "MSH6", "PMS1", "PMS2", "PCNA", "RFC1", "POLD1", "EXO1", "RFC2"],
    },
    KeggPathway {
        id: "hsa03440",
        name: "Homologous recombination",
        genes: &["BRCA1", "BRCA2", "RAD51", "RAD51B", "RAD51C", "RAD51D", "XRCC2", "XRCC3", "PALB2", "NBN", "MRE11", "RAD50", "ATM"],
    },
    KeggPathway {
        id: "hsa03450",
        name: "Non-homologous end-joining",
        genes: &["XRCC6", "XRCC5", "PRKDC", "LIG4", "XRCC4", "DCLRE1C", "NHEJ1", "MRE11", "RAD50", "NBN", "TP53BP1"],
    },
    KeggPathway {
        id: "hsa04144",
        name: "Endocytosis",
        genes: &["EGFR", "ERBB2", "HGF", "MET", "EPS15", "EEA1", "RAB5A", "RAB5B", "RAB7A", "RAB11A", "VPS34", "CLTC", "AP2A1", "LDLR"],
    },
    KeggPathway {
        id: "hsa04145",
        name: "Phagosome",
        genes: &["ITGAM", "ITGB2", "CR3", "FCGR1A", "FCGR2A", "FCGR3A", "SYK", "PI3K", "RAC1", "CDC42", "LAMP1", "LAMP2", "CTSD", "CTSL"],
    },
    KeggPathway {
        id: "hsa04721",
        name: "Synaptic vesicle cycle",
        genes: &["SYP", "SYN1", "SYN2", "SNAP25", "VAMP1", "VAMP2", "STX1A", "NSF", "ATP6V0A1", "SLC6A1", "SLC32A1", "RAB3A", "CPLX1"],
    },
    KeggPathway {
        id: "hsa05200",
        name: "Pathways in cancer",
        genes: &["TP53", "KRAS", "EGFR", "PTEN", "PIK3CA", "AKT1", "MYC", "RB1", "CDKN2A", "CTNNB1", "VHL", "NF1", "NF2", "BRCA1", "BRCA2", "CDH1", "MLH1"],
    },
    KeggPathway {
        id: "hsa05225",
        name: "Hepatocellular carcinoma",
        genes: &["TP53", "CTNNB1", "AXIN1", "ARID1A", "ARID2", "TSC1", "TSC2", "PIK3CA", "PTEN", "CDKN2A", "MET", "FGF19", "VEGFA"],
    },
    KeggPathway {
        id: "hsa05226",
        name: "Gastric cancer",
        genes: &["TP53", "CDH1", "ARID1A", "PIK3CA", "ERBB2", "ERBB3", "KRAS", "FGFR2", "CCND1", "CCNE1", "MYC", "MLH1", "MSH2"],
    },
    KeggPathway {
        id: "hsa05230",
        name: "Central carbon metabolism in cancer",
        genes: &["HK2", "PKM", "LDHA", "PFKM", "G6PD", "GLS", "FASN", "IDH1", "IDH2", "SDHA", "FH", "EGFR", "MYC", "HIF1A", "KRAS"],
    },
];

/// Result of pathway enrichment analysis for a single pathway.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EnrichmentResult {
    pub pathway_id: String,
    pub pathway_name: String,
    /// Number of query genes found in this pathway.
    pub overlap: usize,
    pub pathway_size: usize,
    pub query_size: usize,
    /// Jaccard-like enrichment score: overlap / sqrt(pathway_size × query_size).
    pub score: f64,
}

/// Run enrichment of `query_genes` against the static KEGG table.
///
/// Only returns pathways with `overlap >= min_overlap`.
/// Results are sorted by score, descending.
pub fn enrichment_analysis(query_genes: &[String], min_overlap: usize) -> Vec<EnrichmentResult> {
    let query_set: std::collections::HashSet<String> = query_genes
        .iter()
        .map(|g| g.to_uppercase())
        .collect();

    if query_set.is_empty() {
        return Vec::new();
    }

    let mut results: Vec<EnrichmentResult> = KEGG_PATHWAYS
        .iter()
        .filter_map(|pathway| {
            let pathway_genes: std::collections::HashSet<String> = pathway
                .genes
                .iter()
                .map(|g| g.to_uppercase())
                .collect();

            let overlap = query_set.intersection(&pathway_genes).count();
            if overlap < min_overlap {
                return None;
            }

            let score = overlap as f64
                / ((pathway.genes.len() as f64) * (query_set.len() as f64)).sqrt();

            Some(EnrichmentResult {
                pathway_id: pathway.id.to_string(),
                pathway_name: pathway.name.to_string(),
                overlap,
                pathway_size: pathway.genes.len(),
                query_size: query_set.len(),
                score,
            })
        })
        .collect();

    results.sort_unstable_by(|a, b| {
        b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal)
    });

    results
}
