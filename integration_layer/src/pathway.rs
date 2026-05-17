use ahash::{AHashMap, AHashSet};
use serde::{Deserialize, Serialize};
use std::sync::OnceLock;

use biomics_core::statistics::{benjamini_hochberg, hypergeometric_pvalue};

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
        genes: &[
            "APC", "CTNNB1", "TP53", "KRAS", "SMAD4", "PIK3CA", "BRAF", "MLH1", "MSH2", "MSH6",
            "AXIN1", "AXIN2", "TCF7L2",
        ],
    },
    KeggPathway {
        id: "hsa05212",
        name: "Pancreatic cancer",
        genes: &[
            "KRAS", "TP53", "CDKN2A", "SMAD4", "BRCA2", "PALB2", "ATM", "BRAF", "MKK4", "RB1CC1",
        ],
    },
    KeggPathway {
        id: "hsa05215",
        name: "Prostate cancer",
        genes: &[
            "AR", "PTEN", "TP53", "RB1", "NKX3-1", "CDKN1B", "PIK3CA", "AKT1", "BRCA2", "CDK1",
        ],
    },
    KeggPathway {
        id: "hsa05219",
        name: "Bladder cancer",
        genes: &[
            "TP53", "RB1", "FGFR3", "HRAS", "CDKN2A", "PIK3CA", "EGFR", "VEGFA", "MMP1", "MMP2",
        ],
    },
    KeggPathway {
        id: "hsa05220",
        name: "Chronic myeloid leukemia",
        genes: &[
            "BCR", "ABL1", "GRB2", "SOS1", "KRAS", "RAF1", "MAP2K1", "MAPK1", "MYC", "STAT5A",
            "STAT5B",
        ],
    },
    KeggPathway {
        id: "hsa04110",
        name: "Cell cycle",
        genes: &[
            "CDK1", "CDK2", "CDK4", "CDK6", "CCND1", "CCNE1", "CCNA2", "CCNB1", "RB1", "E2F1",
            "TP53", "CDKN1A", "CDKN2A", "ATM", "CHEK1", "CHEK2",
        ],
    },
    KeggPathway {
        id: "hsa04151",
        name: "PI3K-Akt signaling",
        genes: &[
            "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PTEN", "AKT1", "AKT2", "AKT3", "MTOR", "TSC1",
            "TSC2", "GSK3B", "FOXO1", "MDM2", "TP53", "CCND1",
        ],
    },
    KeggPathway {
        id: "hsa04010",
        name: "MAPK signaling",
        genes: &[
            "KRAS", "BRAF", "RAF1", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "MAPK8", "MAPK9",
            "MAPK14", "JUN", "FOS", "MYC", "HRAS", "NRAS",
        ],
    },
    KeggPathway {
        id: "hsa04310",
        name: "Wnt signaling",
        genes: &[
            "WNT1", "WNT2", "WNT3A", "FZD1", "LRP5", "LRP6", "AXIN1", "APC", "CTNNB1", "TCF7L2",
            "MYC", "CCND1", "DVL1", "GSK3B",
        ],
    },
    KeggPathway {
        id: "hsa04115",
        name: "p53 signaling",
        genes: &[
            "TP53", "MDM2", "CDKN1A", "BBC3", "BAX", "CASP3", "CASP9", "APAF1", "CHEK1", "CHEK2",
            "ATM", "ATR", "GADD45A", "PUMA",
        ],
    },
    KeggPathway {
        id: "hsa04210",
        name: "Apoptosis",
        genes: &[
            "BCL2",
            "BCL2L1",
            "BAX",
            "BAK1",
            "CASP3",
            "CASP8",
            "CASP9",
            "FADD",
            "FAS",
            "TNFRSF10A",
            "TNFRSF10B",
            "APAF1",
            "CYCS",
            "BID",
        ],
    },
    KeggPathway {
        id: "hsa04350",
        name: "TGF-beta signaling",
        genes: &[
            "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4", "SMAD7",
            "ID1", "ID2", "ID3", "MYC",
        ],
    },
    KeggPathway {
        id: "hsa04012",
        name: "ErbB signaling",
        genes: &[
            "EGFR", "ERBB2", "ERBB3", "ERBB4", "KRAS", "PIK3CA", "AKT1", "MAP2K1", "MAPK1", "SRC",
            "STAT3", "CCND1", "MYC",
        ],
    },
    KeggPathway {
        id: "hsa04020",
        name: "Calcium signaling",
        genes: &[
            "CALM1", "CALM2", "CAMK2A", "CAMK4", "PRKCA", "PRKCB", "PLCB1", "PLCG1", "ATP2A1",
            "ATP2A2", "RYR1", "RYR2", "IP3R1",
        ],
    },
    KeggPathway {
        id: "hsa04630",
        name: "JAK-STAT signaling",
        genes: &[
            "JAK1", "JAK2", "JAK3", "TYK2", "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B",
            "STAT6", "IL6", "IL10", "IFNG", "PIK3CA",
        ],
    },
    KeggPathway {
        id: "hsa04620",
        name: "Toll-like receptor signaling",
        genes: &[
            "TLR1", "TLR2", "TLR3", "TLR4", "TLR7", "TLR9", "MYD88", "TRIF", "IRF3", "NF-KB1",
            "RELA", "TNF", "IL6", "IL12A", "IFNB1",
        ],
    },
    KeggPathway {
        id: "hsa04660",
        name: "T cell receptor signaling",
        genes: &[
            "CD3D", "CD3E", "CD3G", "ZAP70", "LCK", "FYN", "LAT", "PLCG1", "PIK3CA", "AKT1",
            "NFATC1", "JUN", "FOS", "RELA", "IL2",
        ],
    },
    KeggPathway {
        id: "hsa04662",
        name: "B cell receptor signaling",
        genes: &[
            "CD19", "CD79A", "CD79B", "SYK", "BTK", "BLNK", "PLCG2", "PIK3CA", "AKT1", "NFKB1",
            "RELA", "JUN", "FOS",
        ],
    },
    KeggPathway {
        id: "hsa04668",
        name: "TNF signaling",
        genes: &[
            "TNF", "TNFRSF1A", "TNFRSF1B", "TRADD", "TRAF2", "RIPK1", "MAP3K7", "NFKB1", "RELA",
            "JUN", "FOS", "CASP8", "FADD", "BID",
        ],
    },
    KeggPathway {
        id: "hsa04370",
        name: "VEGF signaling",
        genes: &[
            "VEGFA", "VEGFB", "VEGFC", "VEGFD", "FLT1", "KDR", "FLT4", "PLCG1", "PIK3CA", "AKT1",
            "SRC", "FAK1", "MAPK1", "NOS3",
        ],
    },
    KeggPathway {
        id: "hsa04520",
        name: "Adherens junction",
        genes: &[
            "CDH1", "CDH2", "CTNNB1", "CSNK1A1", "APC", "AXIN1", "PTPN11", "FYN", "SRC", "EGFR",
            "FGFR1", "MET", "SNAI1",
        ],
    },
    KeggPathway {
        id: "hsa04540",
        name: "Gap junction",
        genes: &[
            "GJA1", "GJA4", "GJA5", "GJB1", "GJB2", "PRKG1", "PRKG2", "PRKCA", "MAPK1", "TUBA1A",
            "TUBB",
        ],
    },
    KeggPathway {
        id: "hsa04810",
        name: "Regulation of actin cytoskeleton",
        genes: &[
            "ACTB", "ACTA1", "ARPC2", "ARPC3", "CDC42", "RAC1", "RHOA", "PFN1", "VCL", "TLN1",
            "ITGA1", "ITGB1", "PTK2", "SRC",
        ],
    },
    KeggPathway {
        id: "hsa04920",
        name: "Adipocytokine signaling",
        genes: &[
            "ADIPOQ", "ADIPOR1", "ADIPOR2", "PPARA", "PPARG", "LEPR", "LEP", "SOCS3", "IRS1",
            "PIK3CA", "AKT1", "AMPK", "PRKAA1",
        ],
    },
    KeggPathway {
        id: "hsa04910",
        name: "Insulin signaling",
        genes: &[
            "INS", "INSR", "IRS1", "IRS2", "PIK3CA", "AKT1", "FOXO1", "GSK3B", "PHKG1", "GYS1",
            "PYGL", "SLC2A4", "KRAS", "RAF1",
        ],
    },
    KeggPathway {
        id: "hsa00190",
        name: "Oxidative phosphorylation",
        genes: &[
            "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND5", "MT-CO1", "MT-CO2", "MT-CO3",
            "MT-ATP6", "MT-ATP8", "NDUFS1", "SDHA", "UQCRC1", "COX4I1", "ATP5A1",
        ],
    },
    KeggPathway {
        id: "hsa00010",
        name: "Glycolysis / Gluconeogenesis",
        genes: &[
            "HK1", "HK2", "GPI", "PFKM", "ALDOA", "GAPDH", "PGK1", "PGAM1", "ENO1", "PKM", "LDHA",
            "G6PC", "PCK1", "FBP1",
        ],
    },
    KeggPathway {
        id: "hsa00020",
        name: "Citrate cycle (TCA cycle)",
        genes: &[
            "CS", "ACO2", "IDH1", "IDH2", "OGDH", "SUCLA2", "SUCLG1", "SDHA", "FH", "MDH2", "PC",
            "PDHB", "PDHA1",
        ],
    },
    KeggPathway {
        id: "hsa00071",
        name: "Fatty acid degradation",
        genes: &[
            "ACSL1", "ACSL4", "CPT1A", "CPT2", "HADHA", "HADHB", "ACADM", "ACADL", "ACADVL",
            "ECHS1", "HADH", "ACAA2",
        ],
    },
    KeggPathway {
        id: "hsa00240",
        name: "Pyrimidine metabolism",
        genes: &[
            "UMPS", "UPP1", "UPP2", "TYMP", "TK1", "TK2", "DCTD", "CDA", "DTYMK", "NME1", "NME2",
            "RRM1", "RRM2",
        ],
    },
    KeggPathway {
        id: "hsa04062",
        name: "Chemokine signaling",
        genes: &[
            "CCL2", "CCL5", "CXCL8", "CXCL12", "CCR2", "CCR5", "CXCR3", "CXCR4", "GNB1", "GNG2",
            "PIK3CA", "AKT1", "KRAS", "RAF1",
        ],
    },
    KeggPathway {
        id: "hsa04150",
        name: "mTOR signaling",
        genes: &[
            "MTOR", "RPTOR", "RICTOR", "TSC1", "TSC2", "RHEB", "AKT1", "PTEN", "PIK3CA", "RPS6KB1",
            "EIF4EBP1", "DEPTOR", "MLST8",
        ],
    },
    KeggPathway {
        id: "hsa04330",
        name: "Notch signaling",
        genes: &[
            "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "DLL1", "DLL3", "DLL4", "JAG1", "JAG2", "RBPJ",
            "HES1", "HES5", "MYC", "CCND1",
        ],
    },
    KeggPathway {
        id: "hsa04340",
        name: "Hedgehog signaling",
        genes: &[
            "SHH", "IHH", "DHH", "PTCH1", "PTCH2", "SMO", "GLI1", "GLI2", "GLI3", "SUFU", "HHIP",
            "CDON", "BOC",
        ],
    },
    KeggPathway {
        id: "hsa04390",
        name: "Hippo signaling",
        genes: &[
            "MST1", "MST2", "LATS1", "LATS2", "MOB1A", "SAV1", "YAP1", "TEAD1", "TEAD2", "TEAD3",
            "TEAD4", "NF2", "CTGF", "CYR61",
        ],
    },
    KeggPathway {
        id: "hsa04064",
        name: "NF-kB signaling",
        genes: &[
            "NFKB1", "NFKB2", "RELA", "RELB", "REL", "IKBKA", "IKBKB", "IKBKG", "TNFRSF1A", "IL1B",
            "TNF", "LTA", "TRAF2", "TRAF3", "TRAF5", "TRAF6",
        ],
    },
    KeggPathway {
        id: "hsa04217",
        name: "Necroptosis",
        genes: &[
            "RIPK1", "RIPK3", "MLKL", "CASP8", "FADD", "TNFRSF1A", "TNF", "TLR3", "TLR4", "DAI",
            "TRIF", "PGAM5",
        ],
    },
    KeggPathway {
        id: "hsa04141",
        name: "Protein processing in endoplasmic reticulum",
        genes: &[
            "HSP90B1", "DNAJB11", "CALR", "CANX", "PDIA3", "ERP29", "ATF6", "IRE1", "PERK", "XBP1",
            "DDIT3", "BIP", "GRP78",
        ],
    },
    KeggPathway {
        id: "hsa03030",
        name: "DNA replication",
        genes: &[
            "POLA1", "POLA2", "POLB", "POLD1", "POLE", "RFC1", "RFC2", "RFC3", "RFC4", "RFC5",
            "PCNA", "MCMC2", "CDC6", "CDT1", "ORC1",
        ],
    },
    KeggPathway {
        id: "hsa03430",
        name: "Mismatch repair",
        genes: &[
            "MLH1", "MLH3", "MSH2", "MSH3", "MSH6", "PMS1", "PMS2", "PCNA", "RFC1", "POLD1",
            "EXO1", "RFC2",
        ],
    },
    KeggPathway {
        id: "hsa03440",
        name: "Homologous recombination",
        genes: &[
            "BRCA1", "BRCA2", "RAD51", "RAD51B", "RAD51C", "RAD51D", "XRCC2", "XRCC3", "PALB2",
            "NBN", "MRE11", "RAD50", "ATM",
        ],
    },
    KeggPathway {
        id: "hsa03450",
        name: "Non-homologous end-joining",
        genes: &[
            "XRCC6", "XRCC5", "PRKDC", "LIG4", "XRCC4", "DCLRE1C", "NHEJ1", "MRE11", "RAD50",
            "NBN", "TP53BP1",
        ],
    },
    KeggPathway {
        id: "hsa04144",
        name: "Endocytosis",
        genes: &[
            "EGFR", "ERBB2", "HGF", "MET", "EPS15", "EEA1", "RAB5A", "RAB5B", "RAB7A", "RAB11A",
            "VPS34", "CLTC", "AP2A1", "LDLR",
        ],
    },
    KeggPathway {
        id: "hsa04145",
        name: "Phagosome",
        genes: &[
            "ITGAM", "ITGB2", "CR3", "FCGR1A", "FCGR2A", "FCGR3A", "SYK", "PI3K", "RAC1", "CDC42",
            "LAMP1", "LAMP2", "CTSD", "CTSL",
        ],
    },
    KeggPathway {
        id: "hsa04721",
        name: "Synaptic vesicle cycle",
        genes: &[
            "SYP", "SYN1", "SYN2", "SNAP25", "VAMP1", "VAMP2", "STX1A", "NSF", "ATP6V0A1",
            "SLC6A1", "SLC32A1", "RAB3A", "CPLX1",
        ],
    },
    KeggPathway {
        id: "hsa05200",
        name: "Pathways in cancer",
        genes: &[
            "TP53", "KRAS", "EGFR", "PTEN", "PIK3CA", "AKT1", "MYC", "RB1", "CDKN2A", "CTNNB1",
            "VHL", "NF1", "NF2", "BRCA1", "BRCA2", "CDH1", "MLH1",
        ],
    },
    KeggPathway {
        id: "hsa05225",
        name: "Hepatocellular carcinoma",
        genes: &[
            "TP53", "CTNNB1", "AXIN1", "ARID1A", "ARID2", "TSC1", "TSC2", "PIK3CA", "PTEN",
            "CDKN2A", "MET", "FGF19", "VEGFA",
        ],
    },
    KeggPathway {
        id: "hsa05226",
        name: "Gastric cancer",
        genes: &[
            "TP53", "CDH1", "ARID1A", "PIK3CA", "ERBB2", "ERBB3", "KRAS", "FGFR2", "CCND1",
            "CCNE1", "MYC", "MLH1", "MSH2",
        ],
    },
    KeggPathway {
        id: "hsa05230",
        name: "Central carbon metabolism in cancer",
        genes: &[
            "HK2", "PKM", "LDHA", "PFKM", "G6PD", "GLS", "FASN", "IDH1", "IDH2", "SDHA", "FH",
            "EGFR", "MYC", "HIF1A", "KRAS",
        ],
    },
    // ── Epigenetic regulation ─────────────────────────────────────────────────
    KeggPathway {
        id: "epi_dnmt",
        name: "DNA methylation (DNMT)",
        genes: &[
            "DNMT1", "DNMT3A", "DNMT3B", "DNMT3L", "TET1", "TET2", "TET3", "UHRF1", "MBD1", "MBD2",
            "MBD4", "MECP2",
        ],
    },
    KeggPathway {
        id: "epi_histone_me",
        name: "Histone methylation",
        genes: &[
            "EZH2", "EZH1", "SUZ12", "EED", "SMYD2", "NSD1", "NSD2", "NSD3", "KMT2A", "KMT2B",
            "KMT2C", "KMT2D", "SETD2", "KDM1A", "KDM5C", "KDM6A",
        ],
    },
    KeggPathway {
        id: "epi_histone_ac",
        name: "Histone acetylation (HAT/HDAC)",
        genes: &[
            "HDAC1", "HDAC2", "HDAC3", "HDAC4", "HDAC5", "HDAC6", "HDAC7", "HDAC8", "EP300",
            "CREBBP", "KAT2A", "KAT2B", "SIRT1", "SIRT2", "SIRT3",
        ],
    },
    KeggPathway {
        id: "epi_swi_snf",
        name: "SWI/SNF chromatin remodeling",
        genes: &[
            "SMARCA4", "SMARCA2", "SMARCB1", "SMARCC1", "SMARCC2", "SMARCD1", "SMARCE1", "ARID1A",
            "ARID1B", "ARID2", "PBRM1", "BRD7", "BRD9",
        ],
    },
    KeggPathway {
        id: "epi_polycomb",
        name: "Polycomb repressive complex",
        genes: &[
            "EZH2", "EED", "SUZ12", "BMI1", "RING1", "RNF2", "RYBP", "CBX2", "CBX4", "CBX6",
            "CBX7", "CBX8", "JARID2", "PHF1", "PHF19",
        ],
    },
    // ── Immune checkpoints ───────────────────────────────────────────────────
    KeggPathway {
        id: "imm_checkpoint",
        name: "Immune checkpoint signaling",
        genes: &[
            "CD274", "PDCD1", "PDCD1LG2", "CTLA4", "CD80", "CD86", "LAG3", "TIM3", "HAVCR2",
            "TIGIT", "BTLA", "VSIR", "IDO1", "CD47", "SIRPA",
        ],
    },
    KeggPathway {
        id: "imm_tcell_exhaust",
        name: "T cell exhaustion",
        genes: &[
            "PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "TOX", "NR4A1", "NR4A2", "NR4A3", "EOMES",
            "TBX21", "PRDM1", "BATF", "IRF4",
        ],
    },
    KeggPathway {
        id: "imm_innate",
        name: "Innate immune sensing",
        genes: &[
            "STING1", "CGAS", "IRF3", "IRF7", "TBK1", "MAVS", "DDX58", "IFIH1", "NLRP3", "PYCARD",
            "CASP1", "IL1B", "IL18", "HMGB1",
        ],
    },
    KeggPathway {
        id: "imm_nk_cell",
        name: "NK cell cytotoxicity",
        genes: &[
            "NCR1", "NCR2", "NCR3", "KLRK1", "KIR2DL1", "KIR3DL1", "FCGR3A", "PRF1", "GZMA",
            "GZMB", "GZMK", "IFNG", "TNF", "FASL",
        ],
    },
    // ── DNA damage response ──────────────────────────────────────────────────
    KeggPathway {
        id: "ddr_atr",
        name: "ATR-mediated DNA damage response",
        genes: &[
            "ATR", "ATRIP", "TOPBP1", "CLSPN", "CHEK1", "RAD9A", "RAD1", "HUS1", "RPA1", "RPA2",
            "RPA3", "RAD17", "RHINO",
        ],
    },
    KeggPathway {
        id: "ddr_base_excision",
        name: "Base excision repair",
        genes: &[
            "OGG1", "MUTYH", "UNG", "SMUG1", "MBD4", "APEX1", "POLB", "XRCC1", "LIG3", "PARP1",
            "PARP2", "NEIL1", "NEIL2", "NEIL3",
        ],
    },
    KeggPathway {
        id: "ddr_nucleotide_excision",
        name: "Nucleotide excision repair",
        genes: &[
            "XPC", "RAD23B", "CETN2", "DDB1", "DDB2", "XPA", "RPA1", "ERCC1", "ERCC2", "ERCC3",
            "ERCC4", "ERCC5", "PCNA", "RFC1", "LIG1",
        ],
    },
    KeggPathway {
        id: "ddr_fanconi",
        name: "Fanconi anemia / BRCA pathway",
        genes: &[
            "FANCA", "FANCC", "FANCD2", "FANCE", "FANCF", "FANCG", "FANCI", "FANCJ", "FANCL",
            "FANCM", "BRCA1", "BRCA2", "PALB2", "RAD51", "USP1",
        ],
    },
    // ── Metabolic reprogramming ──────────────────────────────────────────────
    KeggPathway {
        id: "met_warburg",
        name: "Warburg effect / aerobic glycolysis",
        genes: &[
            "HK1", "HK2", "PFKFB3", "PFKFB4", "PKM", "LDHA", "LDHB", "PDK1", "PDK2", "HIF1A",
            "MYC", "SLC2A1", "SLC2A3", "MCT4",
        ],
    },
    KeggPathway {
        id: "met_glutamine",
        name: "Glutamine metabolism",
        genes: &[
            "GLS", "GLS2", "GLUD1", "GOT1", "GOT2", "ASNS", "ASPA", "SLC1A5", "SLC38A1", "SLC38A2",
            "PPAT", "MYC", "MTOR",
        ],
    },
    KeggPathway {
        id: "met_lipid",
        name: "Lipid biosynthesis",
        genes: &[
            "FASN", "ACACA", "ACACB", "ACLY", "ACSS2", "SREBF1", "SREBF2", "HMGCR", "SQLE",
            "DHCR24", "LDLR", "PPARG", "PPARA",
        ],
    },
    KeggPathway {
        id: "met_one_carbon",
        name: "One-carbon / folate metabolism",
        genes: &[
            "MTHFR", "MTR", "MTRR", "DHFR", "TYMS", "SHMT1", "SHMT2", "MTHFD1", "MTHFD2",
            "ALDH1L1", "ALDH1L2", "MAT1A", "MAT2A", "DNMT1",
        ],
    },
    // ── RNA biology / splicing ───────────────────────────────────────────────
    KeggPathway {
        id: "rna_splicing",
        name: "RNA splicing (spliceosome mutations)",
        genes: &[
            "SF3B1", "U2AF1", "SRSF2", "ZRSR2", "SF3A1", "SF3A2", "SF3A3", "PRPF8", "PRPF3",
            "PRPF31", "U2AF2", "RBM10", "RBM5", "HNRNPK",
        ],
    },
    KeggPathway {
        id: "rna_m6a",
        name: "m6A RNA methylation",
        genes: &[
            "METTL3", "METTL14", "WTAP", "ALKBH5", "FTO", "YTHDF1", "YTHDF2", "YTHDF3", "YTHDC1",
            "YTHDC2", "IGF2BP1", "IGF2BP2", "IGF2BP3",
        ],
    },
    // ── Developmental / stem cell ────────────────────────────────────────────
    KeggPathway {
        id: "stem_wnt_stem",
        name: "Wnt-driven stem cell self-renewal",
        genes: &[
            "CTNNB1", "LGR5", "AXIN2", "CD44", "PROM1", "ALDH1A1", "SOX2", "OCT4", "NANOG", "KLF4",
            "MYC", "SOX9", "ASCL2", "EphB2",
        ],
    },
    KeggPathway {
        id: "stem_emt",
        name: "Epithelial-mesenchymal transition (EMT)",
        genes: &[
            "CDH1", "VIM", "FN1", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1", "TWIST2", "MMP2",
            "MMP9", "TGFB1", "SMAD2", "SMAD3", "WNT5A",
        ],
    },
    // ── Signaling: additional ─────────────────────────────────────────────────
    KeggPathway {
        id: "sig_ros",
        name: "Reactive oxygen species / oxidative stress",
        genes: &[
            "SOD1", "SOD2", "CAT", "GPX1", "GPX4", "PRDX1", "PRDX2", "TXNRD1", "NQO1", "NFE2L2",
            "KEAP1", "HMOX1", "SRXN1", "PARK7",
        ],
    },
    KeggPathway {
        id: "sig_hippo_cancer",
        name: "Hippo pathway in cancer",
        genes: &[
            "YAP1", "WWTR1", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "LATS1", "LATS2", "MST1", "MST2",
            "NF2", "FAT1", "FAT4", "RASSF1", "RASSF2",
        ],
    },
    KeggPathway {
        id: "sig_erbb_resistance",
        name: "EGFR/ErbB resistance mechanisms",
        genes: &[
            "EGFR", "KRAS", "BRAF", "MEK1", "PIK3CA", "PTEN", "AKT1", "MET", "ERBB2", "ERBB3",
            "HGF", "IGF1R", "AXL", "FGFR1",
        ],
    },
    KeggPathway {
        id: "sig_tert",
        name: "Telomere maintenance / TERT",
        genes: &[
            "TERT", "TERC", "DKC1", "NOP10", "NHP2", "GAR1", "TINF2", "ACD", "TERF1", "TERF2",
            "POT1", "TPP1", "RAP1", "RTEL1",
        ],
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
    /// One-sided hypergeometric p-value (Fisher's exact test upper tail).
    pub p_value: f64,
    /// Benjamini-Hochberg FDR-adjusted p-value.
    pub padj: f64,
}

// ── Static inverted gene → pathway index (built once from KEGG_PATHWAYS) ─────
// gene (uppercase) → list of KEGG_PATHWAYS indices.
// Eliminates O(P × G) AHashSet rebuilds per enrichment call; lookup is O(1).

type GeneIndex = AHashMap<&'static str, Vec<u8>>; // u8 sufficient for ≤255 pathways

fn kegg_gene_index() -> &'static (GeneIndex, usize) {
    static IDX: OnceLock<(GeneIndex, usize)> = OnceLock::new();
    IDX.get_or_init(|| {
        let mut map: GeneIndex = AHashMap::new();
        let mut bg: AHashSet<&str> = AHashSet::new();
        for (i, pw) in KEGG_PATHWAYS.iter().enumerate() {
            for &gene in pw.genes {
                map.entry(gene).or_default().push(i as u8);
                bg.insert(gene);
            }
        }
        (map, bg.len())
    })
}

/// Fisher's exact test (one-sided hypergeometric) with Benjamini-Hochberg FDR
/// correction. Background gene universe = union of all genes in KEGG_PATHWAYS.
/// Only pathways with `overlap >= min_overlap` are returned.
/// Results are sorted by `padj` ascending.
pub fn enrichment_analysis(query_genes: &[String], min_overlap: usize) -> Vec<EnrichmentResult> {
    let query_set: AHashSet<String> = query_genes.iter().map(|g| g.to_uppercase()).collect();

    if query_set.is_empty() {
        return Vec::new();
    }

    let (gene_idx, bg_size) = kegg_gene_index();
    let query_size = query_set.len();

    // Count overlaps via inverted index: O(Q) instead of O(P × G)
    let mut overlap_counts = vec![0usize; KEGG_PATHWAYS.len()];
    for gene in &query_set {
        if let Some(idxs) = gene_idx.get(gene.as_str()) {
            for &i in idxs {
                overlap_counts[i as usize] += 1;
            }
        }
    }

    let mut results: Vec<EnrichmentResult> = KEGG_PATHWAYS
        .iter()
        .enumerate()
        .filter_map(|(i, pathway)| {
            let overlap = overlap_counts[i];
            if overlap < min_overlap {
                return None;
            }
            let score =
                overlap as f64 / ((pathway.genes.len() as f64) * (query_size as f64)).sqrt();
            let p_value = hypergeometric_pvalue(overlap, query_size, pathway.genes.len(), *bg_size);
            Some(EnrichmentResult {
                pathway_id: pathway.id.to_string(),
                pathway_name: pathway.name.to_string(),
                overlap,
                pathway_size: pathway.genes.len(),
                query_size,
                score,
                p_value,
                padj: f64::NAN,
            })
        })
        .collect();

    if results.is_empty() {
        return Vec::new();
    }

    // Apply Benjamini-Hochberg FDR correction across all results
    let pvals: Vec<f64> = results.iter().map(|r| r.p_value).collect();
    let padj_vals = benjamini_hochberg(&pvals);
    for (r, padj) in results.iter_mut().zip(padj_vals) {
        r.padj = padj;
    }

    // Sort by padj ascending (most significant first)
    results.sort_unstable_by(|a, b| {
        a.padj
            .partial_cmp(&b.padj)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    results
}
