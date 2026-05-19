//! Cell-cell communication analysis.
//!
//! Implements a CellChat-style ligand-receptor scoring approach:
//! for each LR pair and each (sender, receiver) cluster pair, compute
//! the geometric mean of mean ligand expression in the sender cluster and
//! mean receptor expression in the receiver cluster, with permutation-based
//! significance testing.
//!
//! # References
//! - Jin S, et al. (2021) "Inference and analysis of cell-cell communication using CellChat"
//!   Nature Communications 12:1088.
//! - Browaeys R, et al. (2020) "NicheNet: modeling intercellular communication by linking
//!   ligands to target genes" Nature Methods 17:159–162.

use ahash::AHashMap;
use ndarray::Array2;
use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// A ligand-receptor pair.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LRPair {
    pub ligand: String,
    pub receptor: String,
    /// Signaling pathway name, e.g. "WNT", "NOTCH", "EGF".
    pub pathway: String,
    /// Annotation, e.g. "Secreted Signaling".
    pub annotation: String,
}

/// Communication score between two cell clusters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CommScore {
    pub sender_cluster: u32,
    pub receiver_cluster: u32,
    pub ligand: String,
    pub receptor: String,
    /// Geometric mean of ligand expr in sender × receptor expr in receiver.
    pub score: f64,
    /// Fraction of permutations with score >= observed (100 permutations).
    pub p_value: f64,
}

/// Summary of cell-cell communication.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CommSummary {
    pub scores: Vec<CommScore>,
    /// n_clusters × n_clusters interaction count matrix (significant pairs per cluster pair).
    pub interaction_counts: Vec<Vec<u32>>,
    /// Top ligand-receptor pairs by mean score: (ligand, receptor, mean_score).
    pub top_lr_pairs: Vec<(String, String, f64)>,
}

// ---------------------------------------------------------------------------
// Built-in LR database (~30 canonical pairs)
// ---------------------------------------------------------------------------

/// Return the built-in ligand-receptor database with ~30 canonical pairs
/// spanning key signaling pathways.
pub fn builtin_lr_database() -> Vec<LRPair> {
    let raw: &[(&str, &str, &str, &str)] = &[
        // NOTCH
        ("DLL1", "NOTCH1", "NOTCH", "Membrane-bound Signaling"),
        ("DLL4", "NOTCH1", "NOTCH", "Membrane-bound Signaling"),
        ("JAG1", "NOTCH2", "NOTCH", "Membrane-bound Signaling"),
        // WNT
        ("WNT5A", "FZD1", "WNT", "Secreted Signaling"),
        ("WNT3A", "FZD4", "WNT", "Secreted Signaling"),
        // TGFb
        ("TGFB1", "TGFBR1", "TGFb", "Secreted Signaling"),
        ("BMP4", "BMPR1A", "TGFb", "Secreted Signaling"),
        // EGF
        ("EGF", "EGFR", "EGF", "Secreted Signaling"),
        ("EREG", "EGFR", "EGF", "Secreted Signaling"),
        // FGF
        ("FGF2", "FGFR1", "FGF", "Secreted Signaling"),
        ("FGF7", "FGFR2", "FGF", "Secreted Signaling"),
        // VEGF
        ("VEGFA", "KDR", "VEGF", "Secreted Signaling"),
        ("VEGFC", "FLT4", "VEGF", "Secreted Signaling"),
        // Immune checkpoint
        ("CD274", "PDCD1", "Checkpoint", "Membrane-bound Signaling"),
        ("CD80", "CTLA4", "Checkpoint", "Membrane-bound Signaling"),
        // Cytokines
        ("IL6", "IL6R", "Cytokine", "Secreted Signaling"),
        ("TNF", "TNFRSF1A", "Cytokine", "Secreted Signaling"),
        ("CXCL12", "CXCR4", "Cytokine", "Secreted Signaling"),
        ("CCL2", "CCR2", "Cytokine", "Secreted Signaling"),
        ("IL1B", "IL1R1", "Cytokine", "Secreted Signaling"),
        // Growth factors
        ("IGF1", "IGF1R", "IGF", "Secreted Signaling"),
        ("HGF", "MET", "HGF", "Secreted Signaling"),
        ("PDGFA", "PDGFRA", "PDGF", "Secreted Signaling"),
        // Cell adhesion
        ("CADM1", "CADM1", "Adhesion", "Cell-Cell Contact"),
        ("ICAM1", "ITGAL", "Adhesion", "Cell-Cell Contact"),
        ("FN1", "ITGA5", "Adhesion", "ECM-Receptor"),
    ];

    raw.iter()
        .map(|&(lig, rec, path, ann)| LRPair {
            ligand: lig.to_owned(),
            receptor: rec.to_owned(),
            pathway: path.to_owned(),
            annotation: ann.to_owned(),
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Core computation
// ---------------------------------------------------------------------------

/// Compute cell-cell communication scores.
///
/// For each LR pair and each (sender, receiver) cluster pair:
///   `score = mean_expr(ligand, sender) × mean_expr(receptor, receiver)`
///   (geometric mean, i.e. `sqrt(mean_lig * mean_rec)`)
///   `p_value = fraction of 100 permutations where score >= observed`
///
/// # Arguments
/// * `count_matrix` – n_cells × n_genes normalized log counts.
/// * `gene_ids` – gene names (length n_genes).
/// * `cluster_labels` – cluster label per cell.
/// * `lr_pairs` – ligand-receptor pairs to test.
/// * `min_pct_expr` – minimum fraction of cells in a cluster that must express
///   a gene (> 0) for the gene to be considered expressed in that cluster.
/// * `seed` – RNG seed for the permutation test.
pub fn compute_communication(
    count_matrix: &Array2<f32>,
    gene_ids: &[String],
    cluster_labels: &[u32],
    lr_pairs: &[LRPair],
    min_pct_expr: f64,
    seed: u64,
) -> CommSummary {
    let n_cells = count_matrix.nrows();
    let n_genes = count_matrix.ncols();

    if n_cells == 0 || n_genes == 0 || lr_pairs.is_empty() {
        return CommSummary {
            scores: Vec::new(),
            interaction_counts: Vec::new(),
            top_lr_pairs: Vec::new(),
        };
    }

    // Build gene index map
    let gene_idx: AHashMap<&str, usize> = gene_ids
        .iter()
        .enumerate()
        .map(|(i, g)| (g.as_str(), i))
        .collect();

    // Identify unique clusters and their cell indices
    let n_clusters = cluster_labels
        .iter()
        .copied()
        .max()
        .map_or(0, |m| m as usize + 1);
    if n_clusters == 0 {
        return CommSummary {
            scores: Vec::new(),
            interaction_counts: Vec::new(),
            top_lr_pairs: Vec::new(),
        };
    }

    // cluster_cells[c] = sorted list of cell indices in cluster c
    let mut cluster_cells: Vec<Vec<usize>> = vec![Vec::new(); n_clusters];
    for (cell_idx, &label) in cluster_labels.iter().enumerate() {
        if (label as usize) < n_clusters {
            cluster_cells[label as usize].push(cell_idx);
        }
    }

    // Precompute per-cluster, per-gene mean expression and pct expressed
    // means[c][g], pcts[c][g]
    let mut means = vec![vec![0.0f64; n_genes]; n_clusters];
    let mut pcts = vec![vec![0.0f64; n_genes]; n_clusters];
    for c in 0..n_clusters {
        let cells = &cluster_cells[c];
        if cells.is_empty() {
            continue;
        }
        let n_c = cells.len() as f64;
        for g in 0..n_genes {
            let mut sum = 0.0f64;
            let mut expressed = 0usize;
            for &ci in cells {
                let v = count_matrix[[ci, g]] as f64;
                sum += v;
                if v > 0.0 {
                    expressed += 1;
                }
            }
            means[c][g] = sum / n_c;
            pcts[c][g] = expressed as f64 / n_c;
        }
    }

    // Permutation RNG (xorshift64)
    let mut rng_state = seed ^ 0x9e3779b97f4a7c15;
    let mut xorshift = move || -> u64 {
        rng_state ^= rng_state << 13;
        rng_state ^= rng_state >> 7;
        rng_state ^= rng_state << 17;
        rng_state
    };

    const N_PERM: usize = 100;

    let mut scores: Vec<CommScore> = Vec::new();
    // interaction_counts[sender][receiver] = # significant LR pairs
    let mut interaction_counts: Vec<Vec<u32>> = vec![vec![0u32; n_clusters]; n_clusters];

    // Collect LR scores per (ligand, receptor) key for top_lr_pairs
    let mut lr_score_sums: AHashMap<(String, String), (f64, usize)> = AHashMap::new();

    for lr in lr_pairs {
        let lig_idx = match gene_idx.get(lr.ligand.as_str()) {
            Some(&i) => i,
            None => continue,
        };
        let rec_idx = match gene_idx.get(lr.receptor.as_str()) {
            Some(&i) => i,
            None => continue,
        };

        for sender in 0..n_clusters {
            let sender_cells = &cluster_cells[sender];
            if sender_cells.is_empty() {
                continue;
            }
            // Ligand must be expressed in > min_pct_expr of sender cells
            if pcts[sender][lig_idx] < min_pct_expr {
                continue;
            }
            let mean_lig = means[sender][lig_idx];

            for receiver in 0..n_clusters {
                let receiver_cells = &cluster_cells[receiver];
                if receiver_cells.is_empty() {
                    continue;
                }
                // Receptor must be expressed in > min_pct_expr of receiver cells
                if pcts[receiver][rec_idx] < min_pct_expr {
                    continue;
                }
                let mean_rec = means[receiver][rec_idx];

                // Geometric mean score
                let score = (mean_lig * mean_rec).sqrt();

                // Permutation test: permute cluster labels 100 times
                let n_perm_ge = run_permutation_test(
                    count_matrix,
                    cluster_labels,
                    n_cells,
                    n_clusters,
                    lig_idx,
                    rec_idx,
                    sender as u32,
                    receiver as u32,
                    score,
                    N_PERM,
                    &mut xorshift,
                );

                let p_value = n_perm_ge as f64 / N_PERM as f64;

                // Track interaction counts (score > 0 counts as a communication)
                if score > 0.0 {
                    interaction_counts[sender][receiver] += 1;
                }

                // Accumulate per-LR pair mean score
                let entry = lr_score_sums
                    .entry((lr.ligand.clone(), lr.receptor.clone()))
                    .or_insert((0.0, 0));
                entry.0 += score;
                entry.1 += 1;

                scores.push(CommScore {
                    sender_cluster: sender as u32,
                    receiver_cluster: receiver as u32,
                    ligand: lr.ligand.clone(),
                    receptor: lr.receptor.clone(),
                    score,
                    p_value,
                });
            }
        }
    }

    // Top LR pairs by mean score (descending)
    let mut top_lr_pairs: Vec<(String, String, f64)> = lr_score_sums
        .into_iter()
        .map(|((lig, rec), (sum, cnt))| {
            let mean = if cnt > 0 { sum / cnt as f64 } else { 0.0 };
            (lig, rec, mean)
        })
        .collect();
    top_lr_pairs.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));

    CommSummary {
        scores,
        interaction_counts,
        top_lr_pairs,
    }
}

/// Run a permutation test for a single (sender, receiver, ligand, receptor) tuple.
///
/// Returns the number of permutations where the permuted score >= `observed_score`.
#[allow(clippy::too_many_arguments)]
fn run_permutation_test(
    count_matrix: &Array2<f32>,
    cluster_labels: &[u32],
    n_cells: usize,
    n_clusters: usize,
    lig_idx: usize,
    rec_idx: usize,
    sender: u32,
    receiver: u32,
    observed_score: f64,
    n_perm: usize,
    xorshift: &mut impl FnMut() -> u64,
) -> usize {
    let mut ge_count = 0usize;

    // Build a mutable permuted label array
    let mut perm_labels: Vec<u32> = cluster_labels.to_vec();

    for _ in 0..n_perm {
        // Fisher-Yates shuffle of perm_labels
        for i in (1..n_cells).rev() {
            let j = (xorshift() as usize) % (i + 1);
            perm_labels.swap(i, j);
        }

        // Compute mean ligand in permuted sender
        let (mean_lig_p, pct_lig_p) = cluster_mean_pct(
            count_matrix,
            &perm_labels,
            sender,
            lig_idx,
            n_cells,
            n_clusters,
        );
        if pct_lig_p == 0.0 {
            // effectively no expression → score = 0
            if observed_score <= 0.0 {
                ge_count += 1;
            }
            continue;
        }

        let (mean_rec_p, pct_rec_p) = cluster_mean_pct(
            count_matrix,
            &perm_labels,
            receiver,
            rec_idx,
            n_cells,
            n_clusters,
        );
        if pct_rec_p == 0.0 {
            if observed_score <= 0.0 {
                ge_count += 1;
            }
            continue;
        }

        let perm_score = (mean_lig_p * mean_rec_p).sqrt();
        if perm_score >= observed_score {
            ge_count += 1;
        }
    }

    ge_count
}

/// Compute mean expression and percent expressed for a given gene in a cluster,
/// given a (possibly permuted) label assignment.
#[inline]
fn cluster_mean_pct(
    count_matrix: &Array2<f32>,
    labels: &[u32],
    cluster: u32,
    gene_idx: usize,
    n_cells: usize,
    _n_clusters: usize,
) -> (f64, f64) {
    let mut sum = 0.0f64;
    let mut expressed = 0usize;
    let mut count = 0usize;

    for i in 0..n_cells {
        if labels[i] == cluster {
            let v = count_matrix[[i, gene_idx]] as f64;
            sum += v;
            if v > 0.0 {
                expressed += 1;
            }
            count += 1;
        }
    }

    if count == 0 {
        return (0.0, 0.0);
    }
    (sum / count as f64, expressed as f64 / count as f64)
}

/// Filter `CommSummary` to only significant interactions (p_value < p_threshold).
pub fn filter_significant(summary: &CommSummary, p_threshold: f64) -> CommSummary {
    let scores: Vec<CommScore> = summary
        .scores
        .iter()
        .filter(|s| s.p_value < p_threshold)
        .cloned()
        .collect();

    // Recompute interaction_counts from filtered scores
    let n_clusters = summary.interaction_counts.len();
    let mut interaction_counts: Vec<Vec<u32>> = vec![vec![0u32; n_clusters]; n_clusters];
    for s in &scores {
        let src = s.sender_cluster as usize;
        let dst = s.receiver_cluster as usize;
        if src < n_clusters && dst < n_clusters {
            interaction_counts[src][dst] += 1;
        }
    }

    // Recompute top LR pairs
    let mut lr_score_sums: AHashMap<(&str, &str), (f64, usize)> = AHashMap::new();
    for s in &scores {
        let entry = lr_score_sums
            .entry((s.ligand.as_str(), s.receptor.as_str()))
            .or_insert((0.0, 0));
        entry.0 += s.score;
        entry.1 += 1;
    }
    let mut top_lr_pairs: Vec<(String, String, f64)> = lr_score_sums
        .into_iter()
        .map(|((lig, rec), (sum, cnt))| {
            let mean = if cnt > 0 { sum / cnt as f64 } else { 0.0 };
            (lig.to_owned(), rec.to_owned(), mean)
        })
        .collect();
    top_lr_pairs.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));

    CommSummary {
        scores,
        interaction_counts,
        top_lr_pairs,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    #[test]
    fn lr_database_not_empty() {
        let db = builtin_lr_database();
        assert!(db.len() > 20, "expected > 20 LR pairs, got {}", db.len());
    }

    /// With ligand expressed in cluster 0 and receptor expressed in cluster 1,
    /// the communication score should be > 0.
    #[test]
    fn communication_score_positive() {
        let n_cells = 10;
        let n_genes = 2; // gene 0 = ligand, gene 1 = receptor

        // Cells 0-4 are cluster 0, cells 5-9 are cluster 1.
        // Cluster 0 expresses gene 0 (ligand), cluster 1 expresses gene 1 (receptor).
        let counts = Array2::<f32>::from_shape_fn((n_cells, n_genes), |(i, g)| {
            if i < 5 && g == 0 {
                2.0 // ligand in cluster 0
            } else if i >= 5 && g == 1 {
                3.0 // receptor in cluster 1
            } else {
                0.0
            }
        });

        let cluster_labels: Vec<u32> = (0..n_cells).map(|i| if i < 5 { 0 } else { 1 }).collect();
        let gene_ids = vec!["LIG1".to_owned(), "REC1".to_owned()];

        let lr_pairs = vec![LRPair {
            ligand: "LIG1".to_owned(),
            receptor: "REC1".to_owned(),
            pathway: "TEST".to_owned(),
            annotation: "Secreted Signaling".to_owned(),
        }];

        let summary =
            compute_communication(&counts, &gene_ids, &cluster_labels, &lr_pairs, 0.1, 42);

        // Find score for sender=0, receiver=1
        let comm = summary
            .scores
            .iter()
            .find(|s| s.sender_cluster == 0 && s.receiver_cluster == 1);
        assert!(comm.is_some(), "no CommScore for (0→1)");
        assert!(
            comm.unwrap().score > 0.0,
            "expected score > 0, got {}",
            comm.unwrap().score
        );
    }

    /// interaction_counts should be n_clusters × n_clusters.
    #[test]
    fn comm_summary_shape() {
        let n_cells = 12;
        let n_genes = 2;
        let n_clusters = 3usize;

        let counts = Array2::<f32>::from_shape_fn((n_cells, n_genes), |(i, g)| {
            if g == 0 {
                i as f32 + 1.0
            } else {
                (n_cells - i) as f32
            }
        });
        let cluster_labels: Vec<u32> = (0..n_cells).map(|i| (i % n_clusters) as u32).collect();
        let gene_ids = vec!["G0".to_owned(), "G1".to_owned()];

        let lr_pairs = vec![LRPair {
            ligand: "G0".to_owned(),
            receptor: "G1".to_owned(),
            pathway: "TEST".to_owned(),
            annotation: "Secreted Signaling".to_owned(),
        }];

        let summary = compute_communication(&counts, &gene_ids, &cluster_labels, &lr_pairs, 0.0, 1);

        assert_eq!(
            summary.interaction_counts.len(),
            n_clusters,
            "interaction_counts rows mismatch"
        );
        for row in &summary.interaction_counts {
            assert_eq!(row.len(), n_clusters, "interaction_counts cols mismatch");
        }
    }
}
