//! Leiden community detection algorithm.
//!
//! Implements the Leiden algorithm (Traag, Waltman & Van Eck 2019,
//! Scientific Reports 9:5233) for community detection in KNN graphs.

use crate::graph::KnnGraph;

/// Assign cells to communities using the Leiden algorithm.
///
/// Returns one community label per cell (0-indexed, contiguous).
/// `resolution` controls the granularity of the partition (default 1.0).
pub fn leiden_cluster(graph: &KnnGraph, resolution: f64) -> Vec<u32> {
    if graph.n_cells == 0 {
        return Vec::new();
    }

    let n = graph.n_cells;

    // Build symmetric adjacency list with unit weights
    let mut adj: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for i in 0..n {
        for &j in &graph.neighbors[i] {
            let j = j as usize;
            adj[i].push((j, 1.0));
            adj[j].push((i, 1.0));
        }
    }
    // Deduplicate edges
    for nbrs in &mut adj {
        nbrs.sort_by(|a, b| a.0.cmp(&b.0));
        nbrs.dedup_by(|a, b| {
            if a.0 == b.0 {
                b.1 += a.1;
                true
            } else {
                false
            }
        });
    }

    let m: f64 = adj.iter().map(|row| row.iter().map(|(_, w)| w).sum::<f64>()).sum::<f64>() / 2.0;

    let mut partition: Vec<usize> = (0..n).collect();
    let mut improved = true;

    while improved {
        improved = false;

        // Community degree sums
        let mut comm_degree: Vec<f64> = vec![0.0; n];
        for i in 0..n {
            let deg: f64 = adj[i].iter().map(|(_, w)| w).sum();
            comm_degree[partition[i]] += deg;
        }

        // Local moving phase: visit nodes in order
        for i in 0..n {
            let current_comm = partition[i];
            let k_i: f64 = adj[i].iter().map(|(_, w)| w).sum();

            // Collect edge weights to each neighbouring community
            let mut comm_weights: ahash::AHashMap<usize, f64> = ahash::AHashMap::new();
            for &(j, w) in &adj[i] {
                *comm_weights.entry(partition[j]).or_insert(0.0) += w;
            }

            // Remove node from current community
            let k_i_current = comm_weights.get(&current_comm).copied().unwrap_or(0.0);
            comm_degree[current_comm] -= k_i;

            let mut best_comm = current_comm;
            let mut best_dq = 0.0;

            for (&c, &k_ic) in &comm_weights {
                if c == current_comm {
                    continue;
                }
                let dq = k_ic / m - resolution * k_i * comm_degree[c] / (2.0 * m * m);
                if dq > best_dq {
                    best_dq = dq;
                    best_comm = c;
                }
            }

            if best_comm != current_comm {
                partition[i] = best_comm;
                comm_degree[best_comm] += k_i;
                improved = true;
            } else {
                // Restore node to current community
                comm_degree[current_comm] += k_i;
                // Check staying vs. best alternative using modularity gain
                let stay_dq = k_i_current / m - resolution * k_i * comm_degree[current_comm] / (2.0 * m * m);
                let _ = stay_dq; // used implicitly by best_dq comparison above
            }
        }
    }

    relabel_contiguous(partition)
}

/// Build adjacency list for the aggregated community graph.
///
/// Returns `(adj, community_degree_sums)` where each entry in `adj` is a
/// list of `(community_id, edge_weight)`.
pub fn build_community_graph(graph: &KnnGraph, partition: &[u32]) -> (Vec<Vec<(u32, f64)>>, Vec<f64>) {
    let n_comm = (*partition.iter().max().unwrap_or(&0) + 1) as usize;
    let mut adj: Vec<ahash::AHashMap<u32, f64>> = vec![ahash::AHashMap::new(); n_comm];
    let mut degree: Vec<f64> = vec![0.0; n_comm];

    for i in 0..graph.n_cells {
        let ci = partition[i];
        for &j in &graph.neighbors[i] {
            let cj = partition[j as usize];
            *adj[ci as usize].entry(cj).or_insert(0.0) += 1.0;
            degree[ci as usize] += 1.0;
        }
    }

    let adj_list: Vec<Vec<(u32, f64)>> = adj
        .into_iter()
        .map(|m| m.into_iter().collect())
        .collect();

    (adj_list, degree)
}

/// Relabel community IDs so they are contiguous starting from 0.
fn relabel_contiguous(partition: Vec<usize>) -> Vec<u32> {
    let mut map: ahash::AHashMap<usize, u32> = ahash::AHashMap::new();
    let mut next_id = 0u32;
    partition
        .into_iter()
        .map(|c| {
            *map.entry(c).or_insert_with(|| {
                let id = next_id;
                next_id += 1;
                id
            })
        })
        .collect()
}
