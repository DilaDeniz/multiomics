use std::collections::HashMap;

use biomics_core::BatchAccum;

use crate::types::{GeneRecord, GeneStats, TranscriptomicsSummary};

/// Lock-free accumulator for gene expression statistics.
///
/// During the parallel phase each worker accumulates running sums and
/// sum-of-squares per gene. `merge` combines maps by addition. `finalize`
/// derives mean and standard deviation from the combined sums.
pub struct TranscriptomicsAccum {
    /// gene_id → per-sample running sum of TPM values.
    pub gene_sums: HashMap<String, Vec<f64>>,
    /// gene_id → per-sample running sum of squared TPM values (for std).
    pub gene_sq_sums: HashMap<String, Vec<f64>>,
    /// gene_id → number of records contributing to the sums.
    pub gene_counts: HashMap<String, u64>,
    pub total_genes_seen: u64,
    pub sample_count: usize,
}

impl Default for TranscriptomicsAccum {
    fn default() -> Self {
        Self {
            gene_sums: HashMap::new(),
            gene_sq_sums: HashMap::new(),
            gene_counts: HashMap::new(),
            total_genes_seen: 0,
            sample_count: 0,
        }
    }
}

impl BatchAccum for TranscriptomicsAccum {
    type Record = GeneRecord;
    type Summary = TranscriptomicsSummary;

    fn process(&mut self, r: &GeneRecord) -> anyhow::Result<()> {
        self.total_genes_seen += 1;

        if self.sample_count == 0 {
            self.sample_count = r.samples.len();
        }

        let sums = self
            .gene_sums
            .entry(r.gene_id.clone())
            .or_insert_with(|| vec![0.0; r.samples.len()]);
        let sq_sums = self
            .gene_sq_sums
            .entry(r.gene_id.clone())
            .or_insert_with(|| vec![0.0; r.samples.len()]);
        let count = self.gene_counts.entry(r.gene_id.clone()).or_insert(0);

        for (i, &val) in r.samples.iter().enumerate() {
            if i < sums.len() {
                sums[i] += val;
                sq_sums[i] += val * val;
            }
        }
        *count += 1;

        Ok(())
    }

    fn merge(&mut self, other: Self) {
        self.total_genes_seen += other.total_genes_seen;
        if self.sample_count == 0 {
            self.sample_count = other.sample_count;
        }

        for (gene, sums) in other.gene_sums {
            let entry = self
                .gene_sums
                .entry(gene.clone())
                .or_insert_with(|| vec![0.0; sums.len()]);
            for (i, v) in sums.iter().enumerate() {
                if i < entry.len() {
                    entry[i] += v;
                }
            }
        }

        for (gene, sq_sums) in other.gene_sq_sums {
            let entry = self
                .gene_sq_sums
                .entry(gene.clone())
                .or_insert_with(|| vec![0.0; sq_sums.len()]);
            for (i, v) in sq_sums.iter().enumerate() {
                if i < entry.len() {
                    entry[i] += v;
                }
            }
        }

        for (gene, count) in other.gene_counts {
            *self.gene_counts.entry(gene).or_insert(0) += count;
        }
    }

    fn finalize(self) -> anyhow::Result<TranscriptomicsSummary> {
        let mut gene_stats = std::collections::HashMap::new();
        let mut expressed_genes = 0u64;
        let mut low_expression_genes = Vec::new();
        let mut mean_tpm_by_gene: Vec<(String, f64)> = Vec::new();

        for (gene, sums) in &self.gene_sums {
            let count = *self.gene_counts.get(gene.as_str()).unwrap_or(&1) as f64;
            let sq_sums = self.gene_sq_sums.get(gene.as_str());

            let means: Vec<f64> = sums.iter().map(|&s| s / count).collect();
            let overall_mean = if means.is_empty() {
                0.0
            } else {
                means.iter().sum::<f64>() / means.len() as f64
            };

            let std = if let Some(sq) = sq_sums {
                let variances: Vec<f64> = sq
                    .iter()
                    .zip(sums.iter())
                    .map(|(&sq_s, &s)| {
                        let mean_val = s / count;
                        (sq_s / count - mean_val * mean_val).max(0.0)
                    })
                    .collect();
                let overall_var = if variances.is_empty() {
                    0.0
                } else {
                    variances.iter().sum::<f64>() / variances.len() as f64
                };
                overall_var.sqrt()
            } else {
                0.0
            };

            let max = means.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

            gene_stats.insert(
                gene.clone(),
                GeneStats {
                    mean: overall_mean,
                    std,
                    max,
                },
            );

            if overall_mean >= 1.0 {
                expressed_genes += 1;
            } else {
                low_expression_genes.push(gene.clone());
            }

            mean_tpm_by_gene.push((gene.clone(), overall_mean));
        }

        mean_tpm_by_gene
            .sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

        let top_100_expressed = mean_tpm_by_gene.into_iter().take(100).collect();

        Ok(TranscriptomicsSummary {
            total_genes: self.gene_sums.len() as u64,
            expressed_genes,
            low_expression_genes,
            gene_stats,
            top_100_expressed,
            diff_expr: None, // filled by finalize in diffexpr.rs after fold
            sample_count: self.sample_count,
            sample_names: Vec::new(), // filled by tsv.rs
        })
    }
}
