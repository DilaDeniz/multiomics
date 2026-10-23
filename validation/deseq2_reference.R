#!/usr/bin/env Rscript
# Regenerate the R/DESeq2 reference used by
# transcriptomics_core/tests/deseq2_validation.rs.
#
# Requires R with Bioconductor DESeq2 and jsonlite (see validation/environment.yml).
# Usage: Rscript validation/deseq2_reference.R   (run from the repository root)

suppressPackageStartupMessages({
  library(DESeq2)
  library(jsonlite)
})

out_dir <- "tests/golden/deseq2"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create("test_data/validation", showWarnings = FALSE, recursive = TRUE)

run_reference <- function(counts, condition, name, batch = NULL, test = "Wald") {
  coldata <- data.frame(condition = factor(condition, levels = c("A", "B")),
                        row.names = colnames(counts))
  design <- ~ condition
  if (!is.null(batch)) {
    coldata$batch <- factor(batch)
    design <- ~ batch + condition
  }
  dds <- DESeqDataSetFromMatrix(counts, coldata, design = design)
  if (test == "LRT") {
    dds <- DESeq(dds, test = "LRT", reduced = ~ 1, quiet = TRUE)
    res <- results(dds)
  } else {
    dds <- DESeq(dds, quiet = TRUE)
    res <- results(dds, contrast = c("condition", "B", "A"))
  }
  fit <- attr(dispersionFunction(dds), "coefficients")
  num <- function(x) ifelse(is.na(x), NA, x)
  genes <- lapply(seq_len(nrow(dds)), function(i) list(
    gene = rownames(dds)[i],
    base_mean = num(res$baseMean[i]),
    disp_gene_est = num(mcols(dds)$dispGeneEst[i]),
    disp_fit = num(mcols(dds)$dispFit[i]),
    dispersion = num(mcols(dds)$dispersion[i]),
    log2_fold_change = num(res$log2FoldChange[i]),
    lfc_se = num(res$lfcSE[i]),
    stat = num(res$stat[i]),
    pvalue = num(res$pvalue[i]),
    padj = num(res$padj[i]),
    beta_iter = num(mcols(dds)$betaIter[i])
  ))
  ref <- list(
    deseq2_version = as.character(packageVersion("DESeq2")),
    r_version = R.version.string,
    design = paste(deparse(design), if (test == "LRT") "LRT vs ~ 1" else "Wald, condition B vs A"),
    samples = colnames(counts),
    condition = condition,
    batch = if (is.null(batch)) NULL else as.character(batch),
    size_factors = unname(sizeFactors(dds)),
    disp_fit_coefficients = list(asympt_disp = unname(fit["asymptDisp"]),
                                 extra_pois = unname(fit["extraPois"])),
    disp_prior_var = attr(dispersionFunction(dds), "dispPriorVar"),
    genes = genes
  )
  write_json(ref, file.path(out_dir, paste0(name, ".json")),
             auto_unbox = TRUE, digits = NA, na = "null", pretty = TRUE)
}

# 1. Small hand-made matrix (formerly a self-referential "golden" file).
small <- matrix(c(100, 120,  95,  200,  210,  195,
                  500, 480, 510, 1000,  990, 1010,
                   10,   8,  12,   10,    9,   11,
                  300, 320, 280,  150,  160,  140,
                   50,  55,  45,  100,   95,  105),
                nrow = 5, byrow = TRUE,
                dimnames = list(paste0("GENE_", LETTERS[1:5]), paste0("S", 1:6)))
run_reference(small, rep(c("A", "B"), each = 3), "small_5x6")

# 2. DESeq2's own simulator: 1 000 genes, 3 vs 3, including low counts and zeros.
set.seed(20261023)
sim <- makeExampleDESeqDataSet(n = 1000, m = 6, betaSD = 1)
sim_counts <- counts(sim)
colnames(sim_counts) <- paste0("S", 1:6)
rownames(sim_counts) <- sprintf("gene%04d", seq_len(nrow(sim_counts)))
write.table(data.frame(gene_id = rownames(sim_counts), sim_counts, check.names = FALSE),
            "test_data/validation/deseq2_sim_1000x6.counts.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
run_reference(sim_counts, rep(c("A", "B"), each = 3), "sim_1000x6")
run_reference(sim_counts, rep(c("A", "B"), each = 3), "sim_1000x6_lrt", test = "LRT")

# 3. Two-factor design (~ batch + condition): four distinct design rows for three
#    coefficients, so fitted means come from the NB GLM rather than group means.
set.seed(20261024)
sim8 <- makeExampleDESeqDataSet(n = 1000, m = 8, betaSD = 1)
sim8_counts <- counts(sim8)
colnames(sim8_counts) <- paste0("S", 1:8)
rownames(sim8_counts) <- sprintf("gene%04d", seq_len(nrow(sim8_counts)))
write.table(data.frame(gene_id = rownames(sim8_counts), sim8_counts, check.names = FALSE),
            "test_data/validation/deseq2_sim_1000x8.counts.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
run_reference(sim8_counts, rep(c("A", "B"), each = 4), "sim_1000x8_batch",
              batch = c(1, 2, 1, 2, 1, 2, 2, 1))
