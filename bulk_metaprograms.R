# Bulk metaprogram discovery via NMF
# Uses TMM-normalized CPM (non-log) as input to NMF (nonnegative)

library(edgeR)
library(NMF)
library(pheatmap)
library(ggplot2)

cat("╔═══════════════════════════════════════════════════════════════════╗\n")
cat("║        BULK MALIGNANT METAPROGRAM DISCOVERY (NMF)                 ║\n")
cat("╚═══════════════════════════════════════════════════════════════════╝\n\n")

# 1. Load cleaned counts -----------------------------------------------------

# We start from the cleaned counts produced by your cleaning pipeline
bulk_counts_filtered <- read.delim("Bulk_rna_raw_counts_CLEANED.tsv",
                                   row.names = 1,
                                   check.names = FALSE)

cat("Cleaned bulk counts dimensions: ",
    nrow(bulk_counts_filtered), "genes x",
    ncol(bulk_counts_filtered), "cell lines\n")

# 2. Compute TMM-normalized CPM (non-log, all >= 0)

dge <- DGEList(counts = bulk_counts_filtered)
dge <- calcNormFactors(dge)

bulk_cpm <- cpm(dge, log = FALSE)  # <-- non-log CPM, all non-negative

cat("Bulk CPM matrix dimensions: ",
    nrow(bulk_cpm), "genes x",
    ncol(bulk_cpm), "cell lines\n")

# 3. Optional: filter low-variance genes

gene_var <- apply(bulk_cpm, 1, var)
bulk_cpm_filt <- bulk_cpm[gene_var > quantile(gene_var, 0.25), ]  # top 75% variable

cat("After variance filter:", nrow(bulk_cpm_filt), "genes remain\n")

# NMF requires a numeric matrix
bulk_nmf_input <- as.matrix(bulk_cpm_filt)

# Sanity: ensure no negative values
cat("Min value in NMF input:", min(bulk_nmf_input), "\n")
if (min(bulk_nmf_input) < 0) {
  stop("NMF input still has negative values, something is wrong.")
}

# 4. Choose NMF rank (k)

set.seed(123)

ranks_to_try <- 4:10

nmf_res <- nmf(
  bulk_nmf_input,
  rank  = ranks_to_try,
  method = "brunet",
  nrun   = 20,
  seed   = 123
)

pdf("Bulk_NMF_rank_diagnostics.pdf", width = 10, height = 4)
plot(nmf_res)
dev.off()

cat("NMF rank diagnostics saved to Bulk_NMF_rank_diagnostics.pdf\n")
cat("Inspect this plot and pick a rank (k) where cophenetic is high\n\n")

# 5. Final NMF run with chosen rank
# You can change k_final after inspecting the diagnostic plot

k_final <- 6  # adjust if needed

set.seed(123)
nmf_k <- nmf(
  bulk_nmf_input,
  rank   = k_final,
  method = "brunet",
  nrun   = 50
)

cat("Final NMF run completed with rank =", k_final, "\n")

W <- basis(nmf_k)  # genes x k
H <- coef(nmf_k)   # k x samples

# 6. Inspect and save W and H

pheatmap(H,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "NMF metaprogram activation (H matrix)")

saveRDS(W, file = "Bulk_metaprograms_W.rds")
saveRDS(H, file = "Bulk_metaprograms_H.rds")
cat("Saved W (genes x k) and H (k x samples).\n")

# 7. Extract top genes for each metaprogram

top_n <- 50
metaprogram_genes <- list()

for (i in 1:k_final) {
  loadings <- W[, i]
  top_genes <- names(sort(loadings, decreasing = TRUE))[1:top_n]
  metaprogram_genes[[paste0("Prog", i)]] <- top_genes
}

saveRDS(metaprogram_genes, file = "Bulk_metaprogram_genes_top50.rds")

max_len <- max(lengths(metaprogram_genes))
mp_df <- do.call(cbind, lapply(metaprogram_genes, function(x) {
  c(x, rep(NA, max_len - length(x)))
}))
write.table(mp_df,
            file = "Bulk_metaprogram_genes_top50.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Saved metaprogram gene sets to Bulk_metaprogram_genes_top50.rds and .tsv\n")

# 8. Check known markers

known_markers <- c("MITF", "PMEL", "TYR",
                   "AXL", "NGFR", "WNT5A",
                   "SOX10", "EDNRB",
                   "MKI67", "CDK1", "TOP2A",
                   "HSPA1A", "HSP90AA1",
                   "TP53", "CDKN1A")

marker_loadings <- W[rownames(W) %in% known_markers, , drop = FALSE]
print(marker_loadings)

cat("\nDone. You can now load Bulk_metaprogram_genes_top50.rds in tme_mal.R\n")
