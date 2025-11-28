# Bulk RNA-Seq Data Cleaning 
clean_bulk_melanoma <- function(file_path = "Bulk_rna_raw_counts.tsv") {
  
  cat("╔═══════════════════════════════════════════════════════════════════╗\n")
  cat("║     CLEANING BULK RNA-SEQ DATA (MELANOMA CELL LINES)             ║\n")
  cat("╚═══════════════════════════════════════════════════════════════════╝\n\n")
  
  # 1. LOAD DATA
  cat("Loading data...\n")
  bulk_counts <- read.delim(file_path, row.names = 1, check.names = FALSE)
  
  cat("\n ORIGINAL DATA:\n")
  cat("  • Genes:", nrow(bulk_counts), "\n")
  cat("  • Cell lines:", ncol(bulk_counts), "\n")
  cat("  • Sample names:", paste(head(colnames(bulk_counts), 3), collapse=", "), "...\n\n")
  
  # 2. REMOVE ALL-ZERO GENES
  cat("Step 1: Removing genes with zero counts across all samples...\n")
  non_zero_genes <- rowSums(bulk_counts) > 0
  bulk_counts <- bulk_counts[non_zero_genes, ]
  cat("  ✓ Genes remaining:", nrow(bulk_counts), "\n\n")
  
  # 3. FILTER LOW-EXPRESSION GENES
  cat("Step 2: Filtering lowly expressed genes...\n")
  cat("  • Keeping genes expressed (count > 0) in at least 3 samples\n")
  expressed_samples <- rowSums(bulk_counts > 0)
  keep_genes <- expressed_samples >= 3
  bulk_counts <- bulk_counts[keep_genes, ]
  cat("  ✓ Genes remaining:", nrow(bulk_counts), "\n\n")
  
  # 4. CPM FILTERING
  cat("Step 3: CPM (Counts Per Million) filtering...\n")
  cat("  • Keeping genes with CPM > 1 in at least 3 samples\n")
  cpm_matrix <- edgeR::cpm(bulk_counts)
  keep_cpm <- rowSums(cpm_matrix > 1) >= 3
  bulk_counts_filtered <- bulk_counts[keep_cpm, ]
  cat("  ✓ Genes remaining:", nrow(bulk_counts_filtered), "\n\n")
  
  # 5. CHECK FOR DUPLICATES
  cat("Step 4: Checking for duplicate gene names...\n")
  if(any(duplicated(rownames(bulk_counts_filtered)))) {
    n_dups <- sum(duplicated(rownames(bulk_counts_filtered)))
    cat("  ⚠ WARNING:", n_dups, "duplicate gene names found\n")
    cat("  • Keeping first occurrence of each duplicate\n")
    bulk_counts_filtered <- bulk_counts_filtered[!duplicated(rownames(bulk_counts_filtered)), ]
  } else {
    cat("  ✓ No duplicates found\n")
  }
  cat("\n")
  
  # 6. LIBRARY SIZE QC
  cat("Step 5: Library size quality control...\n")
  lib_sizes <- colSums(bulk_counts_filtered)
  lib_sizes_millions <- lib_sizes / 1e6
  
  cat("\n  Library sizes (millions of reads):\n")
  cat("    • Min:", round(min(lib_sizes_millions), 2), "\n")
  cat("    • Median:", round(median(lib_sizes_millions), 2), "\n")
  cat("    • Max:", round(max(lib_sizes_millions), 2), "\n")
  cat("    • Mean:", round(mean(lib_sizes_millions), 2), "\n")
  
  # Identify outliers (>2x or <0.5x median)
  median_lib <- median(lib_sizes)
  outlier_samples <- names(lib_sizes[lib_sizes < median_lib * 0.5 | 
                                       lib_sizes > median_lib * 2])
  
  if(length(outlier_samples) > 0) {
    cat("\n  ⚠ POTENTIAL OUTLIER SAMPLES (>2x or <0.5x median):\n")
    for(sample in outlier_samples) {
      cat("    •", sample, ":", round(lib_sizes[sample]/1e6, 2), "M reads\n")
    }
    cat("\n  Note: Review these samples. Consider removing if they're technical outliers.\n")
  } else {
    cat("\n  ✓ No outlier samples detected\n")
  }
  cat("\n")
  
  # 7. CREATE QC PLOTS
  cat("Step 6: Creating QC plots...\n")
  pdf("Bulk_RNA_QC_plots.pdf", width = 12, height = 8)
  
  par(mfrow = c(2, 2))
  barplot(lib_sizes_millions, 
          main = "Library Sizes per Sample",
          ylab = "Millions of reads",
          las = 2, cex.names = 0.7,
          col = "steelblue")
  abline(h = median(lib_sizes_millions), col = "red", lty = 2, lwd = 2)
  
  boxplot(log10(bulk_counts_filtered + 1),
          main = "Distribution of Gene Counts (log10)",
          ylab = "log10(counts + 1)",
          las = 2, cex.axis = 0.7,
          col = "lightblue")
  
  genes_detected <- colSums(bulk_counts_filtered > 0)
  barplot(genes_detected,
          main = "Number of Genes Detected per Sample",
          ylab = "Number of genes",
          las = 2, cex.names = 0.7,
          col = "coral")
  
  cpm_for_pca <- log2(edgeR::cpm(bulk_counts_filtered) + 1)
  pca_result <- prcomp(t(cpm_for_pca), scale. = TRUE)
  plot(pca_result$x[,1], pca_result$x[,2],
       main = "PCA of Samples",
       xlab = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
       ylab = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)"),
       pch = 19, col = "darkblue", cex = 1.5)
  text(pca_result$x[,1], pca_result$x[,2], 
       labels = rownames(pca_result$x), 
       pos = 3, cex = 0.6)
  
  dev.off()
  cat("  ✓ QC plots saved to: Bulk_RNA_QC_plots.pdf\n\n")
  
  # 8. SAVE CLEANED DATA
  cat("Step 7: Saving cleaned data...\n")
  write.table(bulk_counts_filtered, 
              "Bulk_rna_raw_counts_CLEANED.tsv",
              sep = "\t", 
              quote = FALSE,
              col.names = NA)
  cat("  ✓ Cleaned data saved to: Bulk_rna_raw_counts_CLEANED.tsv\n\n")

  

  # 9. Normalize cleaned bulk RNA (TMM + logCPM)
  
  cat("Step 8: Normalizing cleaned bulk RNA-seq using TMM (edgeR)...\n")
  library(edgeR)
  
  # Create DGEList object
  dge <- DGEList(counts = bulk_counts_filtered)
  
  # Compute TMM normalization factors
  dge <- calcNormFactors(dge)
  
  # Compute log2 CPM values (prior.count=1 avoids log(0))
  bulk_logCPM <- cpm(dge, log = TRUE, prior.count = 1)
  
  # Save normalized matrix
  write.table(bulk_logCPM,
              file = "Bulk_rna_logCPM_normalized.tsv",
              sep = "\t",
              quote = FALSE,
              col.names = NA)
  
  cat("  ✓ Normalized logCPM matrix saved to: Bulk_rna_logCPM_normalized.tsv\n\n")
  return(bulk_counts_filtered)
}

