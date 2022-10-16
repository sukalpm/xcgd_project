#Calculate DE for all clusters (parallelized).
#DE tables obtained are then corrected for outliers sometimes found in output
#when a gene is only expressed at a low level in one condition -
#strategy used is to exclude all genes which have an absolute fold
#change of more than mean fold change +/- 40 x s.d. of the fold change.

cluster_DE_results <- list()

x <- mclapply(X = used_cell_types, FUN = function(i){
  cluster_DE_results[[i]] <<- analyze_DE_shrunk_prob_con(seurat_obj = pbmc.combined.no.carriers, cluster = i)
}, mc.cores = 8)

cluster_DE_results <- x
names(cluster_DE_results) <- used_cell_types

cluster_DE_corr <- list()
x <- lapply(X = names(cluster_DE_results), function(i){
  l2fc_table <- cluster_DE_results[[i]]$prob_con
  fcs <- l2fc_table %>% filter(!is.na(padj)) %>% dplyr::select(log2FoldChange) %>% "[["(1)
  mean_fc <- mean(fcs)
  sd_fc <- sd(fcs)
  to_corr <- l2fc_table %>%
    filter(!is.na(padj)) %>%
    filter((abs(log2FoldChange) > mean_fc + (40 * sd_fc)) & abs(log2FoldChange) > 5) %>% dplyr::select(gene) %>% "[["(1)

  fcs <- l2fc_table %>% filter(!is.na(padj) & (gene %notin% to_corr)) %>% dplyr::select(log2FoldChange) %>% "[["(1)

  print(paste("correcting ", length(to_corr), " genes for cluster ", i, sep = ""))

  if (length(to_corr) > 0){
    for(j in to_corr){
      orig_fc <- l2fc_table[l2fc_table$gene == j, "log2FoldChange"]
      corr_fc <- ifelse(orig_fc > 0, max(fcs) + 1e-6, min(fcs) - 1e-6)
      l2fc_table[l2fc_table$gene == j, "log2FoldChange"] <- corr_fc
    }
  }
  cluster_DE_corr[[i]] <<- l2fc_table
})

x <- lapply(names(cluster_DE_results), function(i){
  cluster_DE_results[[i]]$prob_con <<- cluster_DE_corr[[i]]
})

used_cell_types_mono <- as.character(unique(Idents(mono.combined.no.carriers)))

cluster_DE_results_mono <- list()
x <- mclapply(X = used_cell_types_mono, FUN = function(i){
  cluster_DE_results_mono[[i]] <<- analyze_DE_shrunk_prob_con(seurat_obj = mono.combined.no.carriers, cluster = i)
}, mc.cores = 8)

cluster_DE_results_mono <- x
names(cluster_DE_results_mono) <- used_cell_types_mono

#DE for isolated monocytes, followed by similar correction for outliers as above
cluster_DE_corr_mono <- list()
x <- lapply(X = names(cluster_DE_results_mono), function(i){
  l2fc_table <- cluster_DE_results_mono[[i]]$prob_con
  fcs <- l2fc_table %>% filter(!is.na(padj)) %>% dplyr::select(log2FoldChange) %>% "[["(1)
  mean_fc <- mean(fcs)
  sd_fc <- sd(fcs)
  to_corr <- l2fc_table %>%
    filter(!is.na(padj)) %>%
    filter((abs(log2FoldChange) > mean_fc + (40 * sd_fc)) & abs(log2FoldChange) > 5) %>% dplyr::select(gene) %>% "[["(1)

  fcs <- l2fc_table %>% filter(!is.na(padj) & (gene %notin% to_corr)) %>% dplyr::select(log2FoldChange) %>% "[["(1)

  print(paste("correcting ", length(to_corr), " genes for cluster ", i, sep = ""))

  if (length(to_corr) > 0){
    for(j in to_corr){
      orig_fc <- l2fc_table[l2fc_table$gene == j, "log2FoldChange"]
      corr_fc <- ifelse(orig_fc > 0, max(fcs) + 1e-6, min(fcs) - 1e-6)
      l2fc_table[l2fc_table$gene == j, "log2FoldChange"] <- corr_fc
    }
  }
  cluster_DE_corr_mono[[i]] <<- l2fc_table
})

x <- lapply(names(cluster_DE_results_mono), function(i){
  cluster_DE_results_mono[[i]]$prob_con <<- cluster_DE_corr_mono[[i]]
})

#Pool DE tables and used_cell_types across PBMCs and isolated monocytes for ease of plotting

all_DE <- c(cluster_DE_results, cluster_DE_results_mono)
all_used_cell_types <- c(used_cell_types, used_cell_types_mono)
