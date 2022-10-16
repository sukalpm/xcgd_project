options(future.globals.maxSize= 53687063712)
`%notin%` <- Negate(`%in%`)

#Global comparison/group lists

ALL_SMPL_TYPE_COMPARISONS <- list(c("Control", "Carrier"), c("Carrier", "Proband"), c("Control", "Proband"))
SMPL_TYPES <- c("Control", "Carrier", "Proband")

DEFAULT_COLORS <- c("steelblue", "springgreen4", "salmon")

#Single-cell quality constants

MIN_UMI <- 500
MIN_GENES <- 250
MIN_LOG10_GENES_PER_UMI <- 0.8
MAX_PCT_MT_RNA <- 20


#Plotting constants

X_AXIS_TITLE_SIZE <- 12
Y_AXIS_TITLE_SIZE <- 12
X_AXIS_TEXT_SIZE <- 12
Y_AXIS_TEXT_SIZE <- 12
PLOT_TITLE_SIZE <- 14
STRIP_TEXT_SIZE <- 12

#Set default ggplot theme - for most plots - some individual plots may have
#modifications within the code itself

theme_set(theme_bw() +  theme(legend.position = "none",
                              axis.title.x = element_text(size = X_AXIS_TITLE_SIZE),
                              axis.title.y = element_text(size = Y_AXIS_TITLE_SIZE),
                              axis.text.x = element_text(size = X_AXIS_TEXT_SIZE, angle = 90),
                              axis.text.y = element_text(size = Y_AXIS_TEXT_SIZE),
                              strip.background = element_blank(),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              strip.text.x = element_text(size = STRIP_TEXT_SIZE, margin = margin(.2, 0, .1, 0, "cm")),
                              strip.text.y = element_text(size = STRIP_TEXT_SIZE, margin = margin(.2, 0, .1, 0, "cm")),
                              plot.title = element_text(size = PLOT_TITLE_SIZE)))


# Define custom function to calculate differential expression in probands vs controls using DESeq2
#' analyze_DE_shrunk_prob_con
#'
#'Perform differential expression analysis given a Seurat object, the slot to pull counts from,
#' and the assay to use
#' @param seurat_obj The Seurat object containing the data for the DE assay
#' @param slot The slot in the Seurat object to be used (default "counts") (optional)
#' @param assay The assay to draw the data from for the DE assay (default "RNA") (optional)
#' @param min_size_gene_set The minimum size of the geneset for enrichment analysis (default = 10)
#' @param max_size_gene_set The maximum size of the geneset for enrichment analysis (default = 101)
#'
#' @return
#' List object with
#' DESeqDataSet ( _dds_ ),
#' vst-transformed object ( _vst_ ),
#' normalized counts matrix ( _counts_ ) ,
#' and DE table ( _prob\_\con_ )
#' @export
#'
#' @examples
#' analyze_DE_shrunk_prob_con(seurat_obj, slot = "counts", assay = "RNA")

analyze_DE_shrunk_prob_con <- function(seurat_obj, cluster, slot = "counts", assay = "RNA"){
  seurat_obj <- subset(seurat_obj, ident = cluster)
  expr <- AggregateExpression(seurat_obj, assays = assay, slot = slot, group.by = c("smpl"))
  counts <- round(x = expr$RNA, digits = 0)
  meta.data <- seurat_obj@meta.data
  rownames(meta.data) <- NULL
  meta.data <- meta.data %>%
    dplyr::select(nUMI, nGene, smpl, src, batch, pMT_RNA, smpl_type, log10GenesPerUMI, nCount_RNA, nFeature_RNA, Sex) %>%
    group_by(smpl) %>%
    mutate(nUMI = mean(nUMI),
           nGene = mean(nGene),
           pMT_RNA = mean(pMT_RNA),
           log10GenesPerUMI = mean(log10GenesPerUMI),
           nCount_RNA = mean(nCount_RNA),
           nFeature_RNA = mean(nFeature_RNA)) %>% distinct %>% column_to_rownames("smpl")

  dds <- DESeqDataSetFromMatrix(countData = counts[, rownames(meta.data)], colData = meta.data, design = ~ smpl_type)

  vsd <- vst(dds)

  dds <- DESeq(dds)

  proband_control <- results(dds, contrast = c("smpl_type", "proband", "control"))

  proband_control <- lfcShrink(dds, contrast = c("smpl_type", "proband", "control"), type = "ashr")

  proband_control <- data.frame(proband_control) %>% rownames_to_column("gene")

  results <-
    list("dds" = dds,
         "vsd" = vsd,
         "counts" = counts(dds, normalized = TRUE),
         "prob_con" = proband_control)

}
