#__3.Integrate all samples across batch using Seurat RPCA pipeline__


#Integration section - very resource (esp. memory) intensive - for ~400,000 cells, requires ~250GB RAM

#Integration of PBMC datasets - adapted from Seurat RPCA integration vignette

gc()

common.filtered.genes <- intersect(rownames(filtered_subset), rownames(filtered_subset_monos))
common.filtered.genes <- intersect(rownames(filtered_subset), rownames(filtered_subset_monos))

filtered_subset <- subset(filtered_subset, features = common.filtered.genes)
filtered_subset_monos <- subset(filtered_subset_monos, features = common.filtered.genes)

split_seurat <- SplitObject(filtered_subset, split.by = "smpl")

plan("sequential")

split_seurat <- lapply(X = split_seurat, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = split_seurat)
split_seurat <- lapply(X = split_seurat, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

gc()

plan("multiprocess", workers = 20)

anchors <- FindIntegrationAnchors(object.list = split_seurat, anchor.features = features, reduction = "rpca")

gc()

plan("sequential")

pbmc.combined <- IntegrateData(anchorset = anchors)

gc()

plan("multiprocess", workers = 20)
DefaultAssay(pbmc.combined) <- "integrated"
pbmc.combined <- ScaleData(pbmc.combined)
pbmc.combined <- RunPCA(pbmc.combined, seed.use = 42)
ElbowPlot(pbmc.combined)
pbmc.combined <- RunUMAP(pbmc.combined, reduction = "pca", dims = 1:10, seed.use = 42)
pbmc.combined <- FindNeighbors(pbmc.combined, reduction = "pca", dims = 1:10)
pbmc.combined <- FindClusters(pbmc.combined, resolution = seq(0.1, 0.5, 0.1), random.seed = 42)
plan("sequential")

gc()


#Integration of Monocyte datasets - adapted from Seurat RPCA integration vignette

split_seurat_monos <- SplitObject(filtered_subset_monos, split.by = "smpl")

length(rownames(pbmc.combined@assays$RNA@counts))
length(rownames(filtered_subset))
length(rownames(filtered_subset_monos))
common.filtered.genes <- intersect(rownames(filtered_subset), rownames(filtered_subset_monos))
sum(rownames(filtered_subset_monos) != common.filtered.genes)

split_seurat_monos <- lapply(X = split_seurat_monos, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

plan("sequential")

features_monos <- SelectIntegrationFeatures(object.list = split_seurat_monos)
split_seurat_monos <- lapply(X = split_seurat_monos, FUN = function(x) {
  x <- ScaleData(x, features = features_monos, verbose = FALSE)
  x <- RunPCA(x, features = features_monos, verbose = FALSE)
})

gc()

plan("multiprocess", workers = 10)
anchors_monos <- FindIntegrationAnchors(object.list = split_seurat_monos, anchor.features = features_monos, reduction = "rpca")

gc()

plan("sequential")

mono.combined <- IntegrateData(anchorset = anchors_monos)

gc()

plan("multiprocess", workers = 10)
DefaultAssay(mono.combined) <- "integrated"
mono.combined <- ScaleData(mono.combined)
mono.combined <- RunPCA(mono.combined, seed.use = 42)
mono.combined <- RunUMAP(mono.combined, reduction = "pca", dims = 1:10, seed.use = 42)
ElbowPlot(mono.combined)
mono.combined <- FindNeighbors(mono.combined, reduction = "pca", dims = 1:10)
mono.combined <- FindClusters(mono.combined, resolution = seq(0.1, 2.0, 0.1), random.seed = 42)
plan("sequential")
