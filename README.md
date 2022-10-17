# Data analysis of single-cell sequencing of X-linked chronic granulomatous disease
Codebase for X-CGD single-cell expression data analysis

__Rmd/analysis_notebook.Rmd__ contains all code necessary to run the analyses as described in
__"Autoimmunity-related disease signatures shared across cell type dominate the transcriptional landscape of X-linked CGD"__

Code for running components of the analysis are in R/ and are referenced at the appropriate location in __Rmd/analysis_notebook.Rmd__ - 
Global variables and supporting function - globals.R
Integration pipeline - integration.R
Differential expression - diff_exp.R
Co-expression analysis - coexp.R


## Project background

Data sources - 10X single-cell sequencing of PBMCs and isolated monocytes from peripheral blood of probands and carriers with X-linked chronic granulomatous disease (X-CGD), as well as age- and sex-matched controls (14 male probands, 10 female carriers, and 15 (8 male, and 7 female) controls).

## Outline of analysis

1. Load all required data into environment
2. Initial quality control - filter out poor quality cells
3. Integrate all samples across batch using Seurat RPCA pipeline
4. Clustering of integrated output, followed by cluster annotation
5. Differential abundance analysis
6. Differential expression analysis - probands vs controls using DESeq2
7. Geneset enrichment analysis of differentially expressed genes in X-CGD probands vs controls
8. Definition of cluster-specific signatures in X-CGD probands
9. Evaluate expression of cluster-specific signatures in X-CGD carriers (incl. correlation with DHR/NOI scores)
10. Evaluate whether expression of cluster-specific signatures can classify carriers from controls (+ calculate null expectation/p-value)
11. Load in reference CoCoCoNet co-expression networks for human and other species (+ subset and re-rank networks)
12. Evaluate co-expression of cluster-specific disease signatures in reference co-expression networks
13. Evaluate conservation of co-expression using 1-1 orthologs across species in reference co-expression networks
14. Evaluate co-expression partners of CYBB with relationship to the type 1 interferon pathway
