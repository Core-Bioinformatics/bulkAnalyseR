## ----options, include = FALSE-------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "##>"
)

## ----workflow, echo = FALSE, out.width = "80%"--------------------------------
knitr::include_graphics("figures/workflow.png") 

## ----cran_install, eval = FALSE-----------------------------------------------
#  packages.cran <- c("ggplot2",
#                     "shiny",
#                     "shinythemes",
#                     "gprofiler2",
#                     "stats",
#                     "ggrepel",
#                     "utils",
#                     "RColorBrewer",
#                     "circlize",
#                     "shinyWidgets",
#                     "shinyjqui",
#                     "dplyr",
#                     "magrittr",
#                     "ggforce",
#                     "rlang",
#                     "glue",
#                     "matrixStats",
#                     "noisyr",
#                     "tibble",
#                     "ggnewscale",
#                     "ggrastr",
#                     "visNetwork",
#                     "shinyLP",
#                     "grid",
#                     "DT",
#                     "scales",
#                     "shinyjs",
#                     "tidyr",
#                     "UpSetR")
#  new.packages.cran <- packages.cran[!(packages.cran %in% installed.packages()[, "Package"])]
#  if(length(new.packages.cran))
#    install.packages(new.packages.cran)

## ----bioc_install, eval = FALSE-----------------------------------------------
#  packages.bioc <- c("edgeR",
#                     "DESeq2",
#                     "preprocessCore",
#                     "GENIE3",
#                     "ComplexHeatmap")
#  
#  new.packages.bioc <- packages.bioc[!(packages.bioc %in% installed.packages()[,"Package"])]
#  if(length(new.packages.bioc)){
#    if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#    BiocManager::install(new.packages.bioc)
#  }

## ----github_install, eval = FALSE---------------------------------------------
#  install.packages("bulkAnalyseR")
#  
#  ### OR
#  
#  if (!requireNamespace("devtools", quietly = TRUE))
#    install.packages("devtools")
#  
#  devtools::install_github("Core-Bioinformatics/bulkAnalyseR")

## ----setup--------------------------------------------------------------------
library(bulkAnalyseR)
library(dplyr)
library(ggplot2)

## ----read---------------------------------------------------------------------
counts.in <- system.file("extdata", "expression_matrix.csv", package = "bulkAnalyseR")
exp <- as.matrix(read.csv(counts.in, row.names = 1))
head(exp)

## ----metadata-----------------------------------------------------------------
meta <- data.frame(
  srr = colnames(exp), 
  timepoint = rep(c("0h", "12h", "36h"), each = 2)
)
meta

## ----preprocess,fig.width=7, fig.height=5-------------------------------------
exp.proc <- preprocessExpressionMatrix(exp, output.plot = TRUE)

## ----bioconductor dbs, eval = FALSE-------------------------------------------
#  BiocManager::available("^org\\.")

## ----generate app, eval=FALSE-------------------------------------------------
#  generateShinyApp(
#    shiny.dir = "shiny_Yang2019",
#    app.title = "Shiny app for visualisation of three timepoints from the Yang 2019 data",
#    modality = "RNA",
#    expression.matrix = exp.proc,
#    metadata = meta,
#    organism = "mmusculus",
#    org.db = "org.Mm.eg.db"
#  )

## ----only QC and DE panels, eval = FALSE--------------------------------------
#  generateShinyApp(
#    shiny.dir = "shiny_Yang2019_onlyQC_DE",
#    app.title = "Shiny app for visualisation of three timepoints from the Yang 2019 data",
#    modality = "RNA",
#    expression.matrix = exp.proc,
#    metadata = meta,
#    organism = "mmusculus",
#    org.db = "org.Mm.eg.db",
#    panels.default = c('QC','DE')
#  )

## ----two modalities, eval = FALSE---------------------------------------------
#  generateShinyApp(
#    shiny.dir = "shiny_Yang2019_two_modalities",
#    app.title = "Shiny app for visualisation of three timepoints from the Yang 2019 data",
#    modality = c("RNA 1", "RNA 2"),
#    expression.matrix = exp.proc,
#    metadata = meta,
#    organism = "mmusculus",
#    org.db = "org.Mm.eg.db",
#    panels.default = list(
#      c("Landing", "SampleSelect", "QC", "DE", "DEplot", "DEsummary",
#        "Enrichment", "Patterns", "Cross", "GRN"),
#      c('QC','DE')
#    )
#  )

## ----sample select, echo = FALSE, out.width = "80%"---------------------------
knitr::include_graphics("figures/SampleSelect.png")

## ----JSI, echo = FALSE, out.width = "80%"-------------------------------------
knitr::include_graphics("figures/JSI.png") 

## ----JSI_code, eval = FALSE---------------------------------------------------
#  jaccard_heatmap(
#    expression.matrix = exp.proc,
#    metadata = meta,
#    top.annotation.ids = 2,
#    n.abundant = 500,
#    show.values = FALSE,
#    show.row.column.names = (nrow(meta) <= 20)
#  )

## ----PCA, echo = FALSE, out.width = "80%"-------------------------------------
knitr::include_graphics("figures/PCA.png") 

## ----PCA_code, eval = FALSE---------------------------------------------------
#  plot_pca(
#    expression.matrix = exp.proc,
#    metadata = meta,
#    annotation.id = 2,
#    n.abundant = 500,
#    show.labels = TRUE,
#    show.ellipses = TRUE
#  )

## ----MA QC, echo = FALSE, out.width = "80%"-----------------------------------
knitr::include_graphics("figures/MA_QC.png") 

## ----DE, echo = FALSE, out.width = "80%"--------------------------------------
knitr::include_graphics("figures/DE.png") 

## ----DE_code, eval = FALSE----------------------------------------------------
#  # first create annotation table
#  anno <- AnnotationDbi::select(
#    org.Mm.eg.db::org.Mm.eg.db,
#    keys = rownames(exp.proc),
#    keytype = 'ENSEMBL',
#    columns = 'SYMBOL'
#  ) %>%
#    distinct(ENSEMBL, .keep_all = TRUE) %>%
#    mutate(NAME = ifelse(is.na(SYMBOL), ENSEMBL, SYMBOL))
#  
#  DEtable <- DEanalysis_edger(
#    expression.matrix = exp.proc[, 1:4],
#    condition = meta$timepoint[1:4],
#    var1 = "0h",
#    var2 = "12h",
#    anno = anno
#  )
#  DEtable_deseq <- DEanalysis_deseq2(
#    expression.matrix = exp.proc[, 1:4],
#    condition = meta$timepoint[1:4],
#    var1 = "0h",
#    var2 = "12h",
#    anno = anno
#  )

## ----volcano, echo = FALSE, out.width = "80%"---------------------------------
knitr::include_graphics("figures/VolcanoPlot.png") 

## ----volcano_code, eval = FALSE-----------------------------------------------
#  volcano_plot(
#    genes.de.results = DEtable,
#    pval.threshold = 0.05,
#    lfc.threshold = 1,
#    log10pval.cap = TRUE
#  )

## ----MA_DE, echo = FALSE, out.width = "80%"-----------------------------------
knitr::include_graphics("figures/MAPlot.png") 

## ----MA_DE_code, eval = FALSE-------------------------------------------------
#  ma_plot(
#    genes.de.results = DEtable,
#    pval.threshold = 0.05,
#    lfc.threshold = 1
#  )

## ----heatmap, echo = FALSE, out.width = "80%"---------------------------------
knitr::include_graphics("figures/ExpressionHeatmap.png") 

## ----heatmap_code, eval = FALSE-----------------------------------------------
#  genes_heatmap <- filter(arrange(DEtable, desc(abs(log2FC))), pvalAdj < 0.05)$gene_id
#  expression_heatmap(
#    expression.matrix.subset = exp.proc[genes_heatmap, ],
#    top.annotation.ids = 2,
#    metadata = meta,
#    type = "Z-score",
#    show.column.names = (nrow(meta) <= 20)
#  )

## ----PCA_DE, echo = FALSE, out.width = "80%"----------------------------------
knitr::include_graphics("figures/PCA_DE.png") 

## ----PCA_DE_code, eval = FALSE------------------------------------------------
#  genes_de <- filter(DEtable, pvalAdj < 0.05, abs(log2FC) > 1)$gene_id
#  plot_pca(
#    expression.matrix = exp.proc[genes_de, ],
#    metadata = meta,
#    annotation.id = 2,
#    n.abundant = NULL,
#    show.labels = TRUE,
#    show.ellipses = TRUE
#  )

## ----enrichment, echo = FALSE, out.width = "80%"------------------------------
knitr::include_graphics("figures/Enrichment.png") 

## ----enrichment_code, eval = FALSE--------------------------------------------
#  gostres <- gprofiler2::gost(
#    query = genes_de,
#    organism = "mmusculus",
#    correction_method = 'fdr',
#    custom_bg = DEtable$gene_id,
#    evcodes = TRUE
#  )
#  gostres$result <-  mutate(
#    gostres$result,
#    parents = sapply(.data$parents, toString),
#    intersection_names = sapply(.data$intersection, function(x){
#      ensids <- strsplit(x, split = ",")[[1]]
#      names <- DEtable$gene_name[match(ensids, DEtable$gene_id)]
#      paste(names, collapse = ",")
#    })
#  )

## ----patterns, echo = FALSE, out.width = "30%"--------------------------------
knitr::include_graphics("figures/ExpressionPatterns.png") 

## ----patterns_code, eval = FALSE----------------------------------------------
#  condition <- factor(meta$timepoint)
#  tbl <- calculate_condition_mean_sd_per_gene(exp.proc, condition)
#  patterns <- make_pattern_matrix(tbl, n_sd = 2)[, "pattern"]

## ----patternLine, echo = FALSE, out.width = "80%"-----------------------------
knitr::include_graphics("figures/ExpressionPatternsLinePlot.png") 

## ----patternHeatmap, echo = FALSE, out.width = "80%"--------------------------
knitr::include_graphics("figures/ExpressionPatternsHeatmap.png") 

## ----cross, echo = FALSE, out.width = "80%"-----------------------------------
knitr::include_graphics("figures/CrossPlot.png") 

## ----cross_code, eval = FALSE-------------------------------------------------
#  cross_plot(
#    DEtable1 = DEtable,
#    DEtable2 = DEtable_deseq,
#    DEtable1Subset = filter(DEtable, pvalAdj < 0.05, abs(log2FC) > 1),
#    DEtable2Subset = filter(DEtable_deseq, pvalAdj < 0.05, abs(log2FC) > 1),
#    lfc.threshold = 1
#  )

## ----GRN, echo = FALSE, out.width = "80%"-------------------------------------
knitr::include_graphics("figures/GRN.png") 

## ----GRN_code, eval = FALSE---------------------------------------------------
#  infer_GRN(
#    expression.matrix = exp.proc,
#    metadata = meta,
#    anno = anno,
#    seed = 13,
#    targets = "Trpm1",
#    condition = "timepoint",
#    samples = c("0h", "12h"),
#    inference_method = "GENIE3"
#  )

## ----add extra panel, eval = FALSE--------------------------------------------
#  generateShinyApp(
#    shiny.dir = "shiny_Yang2019_ExtraQC",
#    app.title = "Shiny app for visualisation of three timepoints from the Yang 2019 data - extra QC",
#    modality = "RNA",
#    expression.matrix = exp.proc,
#    metadata = meta,
#    organism = "mmusculus",
#    org.db = "org.Mm.eg.db",
#    panels.extra = tibble::tibble(
#      name = "RNA2",
#      UIfun = "modalityPanelUI",
#      UIvars = "'RNA2', metadata[[1]], NA, 'QC'",
#      serverFun = "modalityPanelServer",
#      serverVars = "'RNA2', expression.matrix[[1]], metadata[[1]], anno[[1]], NA, 'QC'"
#    )
#  )

## ----add extra data, eval = FALSE---------------------------------------------
#  extra.data1 = matrix(rnorm(36),nrow=6)
#  extra.data2 = matrix(rnorm(60),nrow=10)
#  
#  generateShinyApp(
#    shiny.dir = "shiny_Yang2019_ExtraQC",
#    app.title = "Shiny app for visualisation of three timepoints from the Yang 2019 data - extra QC",
#    modality = "RNA",
#    expression.matrix = exp.proc,
#    metadata = meta,
#    organism = "mmusculus",
#    org.db = "org.Mm.eg.db",
#    panels.extra = tibble::tibble(
#      name = "RNA2",
#      UIfun = "modalityPanelUI",
#      UIvars = "'RNA2', metadata[[1]], NA, 'QC'",
#      serverFun = "modalityPanelServer",
#      serverVars = "'RNA2', expression.matrix[[1]], metadata[[1]], anno[[1]], NA, 'QC'"
#    ),
#    data.extra = list(extra.data1, extra.data2),
#    packages.extra = "somePackage",
#  )
#  

