## ----options, include = FALSE-------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "##>"
)

## ----workflow, echo = FALSE, out.width = "80%"--------------------------------
knitr::include_graphics("figures/workflow.png") 

## ----cran_install, eval=FALSE-------------------------------------------------
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
#                     "visNetwork")
#  new.packages.cran <- packages.cran[!(packages.cran %in% installed.packages()[, "Package"])]
#  if(length(new.packages.cran))
#    install.packages(new.packages.cran)

## ----bioc_install, eval = FALSE-----------------------------------------------
#  packages.bioc <- c("edgeR",
#                     "DESeq2",
#                     "preprocessCore",
#                     "AnnotationDbi",
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
#  if (!requireNamespace("devtools", quietly = TRUE))
#    install.packages("devtools")
#  
#  devtools::install_github("Core-Bioinformatics/bulkAnalyseR")

## ----setup--------------------------------------------------------------------
library(bulkAnalyseR)

## ----read---------------------------------------------------------------------
counts.in <- system.file("extdata", "expression_matrix.csv", package = "bulkAnalyseR")
exp <- as.matrix(read.csv(counts.in, row.names = 1))
head(exp)

## ----metadata-----------------------------------------------------------------
meta <- data.frame(
  srr = colnames(exp), 
  timepoint = rep(c("0h", "12h", "36h"), each = 2)
)

## ----convert type, include = FALSE, eval = FALSE------------------------------
#  meta$srr = as.character(meta$srr)
#  meta$timepoint = as.character(meta$timepoint)

## ----preprocess,fig.width=7, fig.height=5-------------------------------------
exp.proc <- preprocessExpressionMatrix(exp, output.plot = TRUE)

## ----bioconductor dbs---------------------------------------------------------
if(requireNamespace("BiocManager", quietly = TRUE))
  BiocManager::available("^org\\.")

## ----generate app, eval=FALSE-------------------------------------------------
#  generateShinyApp(
#    expression.matrix = exp.proc,
#    metadata = meta,
#    shiny.dir = "shiny_Yang2019",
#    app.title = "Shiny app for visualisation of three timepoints from the Yang 2019 data",
#    organism = "mmusculus",
#    org.db = "org.Mm.eg.db"
#  )

## ----only QC and DE panels, eval = FALSE--------------------------------------
#  generateShinyApp(
#    expression.matrix = exp.proc,
#    metadata = meta,
#    shiny.dir = "shiny_Yang2019_onlyQC_DE",
#    app.title = "Shiny app for visualisation of three timepoints from the Yang 2019 data",
#    organism = "mmusculus",
#    org.db = "org.Mm.eg.db",
#    panels.default = c('QC','DE')
#  )

## ----sample select, echo = FALSE, out.width = "80%"---------------------------
knitr::include_graphics("figures/ScreenShot.png") 

## ----JSI, echo = FALSE, out.width = "80%"-------------------------------------
knitr::include_graphics("figures/JSI.png") 

## ----PCA, echo = FALSE, out.width = "80%"-------------------------------------
knitr::include_graphics("figures/PCA.png") 

## ----MA QC, echo = FALSE, out.width = "80%"-----------------------------------
knitr::include_graphics("figures/MA_QC.png") 

## ----DE, echo = FALSE, out.width = "80%"--------------------------------------
knitr::include_graphics("figures/DE.png") 

## ----volcano, echo = FALSE, out.width = "80%"---------------------------------
knitr::include_graphics("figures/VolcanoPlot.png") 

## ----MA DE, echo = FALSE, out.width = "80%"-----------------------------------
knitr::include_graphics("figures/MAPlot.png") 

## ----PCA DE, echo = FALSE, out.width = "80%"----------------------------------
knitr::include_graphics("figures/PCA_DE.png") 

## ----heatmap, echo = FALSE, out.width = "80%"---------------------------------
knitr::include_graphics("figures/ExpressionHeatmap.png") 

## ----enrichment, echo = FALSE, out.width = "80%"------------------------------
knitr::include_graphics("figures/Enrichment.png") 

## ----patterns, echo = FALSE, out.width = "30%"--------------------------------
knitr::include_graphics("figures/ExpressionPatterns.png") 

## ----patternLine, echo = FALSE, out.width = "80%"-----------------------------
knitr::include_graphics("figures/ExpressionPatternsLinePlot.png") 

## ----patternHeatmap, echo = FALSE, out.width = "80%"--------------------------
knitr::include_graphics("figures/ExpressionPatternsHeatmap.png") 

## ----cross, echo = FALSE, out.width = "80%"-----------------------------------
knitr::include_graphics("figures/CrossPlot.png") 

## ----GRN, echo = FALSE, out.width = "80%"-------------------------------------
knitr::include_graphics("figures/GRN.png") 

## ----add extra panel, eval = FALSE--------------------------------------------
#  generateShinyApp(
#    expression.matrix = exp.proc,
#    metadata = meta,
#    shiny.dir = "shiny_Yang2019_ExtraQC",
#    app.title = "Shiny app for visualisation of three timepoints from the Yang 2019 data - extra QC",
#    organism = "mmusculus",
#    org.db = "org.Mm.eg.db",
#    panels.extra = tibble::tibble(
#      UIfun = "QCpanelUI",
#      UIvars = "'QC2', metadata",
#      serverFun = "QCpanelServer",
#      serverVars = "'QC2', expression.matrix, metadata"
#    )
#  )

## ----add extra data, eval = FALSE---------------------------------------------
#  
#  extra.data1 = matrix(rnorm(36),nrow=6)
#  extra.data2 = matrix(rnorm(60),nrow=10)
#  
#  generateShinyApp(
#    expression.matrix = exp.proc,
#    metadata = meta,
#    shiny.dir = "shiny_Yang2019_ExtraData",
#    app.title = "Shiny app for visualisation of three timepoints from the Yang 2019 data - extra QC",
#    organism = "mmusculus",
#    org.db = "org.Mm.eg.db",
#    panels.extra = tibble::tibble(
#      UIfun = "QCpanelUI",
#      UIvars = "'QC2', metadata",
#      serverFun = "QCpanelServer",
#      serverVars = "'QC2', expression.matrix, metadata"
#    ),
#    data.extra = c("extra.data1", "extra.data2"),
#    packages.extra = "somePackage",
#  )
#  

