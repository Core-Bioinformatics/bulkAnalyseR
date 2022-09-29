## ----options, include = FALSE-------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "##>"
)
Sys.setenv("VROOM_CONNECTION_SIZE" = 1e6)

## ----cran_install, eval = FALSE-----------------------------------------------
#  packages.cran <- c(
#    "ggplot2", "shiny", "shinythemes", "gprofiler2", "stats", "ggrepel",
#    "utils", "RColorBrewer", "circlize", "shinyWidgets", "shinyjqui",
#    "dplyr", "magrittr", "ggforce", "rlang", "glue", "matrixStats",
#    "noisyr", "tibble", "ggnewscale", "ggrastr", "visNetwork", "shinyLP",
#    "grid", "DT", "scales", "shinyjs", "tidyr", "UpSetR", "ggVennDiagram"
#  )
#  new.packages.cran <- packages.cran[!(packages.cran %in% installed.packages()[, "Package"])]
#  if(length(new.packages.cran))
#    install.packages(new.packages.cran)
#  
#  packages.bioc <- c(
#    "edgeR", "DESeq2", "preprocessCore", "GENIE3", "ComplexHeatmap"
#  )
#  new.packages.bioc <- packages.bioc[!(packages.bioc %in% installed.packages()[,"Package"])]
#  if(length(new.packages.bioc)){
#    if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#    BiocManager::install(new.packages.bioc)
#  }
#  
#  install.packages("bulkAnalyseR")

## ----read---------------------------------------------------------------------
download_path <- paste0(tempdir(), "expression_matrix.csv.gz")
download.file(
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE178620&format=file&file=GSE178620%5Fraw%5Fabundances%2Ecsv%2Egz", 
  download_path
)
exp <- as.matrix(read.csv(download_path, row.names = 1))[, c(1,2,19,20)]
head(exp)

## ----clean up, include = FALSE------------------------------------------------
file.remove(download_path)

## ----meta---------------------------------------------------------------------
meta <- data.frame(
  name = colnames(exp),
  condition = sapply(colnames(exp), USE.NAMES = FALSE, function(nm){
    strsplit(nm, "_")[[1]][1]
  })
)
meta

## ----preprocess,fig.width=7, fig.height=5-------------------------------------
exp.proc <- bulkAnalyseR::preprocessExpressionMatrix(exp, output.plot = TRUE)

## ----generate app, eval=FALSE-------------------------------------------------
#  bulkAnalyseR::generateShinyApp(
#    shiny.dir = "shiny_GEO",
#    app.title = "Shiny app for visualisation of GEO data",
#    modality = "RNA",
#    expression.matrix = exp.proc,
#    metadata = meta,
#    organism = "hsapiens",
#    org.db = "org.Hs.eg.db"
#  )

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

