generateShinyApp <- function(
  expression.matrix,
  metadata,
  shiny.dir = "shiny_bulkAnalyser",
  app.title = "Visualisation of RNA-Seq data",
  organism = "hsapiens",
  org.db = "org.Hs.eg.db",
  theme = "flatly",
  panels.default = c("QC", "DE", "DEplot", "Enrichment"),
  panels.extra = tibble::tibble(
    UIfun = NULL, 
    UIvars = NULL, 
    serverFun = NULL, 
    serverVars = NULL
  ),
  data.extra = c(),
  packages.extra = c()
){
  validateAppInputs(
    shiny.dir = shiny.dir,
    expression.matrix = expression.matrix,
    metadata = metadata
  )
  generateAppFile(
    shiny.dir = shiny.dir,
    app.title = app.title,
    organism = organism,
    org.db = org.db,
    theme = theme,
    panels.default = panels.default,
    panels.extra = panels.extra,
    packages.extra = packages.extra
  )
  generateDataFiles(
    shiny.dir = shiny.dir,
    expression.matrix = expression.matrix,
    metadata = metadata,
    data.extra = data.extra
  )
  message("App created! To launch, run shiny::runApp('", shiny.dir, "')")
  invisible(shiny.dir)
}

validateAppInputs <- function(
  shiny.dir = shiny.dir,
  expression.matrix = expression.matrix,
  metadata = metadata
){
  if(ncol(expression.matrix) != nrow(metadata)){
    stop("Detected different number of columns in expression.matrix to rows in metadata")
  }else if(!identical(colnames(expression.matrix), metadata[[1]])){
    stop("The first column on metadata must correspond to the column names of expression.matrix")
  }else if(ncol(metadata) < 2){
    stop("metadata must be a data frame with at least 2 columns")
  }
  if(!dir.exists(shiny.dir)) dir.create(shiny.dir)
  if(length(dir(shiny.dir,  all.files = TRUE, include.dirs = TRUE, no.. = TRUE)) > 0){
    stop("Please specify a new or empty directory")
  }
}

generateDataFiles <- function(
  shiny.dir,
  expression.matrix,
  metadata,
  data.extra
){
  if(any(c("expression_matrix", "metadata") %in% data.extra)){
    stop("expression_matrix and metadata are reserved names, please rename your extra objects")
  }
  save(expression.matrix, file = paste0(shiny.dir, "/expression_matrix.rda"))
  save(metadata, file = paste0(shiny.dir, "/metadata.rda"))
  lapply(data.extra, function(name){
    object<-get(name)
    save(object, file = paste0(shiny.dir, "/", name, ".rda"))
  })
}

generateAppFile <- function(
  shiny.dir,
  app.title,
  organism,
  org.db,
  theme,
  panels.default,
  panels.extra,
  packages.extra
){
  lines.out <- c()
  
  packages.to.load <- c("shiny", "dplyr", "ggplot2", "bulkAnalyseR", packages.extra)
  code.load.packages <- paste0("library(", packages.to.load, ")")
  lines.out <- c(lines.out, code.load.packages, "")
  
  shiny.dir <- normalizePath(shiny.dir)
  code.source.objects <- c(
    paste0("r.files <- list.files(path = '", shiny.dir, "', pattern = '\\.R$')"),
    "r.files <- setdiff(r.files, 'app.R')",
    "for(fl in r.files) source(fl)",
    "rda.files <- list.files(pattern = '\\.rda$')",
    "for(fl in rda.files) load(fl)",
    "anno <- AnnotationDbi::select(",
    glue::glue("getExportedValue('{org.db}', '{org.db}'),"),
    "keys = rownames(expression.matrix),",
    "keytype = 'ENSEMBL',",
    "columns = 'SYMBOL'",
    ") %>%",
    "dplyr::distinct(ENSEMBL, .keep_all = TRUE) %>%",
    "dplyr::mutate(NAME = ifelse(is.na(SYMBOL), ENSEMBL, SYMBOL))"
  )
  lines.out <- c(lines.out, code.source.objects, "")
  
  code.ui <- c(
    "ui <- navbarPage(",
    glue::glue("'{app.title}',"), 
    glue::glue("theme = shinythemes::shinytheme('{theme}'),"),
    "tabPanel('RNAseq',",
    "tabsetPanel("
  )
  if("QC" %in% panels.default) code.ui <- c(code.ui, "QCpanelUI('QC', metadata),")
  if("DE" %in% panels.default){
    code.ui <- c(code.ui, "DEpanelUI('DE', metadata),")
    if("DEplot" %in% panels.default) code.ui <- c(code.ui, "DEplotPanelUI('DEplot',anno),")
    if("Enrichment" %in% panels.default) code.ui <- c(code.ui, "enrichmentPanelUI('Enrichment'),")
  }
  for(i in seq_len(nrow(panels.extra))){
    code.ui <- c(code.ui, glue::glue("{panels.extra$UIfun[i]}({panels.extra$UIvars[i]}),"))
  }
  code.ui <- c(code.ui, ")),", ")")
  lines.out <- c(lines.out, code.ui, "")
  
  code.server <- c("server <- function(input, output, session){")
  if("QC" %in% panels.default){
    code.server <- c(code.server, "QCpanelServer('QC', expression.matrix, metadata)")
  }
  if("DE" %in% panels.default){
    code.server <- c(
      code.server, 
      glue::glue("DEresults <- DEpanelServer('DE', expression.matrix, metadata, '{org.db}')")
    )
    if("DEplot" %in% panels.default){
      code.server <- c(code.server, "DEplotPanelServer('DEplot', DEresults,anno)")
    }
    if("Enrichment" %in% panels.default){
      code.server <- c(
        code.server, 
        glue::glue("enrichmentPanelServer('Enrichment', DEresults, organism = '{organism}')")
      )
    }
  }
  for(i in seq_len(nrow(panels.extra))){
    code.server <- c(code.server, glue::glue("{panels.extra$serverFun[i]}({panels.extra$serverVars[i]})"))
  }
  code.server <- c(code.server, "}")
  lines.out <- c(lines.out, code.server, "")
  
  lines.out <- c(lines.out, "shinyApp(ui, server)")
  
  lines.out <- gsub("\\\\", "\\\\\\\\", lines.out)
  
  write(lines.out, paste0(shiny.dir, "/app.R"))
  
}

#generateShinyApp(
#expression.matrix = exp.proc,
#metadata = meta,
#shiny.dir = "shiny_Yang20191",
#app.title = "Shiny app for two timepoints from the Yang 2019 data",
#organism = "mmusculus",
#org = "org.Mm.eg.db",
#data.extra = c("ChIPseqdata","ATACseqdata"),
#panels.extra = tibble::tibble(
#  UIfun = c("peaksPanelUI","peaksPanelUI"),
#  UIvars = c("'chip', 'ChIPseq', c('control BCL11A IP' = 'control_11AIP',
#                                      'control CHD8 IP' = 'control_CHD8IP',
#                                      'BCL11A KD BCL11A IP' = '11AKD_11AIP',
#                                      'BCL11A KD CHD8 IP' = '11AKD_CHD8IP',
#                                      'CHD8 KD BCL11A IP' = 'CHD8KD_11AIP',
#                                      'CHD8 KD CHD8 IP' = 'CHD8KD_CHD8AIP')",
#             "'atac', 'ATACseq', c('control' = 'control',
#                                      'BCL11A' = 'BCL11A',
#                                      'CHD8' = 'CHD8')"),
#  serverFun = c("peaksPanelServer","peaksPanelServer"),
#  serverVars = c("'chip', ChIPseqdata","'atac',ATACseqdata")
#)
#)