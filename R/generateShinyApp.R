#' Generate all files required for an autonomous shiny app
#' @description This function creates an app.R file and all required objects
#' to run the app in .rda format in the target directory. A basic argument 
#' check is performed to avoid input data problems. The app directory
#' is standalone and can be used on another platform, as long as bulkAnalyseR
#' is installed there. It is recommended to run 
#' \code{\link{preprocessExpressionMatrix}} before this function.
#' @param expression.matrix the expression matrix; rows correspond to genes and
#' columns correspond to samples; usually preprocessed by 
#' \code{\link{preprocessExpressionMatrix}}
#' @param metadata a data frame containing metadata for the samples contained
#' in the expression.matrix; must contain at minimum two columns:
#' the first column must contain the column names of the expression.matrix,
#' while the last column is assumed to contain the experimental conditions
#' that will be tested for differential expression
#' @param shiny.dir directory to store the shiny app; if a non-empty
#' directory with that name already exists an error is generated
#' @param app.title title to be displayed within the app
#' @param organism organism name to be passed on to \code{gprofiler2::gost};
#' organism names are constructed by concatenating the first letter of the 
#' name and the family name; default is human - 'hsapiens'
#' @param org.db database for annotations to transform ENSEMBL IDs to
#' gene names; a list of bioconductor packaged databases can be found with 
#' \code{BiocManager::available("^org\\.")}; default is human - 'org.Hs.eg.db'
#' @param theme shiny theme to be used in the app; default is 'flatly'
#' @param panels.default argument to control which of the default panels
#' will be included in the app; default is all; not that the 'DE' panel is
#' required for 'DEplot' and 'Enrichment'
#' @param panels.extra,data.extra,packages.extra functionality to add new
#' user-created panels to the app to extend functionality or change the default
#' behaviour of existing panels; a data frame of the panel UI and server names
#' and default parameters should be passed to panels.extra (see example);
#' the names of any extra data and/or packages required should be passed to
#' the data.extra and packages.extra arguments
#' @return The path to shiny.dir (invisibly).
#' @export
#' @examples
#' expression.matrix <- as.matrix(read.csv(
#'   system.file("extdata", "expression_matrix.csv", package = "bulkAnalyseR"), 
#'   row.names = 1
#' ))
#' expression.matrix.preproc <- preprocessExpressionMatrix(expression.matrix)
#' metadata <- data.frame(
#'   srr = colnames(expression.matrix.preproc), 
#'   timepoint = rep(c("0h", "12h", "36h"), each = 2)
#' )
#' app.dir <- generateShinyApp(
#'   expression.matrix = expression.matrix.preproc,
#'   metadata = metadata,
#'   shiny.dir = paste0(tempdir(), "ssshiny_Yang2019"),
#'   app.title = "Shiny app for the Yang 2019 data",
#'   organism = "mmusculus",
#'   org.db = "org.Mm.eg.db"
#' )
#' # runApp(app.dir)
#' 
#' # Example of an app with a second copy of the QC panel
#' 
#' app.dir.qc2 <- generateShinyApp(
#'   expression.matrix = expression.matrix.preproc,
#'   metadata = metadata,
#'   shiny.dir = paste0(tempdir(), "shiny_Yang2019_QC2"),
#'   app.title = "Shiny app for the Yang 2019 data",
#'   organism = "mmusculus",
#'   org.db = "org.Mm.eg.db",
#'   panels.extra = tibble::tibble(
#'     UIfun = "QCpanelUI", 
#'     UIvars = "'QC2', metadata", 
#'     serverFun = "QCpanelServer", 
#'     serverVars = "'QC2', expression.matrix, metadata"
#'   )
#' )
#' # runApp(app.dir.qc2)
generateShinyApp <- function(
  expression.matrix,
  metadata,
  shiny.dir = "shiny_bulkAnalyseR",
  app.title = "Visualisation of RNA-Seq data",
  organism = "hsapiens",
  org.db = "org.Hs.eg.db",
  theme = "flatly",
  panels.default = c("SampleSelect", "QC", "DE", "DEplot", "DEsummary",
                     "Patterns", "Enrichment", "Cross", "GRN"),
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
    paste0("r.files <- list.files(path = getwd(), pattern = '\\.R$')"),
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
  if("SampleSelect" %in% panels.default) code.ui <- c(code.ui, "sampleSelectPanelUI('SampleSelect'),")
  if("QC" %in% panels.default) code.ui <- c(code.ui, "QCpanelUI('QC', metadata),")
  if("DE" %in% panels.default){
    code.ui <- c(code.ui, "DEpanelUI('DE', metadata),")
    if("DEplot" %in% panels.default) code.ui <- c(code.ui, "DEplotPanelUI('DEplot'),")
    if("DEsummary" %in% panels.default) code.ui <- c(code.ui, "DEsummaryPanelUI('DEsummary', metadata),")
    if("Enrichment" %in% panels.default) code.ui <- c(code.ui, "enrichmentPanelUI('Enrichment'),")
  }
  if("Patterns" %in% panels.default) code.ui <- c(code.ui, "patternPanelUI('Patterns', metadata),")
  if("Cross" %in% panels.default) code.ui <- c(code.ui, "crossPanelUI('Cross', metadata),")
  if("GRN" %in% panels.default) code.ui <- c(code.ui, "GRNpanelUI('GRN'),")
  for(i in seq_len(nrow(panels.extra))){
    code.ui <- c(code.ui, glue::glue("{panels.extra$UIfun[i]}({panels.extra$UIvars[i]}),"))
  }
  code.ui <- c(code.ui, ")),", ")")
  lines.out <- c(lines.out, code.ui, "")
  
  code.server <- c("server <- function(input, output, session){")
  if("SampleSelect" %in% panels.default){
    code.server <- c(
      code.server,
      "filteredInputs <- sampleSelectPanelServer('SampleSelect', expression.matrix, metadata)",
      "expression.matrix <- reactive(filteredInputs()[['expression.matrix']])",
      "metadata <- reactive(filteredInputs()[['metadata']])"
    )
  }else{
    code.server <- c(
      code.server,
      "expression.matrix <- reactiveVal(expression.matrix)",
      "metadata <- reactiveVal(metadata)"
    )
  }
  if("QC" %in% panels.default){
    code.server <- c(code.server, "QCpanelServer('QC', expression.matrix, metadata, anno)")
  }
  if("DE" %in% panels.default){
    code.server <- c(code.server, "DEresults <- DEpanelServer('DE', expression.matrix, metadata, anno)")
    if("DEplot" %in% panels.default){
      code.server <- c(code.server, "DEplotPanelServer('DEplot', DEresults, anno)")
    }
    if("DEsummary" %in% panels.default){
      code.server <- c(code.server, "DEsummaryPanelServer('DEsummary', expression.matrix, metadata, DEresults, anno)")
    }
    if("Enrichment" %in% panels.default){
      code.server <- c(
        code.server, 
        glue::glue("enrichmentPanelServer('Enrichment', DEresults, organism = '{organism}')")
      )
    }
    if("Cross" %in% panels.default){
      code.server <- c(code.server, "crossPanelServer('Cross', expression.matrix, metadata, anno)")
    }
    if("Patterns" %in% panels.default){
      code.server <- c(code.server, "patternPanelServer('Patterns', expression.matrix, metadata, anno)")
    }
    if("GRN" %in% panels.default){
      code.server <- c(code.server, "GRNpanelServer('GRN', expression.matrix, anno)")
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