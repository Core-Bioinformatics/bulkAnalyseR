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
#' name and the family name; default is NA - enrichment is not included
#' to ensure compatibility with datasets that have non-standard gene names
#' @param org.db database for annotations to transform ENSEMBL IDs to
#' gene names; a list of bioconductor packaged databases can be found with 
#' \code{BiocManager::available("^org\\.")};
#' default in NA, in which case the row names of the expression matrix are
#' used directly - it is recommended to provide ENSEMBL IDs if the database
#' for your model organism is available
#' @param theme shiny theme to be used in the app; default is 'flatly'
#' @param panels.default argument to control which of the default panels
#' will be included in the app; default is all, but the enrichment panel
#' will not appear unless organism is also supplied; note that the 'DE' panel 
#' is required for 'DEplot' and 'Enrichment'
#' @param panels.extra,data.extra,packages.extra functionality to add new
#' user-created panels to the app to extend functionality or change the default
#' behaviour of existing panels; a data frame of the panel UI and server names
#' and default parameters should be passed to panels.extra (see example);
#' the names of any extra data and/or packages required should be passed to
#' the data.extra and packages.extra arguments
#' @return The path to shiny.dir (invisibly).
#' @export
#' @import shiny
#' @import ggplot2
#' @importFrom rlang .data
#' @examples
#' expression.matrix.preproc <- as.matrix(read.csv(
#'   system.file("extdata", "expression_matrix_preprocessed.csv", package = "bulkAnalyseR"), 
#'   row.names = 1
#' ))
#' 
#' metadata <- data.frame(
#'   srr = colnames(expression.matrix.preproc), 
#'   timepoint = rep(c("0h", "12h", "36h"), each = 2)
#' )
#' app.dir <- generateShinyApp(
#'   expression.matrix = expression.matrix.preproc,
#'   metadata = metadata,
#'   shiny.dir = paste0(tempdir(), "/shiny_Yang2019"),
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
#'   shiny.dir = paste0(tempdir(), "/shiny_Yang2019_QC2"),
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
#' 
#' # clean up tempdir
#' unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
generateShinyApp <- function(
  shiny.dir = "shiny_bulkAnalyseR",
  app.title = "Visualisation of RNA-Seq data",
  theme = "flatly",
  modality = "RNA",
  expression.matrix,
  metadata,
  organism = NA,
  org.db = NA,
  panels.default = c("Landing", "SampleSelect", "QC", "DE", "DEplot", "DEsummary",
                     "Patterns", "Enrichment", "Cross", "GRN"),
  panels.extra = tibble::tibble(
    modality = NULL,
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
    modality = modality,
    expression.matrix = expression.matrix,
    metadata = metadata,
    organism = organism,
    org.db = org.db,
    panels.default = panels.default,
    data.extra = data.extra
  )
  
  n_modalities <- length(modality)
  if(!is.list(expression.matrix)){
    expression.matrix <- rep(list(expression.matrix), n_modalities)
  }
  if(is.data.frame(metadata)){
    metadata <- rep(list(metadata), n_modalities)
  }
  if(length(organism) == 1){
    organism <- rep(organism, n_modalities)
  }
  if(length(org.db) == 1){
    org.db <- rep(org.db, n_modalities)
  }
  if(!is.list(panels.default)){
    panels.default <- rep(list(panels.default), n_modalities)
  }
  
  generateAppFile(
    shiny.dir = shiny.dir,
    app.title = app.title,
    theme = theme,
    modality = modality,
    organism = organism,
    org.db = org.db,
    panels.default = panels.default,
    panels.extra = panels.extra,
    packages.extra = packages.extra
  )
  
  generateDataFiles(
    shiny.dir = shiny.dir,
    modality = modality,
    expression.matrix = expression.matrix,
    metadata = metadata,
    data.extra = data.extra
  )
  message("App created! To launch, run shiny::runApp('", shiny.dir, "')")
  invisible(shiny.dir)
}

validateAppInputs <- function(
  shiny.dir,
  modality,
  expression.matrix,
  metadata,
  organism,
  org.db,
  panels.default,
  data.extra
){
  if(!dir.exists(shiny.dir)) dir.create(shiny.dir)
  if(length(dir(shiny.dir,  all.files = TRUE, include.dirs = TRUE, no.. = TRUE)) > 0){
    stop("Please specify a new or empty directory")
  }
  
  n_modalities <- length(modality)
  
  validate_matrix_metadata <- function(expression.matrix, metadata){
    if(!is.matrix(expression.matrix)){
      stop("The expression matrix must be a matrix")
    }
    if(ncol(expression.matrix) != nrow(metadata)){
      stop("Detected different number of columns in expression.matrix to rows in metadata")
    }else if(!identical(colnames(expression.matrix), metadata[[1]])){
      stop("The first column on metadata must correspond to the column names of expression.matrix")
    }else if(ncol(metadata) < 2){
      stop("metadata must be a data frame with at least 2 columns")
    }
  }
  if(is.list(expression.matrix) & !is.data.frame(metadata)){
    if(length(expression.matrix) != n_modalities){
      stop("expression.matrix list must have the same length as modality vector")
    }
    if(length(metadata) != n_modalities){
      stop("metadata list must have the same length as modality vector")
    }
    invisible(lapply(X = seq_len(length(expression.matrix)), FUN = function(i){
      validate_matrix_metadata(expression.matrix[[i]], metadata[[i]])
    }))
  }else if(is.list(expression.matrix)){
    if(length(expression.matrix) != n_modalities){
      stop("expression.matrix list must have the same length as modality vector")
    }
    invisible(lapply(X = expression.matrix, FUN = function(exp){
      validate_matrix_metadata(exp, metadata)
    }))
  }else if(!is.data.frame(metadata)){
    if(length(metadata) != n_modalities){
      stop("metadata list must have the same length as modality vector")
    }
    invisible(lapply(X = metadata, FUN = function(meta){
      validate_matrix_metadata(expression.matrix, meta)
    }))
  }else{
    validate_matrix_metadata(expression.matrix, metadata)
  }
  
  if(length(organism) != 1 & length(organism) != n_modalities){
    stop("organism must be length 1 or have the same length as modality vector")
  }
  if(length(org.db) != 1 & length(org.db) != n_modalities){
    stop("org.db must be length 1 or have the same length as modality vector")
  }
  if(is.list(panels.default)){
    if(length(panels.default) != n_modalities){
      stop("panels.default list must have the same length as modality vector")
    }
  }
  
  if(any(c("expression_matrix", "metadata") %in% data.extra)){
    stop("expression_matrix and metadata are reserved names, please rename your extra objects")
  }
}

generateDataFiles <- function(
  shiny.dir,
  modality,
  expression.matrix,
  metadata,
  data.extra
){
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
  theme,
  modality,
  organism,
  org.db,
  panels.default,
  panels.extra,
  packages.extra
){
  lines.out <- c()
  
  packages.to.load <- c("shiny", "dplyr", "ggplot2", "bulkAnalyseR", packages.extra)
  code.load.packages <- paste0("library(", packages.to.load, ")")
  lines.out <- c(lines.out, code.load.packages, "")
  
  code.source.objects <- c(
    paste0("r.files <- list.files(path = getwd(), pattern = '\\.R$')"),
    "r.files <- setdiff(r.files, 'app.R')",
    "for(fl in r.files) source(fl)",
    "rda.files <- list.files(pattern = '\\.rda$')",
    "for(fl in rda.files) load(fl)"
  )
  
  code.source.objects <- c(code.source.objects, "anno <- list()")
  for(i in seq_len(length(org.db))){
    if(is.na(org.db[i])){
      code.source.objects <- c(
        code.source.objects,
        glue::glue("anno[[{i}]] <- data.frame("),
        "ENSEMBL = rownames(expression.matrix),",
        "NAME = rownames(expression.matrix)",
        ")"
      )
    }else{
      code.source.objects <- c(
        code.source.objects,
        glue::glue("anno[[{i}]] <- AnnotationDbi::select("),
        glue::glue("getExportedValue('{org.db}', '{org.db}'),"),
        glue::glue("keys = rownames(expression.matrix[[{i}]]),"),
        "keytype = 'ENSEMBL',",
        "columns = 'SYMBOL'",
        ") %>%",
        "dplyr::distinct(ENSEMBL, .keep_all = TRUE) %>%",
        "dplyr::mutate(NAME = ifelse(is.na(SYMBOL), ENSEMBL, SYMBOL))"
      )
    }
  }
  lines.out <- c(lines.out, code.source.objects, "")
  
  code.ui <- c(
    "ui <- navbarPage(",
    glue::glue("'{app.title}',"), 
    glue::glue("theme = shinythemes::shinytheme('{theme}'),"),
    "header = tags$head(tags$style('body {overflow-y: scroll;}')),"
  )
  for(i in seq_len(length(modality))){
    panels.default.string <- paste0("c('", paste(panels.default[[i]], collapse = "', '"), "')")
    code.ui <- c(
      code.ui, 
      "tabPanel(",
      glue::glue("title = '{modality[i]}',"),
      "modalityPanelUI(",
      glue::glue("id = '{modality[i]}',"),
      glue::glue("metadata = metadata[[{i}]],"),
      glue::glue("organism = '{organism[i]}',"),
      glue::glue("panels.default = {panels.default.string}"),
      ")"
    )
    panels.extra.subset <- dplyr::filter(panels.extra, modality == modality[i])
    for(j in seq_len(nrow(panels.extra.subset))){
      code.ui <- c(
        code.ui,
        "tabsetPanel(",
        glue::glue("{panels.extra.subset$UIfun[j]}({panels.extra.subset$UIvars[j]}),"),
        ")"
      )
    }
    code.ui <- c(code.ui, "),")
  }
  code.ui <- c(code.ui, ")")
  lines.out <- c(lines.out, code.ui, "")
  
  code.server <- c("server <- function(input, output, session){")
  
  for(i in seq_len(length(modality))){
    panels.default.string <- paste0("c('", paste(panels.default[[i]], collapse = "', '"), "')")
    code.server <- c(
      code.server,
      "modalityPanelServer(",
      glue::glue("id = '{modality[i]}',"), 
      glue::glue("expression.matrix = expression.matrix[[{i}]],"),
      glue::glue("metadata = metadata[[{i}]],"),
      glue::glue("anno = anno[[{i}]],"),
      glue::glue("organism = '{organism[i]}',"), 
      glue::glue("panels.default = {panels.default.string}"),
      ")"
    )
    panels.extra.subset <- dplyr::filter(panels.extra, modality == modality[i])
    for(j in seq_len(nrow(panels.extra.subset))){
      code.server <- c(
        code.server, 
        glue::glue("{panels.extra.subset$serverFun[j]}({panels.extra.subset$serverVars[j]})")
      )
    }
  }
  
  code.server <- c(code.server, "}")
  lines.out <- c(lines.out, code.server, "")
  
  lines.out <- c(lines.out, "shinyApp(ui, server)")
  
  lines.out <- gsub("\\\\", "\\\\\\\\", lines.out)
  
  shiny.dir <- normalizePath(shiny.dir)
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