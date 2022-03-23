#' Generate all files required for an autonomous shiny app
#' @description This function creates an app.R file and all required objects
#' to run the app in .rda format in the target directory. A basic argument 
#' check is performed to avoid input data problems. The app directory
#' is standalone and can be used on another platform, as long as bulkAnalyseR
#' is installed there. It is recommended to run 
#' \code{\link{preprocessExpressionMatrix}} before this function.
#' @param shiny.dir directory to store the shiny app; if a non-empty
#' directory with that name already exists an error is generated
#' @param app.title title to be displayed within the app
#' @param theme shiny theme to be used in the app; default is 'flatly'
#' @param modality name of the modality, or a vector of modalities to be
#' included in the app
#' @param expression.matrix the expression matrix; rows correspond to genes and
#' columns correspond to samples; usually preprocessed by 
#' \code{\link{preprocessExpressionMatrix}}; a list  (of the same length as 
#' modality) can be provided if #' \code{length(modality) > 1}
#' @param metadata a data frame containing metadata for the samples contained
#' in the expression.matrix; must contain at minimum two columns:
#' the first column must contain the column names of the expression.matrix,
#' while the last column is assumed to contain the experimental conditions
#' that will be tested for differential expression; a list  (of the same 
#' length as modality) can be provided if #' \code{length(modality) > 1}
#' @param organism organism name to be passed on to \code{gprofiler2::gost};
#' organism names are constructed by concatenating the first letter of the 
#' name and the family name; default is NA - enrichment is not included
#' to ensure compatibility with datasets that have non-standard gene names; 
#' a vector (of the same length as modality) can be provided if 
#' \code{length(modality) > 1}
#' @param org.db database for annotations to transform ENSEMBL IDs to
#' gene names; a list of bioconductor packaged databases can be found with 
#' \code{BiocManager::available("^org\\.")};
#' default in NA, in which case the row names of the expression matrix are
#' used directly - it is recommended to provide ENSEMBL IDs if the database
#' for your model organism is available; 
#' a vector (of the same length as modality) can be provided if 
#' \code{length(modality) > 1}
#' @param panels.default argument to control which of the default panels
#' will be included in the app; default is all, but the enrichment panel
#' will not appear unless organism is also supplied; note that the 'DE' panel 
#' is required for 'DEplot', 'Enrichment', and GRN; a list  (of the same 
#' length as modality) can be provided if \code{length(modality) > 1}
#' @param cis.integration functionality to integrate extra cis-regulatory information into GRN panel. 
#' Tibble containing names of reference expression matrix, tables of coordinates for elements corresponding to rows 
#' of reference expression matrix (reference.coord), tables of coordinates to compare against 
#' reference.coord (comparison.coord) and names for comparison tables.
#' @param trans.integration functionality to integrate extra trans-regulatory information into GRN panel. 
#' Tibble containing names of reference expression matrix, (reference.expression.matrix), comparison expression matrix (comparison.expression.matrix).
#' Organism database names for each expression matrix and names for each table are also required. 
#' @param custom.integration functionality to integrate custom information related to rows of 
#' reference expression matrix. Tibble containining names of reference expression matrix, tables (comparison.table) 
#' with Reference_ID and Reference_Name (matching ENSEMBL and NAME columns of reference organism database) and 
#' Comparison_ID and Comparison_Name. Names for the reference expression matrix and comparison table (comparison.table.name) are also required.
#' @param panels.extra,data.extra,packages.extra functionality to add new
#' user-created panels to the app to extend functionality or change the default
#' behaviour of existing panels; a data frame of the modality, panel UI and 
#' server names and default parameters should be passed to panels.extra 
#' (see example); the names of any packages required 
#' should be passed to the packages.extra argument; extra data should be a
#' single list and passed to the data.extra argument
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
#' metadata <- data.frame(
#'   srr = colnames(expression.matrix.preproc), 
#'   timepoint = rep(c("0h", "12h", "36h"), each = 2)
#' )
#' app.dir <- generateShinyApp(
#'   shiny.dir = paste0(tempdir(), "/shiny_Yang2019"),
#'   app.title = "Shiny app for the Yang 2019 data",
#'   modality = "RNA",
#'   expression.matrix = expression.matrix.preproc,
#'   metadata = metadata,
#'   organism = "mmusculus",
#'   org.db = "org.Mm.eg.db"
#' )
#' # runApp(app.dir)
#' 
#' # Example of an app with a second copy of the QC panel
#' 
#' app.dir.qc2 <- generateShinyApp(
#'   shiny.dir = paste0(tempdir(), "/shiny_Yang2019_QC2"),
#'   app.title = "Shiny app for the Yang 2019 data",
#'   expression.matrix = expression.matrix.preproc,
#'   metadata = metadata,
#'   organism = "mmusculus",
#'   org.db = "org.Mm.eg.db",
#'   panels.extra = tibble::tibble(
#'     modality = "RNA",
#'     UIfun = "sampleSelectPanelUI", 
#'     UIvars = "'SampleSelect2'", 
#'     serverFun = "sampleSelectPanelServer", 
#'     serverVars = "'SampleSelect2', expression.matrix[[1]], metadata[[1]]"
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
                     "Enrichment", "Patterns", "Cross", "GRN"),
  panels.extra = tibble::tibble(
    name = NULL,
    UIfun = NULL, 
    UIvars = NULL, 
    serverFun = NULL, 
    serverVars = NULL
  ),
  data.extra = list(),
  packages.extra = c(),
  cis.integration = tibble::tibble(
    reference.expression.matrix = NULL,
    reference.org.db = NULL,
    reference.coord = NULL,
    comparison.coord = NULL,
    reference.table.name = NULL,
    comparison.table.name = NULL
  ),
  trans.integration = tibble::tibble(
    reference.expression.matrix = NULL,
    reference.org.db = NULL,
    comparison.expression.matrix = NULL,
    comparison.org.db = NULL,
    reference.table.name = NULL,
    comparison.table.name = NULL
  ),
  custom.integration = tibble::tibble(
    reference.expression.matrix = NULL,
    reference.org.db = NULL,
    comparison.table = NULL,
    reference.table.name = NULL,
    comparison.table.name = NULL
  )
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
  validateIntegrationInputs(
    cis.integration = cis.integration,
    trans.integration = trans.integration,
    custom.integration = custom.integration
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
  metadata <- lapply(metadata, as.data.frame)
  
  generateAppFile(
    shiny.dir = shiny.dir,
    app.title = app.title,
    theme = theme,
    modality = modality,
    organism = organism,
    org.db = org.db,
    panels.default = panels.default,
    panels.extra = panels.extra,
    packages.extra = packages.extra,
    cis.integration = cis.integration,
    trans.integration = trans.integration,
    custom.integration = custom.integration
  )
  generateDataFiles(
    shiny.dir = shiny.dir,
    modality = modality,
    expression.matrix = expression.matrix,
    metadata = metadata,
    data.extra = data.extra
  )
  generateIntegrationDataFiles(
    shiny.dir = shiny.dir,
    cis.integration = cis.integration,
    trans.integration = trans.integration,
    custom.integration = custom.integration
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

validateIntegrationInputs <- function(
  cis.integration = cis.integration,
  trans.integration = trans.integration,
  custom.integration = custom.integration
){
  for (i in seq_len(nrow(cis.integration))){
    if ((cis.integration[i,]$reference.expression.matrix)=='expression.matrix'){
      stop("Reference expression matrix for cis integration cannot be named expression.matrix, this is a reserved name")
    }
    cis.integration.row.reference.coord = get(cis.integration[i,]$reference.coord)
    cis.integration.row.comparison.coord = get(cis.integration[i,]$comparison.coord)
    cis.integration.row.reference.expression.matrix = get(cis.integration[i,]$reference.expression.matrix)
    if((!identical(colnames(cis.integration.row.reference.coord), c("ID","Chrom","Start","Stop","Name"))) | 
       (!identical(colnames(cis.integration.row.comparison.coord), c("ID","Chrom","Start","Stop","Name")))){
      stop("Coordinate tables for cis integration should have 5 columns named ID, Chrom, Start, Stop and Name")
    }
    if(length(intersect(rownames(cis.integration.row.reference.expression.matrix), cis.integration.row.reference.coord$ID)) == 0){
      stop("IDs in the reference coordinate table for cis integration should match row names in reference expression matrix")
    }
    if(length(intersect(cis.integration.row.reference.coord$Chrom, cis.integration.row.comparison.coord$Chrom)) == 0){
      stop("Chromosome names for cis integration should match between reference and comparison coordinate tables")
    }
    if((!is.numeric(cis.integration.row.reference.coord$Start) | (!is.numeric(cis.integration.row.reference.coord$Stop)))){
      stop("Start and stop coordinates for cis integration should be numeric")
    }
    if((!is.numeric(cis.integration.row.comparison.coord$Start) | (!is.numeric(cis.integration.row.comparison.coord$Stop)))){
      stop("Start and stop coordinates for cis integration should be numeric")
    }
    if(length(intersect(rownames(cis.integration.row.reference.expression.matrix),cis.integration.row.comparison.coord$ID)) != 0){
      stop("IDs must be unique to either reference or comparison tables for cis integration")
    }
  }
  
  for (i in seq_len(nrow(trans.integration))){
    if ((trans.integration[i,]$reference.expression.matrix)=='expression.matrix'){
      stop("Reference expression matrix for trans integration cannot be named expression.matrix, this is a reserved name")
    }
    if ((trans.integration[i,]$comparison.expression.matrix)=='expression.matrix'){
      stop("Comparison expression matrix for trans integration cannot be named expression.matrix, this is a reserved name")
    }
    
    trans.integration.row.reference.expression.matrix = get(trans.integration[i,]$reference.expression.matrix)
    trans.integration.row.comparison.expression.matrix = get(trans.integration[i,]$comparison.expression.matrix)
    
    if(!is.matrix(trans.integration.row.reference.expression.matrix)){
      stop("The expression matrix for trans integration must be a matrix")
    }
    if(!is.matrix(trans.integration.row.comparison.expression.matrix)){
      stop("The expression matrix for trans integration must be a matrix")
    }
    if(!identical(colnames(trans.integration.row.reference.expression.matrix),colnames(trans.integration.row.comparison.expression.matrix))){
      stop("The columns of the two expression matrices must be identical for trans integration")
    }
    if(length(intersect(rownames(trans.integration.row.reference.expression.matrix),rownames(trans.integration.row.comparison.expression.matrix)))!=0){
      stop("Row names must be unique to either reference or comparison tables for trans integration")
    }
    if(trans.integration$reference.table.name[i]==trans.integration$comparison.table.name[i]){
      stop("Table names for trans integration must be different")
    }
    }
  
  for (i in seq_len(nrow(custom.integration))){
    
    if ((custom.integration[i,]$reference.expression.matrix)=='expression.matrix'){
      stop("Reference expression matrix for custom integration cannot be named expression.matrix, this is a reserved name")
    }
    
    custom.integration.row.comparison.table = get(custom.integration[i,]$comparison.table)
    custom.integration.row.reference.expression.matrix = get(custom.integration[i,]$reference.expression.matrix)
    
    if(!is.matrix(custom.integration.row.reference.expression.matrix)){
      stop("The expression matrix for custom integration must be a matrix")
    }
                                                             
    if(!identical(colnames(custom.integration.row.comparison.table),c('Reference_ID','Reference_Name','Comparison_ID','Comparison_Name'))){
      stop("The columns of comparison.table for custom integration must be Reference_ID, Reference_Name, Comparison_ID, Comparison_Name")
    }
    if(length(intersect(rownames(custom.integration.row.reference.expression.matrix),custom.integration.row.comparison.table$Reference_ID[i]))==0){
      stop("Reference_ID column for custom integration should match row names from reference expression matrix")
    }
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
  if (length(data.extra) != 0){
    save(data.extra, file = paste0(shiny.dir, "/data_extra.rda"))
  }
}

generateIntegrationDataFiles <- function(
  shiny.dir,
  cis.integration,
  trans.integration,
  custom.integration
){
  if (nrow(cis.integration) > 0){
  cis.integration.data <- list()
  for(i in seq_len(nrow(cis.integration))){
    cis.integration.data[[cis.integration[i,]$reference.expression.matrix]]=get(cis.integration[i,]$reference.expression.matrix)
    cis.integration.data[[cis.integration[i,]$reference.coord]]=get(cis.integration[i,]$reference.coord)
    cis.integration.data[[cis.integration[i,]$comparison.coord]]=get(cis.integration[i,]$comparison.coord)
  }
  save(cis.integration.data,file = paste0(shiny.dir, "/", "cis_integration_data.rda"))
  }
  
  if (nrow(trans.integration) > 0){
    trans.integration.data <- list()
    for(i in seq_len(nrow(trans.integration))){
      trans.integration.data[[trans.integration[i,]$reference.expression.matrix]]=get(trans.integration[i,]$reference.expression.matrix)
      trans.integration.data[[trans.integration[i,]$comparison.expression.matrix]]=get(trans.integration[i,]$comparison.expression.matrix)
    }
    save(trans.integration.data, file = paste0(shiny.dir, "/", "trans_integration_data.rda"))
  }
  
  if (nrow(custom.integration) > 0){
    custom.integration.data <- list()
    for(i in seq_len(nrow(custom.integration))){
      custom.integration.data[[custom.integration[i,]$reference.expression.matrix]]=get(custom.integration[i,]$reference.expression.matrix)
      custom.integration.data[[custom.integration[i,]$comparison.table]]=get(custom.integration[i,]$comparison.table)
    }
    save(custom.integration.data,file = paste0(shiny.dir, "/", "custom_integration.data.rda"))
  }
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
  packages.extra,
  cis.integration,
  trans.integration,
  custom.integration
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
    if(is.na(org.db[[i]])){
      code.source.objects <- c(
        code.source.objects,
        glue::glue("anno[[{i}]] <- data.frame("),
        glue::glue("ENSEMBL = rownames(expression.matrix[[{i}]]),"),
        glue::glue("NAME = rownames(expression.matrix[[{i}]])"),
        ")"
      )
    }else{
      code.source.objects <- c(
        code.source.objects,
        glue::glue("anno[[{i}]] <- AnnotationDbi::select("),
        glue::glue("getExportedValue('{org.db[[i]]}', '{org.db[[i]]}'),"),
        glue::glue("keys = rownames(expression.matrix[[{i}]]),"),
        "keytype = 'ENSEMBL',",
        "columns = 'SYMBOL'",
        ") %>%",
        "dplyr::distinct(ENSEMBL, .keep_all = TRUE) %>%",
        "dplyr::mutate(NAME = ifelse(is.na(SYMBOL), ENSEMBL, SYMBOL))"
      )
    }
  }
  
  if (nrow(cis.integration) > 0){
    code.source.objects <- c(code.source.objects, "")
    code.source.objects <- c(code.source.objects, "anno.cis <- list()")
    for(i in seq_len(nrow(cis.integration))){
      if(cis.integration[i,]$reference.org.db=='NULL'){
        code.source.objects <- c(
          code.source.objects,
          glue::glue("anno.cis[[{i}]] <- data.frame("),
          glue::glue("ENSEMBL = rownames(cis.integration.data[['{cis.integration[i,]$reference.expression.matrix}']]),"),
          glue::glue("SYMBOL = rownames(cis.integration.data[['{cis.integration[i,]$reference.expression.matrix}']]),"),
          glue::glue("NAME = rownames(cis.integration.data[['{cis.integration[i,]$reference.expression.matrix}']])"),
          ")"
        )
      }else{
        code.source.objects <- c(
          code.source.objects,
          glue::glue("anno.cis[[{i}]] <- AnnotationDbi::select("),
          glue::glue("getExportedValue('cis.integration[{i},]$reference.org.db', 'cis.integration[{i},]$reference.org.db'),"),
          glue::glue("keys = rownames(cis.integration.data[['{cis.integration[i,]$reference.expression.matrix}']]),"),
          "keytype = 'ENSEMBL',",
          "columns = 'SYMBOL'",
          ") %>%",
          "dplyr::distinct(ENSEMBL, .keep_all = TRUE) %>%",
          "dplyr::mutate(NAME = ifelse(is.na(SYMBOL), ENSEMBL, SYMBOL))"
        )
      }
    }
  }
  
  if (nrow(trans.integration) > 0){
    code.source.objects <- c(code.source.objects, "")
    code.source.objects <- c(code.source.objects, "anno.trans.reference <- list()")
    for(i in seq_len(nrow(trans.integration))){
      if(trans.integration[i,]$reference.org.db=='NULL'){
        code.source.objects <- c(
          code.source.objects,
          glue::glue("anno.trans.reference[[{i}]] <- data.frame("),
          glue::glue("ENSEMBL = rownames(trans.integration.data[['{trans.integration[i,]$reference.expression.matrix}']]),"),
          glue::glue("SYMBOL = rownames(trans.integration.data[['{trans.integration[i,]$reference.expression.matrix}']]),"),
          glue::glue("NAME = rownames(trans.integration.data[['{trans.integration[i,]$reference.expression.matrix}']])"),
          ")"
        )
      }else{
        code.source.objects <- c(
          code.source.objects,
          glue::glue("anno.trans.reference[[{i}]] <- AnnotationDbi::select("),
          glue::glue("getExportedValue('{trans.integration[i,]$reference.org.db}','{trans.integration[i,]$reference.org.db}'),"),
          glue::glue("keys = rownames(trans.integration.data[['{trans.integration[i,]$reference.expression.matrix}']]),"),
          "keytype = 'ENSEMBL',",
          "columns = 'SYMBOL'",
          ") %>%",
          "dplyr::distinct(ENSEMBL, .keep_all = TRUE) %>%",
          "dplyr::mutate(NAME = ifelse(is.na(SYMBOL), ENSEMBL, SYMBOL))"
        )
      }
    }
    
    code.source.objects <- c(code.source.objects, "anno.trans.comparison <- list()")
    for(i in seq_len(nrow(trans.integration))){
      if(trans.integration[i,]$comparison.org.db=='NULL'){
        code.source.objects <- c(
          code.source.objects,
          glue::glue("anno.trans.comparison[[{i}]] <- data.frame("),
          glue::glue("ENSEMBL = rownames(trans.integration.data[['{trans.integration[i,]$comparison.expression.matrix}']]),"),
          glue::glue("SYMBOL = rownames(trans.integration.data[['{trans.integration[i,]$comparison.expression.matrix}']]),"),
          glue::glue("NAME = rownames(trans.integration.data[['{trans.integration[i,]$comparison.expression.matrix}']])"),
          ")"
        )
      }else{
        code.source.objects <- c(
          code.source.objects,
          glue::glue("anno.trans.comparison[[{i}]] <- AnnotationDbi::select("),
          glue::glue("getExportedValue('trans.integration[{i},]$comparison.org.db', 'trans.integration[{i},]$comparison.org.db'),"),
          glue::glue("keys = rownames(trans.integration.data[['{trans.integration[i,]$comparison.expression.matrix}']]),"),
          "keytype = 'ENSEMBL',",
          "columns = 'SYMBOL'",
          ") %>%",
          "dplyr::distinct(ENSEMBL, .keep_all = TRUE) %>%",
          "dplyr::mutate(NAME = ifelse(is.na(SYMBOL), ENSEMBL, SYMBOL))"
        )
      }
    }
  }
  
  if (nrow(custom.integration) > 0){
    code.source.objects <- c(code.source.objects, "")
    code.source.objects <- c(code.source.objects, "anno.custom <- list()")
    for(i in seq_len(nrow(custom.integration))){
      if(custom.integration[i,]$reference.org.db=='NULL'){
        code.source.objects <- c(
          code.source.objects,
          glue::glue("anno.custom[[{i}]] <- data.frame("),
          glue::glue("ENSEMBL = rownames(custom.integration.data[['{custom.integration[i,]$reference.expression.matrix}']]),"),
          glue::glue("SYMBOL = rownames(custom.integration.data[['{custom.integration[i,]$reference.expression.matrix}']]),"),
          glue::glue("NAME = rownames(custom.integration.data[['{custom.integration[i,]$reference.expression.matrix}']])"),
          ")"
        )
      }else{
        code.source.objects <- c(
          code.source.objects,
          glue::glue("anno.custom[[{i}]] <- AnnotationDbi::select("),
          glue::glue("getExportedValue('{org.db}', '{org.db}'),"),
          glue::glue("keys = rownames(custom.integration.data[['{custom.integration[i,]$reference.expression.matrix}']]),"),
          "keytype = 'ENSEMBL',",
          "columns = 'SYMBOL'",
          ") %>%",
          "dplyr::distinct(ENSEMBL, .keep_all = TRUE) %>%",
          "dplyr::mutate(NAME = ifelse(is.na(SYMBOL), ENSEMBL, SYMBOL))"
        )
      }
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
      "),"
    )
    for(j in seq_len(nrow(panels.extra))){
      code.ui <- c(
        code.ui,
        "tabPanel(",
        glue::glue("title = '{panels.extra$name[j]}',"),
        glue::glue("{panels.extra$UIfun[j]}({panels.extra$UIvars[j]}),"),
        "),"
      )
    }
    code.ui <- c(code.ui, "),")
  }
  for(i in seq_len(nrow(cis.integration))){
    code.ui <- c(code.ui, paste0("GRNCisPanelUI('GRNCis_", 
                                 cis.integration[i,]$reference.table.name,
                                 "_vs_", 
                                 cis.integration[i,]$comparison.table.name, 
                                 "','", 
                                 cis.integration[i,]$reference.table.name, 
                                 "','", cis.integration[i,]$comparison.table.name, 
                                 "'),"))
  }
  for(i in seq_len(nrow(trans.integration))){
    code.ui <- c(code.ui, paste0("GRNTransPanelUI('GRNTrans_", 
                                 trans.integration[i,]$reference.table.name,
                                 "_vs_", trans.integration[i,]$comparison.table.name, 
                                 "','", 
                                 trans.integration[i,]$reference.table.name, 
                                 "','", 
                                 trans.integration[i,]$comparison.table.name, 
                                 "'),"))
  }
  for(i in seq_len(nrow(custom.integration))){
    code.ui <- c(code.ui, paste0("GRNCustomPanelUI('GRNCustom_", 
                                 custom.integration[i,]$reference.table.name,
                                 "_vs_", 
                                 custom.integration[i,]$comparison.table.name, 
                                 "','", 
                                 custom.integration[i,]$reference.table.name, 
                                 "','", 
                                 custom.integration[i,]$comparison.table.name, 
                                 "'),"))
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
    
    for(j in seq_len(nrow(panels.extra))){
      code.server <- c(
        code.server, 
        glue::glue("{panels.extra$serverFun[j]}({panels.extra$serverVars[j]})")
      )
    }
  }
  for(i in seq_len(nrow(cis.integration))){
    code.server <- c(code.server, paste0("GRNCisPanelServer('GRNCis_", 
                                         cis.integration[i,]$reference.table.name,
                                         "_vs_", 
                                         cis.integration[i,]$comparison.table.name, 
                                         "', ", 
                                         "cis.integration.data$", 
                                         cis.integration[i,]$reference.expression.matrix, 
                                         ", anno.cis[[",i,"]], ", 
                                         "cis.integration.data$", 
                                         cis.integration[i,]$reference.coord,
                                         ", ", 
                                         "cis.integration.data$", 
                                         cis.integration[i,]$comparison.coord,
                                         ")"))
  }
  for(i in seq_len(nrow(trans.integration))){
    code.server <- c(code.server, paste0("GRNTransPanelServer('GRNTrans_", 
                                         trans.integration[i,]$reference.table.name,
                                         "_vs_", 
                                         trans.integration[i,]$comparison.table.name, 
                                         "', ", 
                                         "trans.integration.data$", 
                                         trans.integration[i,]$reference.expression.matrix, 
                                         ", anno.trans.reference[[",i,"]], anno.trans.comparison[[",i,"]] ,", 
                                         "trans.integration.data$", 
                                         trans.integration[i,]$comparison.expression.matrix,
                                         ", c('", 
                                         trans.integration[i,]$reference.table.name, 
                                         "','", 
                                         trans.integration[i,]$comparison.table.name, 
                                         "'))"))
  }
  for(i in seq_len(nrow(custom.integration))){
    code.server <- c(code.server, paste0("GRNCustomPanelServer('GRNCustom_", 
                                         custom.integration[i,]$reference.table.name,
                                         "_vs_", 
                                         custom.integration[i,]$comparison.table.name, 
                                         "', ", "custom.integration.data$", 
                                         trans.integration[i,]$reference.expression.matrix, 
                                         ", anno.custom[[",i,"]], ", 
                                         "custom.integration.data$", 
                                         custom.integration[i,]$comparison.table, ")"))
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

# metadata = data.frame(id=c(paste0('control_',1:3),paste0('IDD_',1:3)),rep=rep(1:3,2),type=c(rep(c('control','IDD'),each=3)))
# generateShinyApp('mRNA_miRNA_shiny','Li 2021 Trans Regulatory Example',
#                 modality=c('mRNA','miRNA'),
#                 metadata = metadata,
#                 expression.matrix = list(mrna.expression.matrix.preproc,mirna.expression.matrix.preproc),
#                 org.db = c('org.Hs.eg.db',NA),
#                 organism=c('hsapiens',NA),
#                 trans.integration = tibble::tibble(reference.expression.matrix='mrna.expression.matrix.preproc',
#                                                    reference.org.db='org.Hs.eg.db',
#                                                    comparison.expression.matrix='mirna.expression.matrix.preproc',
#                                                    comparison.org.db='NULL',
#                                                    reference.table.name='mRNA',
#                                                    comparison.table.name='miRNA')
#)
