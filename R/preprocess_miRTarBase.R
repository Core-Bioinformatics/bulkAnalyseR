#' Creates a comparison table for miRTarBase to be used for custom integration
#' @description This function downloads the miRTarBase database for the organism 
#' of choice, filters it according to user-specified values and formats ready for 
#' custom integration in \code{\link{generateShinyApp}}.
#' @param download.dir Directory where miRTarBase database will be downloaded.
#' @param download.method Method for downloading miRTarBase file through download.file,
#' see download.file documentation for options for your operating system.
#' @param mirtarbase.file Path to pre-downloaded miRTarBase file for your organism. 
#' If this is left NULL then the file will be downloaded.
#' @param organism.code Three letter code for the organism of choice. See miRTarBase 
#' website for options. For human, enter 'hsa' and for mouse, 'mmu'.
#' @param org.db database for annotations to transform ENSEMBL IDs to
#' gene names; a list of bioconductor packaged databases can be found with 
#' \code{BiocManager::available("^org\\.")}.
#' @param support.type Subset of entries of the 'Support Type' field in miRTarBase. 
#' Only these values will be kept. To find the options available for your organism 
#' of choice, run the function once with \code{print.support.types = TRUE}.
#' @param validation.method Subset of entries of 'Experiments' field in miRTarBase. 
#' Only these values will be kept. To find the options available for your organism 
#' of choice, run the function once with \code{print.validation.methods = TRUE}. 
#' @param reference Should the reference category be mRNA or miRNA? The reference 
#' category chosen here must match the reference category chosen in 
#' \code{custom.integration} in \code{\link{generateShinyApp}}. Default in mRNA.
#' @param print.support.types,print.validation.methods Should options for 
#' Support Type and Experiments be displayed? Default is FALSE.
#' @return A dataframe with Reference_ID/Name and Comparison_ID/Name columns 
#' which can be supplied to \code{custom.integration} in \code{\link{generateShinyApp}}
#' @export
#' @examples
#' comparison.table <- preprocess_miRTarBase(
#'   mirtarbase.file = system.file("extdata", "mmu_MTI_sub.xls", package = "bulkAnalyseR"),
#'   organism.code = "mmu",
#'   org.db = "org.Mm.eg.db",
#'   support.type = "Functional MTI",
#'   validation.method = "Luciferase reporter assay",
#'   reference = "miRNA")
preprocess_miRTarBase <- function(
  download.dir = '.',
  download.method = 'auto',
  mirtarbase.file = NULL,
  organism.code,
  org.db,
  support.type = c(),
  validation.method = c(),
  reference = c('mRNA','miRNA'),
  print.support.types = FALSE,
  print.validation.methods = FALSE
){
  if (!is.null(mirtarbase.file)) {
    if (!file.exists(mirtarbase.file)) {
      stop('miRTarBase file not found at supplied location, please check the path. Otherwise, download the file directly by setting mirtarbase.file = NULL.')
    }
    mirtarbase <- readxl::read_excel(mirtarbase.file)
  }
  else {
    file.suffix = ifelse(organism.code == 'hsa', 'xlsx', 'xls')
    utils::download.file(
      paste0('https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2019/cache/download/8.0/',
             organism.code, '_MTI.', file.suffix), 
      destfile = paste0(download.dir, '/miRTarBase_', organism.code, '.', file.suffix),
      method = download.method
    )
    mirtarbase <- readxl::read_excel(paste0(download.dir,'/miRTarBase_',organism.code,'.',file.suffix))
  }

  if (print.support.types) {
    message('Available support types:')
    message(paste0(unique(mirtarbase$`Support Type`), collapse = ', '))
  }
  if (length(support.type > 0)) {
    if (length(intersect(support.type,mirtarbase$`Support Type`)) == 0) {
      message('Suggested support types not found in database. Proceeding without filtering.')
    } else {
      mirtarbase <- mirtarbase[mirtarbase$`Support Type` %in% support.type,]
    }
  }
  if (length(validation.method > 0) | print.validation.methods) {
    mirtarbase <- tidyr::separate_rows(mirtarbase, .data$Experiments, convert = TRUE, sep = '//')
    mirtarbase <- tidyr::separate_rows(mirtarbase, .data$Experiments, convert = TRUE, sep = ';')
    mirtarbase$Experiments <- trimws(mirtarbase$Experiments)
    if (print.validation.methods){
      message('Available validation methods:')
      message(paste0(unique(mirtarbase$Experiments),collapse = ', '))
    }
    if (length(intersect(validation.method, mirtarbase$Experiments)) == 0) {
      message('Suggested validation methods not found in database. Proceeding without filtering.')
    } else {
      mirtarbase <- mirtarbase[mirtarbase$Experiments %in% validation.method,]
    }
  }
  suppressMessages(
    anno <- AnnotationDbi::select(
      getExportedValue(org.db, org.db),
      keys = unique(mirtarbase$`Target Gene`),
      keytype = 'SYMBOL',
      columns = 'ENSEMBL'
    ) %>%
      dplyr::distinct(.data$SYMBOL, .keep_all = TRUE) %>%
      dplyr::mutate(NAME = ifelse(is.na(.data$SYMBOL), .data$ENSEMBL, .data$SYMBOL))
  )
  
  if (reference[[1]] == 'mRNA'){
    mirtarbase.comparison.table <- data.frame(
      Reference_ID = anno$ENSEMBL[match(mirtarbase$`Target Gene`, anno$NAME)],
      Reference_Name = mirtarbase$`Target Gene`,
      Comparison_ID = mirtarbase$miRNA,
      Comparison_Name = mirtarbase$miRNA,
      Category = 'miRNA')
    mirtarbase.comparison.table <- unique(mirtarbase.comparison.table)
  } else if (reference[[1]] == 'miRNA'){
    mirtarbase.comparison.table <- data.frame(
      Reference_ID = mirtarbase$miRNA,
      Reference_Name = mirtarbase$miRNA,
      Comparison_ID = anno$ENSEMBL[match(mirtarbase$`Target Gene`, anno$NAME)],
      Comparison_Name = mirtarbase$`Target Gene`,
      Category = 'miRNA')
    mirtarbase.comparison.table = unique(mirtarbase.comparison.table)
  } else { stop('Reference should be one of mRNA or miRNA for miRTarBase integration') }
  return(mirtarbase.comparison.table)
}


