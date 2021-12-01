#' Pre-process the expression matrix before generating the shiny app
#' @description This function denoises the expression matrix using the noisyR package
#' and then normalises it. It is recommended to use this function before using
#' \code{\link{generateShinyApp}}.
#' @param expression.matrix the expression matrix; rows correspond to genes and
#' columns correspond to samples
#' @param normalisation.method the normalisation method to be used; default is
#' quantile
#' @param ... optional arguments passed on to \code{noisyr::noisyr_counts()}
#' @return The denoised, normalised expression matrix; some rows (genes)
#' may have been removed by noisyR.
#' @export
#' @examples
#' expression.matrix <- as.matrix(read.csv(
#'   system.file("extdata", "expression_matrix.csv", package = "bulkAnalyseR"), 
#'   row.names = 1
#' ))
#' expression.matrix.preproc <- preprocessExpressionMatrix(expression.matrix)
preprocessExpressionMatrix <- function(
  expression.matrix,
  normalisation.method = c("quantile", "rpm"),
  ...
){
  expression.matrix <- noisyr::noisyr_counts(expression.matrix, ...)
  if(normalisation.method[1] == "quantile"){
    message("Performing ", normalisation.method[1], " normalisation...")
    expression.matrix.normalised <- preprocessCore::normalize.quantiles(expression.matrix)
    rownames(expression.matrix.normalised) <- base::rownames(expression.matrix)
    colnames(expression.matrix.normalised) <- base::colnames(expression.matrix)
    message("Done!")
  }else if(normalisation.method[1] == "rpm"){
    message("Performing ", normalisation.method[1], " normalisation...")
    expression.matrix <- expression.matrix / colSums(expression.matrix) * stats::median(colSums(expression.matrix))
    message("Done!")
  }else{
    warning("No valid normalisation method selected, proceeding with un-normalised data is not recommended")
  }
  expression.matrix
}