#' Pre-process the expression matrix before generating the shiny app
#' @description This function denoises the expression matrix using the noisyR package
#' and then normalises it. It is recommended to use this function before using
#' \code{\link{generateShinyApp}}.
#' @param expression.matrix the expression matrix; rows correspond to genes and
#' columns correspond to samples
#' @param denoise whether to use noisyR to denoise the expression matrix; 
#' proceeding without denoising data is not recommended
#' @param output.plot whether to create an expression-similarity plot for the
#' noise analysis (printed to the console); default is FALSE
#' @param normalisation.method the normalisation method to be used; default is
#' quantile; any unrecognised input will result in no normalisation being 
#' applied, but proceeding with un-normalised data is not recommended;
#' currently supported normalisation methods are:
#' \describe{
#'   \item{quantile}{Quantile normalisation using the \code{normalize.quantiles}
#'   function from the \code{preprocessCore} package}
#'   \item{rpm}{A version of RPM (reads per million) normalisation, where each
#'   sample is scaled by the median expression in the sample divided by the 
#'   total number of reads in that sample}
#'   \item{tmm}{Trimmed Mean of M values normalisation using the
#'   \code{calcNormFactors} function from the \code{edgeR} package}
#'   \item{deseq2}{Size factor normalisation using the 
#'   \code{estimateSizeFactorsForMatrix} function from the \code{DESeq2} package}
#' }
#' @param ... optional arguments passed on to \code{noisyr::noisyr_counts()}
#' @return The denoised, normalised expression matrix; some rows (genes)
#' may have been removed by noisyR.
#' @export
#' @examples
#' expression.matrix <- as.matrix(read.csv(
#'   system.file("extdata", "expression_matrix.csv", package = "bulkAnalyseR"), 
#'   row.names = 1
#' ))[1:500, 1:4]
#' expression.matrix.preproc <- preprocessExpressionMatrix(expression.matrix)
preprocessExpressionMatrix <- function(
  expression.matrix,
  denoise = TRUE,
  output.plot = FALSE,
  normalisation.method = c("quantile", "rpm", "tmm", "deseq2"),
  ...
){
  if(denoise){
    expression.matrix <- noisyr_counts_with_plot(expression.matrix, output.plot = output.plot, ...)
  }else{
    warning("Denoise was set to FALSE, proceeding without denoising data is not recommended")
  }
  message("Performing ", normalisation.method[1], " normalisation...")
  if(normalisation.method[1] == "quantile"){
    expression.matrix.normalised <- preprocessCore::normalize.quantiles(expression.matrix)
    rownames(expression.matrix.normalised) <- base::rownames(expression.matrix)
    colnames(expression.matrix.normalised) <- base::colnames(expression.matrix)
    expression.matrix <- expression.matrix.normalised
  }else if(normalisation.method[1] == "rpm"){
    csums <- colSums(expression.matrix)
    expression.matrix <- sweep(
      x = expression.matrix, 
      MARGIN = 2, 
      STATS = stats::median(csums) / csums,
      FUN = "*"
    )
  }else if(normalisation.method[1] == "tmm"){
    expression.matrix <- sweep(
      x = expression.matrix, 
      MARGIN = 2, 
      STATS = edgeR::calcNormFactors(expression.matrix),
      FUN = "*"
    )
  }else if(normalisation.method[1] == "deseq2"){
    expression.matrix <- sweep(
      x = expression.matrix, 
      MARGIN = 2, 
      STATS = DESeq2::estimateSizeFactorsForMatrix(expression.matrix),
      FUN = "*"
    )
  }else{
    warning("No valid normalisation method selected, proceeding with un-normalised data is not recommended")
  }
  message("Done!")
  expression.matrix
}

#' Apply a modified noisyR counts pipeline printing a plot
#' @description This function is identical to the noisyr::noisyr_counts
#' function, with the addition of the option to print a line plot of the
#' similarity against expression for all samples.
#' @inheritParams preprocessExpressionMatrix
#' @inherit noisyr::noisyr_counts return params
#' @export
#' @examples
#' expression.matrix <- as.matrix(read.csv(
#'   system.file("extdata", "expression_matrix.csv", package = "bulkAnalyseR"), 
#'   row.names = 1
#' ))[1:500, 1:4]
#' expression.matrix.denoised <- noisyr_counts_with_plot(expression.matrix)
noisyr_counts_with_plot <- function(
  expression.matrix, 
  n.elements.per.window = NULL, 
  optimise.window.length.logical = FALSE, 
  similarity.threshold = 0.25, 
  method.chosen = "Boxplot-IQR", 
  ...,
  output.plot = FALSE
){
  base::message(">>> noisyR counts approach pipeline <<<")
  expression.matrix <- noisyr::cast_matrix_to_numeric(expression.matrix)
  if (optimise.window.length.logical) {
    n.elements.per.window <- noisyr::optimise_window_length(expression.matrix, ...)
  }
  expression.summary <- noisyr::calculate_expression_similarity_counts(expression.matrix, ...)
  
  plotlist <- noisyr::plot_expression_similarity(expression.summary = expression.summary)
  plotdf.line <- tibble::tibble()
  for(i in seq_len(ncol(expression.matrix))){
    lineid <- i * 2 - 1
    plotdf.line <- rbind(
      plotdf.line, 
      dplyr::mutate(plotlist[[lineid]]$data, Sample = colnames(expression.matrix)[i]))
  }
  p <- ggplot2::ggplot(plotdf.line) +
    ggplot2::theme_minimal() + 
    ggplot2::geom_line(ggplot2::aes(x = .data$x, y = .data$y, colour = .data$Sample)) +
    ggplot2::geom_smooth(ggplot2::aes(x = .data$x, y = .data$y, colour = .data$Sample), 
                         method = "loess", formula = .data$y ~ .data$x, span = 0.1) +
    ggplot2::ylim(0, 1) +
    ggplot2::xlab("log2(expression)") +
    ggplot2::ylab("Similarity") +
    ggplot2::geom_hline(yintercept = similarity.threshold, color = "black")
  suppressWarnings(print(p))
  
  if (base::length(similarity.threshold) > 1 | base::length(method.chosen) > 1) {
    base::message("Selecting parameters that minimise the coefficient of variation...")
    stats.table <- noisyr::calculate_noise_threshold_method_statistics(
      expression = expression.summary, 
      similarity.threshold.sequence = similarity.threshold, 
      method.chosen.sequence = method.chosen, 
      ... = ...
    )
    row.min.coef.var <- base::which.min(stats.table$noise.threshold.coefficient.of.variation)
    similarity.threshold <- stats.table$similarity.threshold[row.min.coef.var]
    method.chosen <- paste(stats.table$approach[row.min.coef.var], stats.table$method[row.min.coef.var], sep = "-")
  }
  noise.thresholds <- noisyr::calculate_noise_threshold(
    expression = expression.summary, 
    similarity.threshold = similarity.threshold, method.chosen = method.chosen, 
    ...
  )
  expression.matrix.denoised <- noisyr::remove_noise_from_matrix(
    expression.matrix = expression.matrix, 
    noise.thresholds = noise.thresholds, 
    ... = ...
  )
  base::message(">>> Done! <<<")
  expression.matrix.denoised
}