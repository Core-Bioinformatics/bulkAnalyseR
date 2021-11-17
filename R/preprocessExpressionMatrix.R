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