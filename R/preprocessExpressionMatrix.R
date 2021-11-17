preprocessExpressionMatrix <- function(
  expression.matrix,
  normalisation.method = c("quantile", "rpm"),
  ...
){
  
  expression.matrix <- noisyr::noisyr_counts(expression.matrix, ...)
  
  if(normalisation.method[1] == "quantile"){
    expression.matrix.normalised <- preprocessCore::normalize.quantiles(expression.matrix)
    rownames(expression.matrix.normalised) <- base::rownames(expression.matrix)
    colnames(expression.matrix.normalised) <- base::colnames(expression.matrix)
  }else if(normalisation.method[1] == "rpm"){
    expression.matrix <- expression.matrix / colSums(expression.matrix) * stats::median(colSums(expression.matrix))
  }
  
  expression.matrix

}