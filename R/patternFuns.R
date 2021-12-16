#' Calculate statistics for each gene of an expression matrix given a grouping
#' @description This function calculates the mean and standard deviation of
#' the expression of each gene in an expression matrix, grouped by the conditions
#' supplied.
#' @inheritParams generateShinyApp
#' @param condition the condition to group the columns of the expression matrix
#' by; must be a factor of the same length as ncol(expression.matrix)
#' @return A tibble in long format, with the mean and standard deviation of
#' each gene in each condition. The standard deviation is increased to the 
#' minimum value in the expression matrix (the noise threshold) if it is lower,
#' in order to avoid sensitivity to small changes.
#' @export
#' @examples
#' expression.matrix <- as.matrix(read.csv(
#'   system.file("extdata", "expression_matrix.csv", package = "bulkAnalyseR"), 
#'   row.names = 1
#' ))
#' expression.matrix.preproc <- preprocessExpressionMatrix(expression.matrix)
#' condition <- factor(rep(c("0h", "12h", "36h"), each = 2))
#' tbl <- calculate_condition_mean_sd_per_gene(expression.matrix.preproc[1:10, ], condition)
#' tbl
calculate_condition_mean_sd_per_gene <- function(expression.matrix, condition){
  if(length(unique(condition)) < 2) stop("Please select at least two states")
  
  if(!is.factor(condition) | length(condition) != ncol(expression.matrix)){
    stop("condition must be a factor of the same length as ncol(expression.matrix)")
  }
  noise.threshold <- min(expression.matrix)
  
  tbl <- tibble::tibble(
    gene = rep(rownames(expression.matrix), each = length(levels(condition))),
    condition = rep(levels(condition), times = nrow(expression.matrix)),
    mean = NA,
    sd = NA
  ) %>%
    dplyr::mutate(condition = factor(condition, levels = unique(condition)))
  
  pb = utils::txtProgressBar(min = 0, max = nrow(tbl), initial = 0, style = 3) 
  for(i in seq_len(nrow(tbl))){
    vec <- expression.matrix[tbl$gene[i], condition == tbl$condition[i]]
    tbl$mean[i] <- mean(vec)
    tbl$sd[i] <- stats::sd(vec)
    if(i %% 1000 == 0 | i == nrow(tbl)) utils::setTxtProgressBar(pb, i)
  }
  tbl <- tbl %>% dplyr::mutate(sd = pmax(sd, noise.threshold))
}

#' Determine the pattern between two intervals
#' @description This function checks if the two input intervals oferlap and
#' outputs the corresponding pattern (up, down, or straight) based on that.
#' @param min1,max1,min2,max2 the endpoints of the two intervals
#' @return A single character (one of "U", "D", "S") representing the pattern
#' @export
#' @examples
#' determine_uds(10, 20, 15, 25) # overlap
#' determine_uds(10, 20, 25, 35) # no overlap
determine_uds <- function(min1, max1, min2, max2){
  if(max1 < min2) return("U") else if(min1 > max2) return("D") else return("S")
}

#' Create a matrix of the patterns between conditions
#' @description This function determines the patterns between different
#' conditions of each gene. It should be applied to the output of
#' \code{\link{calculate_condition_mean_sd_per_gene}}.
#' @param tbl the output of \code{\link{calculate_condition_mean_sd_per_gene}}
#' @param n_sd number of standard deviations from the mean to use to
#' construct the intervals; default is 2
#' @return A matrix of single character patterns between conditions. The last
#' column is named pattern and is a concatenation of all other columns.
#' @export
#' @examples
#' expression.matrix <- as.matrix(read.csv(
#'   system.file("extdata", "expression_matrix.csv", package = "bulkAnalyseR"), 
#'   row.names = 1
#' ))
#' expression.matrix.preproc <- preprocessExpressionMatrix(expression.matrix)
#' condition <- factor(rep(c("0h", "12h", "36h"), each = 2))
#' tbl <- calculate_condition_mean_sd_per_gene(expression.matrix.preproc[1:10, ], condition)
#' patmat <- make_pattern_matrix(tbl)
#' patmat
make_pattern_matrix <- function(tbl, n_sd = 2){
  tbl <- tbl %>% dplyr::mutate(min = mean - n_sd * sd, max = mean + n_sd * sd)
  
  genes <- unique(tbl$gene)
  conditions <- unique(tbl$condition)
  nc <- length(conditions)
  patmat <- matrix(nrow = length(genes), ncol = nc - 1) %>% as.data.frame()
  rownames(patmat) <- genes
  colnames(patmat) <- sapply(X = seq_len(nc - 1), FUN = function(j){
    paste(conditions[j], "->", conditions[j + 1])
  })
  
  pb = utils::txtProgressBar(min = 0, max = nrow(patmat), initial = 0, style = 3) 
  for(i in seq_len(nrow(patmat))){
    tbl.sub <- tbl[((i - 1) * nc + 1):(i * nc), ]
    for(j in seq_len(nc - 1)){
      patmat[i, j] <- determine_uds(tbl.sub$min[j], tbl.sub$max[j], tbl.sub$min[j + 1], tbl.sub$max[j + 1])
    }
    if(i %% 1000 == 0 | i == nrow(patmat)) utils::setTxtProgressBar(pb, i)
  }
  
  patmat <- patmat %>% 
    dplyr::mutate(pattern = unname(apply(patmat, 1, paste0, collapse = ""))) %>%
    as.matrix()
}

#' Create a matrix of the average expression of each gene in each condition
#' @description This function reshapes the tibble output of 
#' \code{\link{calculate_condition_mean_sd_per_gene}} into a matrix of
#' average expression by condition. Its output can be used by
#' \code{\link{expression_heatmap}}.
#' @inheritParams make_pattern_matrix
#' @param genes gene names to use for the output; if NULL (the default),
#' all genes will be used
#' @return A matrix of averaged expression per gene in each condition.
#' @export
#' @examples
#' expression.matrix <- as.matrix(read.csv(
#'   system.file("extdata", "expression_matrix.csv", package = "bulkAnalyseR"), 
#'   row.names = 1
#' ))
#' expression.matrix.preproc <- preprocessExpressionMatrix(expression.matrix)
#' condition <- factor(rep(c("0h", "12h", "36h"), each = 2))
#' tbl <- calculate_condition_mean_sd_per_gene(expression.matrix.preproc[1:10, ], condition)
#' heatmat <- make_heatmap_matrix(tbl)
#' heatmat
make_heatmap_matrix <- function(tbl, genes = NULL){
  if(is.null(genes)) genes <- unique(tbl$gene)
  tbl %>%
    dplyr::select(-sd) %>%
    dplyr::filter(gene %in% genes) %>%
    tidyr::pivot_wider(names_from = "condition", values_from = "mean") %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()
}

#' Create a line plot of average expression across conditions
#' @description This function creates a line plot of average expression
#' across conditions for a selection of genes, usually to visualise
#' an expression pattern.
#' @inheritParams make_heatmap_matrix
#' @param type whether the expression values should be scaled using their mean
#' (default), log-transformed, or not adjusted for the plot
#' @param show.legend whether to show the gene names in the legend;
#' should be avoided in many genes are plotted
#' @return A matrix of average gene expression per gene in each condition.
#' @export
#' @examples
#' expression.matrix <- as.matrix(read.csv(
#'   system.file("extdata", "expression_matrix.csv", package = "bulkAnalyseR"), 
#'   row.names = 1
#' ))
#' expression.matrix.preproc <- preprocessExpressionMatrix(expression.matrix)
#' condition <- factor(rep(c("0h", "12h", "36h"), each = 2))
#' tbl <- calculate_condition_mean_sd_per_gene(expression.matrix.preproc[1:10, ], condition)
#' plot_line_pattern(tbl)
plot_line_pattern <- function(
  tbl, 
  genes = NULL, 
  type = c("Mean Scaled", 'Log2 Expression', 'Expression'),
  show.legend = FALSE
){
  if(is.null(genes)) genes <- unique(tbl$gene)
  
  tbl <- tbl %>%
    dplyr::filter(gene %in% genes)
  
  type <- type[1]
  if(type == 'Expression'){
    tbl <- dplyr::mutate(tbl, scale = mean)
  }else if(type == 'Log2 Expression'){
    tbl <- dplyr::mutate(tbl, scale = log2(mean + 1))
  }else if(type == "Mean Scaled"){
    tbl$scale <- sapply(X = seq_len(nrow(tbl)), FUN = function(i){
      m <- mean(tbl$mean[tbl$gene == tbl$gene[i]])
      tbl$mean[i] / m
    })
  }
  
  p <- ggplot(tbl) +
    theme_minimal() +
    geom_line(aes(x = condition, y = scale, group = gene, colour = gene))
  
  if(!show.legend) p <- p + theme(legend.position = "none")
  
  p
}

#' Rescale a matrix
#' @description This function rescales the rows of a matrix according to the
#' specified type.
#' @param mat the matrix to rescale
#' @param type type of rescaling; one of "Expression" (defautl, does nothing),
#' "Log2 Expression" (returns log2(x + 1) for every value), "Mean Scaled" (each
#' row is scaled by its average), "Z-score" (each row is centered and scaled
#' to mean = 0 and sd = 1)
#' @return The rescaled matrix.
#' @export
#' @examples
#' mat = matrix(1:10, nrow = 2, ncol = 5)
#' rescale_matrix(mat, type = "Expression")
#' rescale_matrix(mat, type = "Log2 Expression")
#' rescale_matrix(mat, type = "Mean Scaled")
#' rescale_matrix(mat, type = "Z-score")
rescale_matrix <- function(
  mat, 
  type = c('Expression', 'Log2 Expression', 'Mean Scaled', 'Z-score')
){
  type <- type [1]
  if(type == 'Expression'){
    mat <- mat
  }else if(type == 'Log2 Expression'){
    mat <- log2(mat + 1)
  }else if(type == 'Mean Scaled'){
    mat <- mat / rowMeans(mat)
  }else if(type == 'Z-score'){
    mat <- t(scale(t(mat)))
  }
  mat
}
