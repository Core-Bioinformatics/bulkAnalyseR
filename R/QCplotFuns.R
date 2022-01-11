#' Calculate the Jaccard similarity index (JSI) between two vectors
#' @param a,b two vectors 
#' @return The JSI of the two vectors, a single value between 0 and 1.
#' @export
#' @examples
#' jaccard_index(1:4, 2:6)
jaccard_index <- function(a, b){
  if ((length(a) == 0) & (length(b) == 0)) {
    return(1)
  } else {
    u <- length(union(a, b))
    i <- length(intersect(a, b))
    return(i / u)
  }
}

#' Create a heatmap of the Jaccard similarity index (JSI) between samples of
#' an experiment
#' @description This function creates a JSI heatmap between all samples in the
#' expression matrix using the specified number of most abundant genes as
#' input. Metadata columns are used as annotations.
#' @inheritParams generateShinyApp
#' @param top.annotation.ids a vector of column indices denoting which columns
#' of the metadata should become heatmap annotations
#' @param n.abundant number of most abundant genes to use for the JSI calculation
#' @param show.values whether to show the JSI values within the heatmap squares
#' @return The JSI heatmap as detailed in the ComplexHeatmap package.
#' @export
#' @examples
#' expression.matrix <- as.matrix(read.csv(
#'   system.file("extdata", "expression_matrix.csv", package = "bulkAnalyseR"), 
#'   row.names = 1
#' ))
#' expression.matrix.preproc <- preprocessExpressionMatrix(expression.matrix)[1:500, ]
#' metadata <- data.frame(
#'   srr = colnames(expression.matrix.preproc), 
#'   timepoint = rep(c("0h", "12h", "36h"), each = 2)
#' )
#' print(jaccard_heatmap(expression.matrix.preproc, metadata, n.abundant = 100))
jaccard_heatmap <- function(
  expression.matrix,
  metadata,
  top.annotation.ids = NULL, 
  n.abundant = NULL, 
  show.values = TRUE
){
  n.abundant <- min(n.abundant, nrow(expression.matrix))
  n.samples <- ncol(expression.matrix)
  heatmat <- matrix(0, nrow = n.samples, ncol = n.samples)
  for (i in seq_len(n.samples)){
    for (j in seq_len(i)){
      i.gene.indices <- order(expression.matrix[,i], decreasing = TRUE)[1:n.abundant]
      j.gene.indices <- order(expression.matrix[,j], decreasing = TRUE)[1:n.abundant]
      heatmat[i, j] <- heatmat[j, i] <- jaccard_index(i.gene.indices, j.gene.indices)
    }
  }
  rownames(heatmat) <- colnames(heatmat) <- colnames(expression.matrix)
  
  if(!is.null(top.annotation.ids)){
    qual.col.pals = dplyr::filter(RColorBrewer::brewer.pal.info, .data$category == 'qual')
    col.vector = unique(unlist(mapply(RColorBrewer::brewer.pal, 
                                      qual.col.pals$maxcolors, 
                                      rownames(qual.col.pals))))
    
    top.annotation.colour.list=list()
    colind <- 1
    for(annos in seq_len(length(top.annotation.ids))){
      values <- as.character(unique(metadata[, top.annotation.ids[annos]]))
      vec <- vector(mode="character")
      for(i in seq_len(length(values))){
        vec <- c(vec, col.vector[colind])
        names(vec)[i] <- values[i]
        colind <- colind + 1
      }
      top.annotation.colour.list[[colnames(metadata)[top.annotation.ids[annos]]]] <- vec
    }
    
    top.annotation.df <- as.data.frame(metadata[, top.annotation.ids])
    colnames(top.annotation.df) <- colnames(metadata)[top.annotation.ids]
    
    top.annotation <- ComplexHeatmap::HeatmapAnnotation(df = top.annotation.df,
                                                        col = top.annotation.colour.list,
                                                        show_annotation_name = FALSE)
  }else{
    top.annotation <- NULL
  }
  
  if(show.values){
    cell.fun <- function(j, i, x, y, width, height, fill) {
      grid::grid.text(sub("0.", ".", sprintf("%.2f", heatmat[i,j])), x, y, gp = grid::gpar(fontsize = 8))}
  }else{
    cell.fun = NULL
  }
  
  jaccard.plot <- ComplexHeatmap::Heatmap(
    matrix = heatmat,
    name = "Scale",
    col = circlize::colorRamp2(
      breaks = (1:10)/10, 
      colors = c("#FFFFFF", RColorBrewer::brewer.pal(n = 9, name = "YlOrRd"))
    ),
    cell_fun = cell.fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    top_annotation = top.annotation,
    heatmap_legend_param = list(at = (1:10)/10)
  )
}

#' Create a principal component analysis (PCA) plot the samples of an experiment
#' @description This function creates a PCA plot between all samples in the
#' expression matrix using the specified number of most abundant genes as
#' input. A metadata column is used as annotation.
#' @inheritParams jaccard_heatmap
#' @inheritParams volcano_enhance
#' @param annotation.id a column index denoting which column of the metadata
#' should be used to colour the points and draw confidence ellipses
#' @param show.labels whether to label the points with the sample names
#' @param show.ellipses whether to draw confidence ellipses
#' @return The PCA plot as a ggplot object.
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
#' plot_pca(expression.matrix.preproc, metadata, 2)
plot_pca <- function(
  expression.matrix,
  metadata,
  annotation.id, 
  n.abundant = NULL,
  show.labels = FALSE,
  show.ellipses = TRUE,
  label.force = 1
){
  annotation.name <- colnames(metadata)[annotation.id]
  n.abundant <- min(n.abundant, nrow(expression.matrix))
  
  expr.PCA.list <- expression.matrix %>%
    as.data.frame() %>%
    dplyr::filter(seq_len(nrow(expression.matrix)) %in% 
                    utils::tail(order(rowSums(expression.matrix)), n.abundant)) %>%
    t()
  expr.PCA.list <- expr.PCA.list[, apply(expr.PCA.list, 2, function(x) max(x) != min(x))] %>%
    stats::prcomp(center = TRUE, scale = TRUE)
  
  expr.PCA <- dplyr::mutate(
    as.data.frame(expr.PCA.list$x),
    name = factor(metadata[, 1], levels = metadata[, 1]),
    condition = factor(metadata[, annotation.id], levels = unique(metadata[, annotation.id]))
  )
  if(min(table(metadata[, annotation.id])) <= 2){
    expr.PCA.2 <- expr.PCA
    expr.PCA.2$PC1 <- expr.PCA.2$PC1 * 1.001
    expr.PCA.2$PC2 <- expr.PCA.2$PC2 * 1.001
    expr.PCA.full <- rbind(expr.PCA, expr.PCA.2)
  }
  else {expr.PCA.full <- expr.PCA}
  pca.plot <- ggplot(expr.PCA.full, aes(x = .data$PC1, y = .data$PC2, colour = .data$condition)) +
    theme_minimal() +
    geom_point() +
    labs(x = paste0("PC1 (proportion of variance = ", summary(expr.PCA.list)$importance[2, 1] * 100, "%)"),
         y = paste0("PC2 (proportion of variance = ", summary(expr.PCA.list)$importance[2, 2] * 100, "%)"),
         colour = annotation.name)
  if(show.ellipses){
    pca.plot <- pca.plot + 
      ggforce::geom_mark_ellipse(aes(fill = .data$condition, colour = .data$condition), show.legend = FALSE)
  }
  if(show.labels){
    pca.plot <- pca.plot +
      ggrepel::geom_label_repel(
        data = expr.PCA,
        mapping = aes(x = .data$PC1, y = .data$PC2, colour = .data$condition, label = .data$name),
        max.overlaps = nrow(expr.PCA),
        force = label.force,
        point.size = NA
      )
  }
  
  pca.plot
}

qc_ma_plot = function(expression.matrix, metadata, i, j, include.guidelines = TRUE){
  # mask away the zeros
  zero.mask <- !(expression.matrix[, i] == 0 | expression.matrix[,j] == 0)
  l1 <- log2(expression.matrix[zero.mask, i])
  l2 <- log2(expression.matrix[zero.mask, j])
  
  #calculate M and A
  m <- l1 - l2
  a <- 0.5 * (l1 + l2)
  upper.lim <- max(abs(m))
  lower.lim <- -upper.lim
  expr.MA <- data.frame(A = a, M = m)
  #define MA plot
  p <- ggplot(expr.MA, aes(x = .data$A, y = .data$M, color = 'red', fill = 'red')) +
    geom_point(alpha = 0.05) +
    theme_minimal() +
    theme(legend.position = "none") + 
    ylim(lower.lim, upper.lim)
  
  if(include.guidelines){
    p <- p + geom_hline(yintercept = 0.5, colour = 'gray60') +
      geom_hline(yintercept = -0.5, colour = 'gray60') +
      geom_hline(yintercept = 1, colour = 'gray60') +
      geom_hline(yintercept = -1, colour = 'gray60')
  }
  p
}
