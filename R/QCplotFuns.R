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
#' @param show.row.column.names whether to show the row and column names below 
#' the heatmap; default is TRUE
#' @return The JSI heatmap as detailed in the ComplexHeatmap package.
#' @export
#' @examples
#' expression.matrix.preproc <- as.matrix(read.csv(
#'   system.file("extdata", "expression_matrix_preprocessed.csv", package = "bulkAnalyseR"), 
#'   row.names = 1
#' ))[1:500,]
#' 
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
  show.values = TRUE,
  show.row.column.names = TRUE
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
  if(show.row.column.names){
    rownames(heatmat) <- colnames(heatmat) <- colnames(expression.matrix)
  }
  
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
        if(colind > length(col.vector)){colind <- colind %% length(col.vector)}
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
#' expression.matrix.preproc <- as.matrix(read.csv(
#'   system.file("extdata", "expression_matrix_preprocessed.csv", package = "bulkAnalyseR"), 
#'   row.names = 1
#' ))[1:500,]
#' 
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
    condition = if(!is.factor(metadata[,annotation.id])){
      factor(metadata[, annotation.id], levels = unique(metadata[, annotation.id]))}else{metadata[,annotation.id]}

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

qc_ma_plot <- function(expression.matrix, metadata, i, j, include.guidelines = TRUE){
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

#' Create a density plot of log2 expression across samples of an experiment
#' @description This function creates a density plot between all samples in the
#' expression matrix. Metadata columns are used to group samples.
#' @inheritParams generateShinyApp
#' @param annotation.id name of metadata column on which to group samples
#' @return The density plot as a ggplot object.
#' @export
#' @examples
#' expression.matrix.preproc <- as.matrix(read.csv(
#'   system.file("extdata", "expression_matrix_preprocessed.csv", package = "bulkAnalyseR"), 
#'   row.names = 1
#' ))[1:500,]
#' 
#' metadata <- data.frame(
#'   srr = colnames(expression.matrix.preproc), 
#'   timepoint = rep(c("0h", "12h", "36h"), each = 2)
#' )
#' print(qc_density_plot(expression.matrix.preproc, metadata, 'timepoint'))
#' 
qc_density_plot <- function(expression.matrix,
                            metadata,
                            annotation.id){

  log.expression.matrix <- log2(expression.matrix + 1)
  melted.expression.matrix <- tidyr::pivot_longer(as.data.frame(log.expression.matrix),
                                                  cols=colnames(log.expression.matrix))
  melted.expression.matrix$colour <- metadata[,annotation.id][ match(melted.expression.matrix$name,metadata[,1]) ]
  p  <- ggplot(melted.expression.matrix, aes(x = .data$value, 
                                             group = .data$name, 
                                             color = .data$colour)) + 
           geom_density() + 
           theme_minimal() + 
           xlab('log2 expression')
  if (ncol(expression.matrix)>20){
    p <- p + theme(legend.position = "none")
  }
  p
}

#' Create a violin/box plot of expression across samples of an experiment
#' @description This function creates a combined violin and box plot between 
#' all samples in the expression matrix. Metadata columns are used to colour samples.
#' @inheritParams generateShinyApp
#' @param annotation.id name of metadata column on which to group samples
#' @param log.transformation whether expression should be shown on log (default) or 
#' linear scale
#' @return The violin/box plot as a ggplot object.
#' @export
#' @examples
#' expression.matrix.preproc <- as.matrix(read.csv(
#'   system.file("extdata", "expression_matrix_preprocessed.csv", package = "bulkAnalyseR"), 
#'   row.names = 1
#' ))[1:500,]
#' 
#' metadata <- data.frame(
#'   srr = colnames(expression.matrix.preproc), 
#'   timepoint = rep(c("0h", "12h", "36h"), each = 2)
#' )
#' print(qc_violin_plot(expression.matrix.preproc, metadata, 'timepoint'))
#' 
qc_violin_plot = function(expression.matrix,
                          metadata,
                          annotation.id,
                          log.transformation = TRUE){
  # mask away the zeros
  if (log.transformation){
    log.expression.matrix <- log2(expression.matrix + 1)
  } else {log.expression.matrix <- expression.matrix}
  
  melted.expression.matrix <- tidyr::pivot_longer(as.data.frame(log.expression.matrix),
                                                  cols = colnames(log.expression.matrix))
  melted.expression.matrix$colour <- metadata[,annotation.id][ match(melted.expression.matrix$name,metadata[,1]) ]
  
  if (is.factor(metadata[,1])) {
    melted.expression.matrix$name = factor(melted.expression.matrix$name, levels = levels(metadata[,1]))
  }
  p  <- ggplot(melted.expression.matrix, 
               aes(y = .data$value, 
                   x = .data$name, 
                   color = .data$colour)) + 
          geom_violin() + 
          geom_boxplot(width = 0.1) +
          theme_minimal() +
          scale_color_discrete(name = annotation.id) + 
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  p
}

#' Create a bar plot of expression for selected genes across samples in an experiment
#' @description This function creates a clustered bar plot between all samples in the 
#' expression matrix for the selection of genes.
#' @param sub.expression.matrix subset of the expression matrix containing only selected 
#' genes
#' @param log.transformation whether expression should be shown on log (default) or 
#' linear scale
#' @return The bar plot as a ggplot object.
#' @export
#' @examples
#' expression.matrix.preproc <- as.matrix(read.csv(
#'   system.file("extdata", "expression_matrix_preprocessed.csv", package = "bulkAnalyseR"), 
#'   row.names = 1
#' ))[1:500,]
#' 
#' print(genes_barplot(head(expression.matrix.preproc,5)))
#' 
genes_barplot <- function(sub.expression.matrix,
                          log.transformation = TRUE){
  if (log.transformation){
    log.expression.matrix <- data.frame(log2(as.matrix(sub.expression.matrix) + 1))
  } else {
    log.expression.matrix <- data.frame(as.matrix(sub.expression.matrix))
    }
  log.expression.matrix$gene <- rownames(log.expression.matrix)
  melted.expression.matrix <- tidyr::pivot_longer(log.expression.matrix,
                                                  cols = colnames(log.expression.matrix)[1:(ncol(log.expression.matrix)-1)])
  
  p <- ggplot(melted.expression.matrix, aes(x = .data$name, 
                                            y = .data$value, 
                                            fill = .data$gene)) +
                geom_bar(stat='identity',position='dodge') +
                theme_minimal() + 
                ylab(ifelse(log.transformation,'log2 expression','expression')) +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  p
}

#' Create a scatter plot of expression between two samples of an experiment
#' @description This function creates a scatter plot between two samples.
#' @param sub.expression.matrix subset of the expression matrix containing only the two
#' selected samples
#' @param genes.to.highlight vector of gene names to highlight. 
#' These should match entries in the anno NAME column.
#' @param anno annotation data frame containing a match between the row names
#' of the expression.matrix (usually ENSEMBL IDs) and the gene names that
#' should be rendered within the app and in output files; this object is
#' created by \code{\link{generateShinyApp}} using the org.db specified
#' @param log.transformation whether expression should be shown on log (default) or 
#' linear scale
#' @return The scatter plot as a ggplot object.
#' @export
#' @examples
#' expression.matrix.preproc <- as.matrix(read.csv(
#'   system.file("extdata", "expression_matrix_preprocessed.csv", package = "bulkAnalyseR"), 
#'   row.names = 1
#' ))[,1:2]
#' 
#' print(scatter_plot(expression.matrix.preproc, c()))
#' 
scatter_plot <- function(sub.expression.matrix,
                         anno,
                         genes.to.highlight = c(),
                         log.transformation = TRUE){
  
  if (log.transformation){
    log.expression.matrix <- data.frame(log2(as.matrix(sub.expression.matrix) + 1))
  } else {
    log.expression.matrix <- sub.expression.matrix
  }
  
  x.name <- colnames(log.expression.matrix)[1]
  y.name <- colnames(log.expression.matrix)[2]
  
  colnames(log.expression.matrix) <- c('x_var','y_var')
  
  p <- ggplot(log.expression.matrix, aes(x = .data$x_var, 
                                         y = .data$y_var)) +
    geom_point(stat = 'identity',
               position = 'dodge',
               alpha = 0.1,
               aes(fill = 'red', color = 'red')) +
    geom_abline(slope = 1,
                intercept = 0) +
    theme_minimal() + 
    theme(legend.position = "none") +
    ylab(paste0(ifelse(log.transformation,'log2 expression ','expression '), y.name)) +
    xlab(paste0(ifelse(log.transformation,'log2 expression ','expression '), x.name))
  
  if (length(genes.to.highlight)){
    label.expression.matrix <- log.expression.matrix[ anno$ENSEMBL[ match(genes.to.highlight,anno$NAME) ], ]
    label.expression.matrix$name <- genes.to.highlight
    
    p <- p + ggrepel::geom_label_repel(
                data = label.expression.matrix, 
                mapping = aes(x = .data$x_var, 
                              y = .data$y_var, 
                              label = .data$name),
                max.overlaps = Inf,
                force = 1,
                point.size = NA
              )
  }
  
  p
}
