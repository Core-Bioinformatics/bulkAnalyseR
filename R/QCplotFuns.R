jaccard_index <- function(a, b){
  if ((length(a) == 0) & (length(b) == 0)) {
    return(1)
  } else {
    u <- length(union(a,b))
    i <- length(intersect(a,b))
    return(i/u)
  }
}

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

plot_pca <- function(
  expression.matrix,
  metadata,
  annotation.id = NULL, 
  n.abundant = NULL,
  show.labels
){
  annotation.name <- colnames(metadata)[annotation.id]
  n.abundant <- min(n.abundant, nrow(expression.matrix))
  
  expr.PCA.list <- expression.matrix %>%
    dplyr::filter(seq_len(nrow(expression.matrix)) %in% 
                    utils::tail(order(rowSums(expression.matrix)), n.abundant)) %>%
    t() %>%
    stats::prcomp(center = TRUE, scale = TRUE)
  
  expr.PCA <- dplyr::mutate(
    as.data.frame(expr.PCA.list$x),
    name = factor(metadata[, 1], levels = metadata[, 1]),
    condition = factor(metadata[, annotation.id], levels = unique(metadata[, annotation.id]))
  )
  
  pca.plot <- ggplot(expr.PCA, aes(x = .data$PC1, y = .data$PC2, colour = .data$condition)) +
    theme_minimal() +
    geom_point() +
    labs(x = paste0("PC1 (proportion of variance = ", summary(expr.PCA.list)$importance[2, 1] * 100, "%)"),
         y = paste0("PC2 (proportion of variance = ", summary(expr.PCA.list)$importance[2, 2] * 100, "%)"),
         colour = annotation.name) +
    ggforce::geom_mark_ellipse(aes(fill = .data$condition, colour = .data$condition), show.legend = FALSE)
  
  if(show.labels){
    pca.plot <- pca.plot +
      ggrepel::geom_label_repel(aes(label = .data$name))
  }
  
  pca.plot
}
