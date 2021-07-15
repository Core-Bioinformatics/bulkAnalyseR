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
  n.abundant = NULL, 
  top.annotation.df = NULL, 
  top.annotation.colours = NULL,
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
  
  if(show.values){
    cell.fun <- function(j, i, x, y, width, height, fill) {
      grid::grid.text(sub("0.", ".", sprintf("%.2f", heatmat[i,j])), x, y, gp = grid::gpar(fontsize = 8))}
  }else{
    cell.fun = NULL
  }
  
  if(!is.null(top.annotation.df) & !is.null(top.annotation.colours)){
    top.annotation <- ComplexHeatmap::HeatmapAnnotation(df = top.annotation.df,
                                                        col = top.annotation.colours,
                                                        show_annotation_name = FALSE)
  }else{
    top.annotation <- NULL
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

plotPCA <- function(expression.matrix,geom.ind=c("point","text"))
{
  expr.PCA <- prcomp(expression.matrix[,3:ncol(expression.matrix)],
                     center = TRUE, scale = TRUE)
  
  staph.p <- fviz_pca_ind(expr.PCA,
                          geom.ind = geom.ind, # show points only, text only, or both
                          col.ind = expression.matrix[,2], # color by groups
                          # palette = c("#00AFBB", "#FC4E07","#FF7256"),
                          addEllipses = TRUE, # Concentration ellipses
                          ellipse.type = "confidence", # for confidence ellipses
                          legend.title = "Groups") 
  
  toPlot <- ggpubr::ggpar(staph.p,
                          title = "Principal Component Analysis",
                          subtitle = title,
                          # caption = "Source: nanostring",
                          xlab = paste("PC1\n","Proportion of variance = ",summary(expr.PCA)$importance[2,1]*100,"%"), 
                          ylab = paste("PC2\n","Proportion of variance = ",summary(expr.PCA)$importance[2,2]*100,"%"),
                          legend.title = "Condition", legend.position = "top",
                          ggtheme = theme_gray())
  return(toPlot)
}
