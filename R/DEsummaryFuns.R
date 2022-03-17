#' Create heatmap of an expression matrix
#' @description This function creates a heatmap to visualise an expression matrix
#' @inheritParams generateShinyApp
#' @inheritParams rescale_matrix
#' @param expression.matrix.subset a subset of rows from the expression matrix; 
#' rows correspond to genes and columns correspond to samples
#' @param top.annotation.ids a vector of column indices denoting which columns
#' of the metadata should become heatmap annotations
#' @param show.column.names whether to show the column names below the heatmap;
#' default is TRUE
#' @return The heatmap as detailed in the ComplexHeatmap package.
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
#' print(expression_heatmap(head(expression.matrix.preproc), NULL, metadata))
expression_heatmap <- function(
  expression.matrix.subset,
  top.annotation.ids = NULL,
  metadata,
  type = c('Z-score', 'Log2 Expression', 'Expression'),
  show.column.names = TRUE
){
  heatmat <- as.matrix(expression.matrix.subset)
  
  type <- type[1]
  heatmat <- rescale_matrix(heatmat, type)
  if(!show.column.names){colnames(heatmat) <- NULL}
  
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
  if (type != 'Z-score'){
    breaks <- seq(min(heatmat), max(heatmat), (max(heatmat) - min(heatmat)) / 9)
    colours = c("#FFFFFF", RColorBrewer::brewer.pal(n = 9, name = "YlOrRd"))
  }else{
    breaks <- seq(-3, 3, 6 / 9)
    colours = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
    heatmat[heatmat > 3] <- 3
    heatmat[heatmat < (-3)] <- (-3)
  }
  ComplexHeatmap::Heatmap(
    matrix = heatmat,
    name = "Scale",
    col = circlize::colorRamp2(
      breaks = breaks,
      colors = colours
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    top_annotation = top.annotation
  )
}
