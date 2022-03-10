#' Perform GRN inference
#' @description This function performs Gene Regulatory Network inference on
#' a subset of the expression matrix, for a set of potential regulators
#' @inheritParams GRNpanel
#' @param seed the random seed to be set when running GRN inference, to ensure
#' reproducibility of outputs
#' @param regulators the regulators of interest around which the GRN is built;
#' must be row names of the expression matrix
#' @param condition name of the metadata column to select samples from
#' @param samples names of the sample groups to select; must appear in
#' \code{metadata[[condition]]}
#' @param inference_method method used for GRN inference; supported methods
#' are GENIE3 and GNET2.
#' @return The adjacency matrix of the inferred network
#' @export
#' @examples 
#' expression.matrix.preproc <- as.matrix(read.csv(
#'   system.file("extdata", "expression_matrix_preprocessed.csv", package = "bulkAnalyseR"), 
#'   row.names = 1
#' ))[1:500, ]
#' 
#' metadata <- data.frame(
#'   srr = colnames(expression.matrix.preproc), 
#'   timepoint = rep(c("0h", "12h", "36h"), each = 2)
#' )
#' 
#' anno <- AnnotationDbi::select(
#'   getExportedValue('org.Mm.eg.db', 'org.Mm.eg.db'),
#'   keys = rownames(expression.matrix.preproc),
#'   keytype = 'ENSEMBL',
#'   columns = 'SYMBOL'
#' ) %>%
#'   dplyr::distinct(ENSEMBL, .keep_all = TRUE) %>%
#'   dplyr::mutate(NAME = ifelse(is.na(SYMBOL), ENSEMBL, SYMBOL))
#' 
#' res <- infer_GRN(
#'   expression.matrix = expression.matrix.preproc,
#'   metadata = metadata,
#'   anno = anno,
#'   seed = 13,
#'   regulators = c("Hecw2", "Akr1cl"),
#'   condition = "timepoint",
#'   samples = "0h",
#'   inference_method = "GENIE3"
#' )
infer_GRN <- function(expression.matrix, metadata, anno, seed = 13, 
                      regulators, condition, samples, inference_method){
  inference_method <- inference_method[1]
  regulator.ids <- anno$ENSEMBL[match(regulators, anno$NAME)]
  samples <- metadata[[condition]] %in% samples
  set.seed(seed)
  if(inference_method == "GENIE3"){
    res <- GENIE3::GENIE3(expression.matrix[, samples], regulators = regulator.ids)
  }else if(inference_method == "GNET2"){
    gnet_out <- GNET2::gnet(expression.matrix[, samples], reg_names = regulator.ids)
    res <- GNET2::extract_edges(gnet_out)
  }
  res
}

#' Convert the adjacency matrix to network links
#' @description This function converts an adjacency matrix to a data frame
#' of network links, subset to the most important ones.
#' @param weightMat the (weighted) adjacency matrix - regulators in rows,
#' targets in columns
#' @param plotConnections the number of connections to subset to
#' @return A data frame with fields from, to and value, describing the edges
#' of the network
#' @export
#' @examples 
#' weightMat <- matrix(
#'   c(0.1, 0.4, 0.8, 0.3), nrow = 2, ncol = 2,
#'   dimnames = list("regulators" = c("r1", "r2"), "targets" = c("t1", "t2"))
#' )
#' get_link_list_rename(weightMat, 2)
get_link_list_rename <- function(weightMat, plotConnections){
  GENIE3::getLinkList(weightMat, plotConnections) %>%
    dplyr::mutate(from = as.character(.data$regulatoryGene), 
                  to = as.character(.data$targetGene), 
                  value = .data$weight, 
                  regulatoryGene = NULL, 
                  targetGene = NULL,
                  weight = NULL)
}

#' Find recurring regulatory targets
#' @description This function finds regulatory targets that appear as the 
#' same network edge in more than one of the input networks.
#' @inheritParams get_link_list_rename
#' @param weightMatList a list of (weighted) adjacency matrices; 
#' each list element must be an adjacency matrix with regulators in rows,
#' targets in columns
#' @return A vector containing the names of the recurring targets
#' @export
#' @examples 
#' weightMat1 <- matrix(
#'   c(0.1, 0.4, 0.8, 0.3), nrow = 2, ncol = 2,
#'   dimnames = list("regulators" = c("r1", "r2"), "targets" = c("t1", "t2"))
#' )
#' weightMat2 <- matrix(
#'   c(0.1, 0.2, 0.8, 0.3), nrow = 2, ncol = 2,
#'   dimnames = list("regulators" = c("r1", "r2"), "targets" = c("t1", "t2"))
#' )
#' find_targets_with_recurring_edges(list(weightMat1, weightMat2), 2)
find_targets_with_recurring_edges <- function(weightMatList, plotConnections){
  edges <- lapply(weightMatList, function(wm){
    get_link_list_rename(wm, plotConnections)
  }) %>%
    dplyr::bind_rows() 
  edges %>%
    dplyr::group_by(.data$from, .data$to) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "keep") %>%
    dplyr::filter(.data$n > 1) %>%
    dplyr::pull(.data$to) %>%
    dplyr::setdiff(edges$from)
}

#' Plot a GRN
#' @description This function creates a network plot of a GRN.
#' @inheritParams GRNpanel
#' @inheritParams get_link_list_rename
#' @param plot_position_grid,n_networks the position of the plot in 
#' the grid (1-4) and the number of networks shown (1-4); these are
#' solely used for hiding unwanted plots in the shiny app
#' @param recurring_targets targets to be highlighted; usually the
#' result of \code{\link{find_targets_with_recurring_edges}}
#' @return A network plot. See visNetwork package for more details.
#' @export
#' @examples 
#' weightMat1 <- matrix(
#'   c(0.1, 0.4, 0.8, 0.3), nrow = 2, ncol = 2,
#'   dimnames = list("regulators" = c("r1", "r2"), "targets" = c("t1", "t2"))
#' )
#' weightMat2 <- matrix(
#'   c(0.1, 0.2, 0.8, 0.3), nrow = 2, ncol = 2,
#'   dimnames = list("regulators" = c("r1", "r2"), "targets" = c("t1", "t2"))
#' )
#' anno <- tibble::tibble(ENSEMBL = c("r1", "r2", "t1", "t2"), NAME = ENSEMBL)
#' recurring_targets <- find_targets_with_recurring_edges(list(weightMat1, weightMat2), 2)
#' plot_GRN(weightMat1, anno, 2, 1, 1, recurring_targets)
#' plot_GRN(weightMat2, anno, 2, 1, 1, recurring_targets)
plot_GRN <- function(weightMat, anno, plotConnections, 
                     plot_position_grid, n_networks, recurring_targets){
  
  if(n_networks >= plot_position_grid){
    edges <- get_link_list_rename(weightMat, plotConnections)
    nodes <- tibble::tibble(
      id = c(edges$from, edges$to),
      label = anno$NAME[match(.data$id, anno$ENSEMBL)],
      group = rep(c("regulator", "target"), each = nrow(edges)),
      color = NA
    ) %>%
      dplyr::distinct(.data$id, .keep_all = TRUE)
    nodes$group[nodes$id %in% recurring_targets] <- "recurring_target"
    
    color_regulator <- list("background" = '#D2E5FF')
    color_target <- list("background" = '#E0E0E0')
    color_recurring_target <- list("background" = '#ACE9B4')
    colors.list <- list(color_regulator, color_target, color_recurring_target)
    
    for(i in seq_len(nrow(nodes))){
      nodes$color[i] <- colors.list[[match(
        nodes$group[i], c("regulator", "target", "recurring_target")
      )]]
    }
    
    visNetwork::visNetwork(nodes, edges)
  }else{
    NULL
  }
}

#' Visualise the overlap of edges between different networks
#' @description This function creates an UpSet plot of the intersections
#' and specific differences of the edges in the input networks.
#' @inheritParams get_link_list_rename
#' @inheritParams find_targets_with_recurring_edges
#' @return An UpSet plot. See UpSetR package for more details.
#' @export
#' @examples 
#' weightMat1 <- matrix(
#'   c(0.1, 0.4, 0.8, 0.3), nrow = 2, ncol = 2,
#'   dimnames = list("regulators" = c("r1", "r2"), "targets" = c("t1", "t2"))
#' )
#' weightMat2 <- matrix(
#'   c(0.1, 0.2, 0.8, 0.3), nrow = 2, ncol = 2,
#'   dimnames = list("regulators" = c("r1", "r2"), "targets" = c("t1", "t2"))
#' )
#' plot_upset(list(weightMat1, weightMat2), 2)
plot_upset <- function(weightMatList, plotConnections){
  if(length(weightMatList) > 1){
    edge.list <- lapply(weightMatList, function(wm){
      get_link_list_rename(wm, plotConnections) %>%
        dplyr::mutate(connection = paste0(.data$from, "-", .data$to)) %>%
        dplyr::pull(.data$connection)
    })
    names(edge.list) <- paste0("Network", seq_len(length(edge.list)))
    UpSetR::upset(UpSetR::fromList(edge.list), order.by = "degree", keep.order = TRUE)
  }else{
    NULL
  }
}

