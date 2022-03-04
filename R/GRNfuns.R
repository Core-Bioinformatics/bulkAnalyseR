infer_GRN <- function(expression.matrix, metadata, anno, seed, 
                      targetGenes, condition, samples, inference_method){
  set.seed(seed)
  target.genes <- anno$ENSEMBL[match(targetGenes, anno$NAME)]
  samples <- metadata[[condition]] %in% samples
  if(inference_method == "GENIE3"){
    GENIE3::GENIE3(expression.matrix[, samples], targets = target.genes)
  }
}

get_link_list_rename <- function(weightMat, plotConnections){
  GENIE3::getLinkList(weightMat, plotConnections) %>%
    dplyr::mutate(from = as.character(.data$regulatoryGene), 
                  to = as.character(.data$targetGene), 
                  value = .data$weight, 
                  regulatoryGene = NULL, 
                  targetGene = NULL,
                  weight = NULL)
}

find_regulators_with_recurring_edges <- function(weightMatList, plotConnections){
  lapply(weightMatList, function(wm){
    get_link_list_rename(wm, plotConnections)[, 1:2]
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(.data$from, .data$to) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "keep") %>%
    dplyr::filter(.data$n > 1) %>%
    dplyr::pull(.data$from)
}

plot_GRN <- function(weightMat, anno, plotConnections, 
                     plot_position_grid, n_networks, recurring_regulators){
  
  if(n_networks >= plot_position_grid){
    edges <- get_link_list_rename(weightMat, plotConnections)
    nodes <- tibble::tibble(
      id = c(colnames(weightMat), edges$from),
      label = anno$NAME[match(.data$id, anno$ENSEMBL)],
      group = c(rep("target", ncol(weightMat)), 
                rep("regulator", nrow(edges))),
      color = NA
    ) %>%
      dplyr::distinct(.data$id, .keep_all = TRUE)
    nodes$group[nodes$id %in% recurring_regulators] <- "recurring_regulator"
    
    color_target <- list("background" = '#D2E5FF')
    color_regulator <- list("background" = '#E0E0E0')
    color_recurring_regulator <- list("background" = '#ACE9B4')
    colors.list <- list(color_target, color_regulator, color_recurring_regulator)
    
    for(i in seq_len(nrow(nodes))){
      nodes$color[i] <- colors.list[[match(
        nodes$group[i], c("target", "regulator", "recurring_regulator")
      )]]
    }
    # fix colours
    visNetwork::visNetwork(nodes, edges)
  }else{
    NULL
  }
}

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

