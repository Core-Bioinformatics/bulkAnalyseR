calculate_condition_mean_sd_per_gene <- function(expression.matrix, condition){
  if(length(unique(condition)) < 2) stop("Please select at least two states")
  
  if(!is.factor(condition)) condition <- factor(condition, levels = unique(condition))
  noise.threshold <- min(expression.matrix)
  
  tbl <- tibble::tibble(
    gene = rep(rownames(expression.matrix), each = length(levels(condition))),
    condition = rep(levels(condition), times = nrow(expression.matrix)),
    mean = NA,
    sd = NA
  )
  
  pb = utils::txtProgressBar(min = 0, max = nrow(tbl), initial = 0, style = 3) 
  for(i in seq_len(nrow(tbl))){
    vec <- expression.matrix[tbl$gene[i], condition == tbl$condition[i]]
    tbl$mean[i] <- mean(vec)
    tbl$sd[i] <- stats::sd(vec)
    if(i %% 1000 == 0 | i == nrow(tbl)) utils::setTxtProgressBar(pb, i)
  }
  tbl <- tbl %>% dplyr::mutate(sd = pmax(sd, noise.threshold))
}

determine_uds <- function(min1, max1, min2, max2){
  if(max1 < min2) return("U") else if(min1 > max2) return("D") else return("S")
}

make_pattern_matrix <- function(tbl, n_sd = 1){
  
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

make_heatmap_matrix <- function(tbl, genes = NULL){
  if(is.null(genes)) genes <- unique(tbl$gene)
  tbl %>%
    dplyr::select(-sd) %>%
    dplyr::filter(gene %in% genes) %>%
    tidyr::pivot_wider(names_from = "condition", values_from = "mean") %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()
}

plot_line_pattern <- function(
  tbl, 
  genes, 
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

    p <- ggplot(tbl, aes(x = condition, y = scale, group = gene, colour = gene)) +
      theme_minimal() +
      stat_summary(fun = sum, geom = "line")
    
    if(!show.legend) p <- p + theme(legend.position = "none")
    
    p
        
}
