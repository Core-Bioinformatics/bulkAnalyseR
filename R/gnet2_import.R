calc_likelihood_vec <- function(x,group_labels){
  score_total <- 0
  for(i in unique(group_labels)){
    if(sum(group_labels==i) > 1){
      x_i <- x[group_labels==i]
      score_total <- score_total+sum(dnorm(x_i,mean = mean(x_i),sd = sd(x_i),log = TRUE))
    }
  }
  return(score_total)
}

calc_correlation <- function(x){
  # remove any columns with 0 variance
  x <- x[,apply(x,2,var)>0]
  
  if(nrow(x)<=2 | ncol(x)<=2){
    # should not go further divide
    return(1)
  }else{
    x_cor <- cor(x)
    return(mean(abs(x_cor[upper.tri(x_cor)])))
  }
}

get_leaf_group_labels <- function(group_table,format_plot=FALSE){
  if(!format_plot){
    group_table <- group_table[,2:ncol(group_table),drop=FALSE]
  }
  leaf_label <- rep(0,ncol(group_table))
  next_label <- 0
  for(i in seq_len(nrow(group_table))){
    current_label <- group_table[i,]
    for(j in seq_len(2)-1){
      leaf_label[current_label==j] <- next_label
      next_label <- next_label+1
    }
  }
  return(leaf_label)
}

get_group_kmeans <- function(x){
  suppressWarnings({ km <- kmeans(x,centers = 3)})
  grp <- rep(0,length(x))
  o <- order(km$centers)
  grp[km$cluster==o[1]] <- -1
  grp[km$cluster==o[3]] <- 1
  return(grp)
}

build_moduleR_heuristic2 <- function(X,Y,max_depth,cor_cutoff,min_divide_size){
  feature_remaining <- seq_len(ncol(X))
  feature_num <- 0
  subgroup_indicator <- rep(0,nrow(X))
  groups <- c(0)
  group_table <- matrix(-1, nrow = max_depth, ncol = nrow(X)+1)
  while (feature_num < max_depth & length(groups)>0 & length(feature_remaining)>0){
    best_score <- (-(10^6))
    best_feature <- (-1)
    for(i in groups){
      current_split_group <- subgroup_indicator == i
      y_current_split <- Y[current_split_group,,drop=FALSE]
      for (j in feature_remaining){
        feature_vals <- X[current_split_group,j]
        score_nosplit <- calc_likelihood_score(y_current_split,rep(1,length(feature_vals)))
        subgroup_divide <- suppressWarnings({as.numeric(kmeans(feature_vals,centers = 2)$cluster)-1})
        if(mean(feature_vals[subgroup_divide==0]) > mean(feature_vals[subgroup_divide==1])){
          # 0 is always left and 1 is always right
          subgroup_divide <- 1- subgroup_divide
        }
        score <- calc_likelihood_score(y_current_split,subgroup_divide) - score_nosplit
        if (score > best_score){
          best_score <- score
          best_feature <- j
          subgroup_group_labels <- rep(-1,length(current_split_group))
          subgroup_group_labels[current_split_group] <- subgroup_divide
        }
      }
    }
    feature_num <- feature_num+1
    group_table[feature_num,] <- c(best_feature-1,subgroup_group_labels)
    subgroup_indicator_new <- get_leaf_group_labels(group_table)
    subgroup_indicator_new[subgroup_indicator==-1] <- (-1)
    subgroup_indicator <- subgroup_indicator_new
    for (i in seq_len(2)-1){
      new_divide_i <- subgroup_group_labels==i
      if (sum(new_divide_i) <= min_divide_size){
        subgroup_indicator[new_divide_i] <- (-1)
      }else if(calc_correlation(t(Y[new_divide_i,])) >= cor_cutoff){
        subgroup_indicator[new_divide_i] <- (-1)
      }
    }
    groups <- unique(subgroup_indicator[subgroup_indicator!=-1])
    feature_remaining <- feature_remaining[feature_remaining!=best_feature]
  }
  # remove unused rows
  group_table <- group_table[apply(group_table, 1, function(x)sum(x!=-1))>0,]
  return(group_table)
}

assign_regul <- function(regulator_data,gene_data,gene_group_table,min_group_size,
                         max_depth,cor_cutoff,min_divide_size,heuristic = TRUE,split_table = NULL){
  X <- t(regulator_data)
  group_group_labels <- unique(gene_group_table)
  group_group_labels <- group_group_labels[group_group_labels!=-1]
  reg_group_table <- matrix(-1,nrow = max(gene_group_table)+1,ncol = ncol(regulator_data))
  group_table <- matrix(-1,nrow = max_depth*(max(gene_group_table)+1),
                        ncol = ncol(regulator_data)+2)
  i <- 0
  group_table_start <- 1
  for (group_idx in group_group_labels){
    gene_idx <- gene_group_table == group_idx
    if (sum(gene_idx) >= min_group_size){
      Y <- t(gene_data[gene_idx,,drop=FALSE])
      if(heuristic){
        group_table_i <- build_moduleR_heuristic(X,Y,max_depth,cor_cutoff,min_divide_size,split_table)
      }else{
        group_table_i <- build_module(X,Y,max_depth,cor_cutoff,min_divide_size)
      }
      if(length(ncol(group_table_i))!=0){
        group_table_i <- group_table_i[apply(group_table_i, 1, function(x)(sum(x!=0)))!=0,]
        if(length(ncol(group_table_i))!=0){
          reg_group_table[i+1,] <- get_leaf_group_labels(group_table_i)
          group_table_end <- group_table_start+nrow(group_table_i)-1
          group_table[group_table_start:group_table_end,2:ncol(group_table)] <- group_table_i
          group_table[group_table_start:group_table_end,1] <- i
          i <- i+1
          group_table_start <- group_table_end+1
        }
      }
    }
  }
  # remove unused rows
  group_table <- group_table[apply(group_table, 1, function(x)sum(x!=-1))>0,]
  return(list(group_table,reg_group_table))
}

assign_gene <- function(gene_data,reg_group_table){
  gene_group_table <- rep(-1,nrow(gene_data))
  for(gene_idx in seq_len(nrow(gene_data))){
    exp_gene <- as.numeric(gene_data[gene_idx,])
    gene_group_table[gene_idx] <- which.max(
      apply(reg_group_table, 1, function(x)calc_likelihood_vec(exp_gene,x)))-1
  }
  return(gene_group_table)
}

assign_first_cluster <- function(gene_data,regulator_data,max_depth,init_group_num,
                                 init_method='boosting',max_group=5,nthread = 4){
  if(init_method=='boosting'){
    ipt_mat <- matrix(0,nrow = nrow(gene_data),ncol = nrow(regulator_data))
    rownames(ipt_mat) <- rownames(gene_data)
    colnames(ipt_mat) <- rownames(regulator_data)
    for(i in seq_len(nrow(gene_data))){
      dtrain <- xgb.DMatrix(data = t(regulator_data), label=gene_data[i,])
      bst <- xgb.train(data=dtrain, max.depth=max_depth, eta=0.1, nthread = nthread,nrounds =2)
      if(length(xgb.dump(model = bst, with_stats = TRUE))!=4){
        # A constant tree: no predictors were used. 
        importance_matrix <- xgb.importance(model = bst)
        ipt_mat[i,importance_matrix$Feature] <- importance_matrix$Gain
      }
    }
  }else{
    ipt_mat <- gene_data
  }
  avg_cor_list <- c()
  for(i in 2:min(init_group_num,nrow(gene_data)-1)){
    gene_group_table <- suppressWarnings({as.numeric(kmeans(ipt_mat,centers = i)$cluster)-1})
    avg_cor_list <- c(avg_cor_list,mean(get_correlation_list(t(gene_data),gene_group_table)))
  }
  if(length(avg_cor_list)>=3){
    o <- order(avg_cor_list,decreasing = TRUE)
    y <- avg_cor_list[o]
    kn <- kneepointDetection(y)
    groups_keep <- o[seq_len(max(kn,max_group))]
  }else{
    groups_keep <- seq_len(length(avg_cor_list))
  }
  gene_group_table[!(gene_group_table %in% groups_keep)] <- (-1)
  return(gene_group_table)
}

run_gnet <- function(gene_data,regulator_data,init_method = 'boosting',init_group_num = 5,max_depth = 3,
                     cor_cutoff = 0.9,min_divide_size = 3,min_group_size = 2,max_iter = 5,heuristic = TRUE,
                     max_group=5,force_split = 0.5,nthread = 4){
  message('Determining initial group number...')
  gene_group_table <- assign_first_cluster(gene_data,regulator_data,max_depth,
                                           init_group_num,init_method,max_group,nthread = nthread)
  split_table <- build_split_table(t(regulator_data))
  message('Building module networks...')
  for (i in seq_len(max_iter)) {
    message('Iteration ',i)
    assign_reg_names <- assign_regul(regulator_data,gene_data,gene_group_table,min_group_size,
                                     max_depth,cor_cutoff,min_divide_size,heuristic,split_table)
    gene_group_table_new <- assign_gene(gene_data,assign_reg_names[[2]])
    max_group_idx <- as.numeric(names(which.max(table(gene_group_table_new))))
    if(max(table(gene_group_table_new))>round(length(gene_group_table_new)*force_split)){
      gene_group_table_largest <- assign_first_cluster(gene_data[gene_group_table_new == max_group_idx,],
                                                       regulator_data,max_depth,init_group_num,init_method,max_group)
      gene_group_table_largest[gene_group_table_largest==-1] <- (-1-max(gene_group_table_new))
      gene_group_table_new[gene_group_table_new==max_group_idx] <- gene_group_table_largest+max(gene_group_table_new)
      gene_group_table_new <- as.numeric(as.factor(gene_group_table_new))
      group_to_remove <- as.numeric(names(table(gene_group_table_new))[table(gene_group_table_new) <= min_group_size])
      gene_group_table_new[gene_group_table_new %in% group_to_remove] <- -1
      gene_group_table_new <- as.numeric(as.factor(gene_group_table_new))-1
    }
    
    if(all(length(gene_group_table)==length(gene_group_table_new)) && all(gene_group_table==gene_group_table_new)){
      message('Converged.')
      break
    }else{
      gene_group_table <- gene_group_table_new
    }
  }
  message('Generating final network modules...')
  assign_reg_names <- assign_regul(regulator_data,gene_data,gene_group_table,min_group_size,
                                   max_depth,cor_cutoff,min_divide_size,heuristic,split_table)
  gene_group_table <- assign_gene(gene_data,assign_reg_names[[2]])
  message('Done.')
  tree_table_all <- assign_reg_names[[1]]
  colnames(tree_table_all) <- c('group','feature',colnames(gene_data))
  group_order <- order(gene_group_table)
  gene_list_all <- data.frame('gene'  = rownames(gene_data)[group_order],
                              'group' = gene_group_table[group_order])
  return(list('reg_group_table'=tree_table_all,'gene_group_table'=gene_list_all))
}

sum_scores <- function(el_all,el_input){
  scores <- rep(0,nrow(el_all))
  for (i in seq_len(nrow(el_all))) {
    idx <- (el_input[,1] == el_all[i,1] & el_input[,2] == el_all[i,2]) | (el_input[,2] == el_all[i,1] & el_input[,1] == el_all[i,2])
    scores[i] <- sum(abs(el_input[idx,3]))
  }
  return(scores)
}

expand_xy <- function(x,y,score){
  xy <- expand.grid(x,y,stringsAsFactors =FALSE)
  return(cbind.data.frame(xy,score,stringsAsFactors =F))
}

get_sample_cluster <- function(group_table,group_idx){
  group_table <- group_table[group_table[,1]==group_idx,3:ncol(group_table),drop=FALSE]
  leaf_label <- rep(0,ncol(group_table))
  next_label <- 0
  for(i in seq_len(nrow(group_table))){
    current_label <- group_table[i,]
    for(j in c(0,1)){
      leaf_label[current_label==j] <- next_label
      next_label <- next_label+1
    }
  }
  return(as.numeric(factor(leaf_label)))
}

ari <- function (x, y){
  tab <- table(x, y)
  if (all(dim(tab) == c(1, 1))) 
    return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b + 
                                                     a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  return(ARI)
}

similarity_score_unranked <- function(gnet_result,group){
  s_list<- rep(0,gnet_result$modules_count)
  for (i in seq_len(gnet_result$modules_count)) {
    gnet_groups <- get_sample_cluster(gnet_result$reg_group_table,i)
    s_list[i] <- ari(gnet_groups,group)
  }
  return(s_list)
}

compute_kld <- function(gnet_groups,group){
  group_module_dist <- as.matrix(dist(gnet_groups))
  group_dist <- as.matrix(dist(group))
  p <- dnorm(group_dist*group_dist,0,sqrt(max(group_dist)))
  q <- dnorm(group_module_dist*group_module_dist,0,sqrt(max(group_module_dist)))
  
  for (j in seq_len(nrow(group_dist))) {
    p[j,] <- p[j,]/(sum(p[j,])-p[j,])
    q[j,] <- q[j,]/(sum(q[j,])-q[j,])
  }
  
  p1 <- p[upper.tri(p)]
  q1 <- q[upper.tri(q)]
  
  return(sum(p1*log(p1/q1)))
}

similarity_score_ranked <- function(gnet_result,group){
  kld <- rep(0,gnet_result$modules_count)
  for (i in seq_len(gnet_result$modules_count)) {
    gnet_groups <- get_sample_cluster(gnet_result$reg_group_table,i)
    
    kld[i] <- compute_kld(gnet_groups,group)
  }
  return(kld)
}
