gnet_modified <- function (input, reg_names, init_method = "boosting", init_group_num = 4, 
          max_depth = 3, cor_cutoff = 0.9, min_divide_size = 3, min_group_size = 2, 
          max_iter = 5, heuristic = TRUE, max_group = 0, force_split = 0.5, 
          nthread = 4) 
{
  if (is(input, class2 = "SummarizedExperiment")) {
    input <- assay(input)
  }
  if (max_group == 0) {
    max_group <- init_group_num
  }
  input <- input[apply(input, 1, var) > 1e-04, ]
  gene_data <- input[!rownames(input) %in% reg_names, , drop = FALSE]
  regulator_data <- input[reg_names, , drop = FALSE]
  result_all <- run_gnet(gene_data, regulator_data, init_method, 
                         init_group_num, max_depth, cor_cutoff, min_divide_size, 
                         min_group_size, max_iter, heuristic, max_group, force_split, 
                         nthread = nthread)
  reg_group_table <- result_all[[1]]
  gene_group_table <- result_all[[2]]
  group_idx <- unique(reg_group_table[, 1])
  reg_group_table_filtered <- gene_group_table_filtered <- NULL
  avg_cor_list <- c()
  current_group_idx <- 1
  regulators <- target_genes <- list()
  for (i in seq_len(length(group_idx))) {
    if (sum(gene_group_table[, 2] == group_idx[i]) >= min_group_size) {
      current_tree <- reg_group_table[reg_group_table[, 
                                                      1] == group_idx[i], ]
      current_tree[, 1] <- current_group_idx
      current_gene_group <- gene_group_table[gene_group_table$group == 
                                               group_idx[i], ]
      current_gene_group$group <- current_group_idx
      reg_group_table_filtered <- rbind(reg_group_table_filtered, 
                                        current_tree)
      gene_group_table_filtered <- rbind(gene_group_table_filtered, 
                                         current_gene_group)
      cor_m <- cor(t(gene_data[current_gene_group$gene, 
                               , drop = FALSE]))
      avg_cor_list <- c(avg_cor_list, mean(cor_m[upper.tri(cor_m)]))
      regulators[[i]] <- rownames(regulator_data)[
        current_tree[, 2] %% nrow(regulator_data) + 1 # CHANGE HERE
        ]
      target_genes[[i]] <- current_gene_group$gene
      current_group_idx <- current_group_idx + 1
    }
  }
  if (nrow(reg_group_table_filtered) <= 2) 
    warning("Too few modules generated, you may wish to try with higher cor_cutoff.")
  return(list(gene_data = gene_data, regulator_data = regulator_data, 
              group_score = avg_cor_list, reg_group_table = reg_group_table_filtered, 
              gene_group_table = gene_group_table_filtered, modules_count = current_group_idx - 
                1, regulators = regulators, target_genes = target_genes))
}
