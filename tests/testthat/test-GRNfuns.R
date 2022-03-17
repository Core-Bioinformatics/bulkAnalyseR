# expression.matrix <- as.matrix(read.csv(
#   system.file("extdata", "expression_matrix_preprocessed.csv", package = "bulkAnalyseR"),
#   row.names = 1
# ))
# 
# metadata <- data.frame(
#   srr = colnames(expression.matrix),
#   timepoint = rep(c("0h", "12h", "36h"), each = 2)
# )
# 
# suppressMessages(
#   anno <- AnnotationDbi::select(
#     getExportedValue('org.Mm.eg.db', 'org.Mm.eg.db'),
#     keys = rownames(expression.matrix),
#     keytype = 'ENSEMBL',
#     columns = 'SYMBOL'
#   ) %>%
#     dplyr::distinct(ENSEMBL, .keep_all = TRUE) %>%
#     dplyr::mutate(NAME = ifelse(is.na(SYMBOL), ENSEMBL, SYMBOL))
# )
# 
# DEresults <- DEanalysis_edger(
#   expression.matrix = expression.matrix[, 1:4], 
#   condition = metadata$timepoint[1:4], 
#   var1 = "0h", 
#   var2 = "12h", 
#   anno = anno
# ) %>%
#   dplyr::filter(abs(log2FC) > 1, pvalAdj < 0.05) %>%
#   dplyr::arrange(dplyr::desc(abs(log2FC)))
# DEgenes <- DEresults$gene_id
# expression.matrix <- expression.matrix[DEgenes, ]
# 
# weightMat1 <- infer_GRN(expression.matrix, metadata, anno, 13,
#                         c("Trpm1", "Sp5"), "timepoint", c("0h", "12h", "36h"), "GENIE3")
# weightMat2 <- infer_GRN(expression.matrix, metadata, anno, 13,
#                         c("Trpm1", "Sp5"), "timepoint", c("0h", "12h"), "GENIE3")
# weightMat3 <- infer_GRN(expression.matrix, metadata, anno, 13,
#                         c("Trpm1", "Sp5"), "timepoint", c("0h", "36h"), "GENIE3")
# weightMat4 <- infer_GRN(expression.matrix, metadata, anno, 13,
#                         c("Trpm1", "Sp5"), "timepoint", c("12h", "36h"), "GENIE3")
# weightMatGNET <- infer_GRN(expression.matrix, metadata, anno, 13,
#                            c("Trpm1", "Sp5", "Nupr1", "Dnmt3l", "Enox1", "Itga1"), 
#                            "timepoint", c("0h", "12h", "36h"), "GNET2")
# 
# plotConnections <- 5
# 
# recurring_targets_1 <- find_targets_with_recurring_edges(list(weightMat1), plotConnections)
# recurring_targets <- find_targets_with_recurring_edges(
#   list(weightMat1, weightMat2, weightMat3, weightMat4), plotConnections
# )
# 
# # plot_GRN(weightMat1, anno, plotConnections, 1, 4, recurring_targets)
# # plot_GRN(weightMat2, anno, plotConnections, 1, 4, recurring_targets)
# # plot_GRN(weightMat3, anno, plotConnections, 1, 4, recurring_targets)
# # plot_GRN(weightMat4, anno, plotConnections, 1, 4, recurring_targets)
# 
# test_that("recurring regulators calculation works", {
#   expect_equal(recurring_regulators_1, character(0))
#   expect_equal(recurring_regulators, "ENSMUSG00000028766")
# })
# 
# 
