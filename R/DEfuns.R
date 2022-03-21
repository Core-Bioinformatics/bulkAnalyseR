#' Perform differential expression (DE) analysis on an expression matrix
#' @description This function performs DE analysis on an expression using
#' edgeR or DESeq2, given a vector of sample conditions.
#' @inheritParams DEpanel
#' @param condition a vector of the same length as the number of columns of
#' expression.matrix, containing the sample conditions; this is usually the
#' last column of the metadata
#' @param var1,var2 conditions (contained in condition) to perform DE between;
#' note that DESeq2 requires at least two replicates per condition
#' @return A tibble with the differential expression results for all genes.
#' Columns are 
#' * gene_id (usually ENSEMBL ID matching one of the rows of the
#' expression matrix)
#' * gene_name (name matched through the annotation)
#' * log2exp (average log2(expression) of the gene across samples)
#' * log2FC (log2(fold-change) of the gene between conditions)
#' * pval (p-value of the gene being called DE)
#' * pvalAdj (adjusted p-value using the Benjamini Hochberg correction)
#' @examples
#' expression.matrix.preproc <- as.matrix(read.csv(
#'   system.file("extdata", "expression_matrix_preprocessed.csv", package = "bulkAnalyseR"), 
#'   row.names = 1
#' ))[1:100, 1:4]
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
#' edger <- DEanalysis_edger(
#'   expression.matrix = expression.matrix.preproc,
#'   condition = rep(c("0h", "12h"), each = 2),
#'   var1 = "0h",
#'   var2 = "12h",
#'   anno = anno
#' )
#' deseq <- DEanalysis_edger(
#'   expression.matrix = expression.matrix.preproc,
#'   condition = rep(c("0h", "12h"), each = 2),
#'   var1 = "0h",
#'   var2 = "12h",
#'   anno = anno
#' )
#' # DE genes with log2(fold-change) > 1 in both pipelines
#' intersect(
#'   dplyr::filter(edger, abs(log2FC) > 1, pvalAdj < 0.05)$gene_name,
#'   dplyr::filter(deseq, abs(log2FC) > 1, pvalAdj < 0.05)$gene_name
#' )
#' @name DEanalysis
NULL

#' @rdname DEanalysis
#' @export
DEanalysis_edger <- function(
  expression.matrix,
  condition,
  var1,
  var2,
  anno
){
  expression.matrix <- # Remove genes of constant expression
    expression.matrix[matrixStats::rowMins(expression.matrix) != 
                        matrixStats::rowMaxs(expression.matrix), ]
  
  condition <- factor(condition, levels = unique(condition))  
  design <- stats::model.matrix(~ 0 + condition)
  
  edger <- edgeR::DGEList(counts = expression.matrix, group = condition)
  edger <- edgeR::estimateDisp(edger, design)
  if(is.na(edger$common.dispersion)){ # if there are no replicates
    edger.noreps <- edgeR::DGEList(counts = expression.matrix, group = condition)
    edger.noreps <- edgeR::calcNormFactors(edger.noreps)
    edger.noreps <- edgeR::estimateGLMCommonDisp(edger.noreps,method="deviance",robust=TRUE,subset=NULL)
    edger$common.dispersion <- edger.noreps$common.dispersion
    edger$trended.dispersion <- edger.noreps$trended.dispersion
    edger$tagwise.dispersion <- edger.noreps$tagwise.dispersion
  }
  edger.fit = edgeR::glmFit(y = edger, design = design)
  if(var1 < var2) contrast <- c(-1, 1) else contrast <- c(1, -1)
  edger.lrt <- edgeR::glmLRT(edger.fit, contrast = contrast)
  edger.table <- edger.lrt$table
  
  gene_id <- NULL; pval <- NULL
  output = tibble::tibble(
    gene_id = rownames(expression.matrix),
    gene_name = anno$NAME[match(gene_id, anno$ENSEMBL)],
    log2exp = log2(rowMeans(expression.matrix)),
    log2FC = edger.table$logFC,
    pval = edger.table$PValue,
    pvalAdj = stats::p.adjust(pval, method = 'BH'),
  )
}

#' @rdname DEanalysis
#' @export
DEanalysis_deseq2 <- function(
  expression.matrix,
  condition,
  var1,
  var2,
  anno
){
  
  expression.matrix <- round(expression.matrix)
  expression.matrix <- # Remove genes of constant expression
    expression.matrix[matrixStats::rowMins(expression.matrix) != 
                        matrixStats::rowMaxs(expression.matrix), ]
  
  condition <- as.factor(condition)  
  design <- stats::model.matrix(~ 0 + condition)
  
  deseq <- DESeq2::DESeqDataSetFromMatrix(expression.matrix, data.frame(condition), design)
  DESeq2::sizeFactors(deseq) <- stats::setNames(rep(1, ncol(expression.matrix)), 
                                                colnames(expression.matrix))
  deseq <- DESeq2::DESeq(deseq)
  if(var1 < var2) contrast <- c(-1, 1) else contrast <- c(1, -1)
  deseq.res <- DESeq2::results(deseq, contrast = contrast)
  
  gene_id <- NULL; pval <- NULL
  output <- tibble::tibble(
    gene_id = rownames(expression.matrix),
    gene_name = anno$NAME[match(gene_id, anno$ENSEMBL)],
    log2exp = log2(rowMeans(expression.matrix)),
    log2FC = deseq.res$log2FoldChange,
    pval = deseq.res$pvalue,
    pvalAdj = stats::p.adjust(pval, method = 'BH'),
  )
}