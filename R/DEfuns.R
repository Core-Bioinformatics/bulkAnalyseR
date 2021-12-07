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
#' @export
#' @examples
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
  
  condition <- as.factor(condition)  
  design <- stats::model.matrix(~ 0 + condition)
  
  edger <- edgeR::DGEList(counts = expression.matrix, group = condition)
  edger <- edgeR::estimateDisp(edger, design)
  if(is.na(edger$common.dispersion)){ # if there are no replicates
    edger.noreps <- edgeR::DGEList(counts = expression.matrix, group = condition)
    edger.noreps <- edgeR::calcNormFactors(edger.noreps)
    edger.noreps <- edgeR::estimateDisp(edger.noreps, tagwise=FALSE)
    edger$common.dispersion <- edger.noreps$common.dispersion
    edger$trended.dispersion <- edger.noreps$trended.dispersion
    edger$tagwise.dispersion <- edger.noreps$tagwise.dispersion
  }
  edger.fit = edgeR::glmFit(y = edger, design = design)
  if(var1 < var2) contrast <- c(-1, 1) else contrast <- c(1, -1)
  edger.lrt <- edgeR::glmLRT(edger.fit, contrast = contrast)
  edger.table <- edger.lrt$table
  
  output = tibble::tibble(
    gene_id = rownames(expression.matrix),
    gene_name = anno$NAME[match(gene_id, anno$ENSEMBL)],
    log2exp = log2(rowMeans(expression.matrix)),
    log2FC = edger.table$logFC,
    pval = edger.table$PValue,
    pvalAdj = stats::p.adjust(pval, method = 'BH'),
  )}

#' @rdname DEanalysis
#' @export
DEanalysis_deseq <- function(
  expression.matrix,
  condition,
  var1,
  var2,
  anno
)