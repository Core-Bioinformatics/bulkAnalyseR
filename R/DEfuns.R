DEanalysis_edger <- function(
  expression.matrix,
  condition,
  var1,
  var2,
  org
){
  expression.matrix <- # Remove genes of constant expression
    expression.matrix[matrixStats::rowMins(expression.matrix) != 
                        matrixStats::rowMaxs(expression.matrix), ]
  anno <- AnnotationDbi::select(
    getExportedValue(org, org),
    keys = rownames(expression.matrix),
    keytype = "ENSEMBL",
    columns = "SYMBOL"
  ) %>%
    dplyr::distinct(ENSEMBL, .keep_all = TRUE) %>%
    dplyr::mutate(NAME = ifelse(is.na(SYMBOL), ENSEMBL, SYMBOL))
    
  group <- c(rep(0, sum(condition == var1)), rep(1, sum(condition == var2)))
  design <- stats::model.matrix(~ 0 + as.factor(condition))
  
  edger <- edgeR::DGEList(counts = expression.matrix, group = group)
  edger <- edgeR::estimateDisp(edger, design)
  if(is.na(edger$common.dispersion)){ # if there are no replicates
    edger.noreps <- edgeR::DGEList(counts = expression.matrix, group = group)
    edger.noreps <- edgeR::calcNormFactors(edger.noreps)
    edger.noreps <- edgeR::estimateDisp(edger.noreps, tagwise=FALSE)
    edger$common.dispersion <- edger.noreps$common.dispersion
    edger$trended.dispersion <- edger.noreps$trended.dispersion
    edger$tagwise.dispersion <- edger.noreps$tagwise.dispersion
  }
  edger.fit = edgeR::glmFit(y = edger, design = design)
  edger.lrt <- edgeR::glmLRT(edger.fit, contrast = c(-1, 1))
  edger.table <- edger.lrt$table
  
  output = tibble::tibble(
    gene_id = rownames(expression.matrix),
    gene_name = anno$NAME[match(gene_id, anno$ENSEMBL)],
    avgExp = log2(rowMeans(expression.matrix)),
    log2FC = edger.table$logFC,
    pval = edger.table$PValue,
    pvalAdj = stats::p.adjust(pval, method = 'BH'),
  )}