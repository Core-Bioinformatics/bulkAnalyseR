test_that("DESeq2 normalisation removes sample-wide count scaling", {
  expression.matrix <- cbind(
    sample_A = c(10L, 20L, 40L),
    sample_B = c(40L, 80L, 160L)
  )
  rownames(expression.matrix) <- paste0("gene", seq_len(nrow(expression.matrix)))

  expect_warning(
    normalised <- preprocessExpressionMatrix(
      expression.matrix,
      denoise = FALSE,
      normalisation.method = "deseq2"
    ),
    "^Denoise was set to FALSE"
  )

  expected <- cbind(sample_A = c(20, 40, 80), sample_B = c(20, 40, 80))
  rownames(expected) <- rownames(expression.matrix)
  expect_equal(normalised, expected)
})

test_that("DESeq2 normalisation agrees with DESeq2 normalized counts", {
  expression.matrix <- rbind(
    gene1 = c(10L, 20L, 40L, 80L),
    gene2 = c(30L, 45L, 90L, 180L),
    gene3 = c(50L, 100L, 160L, 320L),
    gene4 = c(0L, 5L, 10L, 20L),
    gene5 = c(0L, 0L, 0L, 0L)
  )
  colnames(expression.matrix) <- c("sample_C", "sample_A", "sample_D", "sample_B")
  deseq <- DESeq2::DESeqDataSetFromMatrix(
    countData = expression.matrix,
    colData = data.frame(row.names = colnames(expression.matrix)),
    design = ~ 1
  )
  deseq <- DESeq2::estimateSizeFactors(deseq)

  expect_warning(
    normalised <- preprocessExpressionMatrix(
      expression.matrix,
      denoise = FALSE,
      normalisation.method = "deseq2"
    ),
    "^Denoise was set to FALSE"
  )

  expect_equal(normalised, DESeq2::counts(deseq, normalized = TRUE))
})
