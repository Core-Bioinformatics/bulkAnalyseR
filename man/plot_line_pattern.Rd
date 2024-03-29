% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/patternFuns.R
\name{plot_line_pattern}
\alias{plot_line_pattern}
\title{Create a line plot of average expression across conditions}
\usage{
plot_line_pattern(
  tbl,
  genes = NULL,
  type = c("Mean Scaled", "Log2 Expression", "Expression"),
  show.legend = FALSE
)
}
\arguments{
\item{tbl}{the output of \code{\link{calculate_condition_mean_sd_per_gene}}}

\item{genes}{gene names to use for the output; if NULL (the default),
all genes will be used}

\item{type}{whether the expression values should be scaled using their mean
(default), log-transformed, or not adjusted for the plot}

\item{show.legend}{whether to show the gene names in the legend;
should be avoided in many genes are plotted}
}
\value{
A matrix of average gene expression per gene in each condition.
}
\description{
This function creates a line plot of average expression
across conditions for a selection of genes, usually to visualise
an expression pattern.
}
\examples{
expression.matrix.preproc <- as.matrix(read.csv(
  system.file("extdata", "expression_matrix_preprocessed.csv", package = "bulkAnalyseR"), 
  row.names = 1
))[1:500,]

condition <- factor(rep(c("0h", "12h", "36h"), each = 2))
tbl <- calculate_condition_mean_sd_per_gene(expression.matrix.preproc[1:10, ], condition)
plot_line_pattern(tbl)
}
