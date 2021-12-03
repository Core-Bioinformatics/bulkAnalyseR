% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QCplotFuns.R
\name{plot_pca}
\alias{plot_pca}
\title{Create a principal component analysis (PCA) plot the samples of an experiment}
\usage{
plot_pca(
  expression.matrix,
  metadata,
  annotation.id = NULL,
  n.abundant = NULL,
  show.labels,
  show.ellipses,
  label.force = 1
)
}
\arguments{
\item{expression.matrix}{the expression matrix; rows correspond to genes and
columns correspond to samples; usually preprocessed by
\code{\link{preprocessExpressionMatrix}}}

\item{metadata}{a data frame containing metadata for the samples contained
in the expression.matrix; must contain at minimum two columns:
the first column must contain the column names of the expression.matrix,
while the last column is assumed to contain the experimental conditions
that will be tested for differential expression}

\item{annotation.id}{a column index denoting which column of the metadata
should be used to colour the points and draw confidence ellipses}

\item{n.abundant}{number of most abundant genes to use for the JSI calculation}

\item{show.labels}{whether to label the points with the sample names}

\item{show.ellipses}{whether to draw confidence ellipses}

\item{label.force}{passed to the force argument of ggrepel::geom_label_repel;
higher values make labels overlap less (at the cost of them being further
away from the points they are labelling)}
}
\value{
The PCA plot as a ggplot object.
}
\description{
This function creates a PCA plot between all samples in the
expression matrix using the specified number of most abundant genes as
input. A metadata column is used as annotation.
}