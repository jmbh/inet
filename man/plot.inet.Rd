\name{plot.inet}
\alias{plot.inet}
\title{
Plot point estimates and confidence intervals
}
\description{
Plot point estimates and confidence intervals for models estimated with the \code{lasso_ms}, \code{lasso_dsp}, \code{lasso_dsp_boot} and \code{OLS} functions.
}
\usage{
\method{plot}{inet}(x, labels = NULL, order = FALSE, subset = NULL,
          cex.labels = 0.80, cex.axis = 0.75, ...)
}
\arguments{
  \item{x}{
The output object from either \code{lasso_ms}, \code{lasso_dsp}, \code{lasso_dsp_boot} or \code{OLS}.
}
  \item{labels}{
A p-vector of characters specifying the labels for variables.
}
  \item{order}{
If \code{order = TRUE}, the edges are listed in decreasing order based on the point estimate.
}
  \item{subset}{
Allows to only display a subset of the edges. For example, if \code{subset=1:20} the first 20 edges are displayed. This is especially useful for larger networks, in which all edges are unlikely to fit into a single figure.
}
  \item{cex.labels}{
The font size of the edge labels.
}
  \item{cex.axis}{
The font size of the axes.
}
  \item{\dots}{
Additional arguments.
}
}

\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
}

\examples{

\dontrun{

# Fit GGM to PTSD data
set.seed(1)
out <- lasso_dsp(data = ptsd_data)

# Plot first 20 edges
plot(out, labels = colnames(ptsd_data),
     order=TRUE, subset = 1:20)

}

}
