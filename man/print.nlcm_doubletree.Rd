% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarize_results.R
\name{print.nlcm_doubletree}
\alias{print.nlcm_doubletree}
\title{\code{print.nlcm_doubletree} summarizes the results from \code{\link[=nlcm_doubletree]{nlcm_doubletree()}}.}
\usage{
\method{print}{nlcm_doubletree}(x, ...)
}
\arguments{
\item{x}{Output from \code{\link[=nlcm_doubletree]{nlcm_doubletree()}}.}

\item{...}{Arguments passed to summary and printing methods.}
}
\value{
Summary showing, for each group of leaf nodes discovered by \code{\link[=nlcm_doubletree]{nlcm_doubletree()}},
the class prevalences, 95\% credible intervals, number of leaf nodes in
per group, and number of observations per group.
}
\description{
\code{print.nlcm_doubletree} summarizes the results from \code{\link[=nlcm_doubletree]{nlcm_doubletree()}}.
}
\concept{lcm_tree results}
