% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarize_results.R
\name{compute_params_dt}
\alias{compute_params_dt}
\title{compute model summaries from the model outputs}
\usage{
compute_params_dt(mod, dsgn, ci_level = 0.95)
}
\arguments{
\item{mod}{output from \code{\link[=fit_nlcm_doubletree]{fit_nlcm_doubletree()}}}

\item{dsgn}{output from \code{\link[=design_doubletree]{design_doubletree()}}}

\item{ci_level}{credible interval level; default to \code{0.95}}
}
\value{
a list with elements: \code{prob_est},\code{prob_est_indiv},\code{pi_list}
}
\description{
This function is used in the wrapper function \code{\link[=nlcm_doubletree]{nlcm_doubletree()}}
}
\seealso{
\code{\link[=get_est_cpp_dt]{get_est_cpp_dt()}}
}
