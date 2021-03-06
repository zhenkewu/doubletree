% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update.R
\name{update_hyperparams_doubletree}
\alias{update_hyperparams_doubletree}
\title{update hyperparameters}
\usage{
update_hyperparams_doubletree(
  h_pau,
  levels,
  v_units,
  X,
  n,
  J,
  p1,
  p2,
  pL1,
  pL2,
  Fg1,
  Fg2,
  ind_obs_i,
  prob1,
  prob2,
  mu_gamma,
  mu_alpha,
  rmat,
  emat,
  dirich_mat,
  sigma_gamma,
  sigma_alpha,
  tau_1_t,
  tau_2_t,
  a1_t,
  b1_t,
  a2_t,
  b2_t,
  E_beta_sq,
  E_eta_sq,
  E_beta,
  E_eta,
  psi,
  g_psi,
  phi,
  g_phi,
  tau_1,
  tau_2,
  a1,
  b1,
  a2,
  b2,
  dmat,
  K,
  LD,
  tau_update_levels,
  update_hyper,
  quiet,
  do_tree1_update
)
}
\arguments{
\item{h_pau}{outputs from \code{\link[=design_doubletree]{design_doubletree()}}. Among them,
\code{subject_id_list} is currently not used in the initialization
function \code{\link[=initialize_nlcm_doubletree]{initialize_nlcm_doubletree()}}, but are useful for updating functions here.}

\item{levels}{outputs from \code{\link[=design_doubletree]{design_doubletree()}}. Among them,
\code{subject_id_list} is currently not used in the initialization
function \code{\link[=initialize_nlcm_doubletree]{initialize_nlcm_doubletree()}}, but are useful for updating functions here.}

\item{v_units}{outputs from \code{\link[=design_doubletree]{design_doubletree()}}. Among them,
\code{subject_id_list} is currently not used in the initialization
function \code{\link[=initialize_nlcm_doubletree]{initialize_nlcm_doubletree()}}, but are useful for updating functions here.}

\item{X}{data and design information computed from the outputs
of \code{\link[=initialize_nlcm_doubletree]{initialize_nlcm_doubletree()}}; we use \code{R.utils::do.Call} to get their values from \code{dsgn}.}

\item{n}{data and design information computed from the outputs
of \code{\link[=initialize_nlcm_doubletree]{initialize_nlcm_doubletree()}}; we use \code{R.utils::do.Call} to get their values from \code{dsgn}.}

\item{J}{data and design information computed from the outputs
of \code{\link[=initialize_nlcm_doubletree]{initialize_nlcm_doubletree()}}; we use \code{R.utils::do.Call} to get their values from \code{dsgn}.}

\item{p1}{data and design information computed from the outputs
of \code{\link[=initialize_nlcm_doubletree]{initialize_nlcm_doubletree()}}; we use \code{R.utils::do.Call} to get their values from \code{dsgn}.}

\item{p2}{data and design information computed from the outputs
of \code{\link[=initialize_nlcm_doubletree]{initialize_nlcm_doubletree()}}; we use \code{R.utils::do.Call} to get their values from \code{dsgn}.}

\item{pL1}{data and design information computed from the outputs
of \code{\link[=initialize_nlcm_doubletree]{initialize_nlcm_doubletree()}}; we use \code{R.utils::do.Call} to get their values from \code{dsgn}.}

\item{pL2}{data and design information computed from the outputs
of \code{\link[=initialize_nlcm_doubletree]{initialize_nlcm_doubletree()}}; we use \code{R.utils::do.Call} to get their values from \code{dsgn}.}

\item{Fg1}{data and design information computed from the outputs
of \code{\link[=initialize_nlcm_doubletree]{initialize_nlcm_doubletree()}}; we use \code{R.utils::do.Call} to get their values from \code{dsgn}.}

\item{Fg2}{data and design information computed from the outputs
of \code{\link[=initialize_nlcm_doubletree]{initialize_nlcm_doubletree()}}; we use \code{R.utils::do.Call} to get their values from \code{dsgn}.}

\item{ind_obs_i}{outputs from \code{\link[=design_doubletree]{design_doubletree()}}. Among them,
\code{subject_id_list} is currently not used in the initialization
function \code{\link[=initialize_nlcm_doubletree]{initialize_nlcm_doubletree()}}, but are useful for updating functions here.}

\item{prob1}{variational parameters}

\item{prob2}{variational parameters}

\item{mu_gamma}{variational parameters}

\item{mu_alpha}{variational parameters}

\item{rmat}{variational parameters}

\item{emat}{variational parameters}

\item{dirich_mat}{variational parameters}

\item{sigma_gamma}{variational parameters}

\item{sigma_alpha}{variational parameters}

\item{tau_1_t}{variational parameters}

\item{tau_2_t}{variational parameters}

\item{a1_t}{variational parameters}

\item{b1_t}{variational parameters}

\item{a2_t}{variational parameters}

\item{b2_t}{variational parameters}

\item{E_beta_sq, E_eta_sq, E_beta, E_eta}{moments computed by \code{\link[=update_vi_params_doubletree]{update_vi_params_doubletree()}}}

\item{psi}{parameters updated by \code{update_hyperparams_doubletree()}}

\item{g_psi}{parameters updated by \code{update_hyperparams_doubletree()}}

\item{phi}{parameters updated by \code{update_hyperparams_doubletree()}}

\item{g_phi}{parameters updated by \code{update_hyperparams_doubletree()}}

\item{tau_1}{parameters updated by \code{update_hyperparams_doubletree()}}

\item{tau_2}{parameters updated by \code{update_hyperparams_doubletree()}}

\item{tau_update_levels}{a numeric vector, specifies which levels of hyperparameters to update}

\item{update_hyper}{Logical, \code{TRUE} or \code{FALSE} to indicate
whether to update \code{tau_1} and \code{tau_2}. This is computed at every iteration
in \code{\link[=fit_nlcm_doubletree]{fit_nlcm_doubletree()}}}

\item{quiet}{default to \code{FALSE}, which prints intermediate updates of hyperparameters}

\item{do_tree1_update}{\code{TRUE} to update \code{mu_gamma},\code{sigma_gamma}, \code{prob1},\code{tau_1}}
}
\value{
a list of updated hyperparameters: tau_1,tau_2,psi,g_psi,phi,g_phi,
along with a new ELBO value.
}
\description{
NB: argument redundancy may exist
}
\seealso{
Other Internal VI functions: 
\code{\link{update_vi_params_doubletree}()}
}
\concept{Internal VI functions}
