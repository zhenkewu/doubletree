if(getRversion() >= "2.15.1") utils::globalVariables(c("i"))

#' wrapper function for fitting and summaries
#'
#' @inheritParams design_doubletree
#' @param ci_level A number between 0 and 1 giving the desired credible interval.
#' For example, `ci_level = 0.95` (the default) returns a 95% credible interval
#' @param get_lcm_by_group If `TRUE`, `doubletree` will also return the maximum likelihood estimates of the
#' coefficients for each leaf_ids group discovered by the model.
#' Default is `TRUE`.
#' @param update_hyper_freq How frequently to update hyperparameters.
#' Default = every 50 iterations.
#' @param print_freq How often to print out iteration number and current value of epsilon
#' (the difference in objective function value for the two most recent iterations).
#' @param quiet default to `FALSE`, which prints empirical class probabilities and updates on
#' tau's
#' @param plot_fig plot figure about `prob` (the probability of each node diffuse from
#' the parent node, i.e., s_u=1 for using the slab component) and response profile (1st node)
#' @param hyper_fixed Fixed values of hyperprior parameters.
#' @param tol Convergence tolerance for the objective function.
#' Default is `1E-8`.
#' @param tol_hyper The convergence tolerance for the objective function
#' between subsequent hyperparameter updates. Typically it is a more generous
#' tolerance than `tol`. Default is `1E-4`.
#' @param max_iter Maximum number of iterations of the VI algorithm.
#' Default is `5000`. NB: check this number before package submission.
#' @param nrestarts Number of random re-starts of the VI algorithm.
#' The restart that gives the highest value of the objective function will
#' be returned. It is recommended to choose `nrestarts > 1`; The default is `3`.
#' @param keep_restarts If `TRUE`, the results from all random restarts
#' will be returned. If `FALSE`, only the restart with the highest objective function is returned. '
#' Default is `TRUE`.
#' @param parallel If `TRUE`, the random restarts will be run in parallel.
#' It is recommended to first set the number of cores using `doParallel::registerDoParallel()`.
#' Otherwise, the default number of cores specified by the `doParallel` package will be used.
#' Default is `TRUE`.
#' @param log_restarts If `TRUE`, when `nrestarts > 1` progress of each random
#' restart will be logged to a text file in `log_dir`. If `FALSE` and `nrestarts > 1`,
#' progress will not be shown.
#' If `nrestarts = 1`, progress will always be printed to the console.
#' Default is `FALSE`.
#' @param log_dir Directory for logging progress of random restarts.
#' Default is the working directory.
#' @param vi_params_init,hyperparams_init Named lists containing initial values for the
#' variational parameters and hyperparameters. Supplying good initial values can be challenging,
#' and `lotR()` provides a way to guess initial values based on transformations
#' of latent class model estimates for each individual leaf_ids (see [initialize_tree_lcm()]).
#' The most common use for `vi_params_init` and `hyperparams_init` is to supply starting
#' values based on previous output from `lotR()`;
#' see the `vignette('lotR')` for examples.
#' The user can provide initial values for all parameters or a subset.
#' When initial values for one or more parameters are not
#' supplied, the missing values will be filled in by [initialize_nlcm_doubletree()].
#' @param random_init
#' If `TRUE`, some random variability will be added to the initial values.
#' The default is `FALSE`, unless `nrestarts > 1`, in which case
#' `random_init` will be set to `TRUE` and a warning message will be printed.
#' The amount of variability is determined by `random_init_vals`.
#' @param random_init_vals If `random_init = TRUE`,
#' this is a list containing the following parameters for randomly permuting
#' the initial values.
#' NB: The following are copied from lotR; so need edits!!!!!!!
#' \describe{
#' \item{`tau_lims`}{a vector of length `2`, where `tau_lims[1]` is between `0` and `1`,
#' and `tau_lims[2] > 1`. The initial values for the hyperparameter `tau` will
#' be chosen uniformly at random in the range `(tau_init * tau_lims[1], tau_init * tau_lims[2])`,
#' where `tau_init` is the initial value for `tau` either supplied in `hyperparams_init`
#' or guessed using [initialize_nlcm_doubletree()].}
#' \item{`psi_sd_frac`}{a value between `0` and `1`. The initial values for the auxiliary parameters
#' `psi` will have a normal random variate added to them with standard deviation equal to
#' `psi_sd_frac` multiplied by the initial value for eta either supplied in `hyperparams_init` or guessed
#' using [initialize_nlcm_doubletree()]. Absolute values are then taken for any
#' values of `psi` that are `< 0`.}
#' \item{`phi_sd_frac`}{same as above}.
#' \item{`mu_gamma_sd_frac`}{a value between 0 and 1. The initial values for
#' `mu` will have a normal random variate added to them with standard deviation equal to
#' `mu_sd_frac` multiplied by the absolute value of the initial value for `mu_gamma_sd_frac` either supplied in
#' `vi_params_init` or guessed using [initialize_nlcm_doubletree()].}
#' \item{`mu_alpha_sd_frac`}{same as above.}
#' \item{`u_sd_frac`}{a value between 0 and 1. The initial value for the node inclusion probabilities
#' will first be transformed to the log odds scale to obtain `u`. A normal random variate will be
#' added to `u` with standard deviation equal to u_sd_frac multiplied by the absolute value of the
#' initial value for `u` either supplied in `vi_params_init` or guessed using `moretrees_init_logistic()`.
#' `u` will then be transformed back to the probability scale.}
#' }
#'
#' @param allow_continue logical, `TRUE` to save results so can continue running the VI
#' updates with the last iteration from the old results.
#'
#' @return a list also of class "nlcm_doubletree"; NB: need to create a simulated example that uses this function!
#'
#' \describe{
#'   res <- make_list(mod,mod_restarts,mytrees,dsgn,prob_est,est_ad_hoc)
#'      class(res) <- c("nlcm_doubletree","list")
#'    }
#'
#' @example
#' /inst/example/example_simulate_doubletree.R
#'
#' @useDynLib doubletree
#' @export
#' @family nlcm_doubletree functions
nlcm_doubletree <- function(Y,leaf_ids,mytrees,# may have unordered nodes.
                            weighted_edges = c(TRUE,TRUE),
                            ci_level = 0.95,
                            get_lcm_by_group = FALSE,
                            update_hyper_freq = 50,
                            print_freq = 10,
                            quiet      = FALSE,
                            plot_fig   = FALSE,
                            hyper_fixed = list(K=2,LD=TRUE),
                            tol        = 1E-8,
                            tol_hyper = 1E-4,
                            max_iter = 5000,
                            nrestarts = 3,
                            keep_restarts = TRUE,
                            parallel = TRUE,
                            log_restarts = FALSE,
                            log_dir = ".",
                            vi_params_init = list(),
                            hyperparams_init = list(),
                            random_init = FALSE,
                            random_init_vals = list(mu_gamma_sd_frac = 0.2,
                                                    mu_alpha_sd_frac = 0.2,
                                                    tau1_lims = c(0.5,1.5),
                                                    tau2_lims = c(0.5,1.5),
                                                    u_sd_frac = 0.2, # this is for logit of prob1.
                                                    psi_sd_frac = 0.2,
                                                    phi_sd_frac = 0.2),
                            allow_continue = FALSE
){
  # logs
  log_dir <- sub("/$", "", log_dir)
  if (log_restarts) message("[doubletree] Algorithm progress for restart i will be printed to ",
                            log_dir, "/restart_i_log.txt\n", sep = "")

  # Fill in some arguments
  if (nrestarts > 1 & !random_init) { # force random starts if nstarts>1.
    message("[doubletree] Setting 'random_init = TRUE' since nrestarts > 1\n")
    random_init <- TRUE
  }
  if (nrestarts == 1) parallel <- FALSE

  # construct designed data; here `design_doubletree` reorders the nodes of the two trees, and the rows of the data.
  dsgn <- design_doubletree(Y,leaf_ids,mytrees,weighted_edges) # root_node,weighted_edge <--- need fixing.

  # Get hyper_fixed if not supplied:
  if (is.null(hyper_fixed$a1) | is.null(hyper_fixed$b1)) {
    L1             <- max(dsgn$levels[[1]])
    hyper_fixed   <- append(hyper_fixed,list(a1 = rep(1, L1)))
    hyper_fixed$b1 <- rep(10, L1)
    warning("[doubletree] No fixed hyperparameters (a1,b1) supplied; we set a*_l=1, b*_l=10 for all levels of hyperparameters in tree1.")
  }

  if (is.null(hyper_fixed$a2) | is.null(hyper_fixed$b2)) {
    L2             <- max(dsgn$levels[[2]])
    pL1            <- length(dsgn$leaf_ids_units[[1]])
    hyper_fixed   <- append(hyper_fixed,list(a2 = matrix(1, nrow=pL1,ncol=L2)))
    hyper_fixed$b2 <- matrix(10, nrow=pL1,ncol=L2)
    warning("[doubletree] No fixed hyperparameters (a2,b2) supplied; we set a_cl=1,b_cl=10 for all levels of hyperparameters in tree2.")
  }
  if (is.null(hyper_fixed$K)) {
    warning("[doubletree] No fixed # of classes supplied;
            supply a named element `K` in the list 'hyper_fixed'.")
  }

  # Setting up parallelization
  if (parallel) {
    `%doRestarts%` <- foreach::`%dopar%`
  } else {
    `%doRestarts%` <- foreach::`%do%`
  }

  # Run algorithm:
  mod_restarts <- foreach::foreach(i = 1:nrestarts) %doRestarts% {
    if (log_restarts) {
      sink(file = paste0(log_dir, "/restart_", i, "_log.txt"))
    }
    cat("\nInitializing restart", i, "...\n\n")
    #cat("Random initialization: ", random_init,"...\n\n")
    mod <- fit_nlcm_doubletree(dsgn = dsgn,
                               vi_params_init = vi_params_init,
                               hyperparams_init = hyperparams_init,
                               random_init = random_init,
                               random_init_vals = random_init_vals,
                               tol = tol,
                               tol_hyper = tol_hyper,
                               max_iter = max_iter,
                               print_freq = print_freq,
                               quiet      = quiet,
                               plot_fig   = plot_fig,
                               update_hyper_freq = update_hyper_freq,
                               hyper_fixed = hyper_fixed)
    cat("\nRestart", i, "complete.\n")
    if (log_restarts) {
      sink()
    }
    mod
  } # END `Run algorithm`.

  # Select random restart that gave the highest ELBO
  ELBO_restarts <- sapply(mod_restarts, FUN = function(mod) mod$ELBO_track[length(mod$ELBO_track)])
  best_restart  <- which.max(ELBO_restarts)
  mod <- mod_restarts[[best_restart]]
  if (keep_restarts) {
    mod_restarts <- mod_restarts[- best_restart]
  } else {
    rm(mod_restarts)
    mod_restarts <- NULL
  }

  #
  # The following needs edits as they either relate to posterior summaries or fit ad hoc models.
  # please see the lotR function for copying additional information.

  mytrees <- dsgn$mytrees # <--- check if this is actually needed - I think lotR likely reordered trees which must be documented before visualization.
  res <- make_list(mod,mod_restarts,mytrees,dsgn)

  class(res) <- c("nlcm_doubletree","list")
  res
}

