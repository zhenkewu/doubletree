#' Fit nested latent class model with double-tree-structured shrinkage over
#' 1) domains and 2) causes
#'
#' @param dsgn a list; data organized according to the trees
#' (see [design_doubletree()])
#' @param vi_params_init,hyperparams_init,random_init,random_init_vals,tol,tol_hyper,max_iter,print_freq,quiet,plot_fig,update_hyper_freq,hyper_fixed
#' initial values and updating protocols. Explained more in the wrapper function
#' [nlcm_doubletree()]
#' @return a list with model updates. Because the variational posterior is
#' comprised of familiar distributional forms that can be determined by the moments,
#' the returned values are these moments:
#' \describe{
#' \item{`vi_params`}{named list of final variational parameter estimates}
#' \item{`hyperparams`}{named list of final hyperparameter estimates}
#' \item{`hyper_fixed`}{named list of fixed hyperparameters}
#' \item{`ELBO_track`}{numeric vector containing the values of the objective function
#' (ELBO) at the end of every iteration.}
#' }
#'
#' @importFrom graphics barplot image abline
#' @family internal VI functions
#' @export
fit_nlcm_doubletree <- function(dsgn,
                          vi_params_init,
                          hyperparams_init,
                          random_init,
                          random_init_vals,
                          tol,
                          tol_hyper,
                          max_iter,
                          print_freq,
                          quiet,
                          plot_fig,
                          update_hyper_freq,
                          hyper_fixed
){
  ## set up the design list: 'dsgn' so this can be passed to initialization or updating functions.
  ## At initialization, we opted for parsimony, so `dsgn` does not
  ## have all the dimension information that can be calculated. We first
  ## add these additional elements to 'dsgn'.

  #-------------------------------BEGIN DESIGN PADDING--------------------------#
  dsgn$X <- 2*dsgn$Y-1 # 1, -1 or missing. The rows have been reordered by design_doubletree().
  # obtain the item ids with NON-missing responses for each subject i:
  dsgn$ind_obs_i <- mapply(FUN = function(v) which(!is.na(v)),
                           split_along_dim(dsgn$X,1),SIMPLIFY=FALSE)
  # obtain the subject ids with NON-missing responses for each item j:
  dsgn$ind_obs_j <- mapply(FUN = function(v) which(!is.na(v)),
                           split_along_dim(dsgn$X,2),SIMPLIFY=FALSE)
  dsgn$n <- nrow(dsgn$X)
  dsgn$J <- ncol(dsgn$X)
  dsgn$p1 <- length(unique(unlist(dsgn$ancestors[[1]])))
  dsgn$p2 <- length(unique(unlist(dsgn$ancestors[[2]])))
  dsgn$pL1 <- length(dsgn$ancestors[[1]])
  dsgn$pL2 <- length(dsgn$ancestors[[2]])
  dsgn$Fg1 <- max(dsgn$levels[[1]])
  dsgn$Fg2 <- max(dsgn$levels[[2]])

  if (is.null(hyper_fixed$K)){stop("[doubletree] # of classes 'K' not specified.")}
  if (!is.null(hyper_fixed$K) && hyper_fixed$K<2){cat("[double] # of classes 'K' is 1; assuming conditional independence.")}
  K  <- hyper_fixed$K
  # number of ancestors for each leaf, in tree1, and tree2:
  dsgn$cardanc1      <- unlist(lapply(dsgn$ancestors[[1]],length))
  dsgn$cardanc2      <- unlist(lapply(dsgn$ancestors[[2]],length))
  #-------------------------------END OF DESIGN PADDING--------------------------#

  if (!quiet){cat("\n [doubletree] working weights (edge lengths): `h_pau`: \n");print(dsgn$h_pau)}
  # initialize: ----------------------
  init <- R.utils::doCall(initialize_nlcm_doubletree,
                          vi_params   = vi_params_init,
                          hyperparams = hyperparams_init,
                          hyper_fixed = hyper_fixed,
                          random_init = random_init,
                          random_init_vals = random_init_vals,
                          args = c(dsgn))
  vi_params   <- init$vi_params
  hyperparams <- init$hyperparams
  dsgn$target_id <- init$target_id
  dsgn$scenario <- init$scenario
  cat("\n|--- Model Initialized.\n")

  # initialize ELBO:
  ELBO_track <- numeric(max_iter)
  line_track <- vector("list",max_iter+1)
  line_track[[1]] <- rep(0,17)

  # run algorithm: ---------------------
  i <- 0
  repeat{
    #if (quiet){pb$tick()} #;Sys.sleep(3 / 100)}
    # iterate i
    i <- i + 1


    # check if max_iter reached:
    if (i > max_iter){
      i <- max_iter
      cat(paste("|--- Iteration", i, "complete. \n"))
      warning("[doubletree] Maximum number of iterations reached! Consider increasing 'max_iter'")
      break
    }

    # update vi params:
    vi_params <- R.utils::doCall(update_vi_params_doubletree,
                                 args = c(dsgn, vi_params, hyperparams,
                                          hyper_fixed))

    # compute ELBO and update psi, phi and hyperparameters (tau_1, tau_2):
    update_hyper <- i %% update_hyper_freq == 0
    hyperparams  <- R.utils::doCall(update_hyperparams_doubletree,
                                    update_hyper = update_hyper,
                                    quiet      = quiet,
                                    args = c(dsgn,vi_params,hyperparams,hyper_fixed))
    ELBO_track[i] <- hyperparams$ELBO
    line_track[[i+1]] <- hyperparams$line_vec

    # print progress:
    if (i %% print_freq ==0){
      #if(ELBO_track[i] - ELBO_track[i-1]<0){
      if (!quiet){
        cat("|--- Iteration", i, "; >>> epsilon = ", ELBO_track[i] - ELBO_track[i-1], "<<<<; ELBO = ", ELBO_track[i],"\n")
        cat("|", i, "; line_vec_delta = \n")
        print(line_track[[i+1]]-line_track[[i]])
        cat("> empirical class probabilities: ", round(colMeans(vi_params$rmat),4),"\n")
        cat("> node_select: ",which(vi_params$prob1>0.5),"\n")
      }
      if (plot_fig){
        barplot(vi_params$prob1)
        abline(h=0.5,col="purple",lty=2)
      }
    }

    # check tolerance
    if (update_hyper & i >= 2 * update_hyper_freq) {
      # if we just updated hyperparameters, check for convergence of hyperparameters
      criterion1 <- abs(ELBO_track[i] - ELBO_track[i - update_hyper_freq]) < tol_hyper
      if (criterion1) {
        # did last VI update reach convergence?
        criterion2 <- abs(ELBO_track[i - 1] - ELBO_track[i - 2]) < tol
        # if yes, both have converged. if not, continue.
        if (criterion2) break else next
      } else next
    } else {
      criterion3 <- (i > 2) && (abs(ELBO_track[i] - ELBO_track[i - 1]) < tol)
      # if criterion 3, fill in results until just before the
      # next hyperparameter update (or max_iter, whichever comes first)
      if (criterion3) { # is the current update good enough?
        # if (i<2*update_hyper_freq){ # if update_hyper, but not yet 2*update_hyper_freq:
        #   # criterion4 <- (abs(ELBO_track[i] - init$hyperparams$ELBO) < tol_hyper) |
        #   #   (abs(ELBO_track[i] - ELBO_track[1]) < tol_hyper)
        #   # if (criterion4)
        #   break
        # }
        i2 <- min(max(ceiling(i / update_hyper_freq) * update_hyper_freq - 1,
                          2 * update_hyper_freq-1),
                  max_iter)
        ELBO_track[(i + 1):i2] <- ELBO_track[i]   # can send this iteration much later; so appears updating more frequent than specified.
        #ELBO_track[(i + 1):i2] <- hyperparams$ELBO  # can send this iteration much later; so appears updating more frequent than specified.
        i <- i2
      }
    }
  } # end 'repeat' loop.

  # return results:
  c(make_list(vi_params, hyperparams, hyper_fixed),
    list(ELBO_track=ELBO_track[1:i]),list(line_track = do.call("rbind",line_track[2:(i+1)])))

}
