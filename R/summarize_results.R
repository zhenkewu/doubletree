#######################################################################
# brief descriptions of functions:
# 1. obtain approximate posterior mean, variance and credible
#    intervals at various levels
# 2. compute some group-specific estimates, where the groups are
#   obtained from the posterior collapsing (in both trees?).
# 3. print the model results: print.nlcm_doubletree(); this can be compact or
#    very detailed.
#######################################################################

#' compute model summaries from the model outputs
#'
#' This function is used in the wrapper function [nlcm_doubletree()]
#'
#' @param mod output from [fit_nlcm_doubletree()]
#' @param dsgn output from [design_doubletree()]
#' @param ci_level credible interval level; default to `0.95`
#'
#' @return a list with elements: `prob_est`,`prob_est_indiv`,`pi_list`
#'
#' @seealso [get_est_cpp_dt()]
#' @importFrom stats qbeta
#' @export
compute_params_dt <- function(mod,dsgn,ci_level=0.95){
  ###########
  ## comment out when building the package:
  # mod <- res$mod
  # dsgn <- dsgn0
  # ci_level <- 0.95
  # B    <- 10000
  ###########

  leaf_ids <- dsgn$leaf_ids
  prob1    <- mod$vi_params$prob1
  prob2    <- mod$vi_params$prob2

  mu_gamma <- mod$vi_params$mu_gamma
  mu_alpha <- mod$vi_params$mu_alpha
  sigma_gamma <- mod$vi_params$sigma_gamma
  sigma_alpha <- mod$vi_params$sigma_alpha

  J <- ncol(dsgn$Y)
  p1 <- length(unique(unlist(dsgn$ancestors[[1]])))
  p2 <- length(unique(unlist(dsgn$ancestors[[2]])))
  pL1 <- length(dsgn$ancestors[[1]])
  pL2 <- length(dsgn$ancestors[[2]])
  Fg1 <- max(dsgn$levels[[1]])
  Fg2 <- max(dsgn$levels[[2]])

  K  <- mod$hyper_fixed$K
  ancestors1 <- dsgn$ancestors[[1]]
  ancestors2 <- dsgn$ancestors[[2]]

  # number of ancestors for each leaf, in tree1, and tree2:
  cardanc1      <- unlist(lapply(dsgn$ancestors[[1]],length))
  cardanc2      <- unlist(lapply(dsgn$ancestors[[2]],length))

  node_select1 <- (prob1 >0.5)+0
  node_select2 <- (as.matrix(do.call("rbind",prob2))>0.5)+0

  z <- stats::qnorm(ci_level+(1-ci_level)/2)

  ## collapsed leaf group estimates: --------------------------------------------
  # for both trees:
  prob_est <- get_est_cpp_dt(
    node_select1,node_select2,
    array(unlist(mu_gamma),c(J,K,p1)),array(unlist(sigma_gamma),c(J,K,p1)),
    array(unlist(mu_alpha),c(pL1,K-1,p2)),array(unlist(sigma_alpha),c(pL1,K-1,p2)),
    ancestors1,ancestors2,
    cardanc1,cardanc2,z)

  prob_est$ci_level <- ci_level

  prob_est$lambda <- split_along_dim(
    array(unlist(lapply(split_along_dim(prob_est$eta_est,3),
                        function(mat){t(apply(mat,1,function(vec)
                          tsb(c(expit(vec),1))))})),c(pL1,K,pL2)),1) # length pL1, each K by pL2.


  ## tree2: (primary interest)
  grp <- matrix(NA,nrow=pL1,ncol=pL2)
  for (v1 in 1:pL1){
    grp[v1,] <- round(prob_est$eta_est[v1,1,],6) # currently ad hoc.
  }
  prob_est$grp <- t(apply(grp,1,function(x) as.integer(factor(x,levels=unique(x)))))
  # this is likely better because it is not on probability scale, which
  # reflects true diffusion on the eta scale.

  prob_est$lambda_collapsed <- vector("list",pL1)
  prob_est$eta_est_collapsed <- vector("list",pL1)
  for (v1 in 1:pL1){
    prob_est$lambda_collapsed[[v1]] <- prob_est$lambda[[v1]][,!duplicated(prob_est$grp[v1,]),drop=FALSE]
    prob_est$eta_est_collapsed[[v1]] <- split_along_dim(prob_est$eta_est,1)[[v1]][,!duplicated(prob_est$grp[v1,]),drop=FALSE]
    prob_est$eta_sd_collapsed[[v1]] <- split_along_dim(prob_est$eta_sd,1)[[v1]][,!duplicated(prob_est$grp[v1,]),drop=FALSE]
  }

  prob_est$members <- vector("list",pL1)
  prob_est$n_obs <- vector("list",pL1)
  names(prob_est$members) <- names(prob_est$n_obs) <- names(dsgn$leaf_ids_units[[1]])[1:pL1]
  for(v1 in 1:pL1){
    n_curr_grps <- length(unique(prob_est$grp[v1,]))
    prob_est$members[[v1]] <- vector("list",n_curr_grps)
    for (g in 1:n_curr_grps){
      prob_est$members[[v1]][[g]] <- names(dsgn$leaf_ids_units[[2]])[prob_est$grp[v1,]==g]
    }
    prob_est$n_obs[[v1]] <- lapply(prob_est$members[[v1]],function(u) sum(names(leaf_ids[[2]])%in%u))
  }

  ## tree1: (secondary interest)
  prob_est$theta <- expit(prob_est$beta_est) # length pL1, each J by K.
  grp_tree1 <- round(prob_est$theta[1,1,],6) # ad hoc.
  prob_est$grp_tree1 <- as.integer(factor(grp_tree1,levels=unique(grp_tree1)))
  n_curr_grps_tree1 <- length(unique(prob_est$grp_tree1))
  prob_est$members_tree1 <- vector("list",n_curr_grps_tree1)
  for (g in 1:n_curr_grps_tree1){
    prob_est$members_tree1[[g]] <-  names(dsgn$leaf_ids_units[[1]])[-(pL1+1)][prob_est$grp_tree1==g]
  }

  prob_est$theta_collapsed <- prob_est$theta[,,!duplicated(prob_est$grp_tree1),drop=FALSE]
  prob_est$beta_est_collapsed <- prob_est$beta_est[,,!duplicated(prob_est$grp_tree1),drop=FALSE]
  prob_est$beta_sd_collapsed <- prob_est$beta_sd[,,!duplicated(prob_est$grp_tree1),drop=FALSE]

  # the following only counts observations with observed leaf_ids in tree1 (excluding NA):
  prob_est$n_obs_tree1 <- sapply(prob_est$members_tree1,function(u) sum(names(leaf_ids[[1]])%in%u))

  # the credible intervals for each group is calculated in summary functions below.

  ## individual leaf estimates: ------------------------------------------------------
  # i.e., no posterior median model selection for either trees:
  prob_est_indiv <- get_est_cpp_dt(
    prob1,as.matrix(do.call("rbind",prob2)),
    array(unlist(mu_gamma),c(J,K,p1)),array(unlist(sigma_gamma),c(J,K,p1)),
    array(unlist(mu_alpha),c(pL1,K-1,p2)),array(unlist(sigma_alpha),c(pL1,K-1,p2)),
    ancestors1,ancestors2,
    cardanc1,cardanc2,z)

  ## tree2: (primary interest)
  prob_est_indiv$lambda <- split_along_dim(
    array(unlist(lapply(split_along_dim(prob_est_indiv$eta_est,3),
                        function(mat){t(apply(mat,1,function(vec)
                          tsb(c(expit(vec),1))))})),c(pL1,K,pL2)),1) # length pL1, each K by pL2.
  ## tree1: (secondary interest)
  prob_est$theta <- expit(prob_est$beta_est) # length pL1, each J by K.


  ## cross-leaf mixture estimation (e.g, cause-specific mortality fractions)---
  pi_list <- list()
  pi_list$pi_est <- sweep(mod$vi_params$dirich_mat,MARGIN=2,
                          colSums(mod$vi_params$dirich_mat),"/")
  pi_list$pi_cil <- apply(mod$vi_params$dirich_mat,2,function(x)
    qbeta(0.025,x,sum(x)-x))
  pi_list$pi_ciu <- apply(mod$vi_params$dirich_mat,2,function(x)
    qbeta(0.975,x,sum(x)-x))

  # return results:
  make_list(prob_est,prob_est_indiv,pi_list)
}


#' `print.nlcm_doubletree` summarizes the results from [nlcm_doubletree()].
#'
#' @param x Output from [nlcm_doubletree()].
#' @param ... Arguments passed to summary and printing methods.
#' @return Summary showing, for each group of leaf nodes discovered by [nlcm_doubletree()],
#' the class prevalences, 95% credible intervals, number of leaf nodes in
#' per group, and number of observations per group.
#'
#' @family lcm_tree results
#' @export
print.nlcm_doubletree <- function(x, ...){
  print(summary(x, ...), ...)
  # Return
  return(invisible(x))
}


#' `summary.nlcm_doubletree` summarizes the results from [nlcm_doubletree()].
#'
#' Will have some tiny randomness associated with the upper and lower bounds,
#' because we simulated the Gaussian variables before converting to probabilities.
#'
#' @param object Output from [nlcm_doubletree()]; of class `nlcm_doubletree`
#' An object of class "lcm_tree".# coeff_type Either "lcm_tree" or "ad_hoc"
#' @param compact If `TRUE`, a more compact summary of results is printed.
#' Only works well when the dimension of the variables is low
#' (say, < 4); otherwise the table takes up too much horizontal space.
#' Default is `FALSE`.
#' @param B Default `10000`: The number of random multivariate Gaussian for logit of V in stick breaking
#' parameterization; used for obtaining the posterior distribution of `lambda^(c,g)_k`
#' @param ... Not used.
#' @return see [print.nlcm_doubletree()]
#'
#' @family nlcm_doubletree results
#' @export
summary.nlcm_doubletree <- function(object,
                                    compact=FALSE,
                                    B=10000,
                                    ...){

  # # <----- temporary during package building.
  # object = res
  # coeff_type = "nlcm_doubletree"
  # compact = FALSE
  # # <---- temporary during package building.

  coeff_type <- "nlcm_doubletree"

  # get estimates:
  # get dimension:
  K <- dim(object$prob_est$theta)[2]
  J <- dim(object$prob_est$theta)[1]
  pL1 <- dim(object$prob_est$theta)[3]
  pL2 <- dim(object$prob_est$eta_est)[3]

  if (coeff_type == "nlcm_doubletree"){ # redundant, used to choose scale.
    coeff_type2 <- "prob"
  } else{
    coeff_type2 <- coeff_type
  }

  est_type <- paste0(coeff_type2,"_est")
  est <- object[est_type][[1]] # the collapse estimates.

  est$n_leaves <- lapply(object$prob_est$members,function(mylist){unlist(lapply(mylist,length))})
  est$leaves   <- lapply(object$prob_est$members,function(mylist){(lapply(mylist,function(v)
  {paste(v,collapse=", ")}))})
  est$n_obs    <- object$prob_est$n_obs
  est$group    <- lapply(object$prob_est$n_obs,function(mylist){1:length(mylist)})


  # pi inference:
  pi_inf <- mapply(function(pi_est,pi_cil,pi_ciu){res = cbind(pi_est,pi_cil,pi_ciu)},pi_est=split_along_dim(object$pi_list$pi_est,2),
                   pi_cil=split_along_dim(object$pi_list$pi_cil,2),
                   pi_ciu=split_along_dim(object$pi_list$pi_ciu,2),SIMPLIFY=FALSE)

  names(pi_inf) <-  names(object$dsgn$leaf_ids_units[[2]])

  rslt <- list()
  if (compact){
    if (coeff_type == "nlcm_doubletree"){
      # simulate the lower and upper intervals in the probability scales:
      for (v1 in 1:pL1){
        curr_grp <- est$group[[v1]]
        tmp_main <- matrix(NA,nrow=length(curr_grp),ncol=3*K)
        for (g in seq_along(curr_grp)){
          tmp_post_simu <- apply(cbind(expit(MASS::mvrnorm(B,
                                                           est$eta_est_collapsed[[v1]][,1,drop=FALSE],
                                                           diag((est$eta_sd_collapsed[[v1]][,1])^2,K-1,K-1))
          ),1),1,tsb)
          ci_mat <- apply(tmp_post_simu,1,stats::quantile,c((1-est$ci_level)/2,est$ci_level+(1-est$ci_level)/2))
          for (k in 1:K){
            tmp_main[g,(3*k-2):(3*k)] <- cbind(est$lambda_collapsed[[v1]][k,g],ci_mat[1,k],ci_mat[2,k])
          }
        }

        rslt[[v1]] <- cbind(names(object$dsgn$leaf_ids_units[[1]])[v1],
                            est$group[[v1]],est$n_leaves[[v1]],unlist(est$n_obs[[v1]]),tmp_main)

        colnames(rslt[[v1]]) <- c("tree1_leaf","group","n_leaves","n_obs",paste0(c("est_lambda", "cil_lambda", "ciu_lambda"),
                                                                                 rep(1:K, each = 3)))
        rslt[[v1]] <- data.frame(rslt[[v1]])
        rslt[[v1]]$tree2_leaves <- unlist(est$leaves[[v1]]) #length equals length of curr_grp.

        rslt[[v1]] <- rslt[[v1]]


      }
    }else{
      stop("[doubletree] did not implement other coeff_type.")
    }
    rslt$est_all <- do.call("rbind",rslt[1:pL1])

    rslt$lambda <- est$lambda_collapsed
    rslt$coeff_type <- coeff_type
    rslt$ci_level <- est$ci_level

    rslt$pi_inf <- pi_inf

    class(rslt) <- "summary.nlcm_doubletree_compact"
    return(rslt)
  }

  # make separate data.frames per group (better for >4 classes):
  grps <- list()
  for (v1 in 1:pL1){
    curr_grp <- est$group[[v1]]
    grps[[v1]] <- list()
    for (g in seq_along(curr_grp)){
      grps[[v1]][[g]] <- list(n_leaves = est$n_leaves[[v1]][g],
                              n_obs    = est$n_obs[[v1]][[g]],
                              leaves   = est$leaves[[v1]][[g]])

      # if (!is.list(est$eta_sd_collapsed)){
      #   est$eta_sd_collapsed <-
      #     split_along_dim(array(est$eta_sd_collapsed,c(1,1,length(est$eta_sd_collapsed))),3)
      #   }
      tmp_post_simu <- apply(cbind(expit(MASS::mvrnorm(B,
                                                       est$eta_est_collapsed[[v1]][,1,drop=FALSE],
                                                       diag((est$eta_sd_collapsed[[v1]][,1])^2,K-1,K-1))
      ),1),1,tsb)
      ci_mat <- apply(tmp_post_simu,1,stats::quantile,c((1-est$ci_level)/2,est$ci_level+(1-est$ci_level)/2))

      grps[[v1]][[g]]$est <- t(rbind(est$lambda_collapsed[[v1]][,g],ci_mat))

      colnames(grps[[v1]][[g]]$est) <- c("est_lambda","cil_lambda","ciu_lambda")
    }
    names(grps[[v1]]) <- paste0("Group", 1:length(grps[[v1]]))
  }

  grps$coeff_type <- coeff_type
  grps$ci_level   <- est$ci_level

  grps$pi_inf <- pi_inf

  class(grps) <- "summary.nlcm_doubletree_long"
  # return if not already (in the 'compact' condition above):
  return(grps)
}


#' Compact printing of [nlcm_doubletree()] model fits
#'
#' `print.summary.nlcm_doubletree_compact` is a print method for class
#' `summary.nlcm_doubletree_compact`.
#'
#' @param x output from `summary.nlcm_doubletree` with `compact = TRUE`.
#' @param print_leaves If `TRUE`, for each discovered group the full list
#' of leaves will be printed. Set this to `FALSE` if these leaf lists
#' make output difficult to read.
#' @param print_pi If `TRUE` print inference about pi; default to `FALSE`.
#' @param digits Number of significant digits to print.
#' @param ... Not used.
#' @return see [print.nlcm_doubletree()]
#'
#' @export
#' @family nlcm_doubletree results
print.summary.nlcm_doubletree_compact <- function(x,
                                                  print_leaves = TRUE,
                                                  print_pi  = FALSE,
                                                  digits = max(3L, getOption("digits") - 3L),
                                                  ...) {
  if (!print_leaves) x$est_all$leaves <- NULL

  if (x$coeff_type == "nlcm_doubletree") {
    cat("Showing latent class model estimates for discovered groups in tree 2; not ad hoc.\n\n")
  }

  x$est_all[,grep("lambda",colnames(x$est_all))] <- round(data.frame(lapply(x$est_all[,grep("lambda",colnames(x$est_all))],
                                                                            as.numeric)),
                                                          digits)

  cat(paste0("Group-specific estimate(s), with credible interval level: ", x$ci_level,"\n"))
  print(knitr::kable(x$est_all))

  if (print_pi){
    cat(paste0("\nPopulation fractions of tree1 leaves, with credible interval level: ", x$ci_level,"\n"))
    for (domain in 1:length(x$pi_inf)){
      cat("\n>---------------- Tree2 Leaf ",domain, ": ", names(x$pi_inf)[domain] ,"-----------------<")
      print(knitr::kable(x$pi_inf[[domain]],digits=digits))
    }

    cat("\n")

  }
  # Return
  return(invisible(x))
}




#' Long-form printing of [nlcm_doubletree()] model fits
#'
#' `print.summary.nlcm_doubletree_long` is a print method for class
#' `summary.nlcm_doubletree_long`.
#'
#' @param x output from `summary.nlcm_doubletree` with `compact = FALSE`.
#' @param print_leaves If `TRUE`, for each discovered group the full list
#' of leaves will be printed. Set this to `FALSE` if these leaf lists
#' make output difficult to read.
#' @param digits Number of significant digits to print.
#' @param ... Not used.
#' @return see [print.nlcm_doubletree()]
#'
#' @export
#' @family nlcm_doubletree results
print.summary.nlcm_doubletree_long <- function(x,
                                               print_leaves = TRUE,
                                               digits = max(3L, getOption("digits") - 3L),
                                               ...) {

  cat("Group-specific lambda estimates for the groups discovered by `doubletree`\n")
  cat(paste0("Credible interval level: ", x$ci_level),"\n")
  if (x$coeff_type == "nlcm_doubletree") {
    cat("Showing lambda estimates for discovered groups; not ad hoc.\n\n")
  }

  pL1 <- length(x)-3
  for (v1 in 1:pL1){
    cat(">---------------- Tree1 Leaf", v1, "-----------------<\n\n")
    for (g in 1:length(x[[v1]])) {

      est <- x[[v1]][[g]]
      cat(">>---------------- Group", g, "----------------<<\n\n")

      cat("Number of tree2 leaves:", est$n_leaves, "\n")
      cat("Number of observations:", est$n_obs, "\n")
      if (print_leaves) {
        cat("List of tree2 leaf nodes:\n")
        cat(est$leaves, "\n\n")
      } else {
        cat("\n\n")
      }

      cat("Class prevalence estimate(s):")
      print(knitr::kable(est$est, digits = digits, row.names = FALSE))

      cat("\n")
    }
  }


  cat(paste0("\nPopulation fractions of tree1 leaves, with credible interval level: ", x$ci_level,"\n"))
  for (domain in 1:length(x$pi_inf)){
    cat("\n>---------------- Tree2 Leaf ",domain, ": ", names(x$pi_inf)[domain] ,"-----------------<")
    print(knitr::kable(x$pi_inf[[domain]],digits=digits))
  }
  cat("\nIf this is hard to read, try print(x, compact = TRUE)\n\n")

  # Return
  return(invisible(x))
}



