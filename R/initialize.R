###############################################################################
# Initialize the variational Bayesian algorithm for the doubletree shrinkage
# approach to domain adaptation.
# 1. variational parameters (in terms of moments):
#    mu_alpha,sigma_alpha
#    mu_gamma,sigma_gamma
#    prob (p_{cu}), prob_star (p^star_u)
#    a'_cl,b'_cl; a''_l, b''_l (for updating rho_cl, rho*_l variational dstns)
# 2. hyperparameters
#    - psi and g_psi,phi and g_phi - not technically hyperparameters;but are local variational
#       tuning parameters that are updated at every iteration
#    - tau, tau_star - these are the variance of the slab component for each node
#       level - they are not updated every iteration, but every few iterations;
#       changing them means changing the prior distributions
# 3. hyperfixed: constants, e.g., number of classes K;
#      - a1,b1: cause-tree: (a*_l,b*_l) for level-specific rho*_l ~beta (a*_l,b*_l).
#      - a2,b2: domain-tree: (a_cl,b_cl) for cause-specific, level-specific rho_cl ~ beta(a_cl,b_cl);
#      - dmat - the Dirichlet prior parameters for modeling CSMFs.
################################################################################

#' Initialize the variational Bayesian algorithm for nested latent class models
#' with double-tree structured shrinkage
#'
#' NB: 1) does this work with missing CODs in all domains (not just the target domain);
#' 2) does this work with partially-known labels in the target domain (not all are
#' missing)
#' Organize intermediate variables,
#' 1. target_id.
#' 2. who are in the target domain
#' 3. who are the people from cause c and domain g
#' 4. check the mapping from integers to leaves; same for all the nodes
#'
#' @param Y matrix of binary observations; preprocessed data from
#' [design_doubletree()] stored in `resB$Y`. `NA` missing data is allowed,
#' which in the nested latent class model will be treated as missing at random.
#' @param A a list of two square binary matrices: the first is `p1` by `p1`,
#' with `1`s in each row indicating the ancestors (the node in that row included)
#' in the first tree (cause tree). The second element is `p2` by `p2`,
#' the counterpart for the second tree (domain tree.) The row names of
#' each square matrix are the orders of the trees.
#' @param leaf_ids_units unit ids nested in each leaf node;
#' a list of two elements (for tree1 and tree2); each element is again a list -
#' each of its element is associated with a leaf in a tree and contains the subject ids in that leaf node.
#' @param leaf_ids a list of length two; each element is a vector of character strings
#' indicating the leaf label an observation is coming from. For the first element,
#' there might be missing values, including unobserved leaf label in tree1.
#' @param leaf_ids_nodes leaf descendants; a list of two elements (for tree1 and tree2);
#' for each element, the leaf descendants for each internal or leaf nodes (a list)
#' @param ancestors a list of two elements (for tree1 and tree2);
#' Each element is a list, with each element containing numeric vector
#' of ancestor nodes - the length is equal to the number of leaves in a tree.
#' @param v_units a list of two elements (for tree1 and tree2) - both are of length `N`;
#' each element is a vector of integers indicating the leaf ids. In the first element, there might
#' be missing entries, indicating the missing leaf label in tree1.
#' @param h_pau edge weights; a list of two elements (for tree1 and tree2);
#' a numeric vector (length = #nodes) the edge length between a node and its parent node.
#' The root node has no parent, we set the edge length toward root node to `1`.
#' @param levels a list of two elements (for tree1 and tree2);
#' For each element, a numeric vector of integers from `1` to `Fg1` for the first tree
#' (or `Fg2` for the second tree), indicating for each node which set of hyperparameters
#' to use. The levels are pre-specified, not estimated. We recommend at least
#' five nodes in any level for learning the slab variance parameters.
#' @param vi_params the list of variational parameters.
#' @param hyperparams the list of hyperparameters
#' @param hyper_fixed a list of fixed hyperparameters: the number of classes `K`;
#' `(a1,b1)` and `(a2,b2)` are the Beta hyperparameters
#' for selecting slab components; `dmat` the hyperparameters in the Dirichlet
#' distributions for each domain (i.e., each leaf of tree2.)
#' @param random_init logical; `TRUE` for adding additional variability to the initial values;
#' This is a must if the algorithm needs multiple random starts to choose
#' the best converged values. Currently the `logit(0.5)` of response probability may produce
#' zero, The sd_frac is a fraction of the logit value, which may not be producing any additional randomness.
#' @param random_init_vals NB: fill out specific elements; this is to be done.
#'
#' @importFrom BayesLCA unMAP
#' @importFrom poLCA poLCA
#' @importFrom stats as.formula rnorm
#' @return a list of two elements:`vi_params`,`hyperparams` with initial values for the
#' nested latent class model with double tree shrinkage. `target_id` and `scenario` are
#' both added (use `?doubletree` to see meanings of scenarios).
#'
initialize_nlcm_doubletree <- function(Y,A,
                                       leaf_ids_units,
                                       leaf_ids,
                                       leaf_ids_nodes,
                                       ancestors,
                                       v_units,
                                       h_pau,
                                       levels,
                                       vi_params,
                                       hyperparams,
                                       hyper_fixed,
                                       random_init,
                                       random_init_vals)
{
  # for LCA in poLCA:
  Y_df <- as.data.frame(Y)
  form <- as.formula(paste0("cbind(",paste(colnames(Y_df),collapse=","),")~1"))

  X <- 2*Y-1
  n <- nrow(X)
  J <- ncol(X)
  p1 <- length(unique(unlist(ancestors[[1]])))
  p2 <- length(unique(unlist(ancestors[[2]])))
  pL1 <- length(ancestors[[1]])
  pL2 <- length(ancestors[[2]])
  Fg1 <- max(levels[[1]])
  Fg2 <- max(levels[[2]])
  cat("\n number of levels of hyperparameters for tree1 and tree2: (L*,L) = (", Fg1,",",Fg2,")... \n")

  K    <- hyper_fixed$K
  dmat <- hyper_fixed$dmat

  ## identify the target domain(s):
  # names of all leaf nodes in tree2 (with proper numbers of repetitions):
  labels2 <- rep(names(leaf_ids_units[[2]]),times=lapply(leaf_ids_units[[2]],length))
  scenario <- NA # ?doubletree for meanings of scenarios.
  if (!is.null(leaf_ids_units[[1]]$NA_tree1)){
    labels2_NA <- labels2[match(leaf_ids_units[[1]]$NA_tree1,unlist(leaf_ids_units[[2]]))] # <---- This would have issue if no leaf label in tree1 is missing.
    labels2_not_NA <- labels2[-match(leaf_ids_units[[1]]$NA_tree1,unlist(leaf_ids_units[[2]]))] # <---- This would have issue if no leaf label in tree1 is missing.
    if (length(unique(labels2_NA)) >1){
      scenario <- "c"
    }else{
      if (length(intersect(unique(labels2_NA),labels2_not_NA))==0){
        scenario <- "b1" # 1 target, all missing
      } else{
        scenario <- "b2"
      }
    }
  } else{
    scenario <- "a"
  }
  cat("[doubletree] scenario `",scenario,"'>>> use `?doubletree` for meaning of the scenario.")

  ## node names and leaf nodes (they keep the order of internal and leaf nodes);
  ## The following is to keep track of the integer and the actual node/leaf names:
  node_nms <- list(rownames(A[[1]]),rownames(A[[2]]))
  leaf_nms <- list(node_nms[[1]][-(1:(p1-pL1))],node_nms[[2]][-(1:(p2-pL2))])
  # the above is okay because after designtree() in lotR, it reorders
  # the nodes by (internal nodes--> leaf nodes - which has pL1 of them for tree1).

  #######################################################################
  ####### initialize: mu_gamma [p1,J,K], mu_alpha [p1,p2,K-1]
  #######################################################################
  if (is.null(vi_params[["mu_gamma"]]) | is.null(vi_params[["mu_alpha"]])){
    beta_map <- array(0,c(p1,J,K))     # transformed parameters for response probs.
    eta_map  <- array(0,c(pL1,p2,K-1)) # transformed parameters for class probs.

    #start_class_member_prob <- matrix(runif(n*K),nrow=n,ncol=K)
    #start_class_member_prob <- start_class_member_prob/rowSums(start_class_member_prob)

    ## fit some simple latent class models to the leaves:
    ## NB: need to check the meaning/ordering of these indices:
    # for (v1 in 1:p1){ # over all nodes in tree1:
    #   u1 <- leaf_ids_nodes[[1]][[v1]] # get leaves in tree1 that are nested under v1.
    #   for (v2 in 1:p2){ # over all nodes in tree2:
    #     u2 <- leaf_ids_nodes[[2]][[v2]] # get leaves in tree2 that are nested under v2.
    #     units <- intersect(unlist(leaf_ids_units[[1]][u1]),unlist(leaf_ids_units[[2]][u2])) # get all subjects that are nested under v1 and v2.
    #     cat(length(units)," ")
    #     if (v2==p2){cat("\n")}
    #
    #     # fit a latent class model, accommodating for potential missing responses:
    #     for (j in 1:ncol(Y_df)){Y_df[,j] <- as.integer(as.factor(Y_df[,j]))}
    #     mod <- poLCA(form,Y_df[units,],nclass=K,na.rm=FALSE) # this is slow because EM - not VB. # <--- import poLCA.
    #
    #     # get the prevalence estimate:
    #     tmp <- sapply(mod$P,function(s) min(max(s,0.01),0.99)) # truncate to reasonable values...
    #     tau <- (tmp/sum(tmp)) # normalize; the above may not sum to 1.
    #     # transform to stick-breaking representation:
    #     eta_map[v1,v2,] <- logit(prob2stick(tau))[-length(tau)]
    #     # get the response profile estimates:
    #     probmat <- t(as.matrix(as.data.frame(lapply(mod$probs,function(x) x[,2])))) # there are some close to zeros; check what are the symptoms.
    #     beta_map[v1,,] <- logit(pmin(pmax(probmat,0.01),0.99)) # logit transformation.
    #     #mod <- blca.vb(Y[units,,drop=FALSE],K,method="vb",verbose=FALSE,start.vals="single")
    #     # The above is fast; but unfortunately here blca.vb cannot deal with missing responses directly...
    #   }
    # }

    # fit a latent class model to ALL data, accommodating for potential missing responses:
    ## The following code first computes a plausible set of parameter values (response probabilities
    ## and prevalence parameters) for each node; then we convert these values back to increments.
    ## Because fitting LCA with missing responses uses poLCA which can be slow,
    ## in the following we only compute a single set of parameter values using all data, which
    ## upon transformation to increments means everything beyond root nodes are all ZEROS.
    for (j in 1:ncol(Y_df)){Y_df[,j] <- as.integer(as.factor(Y_df[,j]))}
    mod <- poLCA(form,Y_df,nclass=K,na.rm=FALSE,verbose=FALSE) # this is slow because EM - not VB. # <--- import poLCA.
    # get the prevalence estimate:
    tmp <- pmin(pmax(mod$P,0.01),0.99) # truncate to reasonable values.
    tau <- (tmp/sum(tmp)) # normalize; the above may not sum to 1 after truncation.
    smart_guess_rmat <- (tmp/sum(tmp)) # normalize; the above may not sum to 1 after truncation.

    # transform to stick-breaking representation:
    for (v in 1:pL1){ # for every cause
      for (v2 in 1:p2){ # for every node in tree2:
        eta_map[v,v2,] <- logit(pmin(pmax(prob2stick(tau),0.001),0.999))[-length(tau)]
      }
    }
    # get the response profile estimates:
    probmat <- t(as.matrix(as.data.frame(lapply(mod$probs,function(x) x[,2]))))
    # there are some close to zeros; check what are the symptoms.
    for (v1 in 1:p1){
      beta_map[v1,,] <- logit(pmin(pmax(probmat,0.01),0.99)) # logit transformation.
    }
    rm("tau")
    # replace infinite values (if any):
    beta_map[is.infinite(beta_map) & beta_map<0] <- - 5
    beta_map[is.infinite(beta_map) & beta_map>0] <-   5
    eta_map[is.infinite(eta_map) & eta_map<0] <- - 5
    eta_map[is.infinite(eta_map) & eta_map>0] <-   5
    # replace NaN values (if any):
    beta_map[is.na(beta_map)] <- 0
    beta_map[is.na(beta_map)] <- 0
    eta_map[is.na(eta_map)]   <- 0
    eta_map[is.na(eta_map)]   <- 0

    # transform from estimates to increments (eta to xi - along tree2; beta to zeta - along tree1):
    A1_inv <- solve(A[[1]])
    A2_inv <- solve(A[[2]])
    # for values on the nodes that represents successive sums (beta), convert them to increments (gamma).
    mu_gamma <- apply(beta_map,c(2,3),function(v) as.matrix(A1_inv%*%matrix(v,ncol=1))) # apply this across nodes for tree1.
    # similar as above, but do this separately for each cause:
    mu_alpha <- aperm(apply(eta_map,c(1,3),function(v) as.matrix(A2_inv%*%matrix(v,ncol=1))),
                      c(2,1,3))
    if (sum(mu_gamma[-1,,])!=0 | sum(mu_alpha[,-1,])!=0){stop("[doubletree] error in overall LCM initialization.")}
  }

  # assign the initialized mu_gamma to `vi_params`:
  if (is.null(vi_params[["mu_gamma"]])){
    #in vi_params:
    # convert to list
    vi_params$mu_gamma <- split_along_dim(mu_gamma,1)
  } else{ # if mu_gamma is specified in the vi_params, we check compatibility:
    check <- is.list(vi_params$mu_gamma) &&
      sum(sapply(vi_params$mu_gamma,function(x) sum(dim(x)==c(J,K))==2))==p1
    if (!check) stop("[doubletree] incompatible dimensions of initial value for `mu_gamma`")
  }
  if (random_init){
    vi_params$mu_gamma <- lapply(vi_params$mu_gamma,
                                 function(mu) matrix(mu+rnorm(J*K,sd = c(pmax(abs(mu),0.1))*random_init_vals$mu_gamma_sd_frac),nrow=J,ncol=K))
    # the pmax part is helpful if mu is zero, which would not produce any additional randomness
    # upon multiplication by the sd_frac.
  }

  # assign the initialized mu_alpha to `vi_params`:
  if (is.null(vi_params[["mu_alpha"]])){
    vi_params$mu_alpha <- split_along_dim(mu_alpha,2)
  } else{
    check <- is.list(vi_params$mu_alpha) &&
      sum(sapply(vi_params$mu_alpha,function(x) sum(dim(x)==c(pL1,K-1))==2))==p2
    if (!check) stop("[doubletree] incompatible dimensions of initial value for `mu_alpha`")
  }
  if (random_init){
    vi_params$mu_alpha <- lapply(vi_params$mu_alpha,
                                 function(mu) matrix(mu+rnorm(pL1*(K-1),sd=c(pmax(abs(mu),0.1))*random_init_vals$mu_alpha_sd_frac),nrow=pL1,ncol=K-1))
    # the pmax part is helpful if mu is zero, which would not produce any additional randomness
    # upon multiplication by the sd_frac.
  }

  ##############################################################################
  ## initialize hyper-parameters
  ##  tau_1: slab component variances for tree1;
  ##  tau_2: slab component variances for tree2.
  ##############################################################################
  # for gammas; cause tree in VA; tree1.
  if (is.null(hyperparams[["tau_1"]])){
    mu_gamma_sq_over_h <- mapply(FUN = function(mat,h){mat^2/h},mat=vi_params$mu_gamma,
                                 h=h_pau[[1]],SIMPLIFY=FALSE)
    hyperparams$tau_1 <- sapply(1:Fg1,function(l)
      mean(unlist(mu_gamma_sq_over_h[levels[[1]]==l])))
    hyperparams$tau_1[hyperparams$tau_1 ==0] <- 0.01
    # currently set to nonzeros, because if zeros, will cause issues in gamma update in the 1/tau_1_t part.
  } else {
    check <- is.numeric(hyperparams$tau_1) &&
      length(hyperparams$tau_1)==Fg1
    if (!check) stop("[doubletree] Incompatible initial values for `tau_1` - logit response probabilities (gamma)")
  }
  if (random_init){
    hyperparams$tau_1 <- sapply(hyperparams$tau_1,function(tau) runif(1,min=tau*random_init_vals$tau1_lims[1],
                                                                      max=tau*random_init_vals$tau2_lims[2]))
  }

  # for alphas; domain tree in VA; tree2.
  if (is.null(hyperparams[["tau_2"]])){
    mu_alpha_sq_over_h <- mapply(FUN = function(mat,h){mat^2/h},mat=vi_params$mu_alpha,
                                 h=h_pau[[2]],SIMPLIFY=FALSE)
    hyperparams$tau_2 <- sapply(1:Fg2,function(l)
      mean(unlist(mu_alpha_sq_over_h[levels[[2]]==l])))
    hyperparams$tau_2[hyperparams$tau_2 ==0] <- 0.01
    # currently set to nonzeros, because if zeros, will cause issues in alpha update in the 1/tau_2_t part.
  }else{
    check <- is.numeric(hyperparams$tau_2) &&
      length(hyperparams$tau_2) == Fg2
    if (!check) stop("[doubletree] Incompatible intial value for `tau_2` - slab component variances - tree2.")
  }
  if (random_init){
    hyperparams$tau_2 <- sapply(hyperparams$tau_2,
                                function(tau) runif(1,min=tau*random_init_vals$tau2_lims[1],
                                                    max=tau*random_init_vals$tau2_lims[2]))
  }

  #############################################################################
  ## initialize prob2_cu (in tree2; for s_cu), prob1_u (in tree1; for s*_u):
  #############################################################################
  if (is.null(vi_params[["prob1"]])){
    vi_params$prob1 <- rep(0.5,p1)
  }else{
    check <- is.numeric(vi_params$prob1) &&
      length(vi_params$prob1)==p1 &&
      sum(vi_params$prob1 >=0) == p1 &&
      sum(vi_params$prob1 <=1) == p1
    if (!check) stop("[doubletree] Incompatible initial values for 'p*_u' (for variational probability of s*_u = 1)")
  }
  if (random_init){
    prob1 <- vi_params$prob1
    prob1 <- pmin(pmax(prob1,0.01),0.99)
    u1 <- log(prob1/(1-prob1))
    u1 <- u1+rnorm(p1)*random_init_vals$u_sd_frac*pmax(abs(u1),0.1)
    vi_params$prob1 <- expit(u1)
  }

  if (is.null(vi_params[["prob2"]])){ # this needs to be for each cause so needs
    ## to be a list of length pL1.
    vi_params$prob2 <- rep(list(rep(0.5,p2)),pL1)
  } else{
    check <- sum(unlist(lapply(vi_params$prob2,is.numeric)))==pL1 &&
      length(vi_params$prob2)==pL1 &&
      sum(unlist(lapply(vi_params$prob2,function(x) sum(x<0))))==0 &&
      sum(unlist(lapply(vi_params$prob2,function(x) sum(x>1))))==0
    if (!check) stop("[doubletree] Incompatible initial values for 'p_cu' (for variational probability of s_cu = 1)")
  }
  if (random_init){
    prob2 <- vi_params$prob2
    prob2 <- lapply(vi_params$prob2,function(v) pmin(pmax(v,0.01),0.99))
    vi_params$prob2 <- lapply(prob2,function(v) {expit(log(v/(1-v))+
                                                         rnorm(p2)*random_init_vals$u_sd_frac*
                                                         pmax(0.1,abs(log(v/(1-v))))) # would not be actually random if all probs are .5.
    })
  }

  #######################################################################
  ## (a'_cl,b'_cl)^t, (a''_l,b''_l)^t - these are variational parameters
  ## when updating the variational distribution for rho_cl and rho*_l:
  #######################################################################
  #----
  if (is.null(vi_params[["a2_t"]])){
    vi_params$a2_t <- matrix(NA,nrow=pL1,ncol=Fg2)
    for (v1 in 1:pL1){
      for (f in 1:Fg2){
        # initialize using VI updates:
        vi_params$a2_t[v1,f] <- hyper_fixed$a2[v1,f]+sum(vi_params$prob2[[v1]][levels[[2]]==f])
      }
    }
  }else{
    check <- is.numeric(vi_params$a2_t) &&
      dim(vi_params$a2_t) == c(pL1,Fg2) &&
      sum(vi_params$a2_t <=0)==0
    if (!check) stop("[doubletree] Incompatible initial value supplied for a2_t.")
  }
  #----
  if (is.null(vi_params[["b2_t"]])){
    vi_params$b2_t <- matrix(NA,nrow=pL1,ncol=Fg2)
    for (v1 in 1:pL1){
      for (f in 1:Fg2){
        # initialize using VI updates:
        vi_params$b2_t[v1,f] <- hyper_fixed$b2[v1,f]+sum(1-vi_params$prob2[[v1]][levels[[2]]==f])
      }
    }
  }else{
    check <- is.numeric(vi_params$b2_t) &&
      dim(vi_params$b2_t) == c(pL1,Fg2) &&
      sum(vi_params$b2_t <=0)==0
    if (!check) stop("[doubletree] Incompatible initial value supplied for b2_t.")
  }

  #----a1_t.
  if (is.null(vi_params[["a1_t"]])){
    vi_params$a1_t <- numeric(Fg1)
    for (f in 1:Fg1){
      # initialize using VI updates:
      vi_params$a1_t[f] <- hyper_fixed$a1[f]+sum(vi_params$prob1[levels[[1]]==f])
    }
  }else{
    check <- is.numeric(vi_params$a1_t) &&
      length(vi_params$a1_t) == Fg1 &&
      sum(vi_params$a1_t <=0)==0
    if (!check) stop("[doubletree] Incompatible initial value supplied for a1_t.")
  }
  #----b1_t.
  if (is.null(vi_params[["b1_t"]])){
    vi_params$b1_t <- numeric(Fg1)
    for (f in 1:Fg1){
      # initialize using VI updates:
      vi_params$b1_t[f] <- hyper_fixed$b1[f]+sum(1-vi_params$prob1[levels[[1]]==f])
    }
  }else{
    check <- is.numeric(vi_params$b1_t) &&
      length(vi_params$b1_t) == Fg1 &&
      sum(vi_params$b1_t <=0)==0
    if (!check) stop("[doubletree] Incompatible initial value supplied for b1_t.")
  }

  #####################################################################
  ## initialize local variational parameters psi, phi; they are
  ## put into hyperparameters because they got updated in the same
  ## function in this package; we also calculate g_psi and g_phi.
  #####################################################################

  #####################################################################
  ## psi: _jk^(c): J,K,C - for response probabilities.
  ## phi: _k^{(c,g)}: K-1, C, pL2 (including target domain - check!) - for class probabilities.
  ## both are only needed at the leaf levels because that is where the response
  ## data is entangled with the parameters.
  #####################################################################

  #
  # psi_jk^c - sqrt of expected beta_squared; only need leaf levels of tree1.
  #
  if (is.null(hyperparams[["psi"]])){
    zeta <- mapply(FUN= function(prob,mu) prob*mu,
                   prob=vi_params$prob1,
                   mu=vi_params$mu_gamma,SIMPLIFY=FALSE) # check if this works - prob is vector
    beta_v <- array(NA,c(pL1,J,K))
    for (v in 1:pL1){beta_v[v,,] <- Reduce('+',zeta[ancestors[[1]][[v]]])}
    hyperparams$psi <- abs(beta_v)# not exactly as in VI update; but close.
  } else{
    check <- is.numeric(hyperparams$psi) &&
      sum(dim(hyperparams$psi)==c(pL1,J,K))==3
    if (!check) stop("[doubletree] Incompatible initial values for 'psi' - local variational parameters.")
  }
  if (random_init){
    hyperparams$psi <- abs(hyperparams$psi*(1+rnorm(pL1*J*K)*random_init_vals$psi_sd_frac))
  }
  hyperparams$g_psi <- g_fun.vec(hyperparams$psi)

  tmp_prob2 <- split_along_dim(do.call("rbind",vi_params$prob2),2)
  # phi_k^{(c,g)}; only need leaf levels of tree1 and tree2:
  if (is.null(hyperparams[["phi"]])){# pL1 by pL2 by K-1.
    xi <- mapply(FUN = function(prob,mu) sweep(mu,MARGIN=1,prob,"*"),
                 prob = tmp_prob2, mu=vi_params$mu_alpha,SIMPLIFY=FALSE)
    hyperparams$phi <- array(0,c(pL1,pL2,K-1))
    for (v in 1:pL2){
      eta_v <- Reduce('+',xi[ancestors[[2]][[v]]])
      hyperparams$phi[,v,] <- abs(eta_v)
    }
  } else{
    check <- is.numeric(hyperparams$phi) &&
      sum(dim(hyperparams$phi)==c(pL1,pL2,K-1))==3
    if (!check) stop("[doubletree] Incompatible initial values for 'phi' - local variational parameters.")
  }
  if (random_init){
    hyperparams$phi <- abs(hyperparams$phi*(1+rnorm(pL1*pL2*(K-1))*random_init_vals$phi_sd_frac))
  }
  hyperparams$g_phi <- g_fun.vec(hyperparams$phi)

  #############################################################################
  ## Other variational parameters: multinomial variational parameters etiology
  ## in the target domain, and class probabilities for any observation
  #############################################################################

  ## multinomial variatoinal parameters: N by pL1 - COD for people in target domain.
  if (is.null(vi_params[["emat"]])){ # CHECK: are any of the pL1 leaves not observed in the leaf_ids[[1]]? If no, then remove that leaf.

    obs_cts <- sapply(1:pL1,function(x){sum(v_units[[1]][!is.na(v_units[[1]])]==x)+1})

    tmp_indicators <- rep(NA,n)

    if (sum(!is.na(v_units[[1]]))>0){ # if there is at least one observed leaf label in tree1.
      tmp_indicators[!is.na(v_units[[1]])] <- v_units[[1]][!is.na(v_units[[1]])]
    }
    if (scenario!="a"){
      tmp_indicators[is.na(v_units[[1]])]  <-
        sample(1:pL1,size=sum(is.na(v_units[[1]])),replace=TRUE)
        # sample(1:pL1,size=sum(is.na(v_units[[1]])),replace=TRUE,prob=obs_cts/sum(obs_cts))
    }
    vi_params$emat <- unMAP(tmp_indicators) # this is for all observation; ignore the ones with observed CODs.
    if (scenario!="a"){
      vi_params$emat[is.na(v_units[[1]]),] <-
        matrix(1/pL1,nrow=sum(is.na(v_units[[1]])),ncol=pL1,byrow=TRUE) # non-zero entries.
        # matrix(obs_cts/sum(obs_cts),nrow=sum(is.na(v_units[[1]])),ncol=pL1,byrow=TRUE) # non-zero entries.
    }
  } # NB: need to selectively ignore the ones with observed COD.

  ## multinomial variational parameters: N by K - class for people in each source and target domain.
  if (is.null(vi_params[["rmat"]])){ # always missing in latent class analysis!
    vi_params$rmat <- matrix(smart_guess_rmat,nrow=n,ncol=K,byrow=TRUE) # non-zero entries.
    #vi_params$rmat <- matrix(rep(1/K,K),nrow=n,ncol=K,byrow=TRUE) # non-zero entries.
    # vi_params$rmat <- cbind(rep(0.9,n),matrix(0.1/(K-1),nrow=n,ncol=K-1))
    #vi_params$rmat <- unMAP(sample(1:K,size=n,replace=TRUE))
  }

  ##############################################################################
  ## initialize the CSMF for all domains (source and target). (check Step 1d of Appendix).
  #############################################################################
  # get the leaf id in tree2 for the unique label in labels2_NA:
  ## needs a mapping from the numbers to the actual domain names.
  if (is.null(vi_params[["dirch_mat"]])){
    vi_params$dirich_mat <- matrix(NA,nrow=pL1,ncol=pL2)
    for (v2 in 1:pL2){
      for (v1 in 1:pL1){
        vi_params$dirich_mat[v1,v2] <- dmat[v1,v2] + sum(vi_params$emat[v_units[[2]]==v2,v1])
      }
    }
  } else{
    check <- is.numeric(vi_params$dirich_mat) &&
      sum(dim(vi_params$dirich_mat)==c(pL1,pL2))==2 &&
      sum(vi_params$dirich_mat<=0)==0
    if (!check) stop("[doubletree] Incompatible initial values for 'dirich_mat'.")
  }

  #######################################################################
  ## initialize the prior variance parameters for alpha and gamma
  ## these are tau_1_t, tau_2_t; just assign tau_1 and tau_2 to respective
  ## nodes according to the levels of the nodes.
  #######################################################################
  if (is.null(vi_params[["tau_1_t"]])){
    vi_params$tau_1_t <- hyperparams$tau_1[levels[[1]]]
  } else{
    check <- is.numeric(vi_params$tau_1_t) &&
      length(vi_params$tau_1_t)==p1 &&
      sum(vi_params$tau_1_t<0)==0
    if(!check) stop("[doubletree] Incompatible initial values for 'tau_1_t'; for alpha")
  }

  if (is.null(vi_params[["tau_2_t"]])){
    vi_params$tau_2_t <- hyperparams$tau_2[levels[[2]]]
  } else{
    check <- is.numeric(vi_params$tau_2_t) &&
      length(vi_params$tau_2_t)==p2 &&
      sum(vi_params$tau_2_t<0)==0
    if(!check) stop("[doubletree] Incompatible initial values for 'tau_2_t'; for gamma.")
  }

  #####################################################################
  ## compute sigma_gamma, sigma_alpha,
  ## these are the slab-component variances in the variational distribution
  ## for gamma and alpha; so they must involve not only the prior variance parameters
  ## above:)
  #####################################################################
  #
  # NB: induced factorization (or verify) that the vector of
  # gamma and alpha has a joint distribution that can factorize
  # in a way that independent, component-specific Gaussians are
  # the optimal variational distribution:
  #
  # initialize the variational variance parameters for the gammas:
  if (is.null(vi_params[["sigma_gamma"]])){# check step 1f in Appendix.
    vi_params$sigma_gamma <- array(NA,c(p1,J,K))
    for (u in 1:p1){
      # get subjects nested under node u:
      leaf_desc_curr <- leaf_ids_nodes[[1]][[u]]
      for (j in 1:J){
        for (k in 1:K){
          vi_params$sigma_gamma[u,j,k] <- 1/(vi_params$tau_1_t[u]*h_pau[[1]][u])
          for (v in leaf_desc_curr){
            vi_params$sigma_gamma[u,j,k] <- vi_params$sigma_gamma[u,j,k]+
              2*hyperparams$g_psi[v,j,k]*sum(vi_params$rmat[,k]*vi_params$emat[,v])
            # the final emat term requires correct assignment of leaf 1 label based on observed leaves in tree1.
          }
        }
      }
    }
    vi_params$sigma_gamma <- split_along_dim(1/vi_params$sigma_gamma,1)
  } else{
    check <- is.list(vi_params$sigma_gamma) &&
      sum(sapply(vi_params$sigma_gamma,function(x) {sum(dim(x) == c(J,K))==2}))==p1
    if (!check) stop("[doubletree] Incompatible initial value for sigma_gamma.")
  }

  # Initialize the variational parameters - the variances for the alpha parameters:
  # This is related to the question of how are they factorized - induced factorization.

  if (is.null(vi_params[["sigma_alpha"]])){ # check step 1e in Appendix:
    vi_params$sigma_alpha <- array(NA,c(p2,pL1,K-1))
    for (u in 1:p2){
      leaf_desc_curr <- leaf_ids_nodes[[2]][[u]]
      leaf_list_tmp <- leaf_ids_units[[2]][leaf_desc_curr]
      units <- unlist(leaf_list_tmp) # get the subject ids.
      v_units_curr <- unlist(mapply(rep,leaf_ids_nodes[[2]][[u]],unlist(lapply(leaf_list_tmp,length))))
      for (v in 1:pL1){
        for (k in 1:(K-1)){
          vi_params$sigma_alpha[u,v,k] <- 1/(vi_params$tau_2_t[u]*h_pau[[2]][u])
          vi_params$sigma_alpha[u,v,k] <- vi_params$sigma_alpha[u,v,k]+
            2*sum(hyperparams$g_phi[v,v_units_curr,k]*vi_params$emat[units,v]*rowSums(vi_params$rmat[units,k:K]))
        }
      }
    }
    vi_params$sigma_alpha <- split_along_dim(1/vi_params$sigma_alpha,1)
  } else{
    check <- is.list(vi_params$sigma_alpha) &&
      sum(sapply(vi_params$sigma_alpha,function(x) {sum(dim(x) == c(pL1,K-1))==2}))==p2
    if (!check) stop("[doubletree] Incompatible initial value for 'sigma_alpha'; variational
                     variance for alpha_k^(c,u) for s_cu=1 slab component")

  }
  ##########################################################################
  ## ELBO:
  #########################################################################
  if (is.null(hyperparams$ELBO)){
    hyperparams$ELBO  <- 1E-16
  } else{
    hyperparams$ELBO <- hyperparams$ELBO[length(hyperparams$ELBO)]
  }
  # return values:
  target_id <- NULL
  if (scenario!="a"){
    target_id <- which(leaf_nms[[2]]%in%unique(labels2_NA))
  }
  make_list(vi_params,hyperparams,scenario,target_id)
}


