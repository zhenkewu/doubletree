rm(list=ls())
library(igraph)

library(doubletree)
library(MASS)
library(poLCA)
library(BayesLCA)
library(matrixStats)

# first tree - over causes:
data("example_cause_edges")
cause_tree <- graph_from_edgelist(example_cause_edges, directed = TRUE)

nodes1  <- names(igraph::V(cause_tree))
leaves1 <- names(igraph::V(cause_tree)[igraph::degree(cause_tree, mode = "out") == 0])
rootnode1 <- names(igraph::V(cause_tree)[igraph::degree(cause_tree, mode = "in") == 0])
pL1 <- length(leaves1)
p1  <- length(nodes1)
# set the levels l*_u for nodes in the cause tree; nodes
# in the same level will share a slab variance multiplier tau_l* (times the edge length
# eminating from the parent).
V(cause_tree)$levels <- rep(1,p1) # this is to set every node to the same level.
E(cause_tree)$weight <- rep(1,length(E(cause_tree))) # set equal edge lengths of 1.

# second tree - over domains:
data("example_domain_edges")
domain_tree <- graph_from_edgelist(example_domain_edges, directed = TRUE)
nodes2  <- names(igraph::V(domain_tree))
leaves2 <- names(igraph::V(domain_tree)[igraph::degree(domain_tree, mode = "out") == 0])
rootnode2 <- names(igraph::V(domain_tree)[igraph::degree(domain_tree, mode = "in") == 0])
pL2 <- length(leaves2)
p2  <- length(nodes2)
V(domain_tree)$levels <- rep(1,p2)
E(domain_tree)$weight <- rep(1,length(E(domain_tree)))

mytrees <- list(tree1 = cause_tree, tree2 = domain_tree)

data("X0") # this is the PHMRC data cleaned by `openVA`
Y <- X0

data("cause_ids") # leaf labels for observations when mapped to the first tree.
data("domain_ids") # leaf labels for observations when mapped to the second tree.


leaf_ids <- list(id1 = cause_ids, id2 = domain_ids)
weighted_edges = c(FALSE,FALSE) # For illustrative purposes, the following ignores
# the input edge lengths.
dsgn  <- design_doubletree(Y,leaf_ids,mytrees,weighted_edges)

table(dsgn$resA2$leaf_ids,useNA="ifany")

resA2 <- dsgn$resA2
resB <- dsgn$resB
sum(resA2$Y-resB$Y,na.rm=TRUE)


#
# leaf labels in tree 1 has missing values; here we take PHMRC data and
# artifically remove ALL the cause-of-death labels (tree1 is cause-of-death tree) for
# a site in India (AP):
#

## The following data creation codes are commented out because it is only needed once
# during package creation:
# cause_ids_except_AP <- cause_ids
# cause_ids_except_AP[which(domain_ids=="AP")] <- NA
# save(cause_ids_except_AP,file="data/cause_ids_except_AP.rda",compress="xz")

data("cause_ids_except_AP") # leaf labels for observations when mapped to the first tree.
leaf_ids2 <- list(id1 = cause_ids_except_AP, id2 = domain_ids)


dsgn2     <- design_doubletree(Y,leaf_ids2,mytrees,weighted_edges)

##table(dsgn2$resA2$leaf_ids,useNA="ifany") # here will list the number of observations
# with missing labels in tree1.

# resA2 <- dsgn2$resA2
# resB <- dsgn2$resB
# sum(resA2$Y-resB$Y,na.rm=TRUE)


# aa = resA$Y[resB$ord[210],] # the row in resA$Y that corresponds to row id 210 in resB$Y
# bb = resB$Y[210,]
#
# setequal(which(is.na(aa)),which(is.na(bb))) & sum(aa[!is.na(aa)]==bb[!is.na(bb)])
#
# aa = resA$Y[7426,] # the row in resA$Y that correspond to a particular row in resB$Y.
# bb = resB$Y[resB$inv_ord[7426],]
# resA$leaf_ids[7426]
# leaf_ids[[2]][resA$ord][7426]
# resB$leaf_ids[resB$inv_ord[7426]]
# setequal(which(is.na(aa)),which(is.na(bb))) & sum(aa[!is.na(aa)]==bb[!is.na(bb)])


#
# FITTING MODELS:
#
# mod <- nlcm_doubletree(Y[,1:20],leaf_ids2,mytrees,weighted_edge = c(FALSE,FALSE),
#                        ci_level = 0.95,
#                        get_lcm_by_group = FALSE,
#                        update_hyper_freq = 10,
#                        print_freq = 1,
#                        quiet      = FALSE,
#                        plot_fig   = FALSE, # <-- used?
#                        tol        = 1E-8,
#                        tol_hyper = 1E-4,
#                        max_iter = 100,
#                        nrestarts = 1,
#                        keep_restarts = TRUE,
#                        parallel = TRUE,
#                        log_restarts = FALSE,
#                        log_dir = ".",
#                        vi_params_init = list(),
#                        hyperparams_init = list(),
#                        random_init = FALSE,
#                        random_init_vals = list(mu_gamma_sd_frac = 0.2,
#                                                mu_alpha_sd_frac = 0.2,
#                                                tau1_lims = c(0.5,1.5),
#                                                tau2_lims = c(0.5,1.5),
#                                                u_sd_frac = 0.2, # this is for logit of prob1.
#                                                psi_sd_frac = 0.2,
#                                                phi_sd_frac = 0.2),
#                        hyper_fixed = list(K=3, # number of latent classes.
#                                           a2=matrix(1,nrow=length(dsgn$ancestors[[1]]),ncol=max(dsgn$levels[[2]])),
#                                           b2=matrix(10,nrow=length(dsgn$ancestors[[1]]),ncol=max(dsgn$levels[[2]])),
#                                           a1 = rep(1,max(dsgn$levels[[1]])),
#                                           b1 = rep(10,max(dsgn$levels[[1]])),
#                                           # both (a1,b1),(a2,b2) encourage shrinkage towards the parent.
#                                           dmat = matrix(1,nrow=length(dsgn$ancestors[[1]]),ncol=length(dsgn$ancestors[[2]])), # (cause,domain).
#                                           s1_u_zeroset = NULL,
#                                           s1_u_oneset = c(1),
#                                           s2_cu_zeroset = NULL,
#                                           s2_cu_oneset = rep(list(1),pL1),
#                                           tau_update_levels = list(1,1)
#                        )
# )

#
#
# weighted_edge = c(FALSE,FALSE)
# ci_level = 0.95
# get_lcm_by_group = FALSE
# update_hyper_freq = 50
# print_freq = 1
# quiet      = FALSE
# plot_fig   = FALSE # <-- used?
# tol        = 1E-8
# tol_hyper = 1E-4
# max_iter = 200
# nrestarts = 1
# keep_restarts = TRUE
# parallel = TRUE
# log_restarts = FALSE
# log_dir = "."
# vi_params_init = list()
# hyperparams_init = list()
# random_init = FALSE
# random_init_vals = list(mu_gamma_sd_frac = 0.2,
#                         mu_alpha_sd_frac = 0.2,
#                         tau1_lims = c(0.5,1.5),
#                         tau2_lims = c(0.5,1.5),
#                         u_sd_frac = 0.2, # this is for logit of prob1.
#                         psi_sd_frac = 0.2,
#                         phi_sd_frac = 0.2)
# hyper_fixed = list(K=3, # number of latent classes.
#                    a2=matrix(1,nrow=length(dsgn$ancestors[[1]]),ncol=max(dsgn$levels[[2]])),
#                    b2=matrix(10,nrow=length(dsgn$ancestors[[1]]),ncol=max(dsgn$levels[[2]])),
#                    a1 = rep(1,max(dsgn$levels[[1]])),
#                    b1 = rep(10,max(dsgn$levels[[1]])),
#                    # both (a1,b1),(a2,b2) encourage shrinkage towards the parent.
#                    dmat = matrix(1,nrow=length(dsgn$ancestors[[1]]),ncol=length(dsgn$ancestors[[2]])), # (cause,domain).
#                    s1_u_zeroset = NULL,
#                    s1_u_oneset = c(1),
#                    s2_cu_zeroset = NULL,
#                    s2_cu_oneset = rep(list(1),pL1),
#                    tau_update_levels = list(1,1)
# )
#
#
# ## design_doubletree:
# leaf_ids <- leaf_ids2 # with missing data.
#
#
# ## nlcm_doubletree:
# # logs
# log_dir <- sub("/$", "", log_dir)
# if (log_restarts) message("[doubletree] Algorithm progress for restart i will be printed to ",
#                           log_dir, "/restart_i_log.txt\n", sep = "")
#
# # Fill in some arguments
# if (nrestarts > 1 & !random_init) { # force random starts if nstarts>1.
#   message("[doubletree] Setting 'random_init = TRUE' since nrestarts > 1\n")
#   random_init <- TRUE
# }
# if (nrestarts == 1) parallel <- FALSE
#
# # construct designed data; here `design_doubletree` reorders the nodes of the two trees, and the rows of the data.
# dsgn <- design_doubletree(Y,leaf_ids,mytrees,weighted_edge) # root_node,weighted_edge <--- need fixing.
#
# # Get hyper_fixed if not supplied:
# if (is.null(hyper_fixed$a1) | is.null(hyper_fixed$b1)) {
#   L1             <- max(dsgn$levels[[1]])
#   hyper_fixed   <- append(hyper_fixed,list(a1 = rep(1, L1)))
#   hyper_fixed$b1 <- rep(10, L1)
#   warning("[doubletree] No fixed hyperparameters (a1,b1) supplied; we set a*_l=1, b*_l=10 for all levels of hyperparameters in tree1.")
# }
#
# if (is.null(hyper_fixed$a2) | is.null(hyper_fixed$b2)) {
#   L2             <- max(dsgn$levels[[2]])
#   pL1            <- length(dsgn$leaf_ids_units[[1]])
#   hyper_fixed   <- append(hyper_fixed,list(a2 = matrix(1, nrow=pL1,ncol=L2)))
#   hyper_fixed$b2 <- matrix(10, nrow=pL1,ncol=L2)
#   warning("[doubletree] No fixed hyperparameters (a2,b2) supplied; we set a_cl=1,b_cl=10 for all levels of hyperparameters in tree2.")
# }
# if (is.null(hyper_fixed$K)) {
#   warning("[doubletree] No fixed # of classes supplied;
#             supply a named element `K` in the list 'hyper_fixed'.")
# }
#
# # Setting up parallelization
# if (parallel) {
#   `%doRestarts%` <- foreach::`%dopar%`
# } else {
#   `%doRestarts%` <- foreach::`%do%`
# }
#
#
#
#
#
#
#
# #---------
#
# ## fit_nlcm_doubletree:
# dsgn = dsgn
# vi_params_init = vi_params_init
# hyperparams_init = hyperparams_init
# random_init = random_init
# random_init_vals = random_init_vals
# tol = tol
# tol_hyper = tol_hyper
# max_iter = max_iter
# print_freq = print_freq
# quiet      = quiet
# plot_fig   = plot_fig
# update_hyper_freq = update_hyper_freq
# hyper_fixed = hyper_fixed
#
# # then it calculates ind_obs_i, ind_obs_j;
# # then it calls initialization function; produces target_id.
#
#
#
# #-------------------------------BEGIN DESIGN PADDING--------------------------#
# dsgn$X <- 2*dsgn$Y-1 # 1, -1 or missing. The rows have been reordered by design_doubletree().
# # obtain the item ids with NON-missing responses for each subject i:
# dsgn$ind_obs_i <- mapply(FUN = function(v) which(!is.na(v)),
#                          split_along_dim(dsgn$X,1),SIMPLIFY=FALSE)
# # obtain the subject ids with NON-missing responses for each item j:
# dsgn$ind_obs_j <- mapply(FUN = function(v) which(!is.na(v)),
#                          split_along_dim(dsgn$X,2),SIMPLIFY=FALSE)
# dsgn$n <- nrow(dsgn$X)
# dsgn$J <- ncol(Y)
# dsgn$p1 <- length(unique(unlist(dsgn$ancestors[[1]])))
# dsgn$p2 <- length(unique(unlist(dsgn$ancestors[[2]])))
# dsgn$pL1 <- length(dsgn$ancestors[[1]])
# dsgn$pL2 <- length(dsgn$ancestors[[2]])
# dsgn$Fg1 <- max(dsgn$levels[[1]])
# dsgn$Fg2 <- max(dsgn$levels[[2]])
#
# if (is.null(hyper_fixed$K)){stop("[doubletree] # of classes 'K' not specified.")}
# if (!is.null(hyper_fixed$K) && hyper_fixed$K<2){stop("[double] # of classes 'K' is 1; theoretically possible, but currently not fully reviewed against K=1.")}
# K  <- hyper_fixed$K
# # number of ancestors for each leaf, in tree1, and tree2:
# dsgn$cardanc1      <- unlist(lapply(dsgn$ancestors[[1]],length))
# dsgn$cardanc2      <- unlist(lapply(dsgn$ancestors[[2]],length))
# #-------------------------------END OF DESIGN PADDING--------------------------#
# #------------------------------------------------------------------------------#
#
# if (!quiet){cat("\n [doubletree] working weights (edge lengths): `h_pau`: \n");print(dsgn$h_pau)}
# # initialize: ----------------------
# init <- R.utils::doCall(initialize_nlcm_doubletree,
#                         vi_params   = vi_params_init,
#                         hyperparams = hyperparams_init,
#                         hyper_fixed = hyper_fixed,
#                         random_init = random_init,
#                         random_init_vals = random_init_vals,
#                         args = c(dsgn))
# vi_params   <- init$vi_params
# hyperparams <- init$hyperparams
# dsgn$target_id <- init$target_id
# dsgn$scenario <- init$scenario
# cat("\n|--- Model Initialized.\n")
#
# # initialize ELBO:
# ELBO_track <- numeric(max_iter)
#
# # update
# Y <- dsgn$Y
# A <- dsgn$A
# leaf_ids_units <- dsgn$leaf_ids_units
# leaf_ids_nodes <- dsgn$leaf_ids_nodes
# ancestors <- dsgn$ancestors
# h_pau <- dsgn$h_pau
# levels <- dsgn$levels
# v_units <- dsgn$v_units
# subject_id_list <- dsgn$subject_id_list
# X <- dsgn$X
# n <- dsgn$n
# J <- dsgn$J
# p1 <- dsgn$p1
# p2 <- dsgn$p2
# pL1 <- dsgn$pL1
# pL2 <- dsgn$pL2
# Fg1 <- dsgn$Fg1
# Fg2 <- dsgn$Fg2
# cardanc1 <- dsgn$cardanc1
# cardanc2 <- dsgn$cardanc2
#
# ind_obs_i = dsgn$ind_obs_i
# ind_obs_j = dsgn$ind_obs_j
#
# target_id = dsgn$target_id
# scenario  = dsgn$scenario
#
# prob1 <- init$vi_params$prob1
# prob2 <- init$vi_params$prob2
#
# mu_gamma <- init$vi_params$mu_gamma
# mu_alpha <- init$vi_params$mu_alpha
# rmat <- init$vi_params$rmat
# emat <- init$vi_params$emat
# dirich_mat <- init$vi_params$dirich_mat
#
# sigma_gamma <- init$vi_params$sigma_gamma
# sigma_alpha <- init$vi_params$sigma_alpha
# tau_1_t <- init$vi_params$tau_1_t
# tau_2_t <- init$vi_params$tau_2_t
# a1_t <- init$vi_params$a1_t
# b1_t <- init$vi_params$b1_t
# a2_t <- init$vi_params$a2_t
# b2_t <- init$vi_params$b2_t
# psi <- init$hyperparams$psi
# g_psi <- init$hyperparams$g_psi
# phi <- init$hyperparams$phi
# g_phi <- init$hyperparams$g_phi
# tau_1 <- init$hyperparams$tau_1
# tau_2 <- init$hyperparams$tau_2
# # local variational parameters and prior variances for the slab components.
# a1 <- hyper_fixed$a1
# b1 <- hyper_fixed$b1
# a2 <- hyper_fixed$a2
# b2 <- hyper_fixed$b2
# dmat <- hyper_fixed$dmat
# K <- hyper_fixed$K
# s1_u_zeroset <- NULL
# s1_u_oneset <- c(1)
# s2_cu_zeroset <- NULL
# s2_cu_oneset <- rep(list(1),pL1) #fixed hyper-parameters; updated.
#
#
#
#
#
#
#
# # run algorithm: ---------------------
# i <- 0
# repeat{
#   #if (quiet){pb$tick()} #;Sys.sleep(3 / 100)}
#   # iterate i
#   i <- i + 1
#
#
#   # check if max_iter reached:
#   if (i > max_iter){
#     i <- max_iter
#     cat(paste("|--- Iteration", i, "complete. \n"))
#     warning("[doubletree] Maximum number of iterations reached! Consider increasing 'max_iter'")
#     break
#   }
#
#
#   ##
#   ##
#   ## The following needs modification:!!!!!!!
#   ##
#   ##
#   # update vi params:
#   vi_params <- R.utils::doCall(update_vi_params_doubletree,
#                                args = c(dsgn, vi_params, hyperparams,
#                                         hyper_fixed))
#
#   # compute ELBO and update psi, phi and hyperparameters (tau_1, tau_2):
#   update_hyper <- i %% update_hyper_freq == 0
#   hyperparams  <- R.utils::doCall(update_hyperparams_doubletree,
#                                   update_hyper = update_hyper,
#                                   quiet      = quiet,
#                                   args = c(dsgn,vi_params,hyperparams,hyper_fixed))
#   ELBO_track[i] <- hyperparams$ELBO
#
#   # print progress:
#   if (i %% print_freq ==0){
#     #if(ELBO_track[i] - ELBO_track[i-1]<0){
#     if (!quiet){
#       cat("|--- Iteration", i, "; epsilon = ", ELBO_track[i] - ELBO_track[i-1], "; ELBO = ", ELBO_track[i],"\n")
#       cat("> empirical class probabilities: ", round(colMeans(vi_params$rmat),4),"\n")
#       cat("> node_select: ",which(vi_params$prob>0.5),"\n")
#     }
#     if (plot_fig){
#       barplot(vi_params$prob)
#       abline(h=0.5,col="purple",lty=2)
#       image(expit(vi_params$mu_gamma[[1]])) # root node.
#     }
#   }
#
#   # check tolerance
#   if (update_hyper & i >= 2 * update_hyper_freq) {
#     # if we just updated hyperparameters, check for convergence of hyperparameters
#     criterion1 <- abs(ELBO_track[i] - ELBO_track[i - update_hyper_freq]) < tol_hyper
#     if (criterion1) {
#       # did last VI update reach convergence?
#       criterion2 <- abs(ELBO_track[i - 1] - ELBO_track[i - 2]) < tol
#       # if yes, both have converged. if not, continue.
#       if (criterion2) break else next
#     } else next
#   } else {
#     criterion3 <- (i > 2) && (abs(ELBO_track[i] - ELBO_track[i - 1]) < tol)
#     # if criterion 3, fill in results until just before the
#     # next hyperparameter update (or max_iter, whichever comes first)
#     if (criterion3) { # is the current update good enough?
#       # if (i<2*update_hyper_freq){ # if update_hyper, but not yet 2*update_hyper_freq:
#       #   # criterion4 <- (abs(ELBO_track[i] - init$hyperparams$ELBO) < tol_hyper) |
#       #   #   (abs(ELBO_track[i] - ELBO_track[1]) < tol_hyper)
#       #   # if (criterion4)
#       #   break
#       # }
#       i2 <- min(ceiling(i / update_hyper_freq) * update_hyper_freq - 1,
#                 max_iter)
#       ELBO_track[(i + 1):i2] <- ELBO_track[i]   # can send this iteration much later; so appears updating more frequent than specified.
#       #ELBO_track[(i + 1):i2] <- hyperparams$ELBO  # can send this iteration much later; so appears updating more frequent than specified.
#       i <- i2
#     }
#   }
# } # end 'repeat' loop.
#
#
#
