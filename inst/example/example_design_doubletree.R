rm(list=ls())
library(igraph)

library(doubletree)
library(MASS)
library(poLCA)
library(BayesLCA)

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
#
# #
# # TRY FITTING MODELS:
# #
#
# mod <- nlcm_doubletree(Y,leaf_ids2,mytrees,weighted_edge = c(FALSE,FALSE),
#                 ci_level = 0.95,
#                 get_lcm_by_group = FALSE,
#                 update_hyper_freq = 50,
#                 print_freq = 1,
#                 quiet      = FALSE,
#                 plot_fig   = FALSE,
#                 tol        = 1E-8,
#                 tol_hyper = 1E-4,
#                 max_iter = 200,
#                 nrestarts = 1,
#                 keep_restarts = TRUE,
#                 parallel = TRUE,
#                 log_restarts = FALSE,
#                 log_dir = ".",
#                 vi_params_init = list(),
#                 hyperparams_init = list(),
#                 random_init = FALSE,
#                 random_init_vals = list(mu_gamma_sd_frac = 0.2,
#                                         mu_alpha_sd_frac = 0.2,
#                                         tau1_lims = c(0.5,1.5),
#                                         tau2_lims = c(0.5,1.5),
#                                         u_sd_frac = 0.2, # this is for logit of prob1.
#                                         psi_sd_frac = 0.2,
#                                         phi_sd_frac = 0.2),
#                 hyper_fixed = list(K=3, # number of latent classes.
# a2=matrix(1,nrow=length(dsgn$ancestors[[1]]),ncol=max(dsgn$levels[[2]])),
# b2=matrix(10,nrow=length(dsgn$ancestors[[1]]),ncol=max(dsgn$levels[[2]])),
# a1 = rep(1,max(dsgn$levels[[1]])),b1 = rep(10,max(dsgn$levels[[1]])),
# # both (a1,b1),(a2,b2) encourage shrinkage towards the parent.
# dmat = matrix(1,nrow=length(dsgn$ancestors[[1]]),ncol=length(dsgn$ancestors[[2]])),
### (cause,domain).
# s1_u_zeroset = NULL,
# s1_u_oneset = c(1),
# s2_cu_zeroset = NULL,
# s2_cu_oneset = rep(list(1),pL1),
# tau_update_levels = list(1,1)
#                 )
# )
