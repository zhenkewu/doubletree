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


#
# FITTING MODELS:
#
mod <- nlcm_doubletree(
  Y[,1:15],leaf_ids2,mytrees,weighted_edges = c(FALSE,FALSE),
  ci_level = 0.95,
  get_lcm_by_group = FALSE,
  update_hyper_freq = 10,
  print_freq = 1,
  quiet      = FALSE,
  plot_fig   = FALSE, # <-- used?
  tol        = 1E-8,
  tol_hyper = 1E-4,
  max_iter = 100,
  nrestarts = 1,
  keep_restarts = TRUE,
  parallel = TRUE,
  log_restarts = FALSE,
  log_dir = ".",
  vi_params_init = list(),
  hyperparams_init = list(),
  random_init = FALSE,
  hyper_fixed = list(K=2, LD=TRUE,# number of latent classes.
                     a1 = rep(20,max(igraph::V(cause_tree)$levels)),
                     b1 = rep(1,max(igraph::V(cause_tree)$levels)),
                     a2=matrix(1,nrow=length(dsgn$ancestors[[1]]),ncol=max(igraph::V(domain_tree)$levels)),
                     # <-- NB: where do we specify levels? in the tree.
                     b2=matrix(10,nrow=length(dsgn$ancestors[[1]]),ncol=max(igraph::V(domain_tree)$levels)),
                     # both (a1,b1),(a2,b2) can encourage shrinkage towards the parent.
                     dmat = matrix(1,nrow=length(dsgn$ancestors[[1]]),ncol=length(dsgn$ancestors[[2]])), # (cause,domain).
                     #s1_u_zeroset = c(2:p1), # force NO diffusion in tree1.
                     s1_u_zeroset = NULL, # not force diffusion in tree1.
                     s1_u_oneset = NULL,#1,    # not force diffusion in tree1.
                     #s1_u_oneset = 1:p1,  # force diffusion in tree1.
                     #s2_cu_zeroset = rep(list(2:p2),pL1), # force NO diffusion in non-roots tree2.
                     s2_cu_zeroset = NULL,            # not force diffusion in tree2.
                     s2_cu_oneset = NULL,#rep(list(1),pL1), # not force diffusion in tree2.
                     #s2_cu_oneset = rep(list(1:p2),pL1), # force diffusion tree2.
                     tau_update_levels = list(1,1)
  )
)

# get the design output; this is needed because the function reorders the data rows:
dsgn0 <- design_doubletree(Y,leaf_ids2,
                           mytrees)

# for each tree1 leaf, look at shrinkage structure across tree2:
par(mfrow=c(ceiling(sqrt(pL1+1)),ceiling(sqrt(pL1+1))),
    mar=c(1,1,1,1))
for (u in 1:pL1){
  plot(mod$mod$vi_params$prob2[[u]],type="h",ylim=c(0,1));abline(h=0.5)
}
# look at shrinkage structure across tree1:
plot(mod$mod$vi_params$prob1,type="h",ylim=c(0,1),col="blue");abline(h=0.5)

do.call("rbind",mod$mod$vi_params$prob2)

heatmap(mod$mod$vi_params$emat[is.na(dsgn0$leaf_ids[[1]]),],Rowv=NA,Colv=NA)
# heatmap(mod$mod$vi_params$emat[!is.na(dsgn0$leaf_ids[[1]]),],Rowv=NA,Colv=NA)
# heatmap(mod$mod$vi_params$emat,Rowv=NA,Colv=NA)

# posterior means of CSMFs:
sweep(mod$mod$vi_params$dirich_mat,MARGIN = 2,colSums(mod$mod$vi_params$dirich_mat),"/")

# visualize tree2 root node class probabilities for each tree1 leaf; can change to
# nodes other than tree2 root node:
heatmap(apply(mod$mod$vi_params$mu_alpha[[1]],1,function(v) tsb(c(expit(v),1))),Rowv=NA,Colv=NA)
apply(mod$mod$vi_params$mu_alpha[[1]],1,function(v) tsb(c(expit(v),1)))

#
# CLASSIFICATION:
#
# MAP cause assignment:
xx <- mod$mod$vi_params$emat[is.na(dsgn0$leaf_ids[[1]]),]
apply(xx,1,which.max)

# true causes:
na_index <- which(is.na(dsgn0$leaf_ids[[1]]))
leaf_ids[[1]][dsgn0$all_ord[na_index]]

# domains of the observations missing tree1 leaf label:
truth <- match(leaf_ids[[1]][dsgn0$all_ord[na_index]],names(dsgn0$leaf_ids_units[[1]]))
map_nlcm <- apply(xx,1,which.max)

table(map_nlcm,truth)
sum(map_nlcm!=truth)/length(map_nlcm)

#
# top k accuracy:
#
k = 3
pred_top <- get_topk_COD(xx,k)
acc_topk(pred_top,truth)

#
# RESPONSE PROBABILITIES:
#
itemprob_list_est <-list()
for (v1 in 1:pL1){
  itemprob_list_est[[v1]] <- t(expit(Reduce("+",mod$mod$vi_params$mu_gamma[dsgn0$ancestors[[1]][[v1]]])))
}

par(mfrow=c(ceiling(sqrt(pL1+1)),ceiling(sqrt(pL1+1))),
    mar=c(1,1,1,1))
for (v1 in 1:pL1){
  image(itemprob_list_est[[v1]],main=v1)
}

#
# LATENT CLASS PROBABILITIES:
#
heatmap(mod$mod$vi_params$rmat,Rowv=NA,Colv=NA)


#
# ELBO trajectory
#
plot(mod$mod$ELBO_track,type="o",main="ELBO trajectory")
plot(mod$mod$ELBO_track[-(1:50)],type="o",main="ELBO trajectory")



# how do we avoid tree1 leaf specific class relabeling. perhaps need to consider the same set of alphas
# no matter which tree1 leaf it is.

## need to confirm the classification performance; here we have very collapsed
## tree for the causes, so the classification performance likely would not be great.







