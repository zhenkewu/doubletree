rm(list=ls())
library(igraph)

library(doubletree)
library(MASS)
library(poLCA)
library(BayesLCA)
library(matrixStats)

################################################################################
######## 1. Prepare the PHMRC Data
################################################################################

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

table(dsgn$leaf_ids,useNA="ifany")

# leaf labels in tree 1 has missing values; here we take PHMRC data and
# artifically remove ALL the cause-of-death labels (tree1 is cause-of-death tree) for
# a site in India (AP):

##
## The following data creation codes are commented out because it is only needed once
## during package creation:
##
# cause_ids_except_AP <- cause_ids
# cause_ids_except_AP[which(domain_ids=="AP")] <- NA
# save(cause_ids_except_AP,file="data/cause_ids_except_AP.rda",compress="xz")

data("cause_ids_except_AP") # leaf labels for observations when mapped to the first tree.
leaf_ids2 <- list(id1 = cause_ids_except_AP, id2 = domain_ids)

dsgn2     <- design_doubletree(Y,leaf_ids2,mytrees,weighted_edges)

table(dsgn2$leaf_ids,useNA="ifany") # here will list the number of observations
# with missing labels in tree1.





#############################################################
### some common fixed global quantities for model fitting
#############################################################

UPDATE_HYPER_FREQ = 10
MAX_ITER          = 1000
TOL               = 1e-6
TOL_HYPER         = 1e-3
K_FIT             = 3
nrestarts         = 1 # the number of random initializations.
CONST_a2          = 1
CONST_b2          = 1
CONST_dirich      = 1


################################################################################
######## 2. Fit the NLCM model to the PHMRC Data
################################################################################
mod <- nlcm_doubletree(
  Y[,1:20],
  leaf_ids2,
  mytrees,
  weighted_edges = c(FALSE,FALSE),
  ci_level = 0.95,
  get_lcm_by_group = FALSE,
  update_hyper_freq = UPDATE_HYPER_FREQ,
  print_freq = 10,
  quiet      = FALSE,
  plot_fig   = FALSE, # <-- used?
  tol        = TOL,
  tol_hyper = TOL_HYPER,
  max_iter = MAX_ITER,
  nrestarts = nrestarts,
  keep_restarts = TRUE,
  parallel = TRUE,
  # log_restarts = TRUE,
  # log_dir = log_dir,
  #log_restarts = FALSE,
  #log_dir = ".",
  # hyperparams_init = list(tau_1=c(1.5^2,1.5^2),
  #                         tau_2=c(1.5^2,1.5^2)),
  random_init = FALSE,
  hyper_fixed = list(
    K=K_FIT, LD=TRUE,# number of latent classes.
    a1 = rep(20,max(igraph::V(cause_tree)$levels)),
    b1 = rep(1,max(igraph::V(cause_tree)$levels)),
    a2=matrix(CONST_a2,nrow=length(leaves1),ncol=max(igraph::V(domain_tree)$levels)),
    # <-- NB: where do we specify levels? in the tree.
    b2=matrix(CONST_b2,nrow=length(leaves1),ncol=max(igraph::V(domain_tree)$levels)),
    # both (a1,b1),(a2,b2) can encourage shrinkage towards the parent.
    dmat = matrix(CONST_dirich,nrow=length(leaves1),ncol=length(leaves2)), # (cause,domain).
    s1_u_zeroset = NULL, # force NO diffusion.
    s1_u_oneset = c(1:p1), # force diffusion.
    #s2_cu_zeroset = rep(list(2:p2),pL1), # force NO diffusion in non-roots tree2.
    # ## option: root diffuse; other nodes data-adaptive:
    s2_cu_zeroset = NULL,
    s2_cu_oneset = rep(list(1),pL1), # no force diffusion tree2.
    # ## option: no grouping of domains- each one by itself:
    # s2_cu_zeroset = rep(list(1:3),pL1),
    # s2_cu_oneset = rep(list(4:9),pL1),
    tau_update_levels = list(c(1,2),c(1,2)))
)



## The following will produce results for each method, and then organize them into data sets
## that will be combined to produce final plots; check the code in Simulation I as well.


mod0=mod


lapply(mod0$mod$vi_params$prob2,plot,type="h")

# NB: need to regenerate example data in doubletree package; with new cause names.

# get the design output; this is needed because the function reorders the data rows:
dsgn0 <- design_doubletree(Y[,1:20],leaf_ids2,
                           mytrees)


# for each tree1 leaf, look at shrinkage structure across tree2:
par(mfrow=c(ceiling(sqrt(pL1+1)),ceiling(sqrt(pL1+1))),
    mar=c(1,1,1,1))
for (u in 1:pL1){
  plot(mod0$mod$vi_params$prob2[[u]],type="h",ylim=c(0,1));abline(h=0.5)
}
# look at shrinkage structure across tree1:
plot(mod0$mod$vi_params$prob1,type="h",ylim=c(0,1),col="blue");abline(h=0.5)

do.call("rbind",mod0$mod$vi_params$prob2)

heatmap(mod0$mod$vi_params$emat[is.na(dsgn0$leaf_ids[[1]]),],Rowv=NA,Colv=NA)
# heatmap(mod$mod$vi_params$emat[!is.na(dsgn0$leaf_ids[[1]]),],Rowv=NA,Colv=NA)
# heatmap(mod$mod$vi_params$emat,Rowv=NA,Colv=NA)










# posterior means of CSMFs:
summary(mod0)$pi_inf

#
# CLASSIFICATION:
#
# MAP cause assignment:
xx <- mod0$mod$vi_params$emat[is.na(dsgn0$leaf_ids[[1]]),]
apply(xx,1,which.max)

# true causes:
na_index <- which(is.na(dsgn0$leaf_ids[[1]]))
simu$truth$true_leaf_ids[dsgn0$all_ord[na_index],1]


# domains of the observations missing tree1 leaf label:
truth <- simu$truth$true_leaf_ids[dsgn0$all_ord[na_index],1]
map_nlcm <- apply(xx,1,which.max)

table(map_nlcm,truth)
sum(map_nlcm!=truth)/length(map_nlcm)

#
# CSMF accuracy (the function is obtained from openVA): <---- metric to assess.
#
A <- simu$truth$pi_mat
B <- mod0$pi_list$pi_est
acc_CSMF <- rep(NA,pL2)
for (v2 in 1:pL2){
  acc_CSMF[v2] <- getCSMFacc(B[,v2],A[,v2])
}
print(acc_CSMF)

#
# RMSE:                                                 <---- metric to assess.
#
rmse_curr <- mapply(rmse_bias_AB,
                    A = split_along_dim(array(A,dim=c(1,dim(A))),3),
                    B = split_along_dim(array(B,dim=c(1,dim(B))),3),
                    SIMPLIFY=FALSE)

#
# top k accuracy:                                       <---- metric to assess.
#
# overall:
topk <- 3
acc_k <- rep(NA,topk)
names(acc_k) <- paste("top",1:topk,sep = "_")

for (j in 1:topk){
  pred_top <- get_topk_COD(xx,j)
  acc_k[j] <- acc_topk(pred_top,truth)
}
acc_k

# by domain:
acc_k_by_domain <- matrix(NA,nrow=length(id_target_domain),ncol=topk)
colnames(acc_k_by_domain) <- paste("top",1:topk,sep = "_")
for (j in 1:topk){
  for (s in seq_along(id_target_domain)){
    g <- id_target_domain[s]
    xx <- mod0$mod$vi_params$emat[is.na(dsgn0$leaf_ids[[1]]) & dsgn0$leaf_ids[[2]]==g,]
    pred_top <- get_topk_COD(xx,j)
    na_index <- which(is.na(dsgn0$leaf_ids[[1]]) & dsgn0$leaf_ids[[2]]==g)
    truth <- simu$truth$true_leaf_ids[dsgn0$all_ord[na_index],1]
    acc_k_by_domain[s,j] <- acc_topk(pred_top,truth)
  }
}

#
# coverage:                                             <---- metric to assess.
#

covered_curr <- mapply(FUN=function(mat,true_vec){(true_vec >= mat[,"pi_cil"]) &
    (true_vec <= mat[,"pi_ciu"])},mat=summary(mod0)$pi_inf,
    true_vec =split_along_dim(simu$truth$pi_mat,2))

print(covered_curr)


#
# RESPONSE PROBABILITIES:
#
par(mfrow=c(ceiling(sqrt(pL1+1)),ceiling(sqrt(pL1+1))),
    mar=c(2,2,2,2))
for (v1 in 1:pL1){
  image(1:J,1:K,(mod0$prob_est$theta[,,v1]),main=paste0("leaf_tree1: ",v1))
}


#
# LATENT CLASS PROBABILITIES:
#

# mixture:
print(mod0,compact=TRUE) # currently unpermuted.

# individual-level:
heatmap(mod0$mod$vi_params$rmat,Rowv=NA,Colv=NA,main="rmat")

#
# ELBO trajectory
#
plot(mod0$mod$ELBO_track,type="o",main="ELBO trajectory")






















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










