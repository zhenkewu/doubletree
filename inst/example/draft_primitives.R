rm(list=ls())

library(igraph)
library(doubletree)
library(MASS)
library(poLCA)
library(BayesLCA)

data("example_data_doubletree")

# second tree - over domains:
data("example_domain_edges")
domain_tree <- graph_from_edgelist(example_domain_edges, directed = TRUE)

# set the levels l*_u for nodes in the cause tree; nodes
# in the same level will share a slab variance multiplier tau_l* (times the edge length
# eminating from the parent).
igraph::V(domain_tree)$levels <- c(1,rep(2,length(igraph::V(domain_tree))-1))
igraph::E(domain_tree)$weight <- rep(1,length(E(domain_tree)))
nodes2  <- names(igraph::V(domain_tree))
leaves2 <- names(igraph::V(domain_tree)[igraph::degree(domain_tree, mode = "out") == 0])
rootnode2 <- names(igraph::V(domain_tree)[igraph::degree(domain_tree, mode = "in") == 0])
pL2 <- length(leaves2)
p2  <- length(nodes2)

# # first tree - over causes:
cause_tree <- domain_tree # for simplicity, make then the same for illustration.

igraph::V(cause_tree)$levels <- c(1,rep(2,length(igraph::V(cause_tree))-1))
igraph::E(cause_tree)$weight <- rep(1,length(E(cause_tree)))

nodes1  <- names(igraph::V(cause_tree))
leaves1 <- names(igraph::V(cause_tree)[igraph::degree(cause_tree, mode = "out") == 0])
rootnode1 <- names(igraph::V(cause_tree)[igraph::degree(cause_tree, mode = "in") == 0])
pL1 <- length(leaves1)
p1  <- length(nodes1)

# create a new doubletree list for potentially modifying the levels for the nodes
# in the two trees:
working_mytrees <- list(tree1 = cause_tree, tree2 = domain_tree)

# get lists of ancestors for each leaf_ids:
d1 <- igraph::diameter(working_mytrees[[1]],weights=NA)
# need to set weight=NA to prevent the use of edge lengths in determining the diameter.
ancestors1 <- igraph::ego(working_mytrees[[1]], order = d1 + 1, nodes = leaves1, mode = "in")
ancestors1 <- sapply(ancestors1, names, simplify = FALSE)
ancestors1 <- sapply(ancestors1, function(a, nodes) which(nodes %in% a),
                     nodes = nodes1, simplify = FALSE)
names(ancestors1) <- leaves1

# get lists of ancestors for each leaf_ids:
d2 <- igraph::diameter(working_mytrees[[2]],weights=NA)
# need to set weight=NA to prevent the use of edge lengths in determining the diameter.
ancestors2 <- igraph::ego(working_mytrees[[2]], order = d2 + 1, nodes = leaves2, mode = "in")
ancestors2 <- sapply(ancestors2, names, simplify = FALSE)
ancestors2 <- sapply(ancestors2, function(a, nodes) which(nodes %in% a),
                     nodes = nodes2, simplify = FALSE)
names(ancestors2) <- leaves2

# print("counts in simulation:")
example_data_doubletree$N_sim_mat

# FITTING MODELS:
working_leaf_ids <- vector("list",2)
working_leaf_ids[[1]] <- leaves1[example_data_doubletree$truth$true_leaf_ids[,1]]
working_leaf_ids[[1]][example_data_doubletree$truth$true_leaf_ids[,2]==1] <- NA
# working_leaf_ids[[1]][example_data_doubletree$truth$true_leaf_ids[,2]==2] <- NA
working_leaf_ids[[2]] <- leaves2[example_data_doubletree$truth$true_leaf_ids[,2]]


# Not run:------------------------------------------------
mod0 <- nlcm_doubletree(
  example_data_doubletree$Y,
  working_leaf_ids,
  working_mytrees,
  weighted_edges = c(FALSE,FALSE),
  ci_level = 0.95,
  get_lcm_by_group = FALSE,
  update_hyper_freq = 20,
  print_freq = 20,
  quiet      = FALSE,
  plot_fig   = FALSE, # <-- used?
  tol        = 1E-8,
  tol_hyper = 1E-4,
  max_iter = 40,
  nrestarts = 1,
  keep_restarts = TRUE,
  parallel = TRUE,
  log_restarts = FALSE,
  log_dir = ".",
  vi_params_init = list(),
  hyperparams_init = list(),
  random_init = FALSE,
  hyper_fixed = list(
    K=2, LD=TRUE,# number of latent classes.
    a1 = rep(20,max(igraph::V(cause_tree)$levels)),
    b1 = rep(1,max(igraph::V(cause_tree)$levels)),
    a2=matrix(1,nrow=length(ancestors1),ncol=max(igraph::V(domain_tree)$levels)),
    # <-- NB: where do we specify levels? in the tree.
    b2=matrix(10,nrow=length(ancestors1),ncol=max(igraph::V(domain_tree)$levels)),
    # both (a1,b1),(a2,b2) can encourage shrinkage towards the parent.
    dmat = matrix(1,nrow=length(ancestors1),ncol=length(ancestors2)), # (cause,domain).
    s1_u_zeroset = NULL,
    #s1_u_oneset = NULL, # not force diffusion.
    s1_u_oneset = 1:p1, # force diffusion.
    #s2_cu_zeroset = rep(list(2:p2),pL1), # force NO diffusion in non-roots tree2.
    s2_cu_zeroset = NULL,
    s2_cu_oneset = rep(list(1),pL1), # no force diffusion tree2.
    tau_update_levels = list(c(1,2),c(1,2)))
)


# get the design output; this is needed because the function reorders the data rows:
dsgn0 <- design_doubletree(example_data_doubletree$Y,working_leaf_ids,
                           working_mytrees)

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

#
# CLASSIFICATION:
#
# MAP cause assignment:
xx <- mod$mod$vi_params$emat[is.na(dsgn0$leaf_ids[[1]]),]
apply(xx,1,which.max)

# true causes:
na_index <- which(is.na(dsgn0$leaf_ids[[1]]))
example_data_doubletree$truth$true_leaf_ids[dsgn0$all_ord[na_index],1]


# domains of the observations missing tree1 leaf label:
truth <- example_data_doubletree$truth$true_leaf_ids[dsgn0$all_ord[na_index],1]
map_nlcm <- apply(xx,1,which.max)

table(map_nlcm,truth)
sum(map_nlcm!=truth)/length(map_nlcm)


#
# CSMF accuracy (the function is obtained from openVA):
# check how to get top1, top3 cause classification accuracy.
#
acc_CSMF <- rep(NA,pL2)
for (g in 1:pL2){
  acc_CSMF[g] <- openVA::getCSMF_accuracy(
    sweep(mod$mod$vi_params$dirich_mat,MARGIN = 2,
          colSums(mod$mod$vi_params$dirich_mat),"/")[,g],
    example_data_doubletree$truth$pi_mat[,g])
}
print(acc_CSMF)

#
# top k accuracy:
#
k = 1
pred_top <- get_topk_COD(xx,1)
acc_topk(pred_top,truth)

#
# RESPONSE PROBABILITIES:
#
itemprob_list_est <-list()
for (v1 in 1:pL1){
  itemprob_list_est[[v1]] <- t(expit(Reduce("+",mod$mod$vi_params$mu_gamma[ancestors1[[v1]]])))
}

par(mfrow=c(ceiling(sqrt(pL1+1)),ceiling(sqrt(pL1+1))),
    mar=c(1,1,1,1))
for (v1 in 1:pL1){
  image(itemprob_list_est[[v1]],main=v1)
}

#
# LATENT CLASS PROBABILITIES:
#

# mixture probabilities:
mixprob_list_est <-list()
for (v2 in 1:pL2){
  tmp <- t(expit(Reduce("+",mod$mod$vi_params$mu_alpha[ancestors2[[v2]]])))
  mixprob_list_est[[v2]] <- apply(tmp,2,function(v){tsb(c(v,1))})
}

par(mfrow=c(ceiling(sqrt(pL2+1)),ceiling(sqrt(pL2+1))),
    mar=c(1,1,1,1))
for (v2 in 1:pL2){
  image(mixprob_list_est[[v2]],main=v2)
}

# individual-level:
heatmap(mod$mod$vi_params$rmat,Rowv=NA,Colv=NA)

#
# ELBO trajectory
#
plot(mod$mod$ELBO_track,type="o",main="ELBO trajectory")

## show the ELBO trajectory:
# library(plotly)
# tmp_df <- data.frame(iteration = 1:length(mod$mod$ELBO_track),
#                      ELBO=mod$mod$ELBO_track)
# plot_ly(tmp_df,x=~iteration,y = ~ELBO,type = 'scatter', mode = 'lines')


# need to write a function to determine the misclassification rates.

# how do we avoid tree1 leaf specific class relabeling.
# perhaps need to consider the same set of alphas
# no matter which tree1 leaf it is.

## need to confirm the classification performance; here we have very collapsed
## tree for the causes, so the classification performance likely would not be great.

# window <- -c(1:400)
# plot(mod$mod$line_track[window,1],type="o")
# for (l in 2:17){
#   plot(mod$mod$line_track[window,l],type="o",pch=l)
# }

