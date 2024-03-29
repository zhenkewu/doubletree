rm(list=ls())
library(igraph)
library(doubletree)
library(MASS)
library(poLCA)
library(BayesLCA)

# first tree - over causes:
data("example_cause_edges")
cause_tree <- graph_from_edgelist(example_cause_edges, directed = TRUE)

igraph::V(cause_tree)$levels <- rep(1,length(igraph::V(cause_tree)))
igraph::E(cause_tree)$weight <- rep(1,length(E(cause_tree)))

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
igraph::V(domain_tree)$levels <- rep(1,length(igraph::V(domain_tree)))
igraph::E(domain_tree)$weight <- rep(1,length(E(domain_tree)))


nodes2  <- names(igraph::V(domain_tree))
leaves2 <- names(igraph::V(domain_tree)[igraph::degree(domain_tree, mode = "out") == 0])
rootnode2 <- names(igraph::V(domain_tree)[igraph::degree(domain_tree, mode = "in") == 0])
pL2 <- length(leaves2)
p2  <- length(nodes2)
V(domain_tree)$levels <- rep(1,p2)
E(domain_tree)$weight <- rep(1,length(E(domain_tree)))

mytrees <- list(tree1 = cause_tree, tree2 = domain_tree)
###############################################################################
## Begin simulating data using the two trees above.
###############################################################################
n     <- 5000
K     <- 2
J     <- 18 # J =168 in the data

fracs_leaves1 <- c(5,rep(1,pL1-1))
pi_mat        <- matrix(fracs_leaves1/sum(fracs_leaves1),nrow=pL1,ncol=pL2)

# create tree1 leaf level class-specific response probabilities:
itemprob0 <- rbind(rep(rep(c(0.9, 0.9), each = 1),9),
                   rep(rep(c(0.1, 0.1), each = 1),9))
gamma_mat_list <- list(logit(itemprob0))

for (u in 1:(p1-1)){ # random increments on random columns
  increment_mat <- matrix(rnorm(J*K,0,1),nrow=K,ncol=J)
  ind_j         <- sample(1:J,floor(J/2))
  increment_mat[,ind_j] <- 0
  gamma_mat_list <- append(gamma_mat_list,list(increment_mat))
}
# matrix(rnorm(J*K,0,0.5),nrow=K,ncol=J),
# matrix(rnorm(J*K,0,0.5),nrow=K,ncol=J)),
# rep(list(matrix(0,nrow=K,ncol=J)),p1-3))


# matrix for each of the pL1 leaves in tree1.
#########################################################################
## diffusion along tree1 for itemprob_list
#########################################################################

# get lists of ancestors for each leaf_ids:
d1 <- igraph::diameter(mytrees[[1]],weights=NA)
# need to set weight=NA to prevent the use of edge lengths in determining the diameter.
ancestors1 <- igraph::ego(mytrees[[1]], order = d1 + 1, nodes = leaves1, mode = "in")
ancestors1 <- sapply(ancestors1, names, simplify = FALSE)
ancestors1 <- sapply(ancestors1, function(a, nodes) which(nodes %in% a),
                     nodes = nodes1, simplify = FALSE)
names(ancestors1) <- leaves1
# the list of item response probability
itemprob_list <-list()
for (v1 in 1:pL1){
  curr_mat <- itemprob0
  ind_j <- sample(1:J,floor(J/2))
  for (j in ind_j){
    curr_mat[,j] <- itemprob0[sample(1:nrow(curr_mat)),j]
  }
  # itemprob_list[[v1]] <- expit(Reduce("+",gamma_mat_list[ancestors1[[v1]]]))
  itemprob_list[[v1]] <- curr_mat
}

for (v1 in 1:pL1){
  image(itemprob_list[[v1]],main=v1)
}

############################################################################
## diffusion along tree2 for alpha_mat_list
############################################################################
# specify the nodes that have non-trivial alpha^cu_k, this was called
# xi^cu_k, because xi^cu_k = s_cu*alpha^cu_k and s_cu = 1 if we set it in the simulation.
alpha0_mat <- myrdirich(pL1,rep(1,K)) # root
alpha_mat_list <- vector("list",pL1)
for (l in 1:pL1){
  # alpha0 <- c(0.6,0.3,0.1) # root
  alpha0 <- alpha0_mat[l,]
  # for each leaf node in tree1:
  alpha_mat = rbind(logit(prob2stick(alpha0)[-K]),
                    0,
                    0,
                    matrix(0,nrow=p2-3,ncol=K-1)
  )
  alpha_mat_list[[l]] <- alpha_mat
}
# repeat for every leaf node in tree1:
# alpha_mat_list <- rep(list(alpha_mat),pL1)

# get lists of ancestors for each leaf_ids:
d2 <- igraph::diameter(mytrees[[2]],weights=NA)
# need to set weight=NA to prevent the use of edge lengths in determining the diameter.
ancestors2 <- igraph::ego(mytrees[[2]], order = d2 + 1, nodes = leaves2, mode = "in")
ancestors2 <- sapply(ancestors2, names, simplify = FALSE)
ancestors2 <- sapply(ancestors2, function(a, nodes) which(nodes %in% a),
                     nodes = nodes2, simplify = FALSE)
names(ancestors2) <- leaves2

# calculate the class probabilities for all leaf nodes in tree2; each leaf node
# should have a K-dim vector that sums to one; Some nodes may share
# the same set of K-dim probability vector, others may differ. There are
# one or more groups of leaf nodes with distinct K-dim probability vectors.
# Note the branch lengths may also be used here.
lambda_mat_list <- list() # will be a list of length pL1, each being K by pL2.
for (v1 in 1:pL1){
  lambda_mat_list[[v1]] <- matrix(NA,nrow=K,ncol=pL2)
  for (v2 in 1:pL2){
    tmp <- colSums(alpha_mat_list[[v1]][ancestors2[[v2]],,drop=FALSE])
    lambda_mat_list[[v1]][,v2] <- tsb(c(expit(tmp),1))
  }
}

# s = c(1, 1,1,0,0, rep(0,pL)) # effective nodes
example_data_doubletree <- simulate_nlcm_doubletree(n,mytrees,itemprob_list,lambda_mat_list,pi_mat)
print("counts in simulation:")
example_data_doubletree$N_sim_mat
# save the simulated data to the R package for illustration:
# save(example_data_doubletree, file = "data/example_data_doubletree.rda", compress = "xz

# FITTING MODELS:
curr_leaf_ids <- vector("list",2)
curr_leaf_ids[[1]] <- leaves1[example_data_doubletree$truth$true_leaf_ids[,1]]
curr_leaf_ids[[1]][example_data_doubletree$truth$true_leaf_ids[,2]==1] <- NA # <--- make this work.
curr_leaf_ids[[2]] <- leaves2[example_data_doubletree$truth$true_leaf_ids[,2]]

mod <- nlcm_doubletree(
  example_data_doubletree$Y,curr_leaf_ids,
  example_data_doubletree$truth$mytrees,weighted_edges = c(FALSE,FALSE),
  ci_level = 0.95,
  get_lcm_by_group = FALSE,
  update_hyper_freq = 20,
  print_freq = 20,
  quiet      = FALSE,
  plot_fig   = FALSE, # <-- used?
  tol        = 1E-8,
  tol_hyper = 1E-4,
  max_iter = 1000,
  nrestarts = 1,
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
                          u_sd_frac = 0.2, # for logit of probs
                          psi_sd_frac = 0.2,
                          phi_sd_frac = 0.2),
  hyper_fixed = list(K=3, # number of latent classes.
                     a1 = rep(20,max(igraph::V(cause_tree)$levels)),
                     b1 = rep(1,max(igraph::V(cause_tree)$levels)),
                     a2=matrix(1,nrow=length(ancestors1),ncol=max(igraph::V(domain_tree)$levels)),
                     # <-- NB: where do we specify levels? in the tree.
                     b2=matrix(10,nrow=length(ancestors1),ncol=max(igraph::V(domain_tree)$levels)),
                     # both (a1,b1),(a2,b2) can encourage shrinkage towards the parent.
                     dmat = matrix(1,nrow=length(ancestors1),ncol=length(ancestors2)), # (cause,domain).
                     s1_u_zeroset = NULL,
                     s1_u_oneset = 1:p2,
                     s2_cu_zeroset = rep(list(2:p2),pL1),
                     # s2_cu_zeroset = NULL,
                     s2_cu_oneset = rep(list(1),pL1),
                     tau_update_levels = list(1,1)
  )
)

# get the design output; this is needed because the function reorders the data rows:
dsgn0 <- design_doubletree(example_data_doubletree$Y,curr_leaf_ids,
                           example_data_doubletree$truth$mytrees)

# root node response profiles:
image(expit(mod$mod$vi_params$mu_gamma[[1]]))

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

# MAP cause assignment:
xx <- mod$mod$vi_params$emat[is.na(dsgn0$leaf_ids[[1]]),]
apply(xx,1,which.max)

# true causes:
na_index <- which(is.na(dsgn0$leaf_ids[[1]]))
example_data_doubletree$truth$true_leaf_ids[dsgn0$all_ord[na_index],1]

table(apply(xx,1,which.max),example_data_doubletree$truth$true_leaf_ids[dsgn0$all_ord[na_index],1])

# posterior means of CSMFs:
sweep(mod$mod$vi_params$dirich_mat,MARGIN = 2,colSums(mod$mod$vi_params$dirich_mat),"/")

# visualize tree2 root node class probabilities for each tree1 leaf; can change to
# nodes other than tree2 root node:
heatmap(apply(mod$mod$vi_params$mu_alpha[[1]],1,function(v) tsb(c(expit(v),1))),Rowv=NA,Colv=NA)

# domains of the observations missing tree1 leaf label:
curr_leaf_ids[[2]][1:20]
truth <- example_data_doubletree$truth$true_leaf_ids[dsgn0$all_ord[na_index],1]
map_nlcm <- apply(xx,1,which.max)
sum(map_nlcm!=truth)/length(map_nlcm)



# need to write a function to determine the misclassification rates.

# how do we avoid tree1 leaf specific class relabeling. perhaps need to consider the same set of alphas
# no matter which tree1 leaf it is.

## need to confirm the classification performance; here we have very collapsed
## tree for the causes, so the classification performance likely would not be great.





