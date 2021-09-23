#' simulate data and subject-specific
#' class indicators from latent class models
#'
#' @param n sample size
#' @param itemprob item probabilities; K by J.
#' @param classprob class probabilities; K.
#'
#' @return a list; `x`: data; `z` a vector of integer indicators
#' of class membership.
#'
#' @examples
#'
#' doubletree_blca(10,matrix(0.5,nrow=3,ncol=20),c(0.98,0.01,0.01))
#'
#' @seealso [BayesLCA::blca()]
#' @importFrom stats rmultinom runif
#' @export
doubletree_blca <- function (n, itemprob = 0.5, classprob = 1)
{
  itemprob <- as.matrix(itemprob)
  G <- nrow(itemprob)
  M <- ncol(itemprob)
  x <- matrix(runif(n * M), nrow = n)
  classvec <- as.vector(rmultinom(1, n, prob = classprob))
  ind <- c(0, cumsum(classvec))
  z <- rep(NA,n)
  for (g in 1:G) {
    if (ind[g]<ind[g+1]){
      x[(ind[g] + 1):ind[g + 1], ] <- t(t(x[(ind[g] +
                                               1):ind[g + 1], ]) < itemprob[g, ]) * 1
      z[(ind[g] + 1):ind[g + 1]] <- g
    }
  }
  make_list(x,z)
}

#' Simulate data from double-tree-structured nested latent class models
#'
#' The observations belong to leaves in tree2 that may further form a few groups, each
#' with its own K-class probabilities. In particular,
#' 1) The `K` by `J` response profile
#' probabilities are at the leaves of tree1 - they may form groups of distinct
#' numbers of leaves in tree1 (e.g., cause groups with distinct response probability
#' profiles).
#' 2) For each leaf in tree1, the K-vector class probability vectors may vary across
#' the leaves in tree2; they may further form a small number of leaf groups in tree2
#' with distinct probability vector values.
#'
#' @param n total sample size
#' @param itemprob_list a list of length `pL1`; each element is a `K` by `J` matrix:
#' `K`-classes and `J` items
#' @param mytrees a list of two trees; see [design_doubletree()]
#' @param lambda_mat_list a list of length `pL1`; each element is a `K` by `pL2` matrix;
#' contains K-dim class probabilities for  `pL2` leaf nodes; so each column sums to 1.
#' @param  pi_mat a `pL1` by `pL2` matrix of probabilities,
#' with each column corresponding to a leaf node in tree2. Each column
#' sums to one, representing the fractions of subjects assigned to each of the `pL1` leaves
#' in tree1; the vectors of fractions do not need to be identical across leaf nodes in
#' tree2, e.g., the cause-specific mortality fractions may vary by domains.
#' @param balanced by default is `TRUE` to uniformly assign observations to the leaf nodes;
#' otherwise set this to `FALSE`.
#' @param ratio only used if `balance=FALSE`;
#' for a pair of leaves in tree2; the sample ratios (larger vs smaller ones);
#' in the event of an odd number of leaves, the smaller leaf in the pair is kept.
#' This ratio is applied once to get total sample sizes in each leaf of tree2.
#'
#' @return a list. In particular, the present function does not make the `true_leaf_ids`
#' missing for a subset of subjects in tree1. So to mimic a situation with missing
#' tree1 leaf labels, additional creation of missing tree1 label indicators are needed!
#' In addition, this function simulates fully observed responses; so to mimic
#' situations where missing responses would happen, additional missing response
#' indicators need to be created.
#' \describe{
#' \item{Y}{observations}
#' \item{truth}{a list that contains the simulation truth:
#' \itemize{
#' \item `true_leaf_ids` all the leaf memberships to both trees
#' \item `Z` true class indicators for all observations
#' }
#' }
#' \item{N_sim_mat}{pL1 by pL2 count matrix}
#' }
#'
#' @example
#' inst/example/example_simulate_doubletree.R
#'
#' @seealso [BayesLCA::blca()]
#' @importFrom stats rmultinom runif
#' @importFrom igraph V
#' @export
simulate_nlcm_doubletree <- function (n,
                                      mytrees,
                                      itemprob_list, # length pL1; K by J
                                      lambda_mat_list,  # length pL1; K by pL2
                                      pi_mat, # pL1 by pL2
                                      balanced=TRUE,ratio=4)
{
  # a few quick calculations:
  K <- nrow(itemprob_list[[1]])
  nodes1  <- names(V(mytrees[[1]]))
  nodes2  <- names(V(mytrees[[2]]))
  leaves1 <- names(V(mytrees[[1]])[degree(mytrees[[1]], mode = "out") == 0])
  leaves2 <- names(V(mytrees[[2]])[degree(mytrees[[2]], mode = "out") == 0])
  pL1 <- length(leaves1)
  pL2 <- length(leaves2)
  p1  <- length(V(mytrees[[1]]))
  p2  <- length(V(mytrees[[2]]))

  # simulate number of observations for each leaf node:
  prob_vec <- rep(1/pL2,pL2) #balanced by default
  if (!balanced){
    prob_vec <- rep(c(1,ratio), c(floor(pL2/2),pL2-floor(pL2/2))) # currently 1:4 ratio.
    prob_vec <- sample(prob_vec/sum(prob_vec))
  }
  N_sim2 <- sample(1:pL2,size=n-2*pL1*pL2,prob=prob_vec,replace=TRUE) # besides the two observations at least in each cell.
  N_sim2 <- as.integer(table(sort(N_sim2)))

  # for each leaf in tree2; get the sample sizes across pL1 leaves in tree1:
  N_sim_list <- mapply(FUN=function(v2,pvec){c(rep(1:pL1,each=2), # two observations at least in each cell.
                                               sample(1:pL1,size=N_sim2[v2],pvec,replace = TRUE))},
                       v2=1:pL2,pvec=split_along_dim(pi_mat,2)) # for each leaf in tree2, assign leaf nodes in tree2 for N_sim2[v2]+2*pL1 subjects.

  N_sim_mat <- mapply(FUN="table",N_sim_list)

  # simulate the multivariate responses
  Y_sim <- vector("list",pL1)
  # simulate the class memberships
  Z_sim <- vector("list",pL1)
  curr_leaves_sim <- list()
  for (v1 in 1:pL1){
    for (v2 in 1:pL2){
      n_cell <- N_sim_mat[v1,v2]
      simu <- doubletree_blca(n_cell, itemprob = itemprob_list[[v1]], classprob = lambda_mat_list[[v1]][,v2])
      if (is.null(Y_sim[[v1]])){
        Y_sim[[v1]] <- simu$x
        Z_sim[[v1]] <- simu$z
        curr_leaves_sim[[v1]] <- rep(v2,n_cell)
      } else{
        Y_sim[[v1]] <- rbind(Y_sim[[v1]],simu$x)
        Z_sim[[v1]] <- c(Z_sim[[v1]],simu$z)
        curr_leaves_sim[[v1]] <- c(curr_leaves_sim[[v1]],rep(v2,n_cell))
      }
    }
  }
  Y <- do.call("rbind",Y_sim)
  true_leaf_ids <-  cbind(rep(1:pL1,unlist(lapply(curr_leaves_sim,length))),unlist(curr_leaves_sim))
  Z <- do.call("c",Z_sim)
  truth <- make_list(true_leaf_ids,Z,mytrees,itemprob_list, lambda_mat_list,pi_mat)
  make_list(Y,truth,N_sim_mat)
}

