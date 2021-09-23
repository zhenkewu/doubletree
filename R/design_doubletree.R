#' Organize the data around two rooted weighted trees
#'
#' Prepare the data for fitting nested latent class models with two separate tree
#' structures (i.e., hierarchies) over two sets of leaf labels. Each observation
#' is associated with a particular combination of label values, one from each set of leaves.
#' For example, we have domain and cause hierarchies in verbal autopsy applications.
#' # NB:
#' need to consider `edge_lengths` and `weighted` under missing labels for tree1.
#'
#' @param Y `N` by `J` binary data matrix; rows for subjects,
#' columns for features; missing entries, if present, are encoded by `NA`s. Note
#' that the rows of `Y` will be reordered twice, once according to leaf nodes/missing
#' in tree1, the second time according to the leaf nodes in tree2.
#'
#' @param leaf_ids A list containing two elements. For example,
#' the first can be a vector of character strings for leaf
#' nodes for each observation, representing the leaf nodes
#' in `tree1`; similarly for the second tree; `NA` represents for missing leaf info.
#' For each observation, the pair of labels indicates
#' the leaf memberships in the two trees contained in `mytrees`, respectively.
#' For example, in verbal autopsy (VA) applications, for data in the source domains,
#' we must have both leaf ids observed; for data in the target domain, we can have the leaf id
#' in the first tree (i.e., cause tree) `NA`, indicating an unobserved cause
#' of death (hence unknown to analysts which leaf in the cause tree should an
#' observation be placed). In this package, we only allow `NA` for leaf labels in tree1;
#' tree2's leaves represent domains, which must be known when doing domain adaptation.
#' Currently the `NA` in tree1, if present, can only contain ALL subjects from a single leaf node
#' in tree2. NB: Extensions to deal with CODs that are partially observed in the
#' target domain need additional work...
#'
#' @param mytrees A list of two elements: `tree1`,`tree2`; both are
#' `igraph` objects. They may contain attributes, such as node, edge, edge lengths.
#' (NB: need refinement)
#'
#' @param weighted_edges a vector of logical values, indicating whether to use weighted
#' edges in the two trees; default to `c(FALSE,FALSE)`, i.e., not using weighted edges
#' and assuming the edges in the trees are unit lengths.
#'
#' @return A list of three elements that contains preprocessed data that
#' are organized around the two trees;
#' Y,A,leaf_ids_units,leaf_ids,leaf_ids_nodes,ancestors,h_pau,levels,
#  v_units,subject_id_list,mytrees
#'
#' Also see the arguments of [initialize_nlcm_doubletree()] for details,
#' because the present function is used to set the stage for initialization.
#'
#' \itemize{
#' \item `resA2`
#' \itemize{
#' \item `A` A matrix of `p1` by `p1`; each column contains some `1`s, indicating
#'       the node in that column is an ancestor of the node represented in the row.
#'      Ancestor of a node include that node itself.
#' \item `A_leaves` A matrix of `pL1` by `p1`; A submatrix of `A` that represents
#' the ancestor matrix but only for leaves
#' \item `leaf_ids` A vector of `N`integers; ordered by the leaves as
#'              specified by `mytrees[[1]]`; `NA` represents missing leaf labels.
#' \item `leaf_ids_units` A list of length `pL1`, each element
#'            is a vector of subject ids belonging to each leaf node;
#'            the list will be of length `pL1+1` if missing leaf labels are
#'            present, in which case the final element `NA_tree1`the ids of the observations
#'            with missing leaf labels. NB: can it accommodate subjects from distinct
#'            leaf nodes in tree2?
#' \item `leaf_ids_nodes` leaf descendants based on the tree1 structure;
#'       a list of length `p1`
#' \item `ancestors` a list of length `pL1`,
#'      each element is the vector of ancestors (between `1` and `p1`; id is among
#'      all nodes); NB: what is the mapping between integers and node meanings?
#' \item `edge_lengths` a list of length `pL1`,
#'       each element is a numeric vector of lengths of edges in the path from the root node
#'       to the leaf. It is computed based on `E(mytrees[[1]])$weight`. It is `NULL`
#'       if `E(mytrees[[1]])$weight` is `NULL`
#' \item `h_pau` a numeric vector of length `p1`; each value is
#'       the edge length from a node`u` to its parent (equals `1` if `u` is a root node).
#'       This vector by default is all `1`s. If `weighted_edge=TRUE`, `h_pau`
#'       is set to `E(mytree[[1]])$weight`, the input edge weights.
#' \item `v_units` (`leaf_ids`) integers; just without the names.
#' \item `subject_id_list` a list of length `p1`; each element is a vector of
#'       subject ids belonging to the leaf descendants of node `u` (internal or leaf node)
#' }
#'
#' \item `resB`
#' \itemize{
#' \item `Y` A matrix of `N` by `J`; binary measurements with rows ordered by
#'       leaf groups  (`leaf_ids[[1]]` and `leaf_ids[[2]]`)
#' \item `A` A matrix of `p2` by `p2`; each column contains some `1`s, indicating
#' the node in that column is an ancestor of the node represented in the row.
#' Ancestor of a node include that node itself.
#' \item `A_leaves` A matrix of `pL2` by `p2`; A submatrix of `A` that represents
#' the ancestor matrix but only for leaves
#' \item `leaf_ids` A vector of `N`integers; ordered by the leaves as
#'              specified by `mytree[[2]]`;
#' \item `leaf_ids_units` A list of length `pL2`, each element
#'            is a vector of subject ids belonging to each leaf node
#' \item `leaf_ids_nodes` leaf descendants based on tree2 structure;
#'            a list of length `p2`, each element
#'            is a vector of integers (between `1` and `pL2`; id
#'            is only for leaf nodes) indicating the leaf nodes
#' \item `ancestors` a list of length `pL2`,
#'      each element is the vector of ancestors (between `1` and `p2`; id is among
#'      all nodes)
#' \item `edge_lengths` a list of length `pL2`,
#'       each element is a numeric vector of edge lengths from the root node
#'       to the leaf. It is computed based on `E(mytrees[[2]])$weight`. It is `NULL`
#'       if `E(mytrees[[2]])$weight` is `NULL`
#' \item `h_pau` a numeric vector of length `p2`; each value is
#'       the edge length from a node`u` to its parent (equals `1` if `u` is a root node).
#'       This vector by default is all `1`s. If `weighted_edge=TRUE`, `h_pau`
#'       is set to `E(mytree[[2]])$weight`, the input edge weights.
#' \item `v_units` (identical to `leaf_ids`) just without the names.
#' \item `subject_id_list` a list of length `p2`; each element is a vector of
#'       subject ids that are in the leaf descendants of node `u` (internal or leaf node)
#' }
#'
#' \item `all_ord` A permuted vector of consecutive integers up to `N`;
#' for each row (i.e., position) in the output `resB$Y`,
#' the corresponding row in the original input data `Y`.
#' }
#'
#' @example
#' inst/example/example_design_doubletree.R
#'
#' @seealso [lotR::design_tree()].
#'
#' @export
#' @import igraph
#' @importFrom lotR design_tree
design_doubletree <- function(Y,leaf_ids,mytrees,
                              weighted_edges = c(FALSE,FALSE)){
  for (j in 1:2){# quick checks:
    if (!is.character(leaf_ids[[j]])) stop(paste0("[doubletree] `leaf_ids[[",j, "]]` is not a character object."))
    if (!igraph::is.igraph(mytrees[[j]])) stop(paste0("[doubletree] 'mytrees[[",j,"]]' is not a graph object."))
    if (!igraph::is.directed(mytrees[[j]])) stop(paste0("[doubletree] 'mytrees[[",j,"]]' is not directed."))
  }

  if(sum(rowSums(!is.na(Y))==0)>=1){stop("[doubletree] `Y` contains at least one row with no observed responses. Delete these rows and retry.")}

  # The importFrom above means we need to either upload lotR to CRAN, or to have
  # lotR available via Github and to doubletree package (need edits in DESCRIPTION)

  # split data into two parts (row blocks): rows with observed causes, and
  # rows with unobserved causes:
  if (sum(is.na(leaf_ids[[1]]))>0){ # if there are missing labels in tree1, e.g., missing
    # cause of death for some subjects:
    cat("[doubletree] leaf_ids[[1]] has missing value;\n counts of missing or not are:")
    tb <- table(is.na(leaf_ids[[1]]),leaf_ids[[2]])
    print(tb) # tabulate NA or not by the leaf nodes in tree2.
    if (sum(tb["TRUE",]>0)>1){
      warning("[doubletree] subjects with missing leaf labels in tree1 are scattered
              across more than 1 leaf node in tree2.")
    }
    if (sum(tb["FALSE",tb["TRUE",]>0]>0)>=1){
      cat(paste0("[doubletree] leaf node(s) in tree2 that contain some subjects with missing tree1 leaf labels
             and others with observed tree1 leaf labels: ",colnames(tb)[tb["TRUE",]>0 & tb["FALSE",]>0],".\n"))
    }
    cat("\n")

    ind_na <- which(is.na(leaf_ids[[1]])) ## subjects with missing tree1 label ids.

    #results for tree1 (e.g., cause tree) related quantities:
    resA_block1 <- design_tree(Y=Y[-ind_na,,drop=FALSE],
                               leaf_ids[[1]][-ind_na],
                               mytree = mytrees[[1]], # NB: do we need  every leaf to have at least an observation?
                               weighted_edge = weighted_edges[[1]])
    resA   <- resA_block1

    # Append observations with missing node info in tree1, e.g., unobserved causes:
    # The following has NA values!
    resA$Y <- rbind(resA_block1$Y,Y[ind_na,,drop=FALSE])

    # leaf_ids[[2]] needs reordering:
    leaf_ids_tree2 <- c((leaf_ids[[2]][-ind_na])[resA_block1$ord],#reorder according to block1 ord.
                        leaf_ids[[2]][ind_na]) #append the leaf_ids in tree2 for observations missing labels in tree1.

    resA$leaf_ids <- c(resA_block1$leaf_ids,rep(NA,length(ind_na))) # mising leaf ids.
    resA$v_units  <- c(resA_block1$v_units,rep(NA,length(ind_na))) # same as above; just integers.

    # ids for those without observed labels in tree1 (an element NA_tree1 is
    # appended at the end of the list):
    resA$leaf_ids_units$NA_tree1  <- (nrow(resA$Y)-length(ind_na)+1):nrow(resA$Y)
    resA$subject_id_list$NA_tree1 <- (nrow(resA$Y)-length(ind_na)+1):nrow(resA$Y)
    # this is for all nodes; the above is for leaves only.

    # result for domain-tree related quantities:
    resB <- design_tree(Y=resA$Y,# already reordered for subjects according to leaves in tree1.
                        leaf_ids_tree2, # reordered for subjects (see resA_block1 and above).
                        mytree = mytrees[[2]], # will use this info to reorder nodes, producing resB$ord.
                        weighted_edge = weighted_edges[[2]])
    # from Y to resB$Y:
    all_ord <-c((1:nrow(Y))[-ind_na][resA_block1$ord],ind_na)[resB$ord]
  } else { # if no missing labels in tree1:
    # results for cause-tree related quantities:
    resA <- design_tree(Y=Y,leaf_ids[[1]],mytree = mytrees[[1]],
                        weighted_edge = weighted_edges[[1]])
    # result for domain-tree related quantities:
    resB <- design_tree(Y=resA$Y,# already reordered for subjects according to leaves in tree1.
                        leaf_ids[[2]][resA$ord], # reordered for subjects (see resA).
                        mytree = mytrees[[2]],   # will use this info to reorder nodes, producing resB$ord.
                        weighted_edge = weighted_edges[[2]])
    # from Y to resB$Y:
    all_ord <- (1:nrow(Y))[resA$ord][resB$ord]
  }
  # up to now, the observed data now has been reordered!
  resB$inv_ord <- order(resB$ord) # inverse mapping of resB$ord; see logic
  # below.

  #
  # need to match to the ids in resB:
  #
  resA2 <- resA
  # need to map id used in resA$Y to the row ids in resB$Y:
  # Y, leaf_ids (same as v_units), leaf_ids_units,
  # v_units (same as leaf_ids; just integers), subject_id_list,
  resA2$Y <- resA$Y[resB$ord,,drop=FALSE]
  resA2$leaf_ids <- resA$leaf_ids[resB$ord]
  resA2$v_units  <- resA$v_units[resB$ord]
  for (l in 1:length(resA$leaf_ids_units)){ # could be of length pL if no missing label in tree1;
    # or could be of length pL+1 if there are observations with missing labels in tree1.
    resA2$leaf_ids_units[[l]] <- resB$inv_ord[resA$leaf_ids_units[[l]]]
    # The logic of the above formula.
    # Each element in resA$leaf_ids_units[[l]] is representing the old position
    # before we reorder the rows by resB$ord (using leaf ids in tree2). In the
    # final output, however the rows in resB$Y is after the ordering (using leaf ids
    # in tree2), the ids are new position, so we need to ask "for an id j in the old
    # position, what is it's new position id?". Let old id be i, new id be j.
    # sigma(i) --> j; then we need inverse_sigma(j) --> i. resB$inv_ord is this
    # inverse mapping.
    #
    # Because l indexes the leaf ids in tree 1, which does not get reordered, so we
    # just keep the same indexing l.
  }

  for (l in 1:length(resA$subject_id_list)){# could be of length p if no missing label in tree1;
    # or could be of length p+1 if there are observations with missing labels in tree1.
    # the subject ids for each of all the nodes (including leaves and internal nodes).
    resA2$subject_id_list[[l]] <- resB$inv_ord[resA$subject_id_list[[l]]]
  }

  ## quick check; should be zero (the ids for subjects missing labels in tree1).
  sum(resA$subject_id_list[[length(resA$subject_id_list)]]-
        resA$leaf_ids_units[[length(resA$leaf_ids_units)]])


  if (sum(leaf_ids[[1]][all_ord]!=names(resA2$leaf_ids),na.rm = TRUE)!=0 |
      sum(leaf_ids[[2]][all_ord]!=names(resB$leaf_ids),na.rm = TRUE)!=0 |
      sum(Y[all_ord,]-resB$Y,na.rm=TRUE)!=0){
    stop("[doubletree] `design_doubletree` error.")
  }

  # remove intermediate useless information:
  resA2$Y    <- resA2$Z_obs <- resA2$ord <- NULL
  resB$Z_obs <- resB$ord <- resB$inv_ord <- NULL

  #make_list(resA2,resB,all_ord) # resA2 accounts for the ordering of subjects by tree2.

  Y <- resB$Y
  A <- list(resA2$A,resB$A)
  leaf_ids_units <- list(resA2$leaf_ids_units,resB$leaf_ids_units)
  leaf_ids <- list(resA2$leaf_ids,resB$leaf_ids)
  leaf_ids_nodes <- list(resA2$leaf_ids_nodes,resB$leaf_ids_nodes)
  ancestors <- list(resA2$ancestors,resB$ancestors)
  h_pau <- list(resA2$h_pau,resB$h_pau)
  levels    <- list(resA2$levels,resB$levels)
  # the following might be redundant:
  v_units <- list(resA2$v_units,resB$v_units)
  subject_id_list <- list(resA2$subject_id_list,resB$subject_id_list)

  ## NB: may need to add D_obs in the returned list.

  mytrees <- list(resA2$mytree,resB$mytree)
  make_list(Y,A,leaf_ids_units,leaf_ids,leaf_ids_nodes,ancestors,h_pau,levels,
            v_units,subject_id_list,mytrees,all_ord) # currently not needed in initialization; but useful for updating VI parameters, hyperaprameters....
  # can consider additional information in resA2 and resB and allord into
  # the final returned list; these might be useful for visualizing results.
}



##In addition, `resB` has an additional element `inv_ord`, which maps
## a unit id for a row of `resA$Y` to a corresponding row in `resB$Y`.
## In particular, `resB$inv_ord=order(resB$ord)`. The following
## elements in `resA2` (the first element in the returned list) have been applied
## with this inverse mapping so that `resA2$Y` is identical to `resB$Y`:
## `Y, leaf_ids, leaf_ids_units,v_units (same as leaf_ids!), subject_id_list`.
