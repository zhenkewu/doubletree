#' doubletree
#'
#' `doubletree` is designed for tree-informed Bayesian domain adaptation
#' based on nested latent class models. The special case that we address
#' is when all the observations in a particular leaf in tree2 do not have
#' observed leaf node label in tree1. More generally, the package can be extended
#' to deal with
#' \describe{
#' \item{Scenario a}{No missing leaf labels in tree1 or tree2, for all observations;
#' So this is reduced to a nested latent class model estimation problem with
#' parameter shrinkage according to the two trees.}
#' \item{Scenario b}{All missing tree1 leaf label occurred for in a single leaf node in tree2 (say v2):}
#' \itemize{
#' \item{Scenario b1: }{No observation with observed tree1 label in leaf v2;}
#' \item{Scenario b2: }{More than 1 observations have observed tree1 label in leaf v2.}
#' }
#' \item{Scenario c}{Missing tree1 leaf labels occurred for observations in 2 or more
#' leaf nodes in tree2, say the leaf node set S:}
#' \itemize{
#' \item {Sub-scenarios: }{0,1,2,... leaf/leaves in S have partially observed leaf label in tree1}
#' }
#' }
#'
#' @seealso
#' \itemize{
#' \item <https://github.com/zhenkewu/doubletree> for the source code
#' and system/software requirements to use `doubletree` for your data.
#' }
#'
#' @section main double wrapper function:
#' [nlcm_doubletree()]
#'
#' @docType package
#' @name doubletree
NULL
#> NULL
