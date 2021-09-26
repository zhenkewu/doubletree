#' Edgelist for constructing an example cause hierarchy.
#'
#' A matrix representing edges
#' from the cause hierarchy
#' See vignette("doubletree") for details of how to construct
#' a tree from this edgelist.
#'
#' @format The edge list is a matrix with 41 rows with 2 columns:
#' the left column
#' represents parent nodes or categories, and the right column
#' represents children or subcategories of the parents.
#' There is one row for every parent-child pair.
"example_cause_edges"


#' Edgelist for constructing an example domain hierarchy.
#'
#' A matrix representing edges
#' from the domain hierarchy
#' See vignette("doubletree") for details of how to construct
#' a tree from this edgelist.
#'
#' @format The edge list is a matrix with 8 rows with 2 columns:
#' the left column
#' represents parent nodes or categories, and the right column
#' represents children or subcategories of the parents.
#' There is one row for every parent-child pair.
"example_domain_edges"

#' example cause ids
#'
#' a vector of character strings
#'
#' @format vector; the length is 7841
"cause_ids"

#' example cause ids with NAs
#'
#' a vector of character strings
#'
#' @format vector; the length is 7841
"cause_ids_except_AP"

#' example cause ids
#'
#' a vector of character strings
#'
#' @format vector; the length is 7841
"domain_ids"

#' example VA data from PHMRC - adult
#'
#' a matrix of data with potential missing data
#' `1`: yes; `0`: no; `NA`: missing.
#' The rows of this data match with `cause_ids`
#' and `domain_ids`
#'
#' @format matrix; the length is 7841 by 168
"X0"


#' simulated data (K=2,J=20,N=1000)
#'
#' a matrix of simulated data
#' `1`: yes; `0`: no
#'
#' @format matrix; 1000 by 20
"example_data_doubletree"
