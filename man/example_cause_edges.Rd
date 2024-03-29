% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doubletree_example_data.R
\docType{data}
\name{example_cause_edges}
\alias{example_cause_edges}
\title{Edgelist for constructing an example cause hierarchy.}
\format{
The edge list is a matrix with 41 rows with 2 columns:
the left column
represents parent nodes or categories, and the right column
represents children or subcategories of the parents.
There is one row for every parent-child pair.
}
\usage{
example_cause_edges
}
\description{
A matrix representing edges
from the cause hierarchy
See vignette("doubletree") for details of how to construct
a tree from this edgelist.
}
\keyword{datasets}
