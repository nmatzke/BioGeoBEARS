\name{get_leftright_nodes_matrix_from_results}
\alias{get_leftright_nodes_matrix_from_results}
\title{Make a table of the Right and Left nodes descending from each node}
\usage{
  get_leftright_nodes_matrix_from_results(tr,
    results_object, nodes)
}
\arguments{
  \item{tr}{An ape phylo object}

  \item{results_object}{The results from a BioGeoBEARS ML
  search.}

  \item{nodes}{A list of internal node numbers for tree
  \code{tr}.}
}
\value{
  leftright_nodes_matrix A table with the Right, the Left,
  nodes
}
\description{
  This table shows the Right, then Left, descendant
  nodenums for each node. This gets used later to plot
  splits at corners.
}
\note{
  Go BEARS!
}
\examples{
test=1
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
  \url{https://code.google.com/p/lagrange/}

  Matzke_2012_IBS

  ReeSmith2008
}
\seealso{
  \code{\link{get_lagrange_nodenums}},
  \code{\link{LGpy_splits_fn_to_table}},
  \code{\link{LGcpp_splits_fn_to_table}}
}

