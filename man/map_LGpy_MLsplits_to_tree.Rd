\name{map_LGpy_MLsplits_to_tree}
\alias{map_LGpy_MLsplits_to_tree}
\title{Take the table of ML splits and node number and map on tree (Python version)}
\usage{
  map_LGpy_MLsplits_to_tree(MLsplits_LGpy, tr,
    tiprange_names)
}
\arguments{
  \item{MLsplits_LGpy}{A data.frame containing the node
  numbers, splits, and split probabilities.}

  \item{tr}{An ape phylo object}

  \item{tiprange_names}{The geographic ranges at the tips
  (i.e. the input data)}
}
\value{
  \code{MLsplits_LGpy} A data.frame containing the node
  numbers, ML splits, and split probabilities; reordered
  for this plot
}
\description{
  Given a table of splits probabilities from
  \code{\link{LGpy_splits_fn_to_table}}, map the splits on
  the tree.
}
\details{
  See \code{\link{get_lagrange_nodenums}} for connecting
  these node numbers to APE node numbers.
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

