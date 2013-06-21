\name{order_LGnodes}
\alias{order_LGnodes}
\title{Order LAGRANGE-numbered nodes so that they can be plotted in R}
\usage{
  order_LGnodes(MLsplits_LGcpp, tr = NULL,
    removechar = NULL, type = "C++", type2 = "splits")
}
\arguments{
  \item{MLsplits_LGcpp}{A data.frame containing the node
  numbers, splits, and split probabilities.}

  \item{tr}{An ape phylo object}

  \item{removechar}{The character to remove, if needed.}

  \item{type}{The type of LAGRANGE input (default C++)}

  \item{type2}{"splits" or "states"}
}
\value{
  \code{MLsplits} The splits table, ordered appropriately.
}
\description{
  What it says.
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

