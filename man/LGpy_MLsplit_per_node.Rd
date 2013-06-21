\name{LGpy_MLsplit_per_node}
\alias{LGpy_MLsplit_per_node}
\title{Get the ML splits per node, from a splits table}
\usage{
  LGpy_MLsplit_per_node(splits)
}
\arguments{
  \item{splits}{A data.frame containing the node numbers,
  splits, and split probabilities.}
}
\value{
  \code{MLsplits} A data.frame containing the node numbers,
  ML splits, and split probabilities.
}
\description{
  Given a table of splits probabilities from either
  \code{\link{LGpy_splits_fn_to_table}} or
  \code{\link{LGcpp_splits_fn_to_table}}, get the ML state
  for each node.
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

