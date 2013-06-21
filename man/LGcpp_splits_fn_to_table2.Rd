\name{LGcpp_splits_fn_to_table2}
\alias{LGcpp_splits_fn_to_table2}
\title{Get the ML splits per node, from Python LAGRANGE output}
\usage{
  LGcpp_splits_fn_to_table2(splits_fn)
}
\arguments{
  \item{splits_fn}{The filename of a Python LAGRANGE output
  file.}
}
\value{
  \code{splits} A data.frame containing the node numbers,
  splits, and split probabilities.
}
\description{
  Python LAGRANGE outputs a list of splits and split
  probabilities for each node. This function converts them
  to a table.
}
\details{
  LAGRANGE outputs just the splits making up the top 95% of
  the probability, or 15 states, whichever comes first.

  See \code{\link{LGpy_MLsplit_per_node}} for choosing the
  single ML split at each node, and see
  \code{\link{get_lagrange_nodenums}} for connecting these
  node numbers to APE node numbers.
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
  \code{\link{LGpy_MLsplit_per_node}}
}

