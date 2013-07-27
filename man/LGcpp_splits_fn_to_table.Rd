\name{LGcpp_splits_fn_to_table}
\alias{LGcpp_splits_fn_to_table}
\title{Get the ML splits per node, from C++ LAGRANGE output}
\usage{
  LGcpp_splits_fn_to_table(splits_fn)
}
\arguments{
  \item{splits_fn}{The filename of a C++ \code{LAGRANGE}
  output file.}
}
\value{
  \code{splits} A data.frame containing the node numbers,
  splits, and split probabilities.
}
\description{
  C++ \code{LAGRANGE} outputs a list of splits and split
  probabilities for each node. This function converts them
  to a table.
}
\details{
  \code{LAGRANGE} outputs just the splits making up the top
  95% of the probability, or 15 states, whichever comes
  first.

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

# splits_fn = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/
# examples/Psychotria_M0/LAGRANGE_C++/Psychotria_M0_lgcpp_out_splits00001.txt"
# LGcpp_splits_fn_to_table(splits_fn)
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

