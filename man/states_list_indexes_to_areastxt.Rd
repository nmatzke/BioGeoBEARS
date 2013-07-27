\name{states_list_indexes_to_areastxt}
\alias{states_list_indexes_to_areastxt}
\title{States (ranges) lists to txt string of the areas}
\usage{
  states_list_indexes_to_areastxt(states_list, areanames,
    counting_base = 0, concat = TRUE, sep = "")
}
\arguments{
  \item{states_list}{A list of states, where each state
  consists of a list of areas.}

  \item{areanames}{A list of areanamess.}

  \item{counting_base}{Does states_list start indexing
  areas from 0 (default) or 1?}

  \item{concat}{If TRUE (default), merge the areas in a
  state into a single string.}

  \item{sep}{Character to merge on, as in
  \code{\link[base]{paste}}. Default "".}
}
\value{
  \code{tiparea} A string.
}
\description{
  This is a utility function.
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

  Matzke_2012_IBS
}
\seealso{
  \code{\link{getname}},
  \code{\link{order_tipranges_by_tree_tips}},
  \code{\link{define_tipranges_object}},
  \code{\link{save_tipranges_to_LagrangePHYLIP}}
}

