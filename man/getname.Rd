\name{getname}
\alias{getname}
\title{Collapse range abbreviations to strings}
\usage{
  getname(TFrow, tiparea_names, concat = TRUE, sep = "")
}
\arguments{
  \item{TFrow}{A list of TRUE and FALSE}

  \item{tiparea_names}{The names of each area}

  \item{concat}{If TRUE (default), merge the areas in a
  state into a single string.}

  \item{sep}{The sep argument for
  \code{\link[base]{paste}}.}
}
\value{
  \code{tiparea} A string.
}
\description{
  This is a utility function used by apply in
  \code{\link{tipranges_to_area_strings}}. It extracts the
  present areas and concatenates the abbreviations for one
  row.
}
\note{
  Go BEARS!
}
\examples{
getname(TFrow=c(FALSE, TRUE, TRUE, FALSE),
tiparea_names=c("K", "O", "M", "H"), sep="")
getname(TFrow=c(FALSE, TRUE, TRUE, FALSE),
tiparea_names=c("K", "O", "M", "H"), sep="_")
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link{states_list_indexes_to_areastxt}},
  \code{\link{order_tipranges_by_tree_tips}},
  \code{\link{define_tipranges_object}},
  \code{\link{save_tipranges_to_LagrangePHYLIP}}
}

