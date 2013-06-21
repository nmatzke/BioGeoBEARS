\name{define_tipranges_object}
\alias{define_tipranges_object}
\title{Define a tipranges class and object}
\usage{
  define_tipranges_object(tmpdf = NULL)
}
\arguments{
  \item{tmpdf}{The user may input a \code{data.frame}
  holding the range data, if they like. Default is
  \code{NULL}, which means the function will produce a
  temporary \code{data.frame} as an example.}
}
\value{
  \code{tipranges_object} The tipranges object, of class
  \code{tipranges}
}
\description{
  Class \code{tipranges} is an extension of the
  \code{\link{data.frame}} class.  It is used for holding
  discrete geographic range data for the tips on a
  phylogeny. Geographic ranges are represented with bit
  encoding (0/1) indicating absence or presence in each
  possible area.
}
\details{
  This is just a data.frame with: rows = taxanames\cr
  columns = area names\cr cells = 0/1 representing
  empty/occupied\cr
}
\note{
  Go BEARS!
}
\examples{
testval=1
tipranges_object = define_tipranges_object()
tipranges_object
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link{getareas_from_tipranges_object}},
  \code{\link[cladoRcpp]{areas_list_to_states_list_old}},
  \code{\link{areas_list_to_states_list_new}},
  \code{\link{tipranges_to_tip_condlikes_of_data_on_each_state}}
}

