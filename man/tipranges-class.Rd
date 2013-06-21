\name{tipranges}
\alias{tipranges}
\title{The tipranges class}
\description{
  This class holds geographic range data for each tip in a
  phylogeny.
}
\details{
  Geographic range data can be read into a tipranges class
  object with BioGeoBEARS functions, e.g.
  \code{define_tipranges_object} or
  \code{getareas_from_tipranges_object}.

  Class \code{tipranges} is an extension of the
  \code{\link{data.frame}} class.  It is used for holding
  discrete geographic range data for the tips on a
  phylogeny. Geographic ranges are represented with bit
  encoding (0/1) indicating absence or presence in each
  possible area.

  This is just a data.frame with: rows = taxanames\cr
  columns = area names\cr cells = 0/1 representing
  empty/occupied\cr
}
\note{
  Go BEARS!
}
\section{Slots}{
  \describe{ \item{\code{df}:}{Data.frame of class
  \code{"numeric"}, containing data from df} }
}
\examples{
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
  \code{\link{define_tipranges_object}},
  \code{\link{getareas_from_tipranges_object}},
  \code{\link[cladoRcpp]{areas_list_to_states_list_old}},
  \code{\link{areas_list_to_states_list_new}},
  \code{\link{tipranges_to_tip_condlikes_of_data_on_each_state}}
}

