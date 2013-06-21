\name{BioGeoBEARS_model}
\alias{BioGeoBEARS_model}
\title{An object of class BioGeoBEARS_model holding the model inputs}
\description{
  An object of class BioGeoBEARS_model holding the model
  inputs
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

