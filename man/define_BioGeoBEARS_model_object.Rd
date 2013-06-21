\name{define_BioGeoBEARS_model_object}
\alias{define_BioGeoBEARS_model_object}
\title{Define a BioGeoBEARS_model class and object}
\usage{
  define_BioGeoBEARS_model_object(minval_anagenesis = 1e-15,
    minval_cladogenesis = 1e-05, maxval = 5)
}
\arguments{
  \item{minval_anagenesis}{Minimum value above zero for d,
  e, a, b parameters.}

  \item{minval_cladogenesis}{Minimum value above zero for
  j, v, etc.}

  \item{maxval}{Maximum value for d, e, a}
}
\value{
  \code{BioGeoBEARS_model_object} The BioGeoBEARS_model
  object, of class \code{BioGeoBEARS_model}
}
\description{
  Class \code{BioGeoBEARS_model} is an extension of the
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
BioGeoBEARS_model_object = define_BioGeoBEARS_model_object()
BioGeoBEARS_model_object
define_BioGeoBEARS_model_object()
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link[cladoRcpp]{areas_list_to_states_list_old}},
  \code{\link{areas_list_to_states_list_new}}
}

