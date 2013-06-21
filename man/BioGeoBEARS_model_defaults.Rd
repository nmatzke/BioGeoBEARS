\name{BioGeoBEARS_model_defaults}
\alias{BioGeoBEARS_model_defaults}
\title{Set up a default BioGeoBEARS model object}
\usage{
  BioGeoBEARS_model_defaults(minval_anagenesis = 1e-15,
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
  \code{param_table} Return the parameter table object
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

  Matzke_2012_IBS
}
\seealso{
  \code{\link[base]{rbind}}
}

