\name{plot_cladogenesis_size_probabilities}
\alias{plot_cladogenesis_size_probabilities}
\title{Graphical display of P(daughter rangesize) for your input or inferred speciation model}
\usage{
  plot_cladogenesis_size_probabilities(BioGeoBEARS_run_object,
    plotwhat = "est", statenames = NULL)
}
\arguments{
  \item{BioGeoBEARS_run_object}{The input run object.}

  \item{plotwhat}{Default is "input", which means plotting
  the starting model.}

  \item{statenames}{State names to pass to
  \code{\link{plot_cladogenesis_size_probabilities}}. If
  \code{NULL} (default), these are auto-generated assuming
  all areas up to the maximum number are allowed.}
}
\value{
  \code{Nothing}
}
\description{
  This function produces a graphical summary of the
  daughter rangesize aspect of the cladogenesis model
  stored in a \code{BioGeoBEARS_run_object}. This could be
  either an input model, or the result of the ML parameter
  search.
}
\details{
  The \code{LAGRANGE} DEC model assumes that at
  cladogenesis events, one daughter species has a range
  size of 1 area, and the other daughter either inherits
  the full ancestral range (sympatric-subset speciation),
  inherits the remainder of the ancestral range
  (vicariance), or as the same range (sympatric-range
  copying, which is the only option when the ancestor range
  is of size 1 area.

  BioGeoBEARS enables numerous additional models. To see
  how these are similar or different from the LAGRANGE DEC
  cladogenesis model, this function can be used.  E.g.,
  comparison of \code{LAGRANGE} DEC to a \code{DIVA}-like
  model is instructive: see examples. DIVA disallows
  sympatric-subset speciation (probability 0 under this
  model), but allows classic vicariance (a species with 4
  areas splitting into 2 daughters, each occupying 2
  areas).  LAGRANGE DEC gives 0 probability to a
  \code{4->(2,2)} history, allowing only \code{4->(3,1)} or
  \code{4->(1,3)} histories.

  Several additional plots relating to the cladogenesis
  model are also produced.  Best used via
  \code{\link{plot_BioGeoBEARS_model}}.
}
\note{
  Go BEARS!
}
\examples{
blah=1
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link{plot_BioGeoBEARS_model}},
  \code{\link{define_BioGeoBEARS_run}},
  \code{\link{define_BioGeoBEARS_model_object}}
}

