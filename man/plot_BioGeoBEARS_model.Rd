\name{plot_BioGeoBEARS_model}
\alias{plot_BioGeoBEARS_model}
\title{Graphical display of your anagenetic and cladogenetic biogeography models}
\usage{
  plot_BioGeoBEARS_model(obj, obj_is_run_or_results = NULL,
    plotwhat = "init", titletxt = "", statenames = NULL)
}
\arguments{
  \item{obj}{The input object, either a
  \code{BioGeoBEARS_run_object} (if so, set
  \code{obj_is_run_or_results="run"} or an output object
  from \code{\link{bears_optim_run}} (if so, specify
  \code{obj_is_run_or_results="results"}.}

  \item{obj_is_run_or_results}{Specify \code{"run"} or
  \code{"results"}, as described above for parameter
  \code{obj}.}

  \item{plotwhat}{Default is \code{"init"}, which means
  plotting the starting model parameters. \code{"est"}
  plots the estimated model parameters.}

  \item{titletxt}{Additional text for the title of the
  plot}

  \item{statenames}{State names to pass to
  \code{\link{plot_cladogenesis_size_probabilities}}. If
  \code{NULL} (default), these are auto-generated assuming
  all areas up to the maximum number are allowed.}
}
\value{
  nada
}
\description{
  This function produces a graphical summary of the model
  stored in a \code{BioGeoBEARS_run_object}. This could be
  either an input model, or the result of the ML parameter
  search.
}
\details{
  Understanding of phylogenetic methoods in historical
  biogeography methods is hampered by the difficulty of
  displaying the models the computer is using.  This
  function is one attempt to improve the situation, by
  plotting the relative weights of the various parameters.
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
  \code{\link{plot_cladogenesis_size_probabilities}},
  \code{\link{define_BioGeoBEARS_run}},
  \code{\link{define_BioGeoBEARS_model_object}}
}

