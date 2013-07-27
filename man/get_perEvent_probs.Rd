\name{get_perEvent_probs}
\alias{get_perEvent_probs}
\title{Get the per-event probabilities at cladogenesis}
\usage{
  get_perEvent_probs(params_table, sumval = 1,
    plotwhat = "est")
}
\arguments{
  \item{params_table}{The \code{params_table} from a
  \code{BioGeoBEARS_model_object}.}

  \item{sumval}{Default=1.}

  \item{plotwhat}{Default "est", use "init" to get the
  initial starting values instead.}
}
\value{
  \code{wts} Return the per-event weights
}
\description{
  At a cladogenesis event, a large number of events are
  possible. The simplest way to compute these is just to
  assign some weight to each event, then sum all the events
  and divide by the sum to get the probabilities.  More
  complex schemes can be imagined, but these are fairly
  pointless as they would all break down once e.g.
  distance-dependence, user-specified connectivities, etc.,
  are imposed.
}
\details{
  In addition, one could imagine trying to assign total
  probabilities to each category of event, but each row of
  the cladogenesis matrix may have a different count of the
  different types of events (one row may have 1 y event and
  2 j events; another row may have 4 j, 2 v, and 2 s, and 0
  y events; etc.).

  One thing that IS meaningful is the per-event weight,
  i.e. the values that the program is using for j, v, y,
  and s.  These ARE meaningful, as long as they are forced
  to sum to some value (default 4).  This ensures that they
  are identifiable (otherwise, j,v,y,s=1 and j,v,y,s=2
  would be the same model).

  This function calculates the per-event weight as a
  proportion of some total weight, e.g. default 1.  If the
  optim result was j=0, s=1, y=1, v=1, the
  \code{get_perEventprobs()} result would be 0, 0.333,
  0.333, 0.333.
}
\note{
  Go BEARS!
}
\examples{
# default DEC+J model
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
params_table = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
params_table

get_perEvent_probs(params_table)


# DEC+J model
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 1
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 1

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

BioGeoBEARS_run_object$BioGeoBEARS_model_object =
calc_linked_params_BioGeoBEARS_model_object(
BioGeoBEARS_model_object=BioGeoBEARS_run_object$BioGeoBEARS_model_object,
update_init=TRUE)


BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
params_table = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

get_perEvent_probs(params_table)
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

