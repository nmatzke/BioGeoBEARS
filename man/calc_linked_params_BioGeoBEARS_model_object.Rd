\name{calc_linked_params_BioGeoBEARS_model_object}
\alias{calc_linked_params_BioGeoBEARS_model_object}
\title{Update parameters that are deterministic functions of free parameters}
\usage{
  calc_linked_params_BioGeoBEARS_model_object(BioGeoBEARS_model_object,
    update_init = FALSE)
}
\arguments{
  \item{BioGeoBEARS_model_object}{The BioGeoBEARS_model
  object, of class \code{BioGeoBEARS_model}}

  \item{update_init}{If \code{TRUE}, put the estimates into
  the initial values in the \code{params_table}. Default:
  \code{FALSE}.}
}
\value{
  \code{BioGeoBEARS_model_object} Updated version of the
  BioGeoBEARS_model object, of class
  \code{BioGeoBEARS_model}.
}
\description{
  This function updates the linked parameters (which are
  listed as neither "fixed" nor "free" in
  \code{params_table$type}; i.e., they are equations which
  are calculated from #' the fixed and free parameters,
  which should have already been set by other functions).
}
\details{
  \code{params_table$type} is typically stored in:
  \code{BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table}.
}
\note{
  Go BEARS!
}
\examples{
# Define a BioGeoBEARS run object
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

# Set 'j' to be free, i.e. as in a DEC+J model (adding jump dispersal
# to the LAGRANGE DEC model)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.25
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.25

# Display result
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

# Update the other parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object =
calc_linked_params_BioGeoBEARS_model_object(
BioGeoBEARS_model_object=BioGeoBEARS_run_object$BioGeoBEARS_model_object)

# Display result
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link[BioGeoBEARS]{define_BioGeoBEARS_model_object}}
  \code{\link[BioGeoBEARS]{define_BioGeoBEARS_run}}
}

