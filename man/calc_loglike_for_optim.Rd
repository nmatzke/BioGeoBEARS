\name{calc_loglike_for_optim}
\alias{calc_loglike_for_optim}
\title{Take model parameters and the data and calculate the log-likelihood}
\usage{
  calc_loglike_for_optim(params, BioGeoBEARS_model_object,
    phy, tip_condlikes_of_data_on_each_state,
    print_optim = TRUE, areas_list = areas_list,
    states_list = states_list, force_sparse = force_sparse,
    cluster_already_open = cluster_already_open)
}
\arguments{
  \item{params}{A vector of parameters for optimization.}

  \item{BioGeoBEARS_model_object}{Object containing the
  model.}

  \item{phy}{An ape tree object}

  \item{tip_condlikes_of_data_on_each_state}{A numeric
  matrix with rows representing tips, and columns
  representing states/geographic ranges.  The cells give
  the likelihood of the observation data under the
  assumption that the tip has that state; typically this
  means that the known geographic range gets a '1' and all
  other states get a 0.}

  \item{force_sparse}{Should sparse matrix exponentiation
  be used?}

  \item{print_optim}{If TRUE (default), print the
  optimization steps as ML estimation progresses.}

  \item{areas_list}{A list of the desired area
  names/abbreviations/letters (?).}

  \item{states_list}{A list of the possible
  states/geographic ranges, in 0-based index form.}

  \item{cluster_already_open}{If the user wants to
  distribute the matrix exponentiation calculations from
  all the branches across a number of processors/nodes on a
  cluster, specify the cluster here.  E.g.
  \code{cluster_already_open =
  makeCluster(rep("localhost",num_cores_to_use), type =
  "SOCK")}.  Note: this will work on most platforms,
  including Macs running R from command line, but will NOT
  work on Macs running the R GUI \code{R.app}, because
  parallel processing functions like \code{MakeCluster}
  from e.g. \code{library(parallel)} for some reason crash
  R.app.  The program runs a check for R.app and will just
  run on 1 node if found.}
}
\value{
  \code{ttl_loglike} The log-likelihood of the data under
  the input model and parameters.
}
\description{
  This function is an input to optim or optimx, the ML
  estimation routines.
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
  \code{\link[stats]{convolve}}
}

