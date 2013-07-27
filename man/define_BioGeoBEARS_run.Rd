\name{define_BioGeoBEARS_run}
\alias{define_BioGeoBEARS_run}
\title{Define a maximum likelihood search, perhaps stratified}
\usage{
  define_BioGeoBEARS_run(abbr = "default",
    description = "defaults",
    BioGeoBEARS_model_object = define_BioGeoBEARS_model_object(),
    trfn = "Psychotria_5.2.newick",
    geogfn = "Psychotria_geog.data", timesfn = NA,
    distsfn = NA, dispersal_multipliers_fn = NA,
    area_of_areas_fn = NA, areas_allowed_fn = NA,
    detects_fn = NA, controls_fn = NA, max_range_size = NA,
    states_list = NULL, force_sparse = FALSE,
    use_detection_model = FALSE, print_optim = TRUE,
    num_cores_to_use = NA, cluster_already_open = FALSE,
    use_optimx = TRUE, return_condlikes_table = FALSE,
    calc_TTL_loglike_from_condlikes_table = TRUE,
    calc_ancprobs = TRUE, fixnode = NULL, fixlikes = NULL,
    speedup = TRUE, tmpwd = getwd())
}
\arguments{
  \item{abbr}{Text abbreviation of run, e.g. "default"}

  \item{description}{Text description of run, e.g.
  "defaults"}

  \item{BioGeoBEARS_model_object}{Default is
  \code{define_BioGeoBEARS_model_object()}}

  \item{trfn}{The filename of the phylogenetic tree, in
  NEWICK format
  (\url{http://evolution.genetics.washington.edu/phylip/newicktree.html}).
  Tipnames should match the names in geogfn.  See
  \code{\link[ape]{read.tree}} in APE for reading in
  phylogenetic trees. Default "Psychotria_5.2.newick"}

  \item{geogfn}{A PHYLIP-style file with geographic range
  data (see \code{\link{getranges_from_LagrangePHYLIP}})
  for each tipname. This is the same format used by C++
  LAGRANGE (\cite{SmithRee2010_CPPversion}). Default
  "Psychotria_geog.data"}

  \item{timesfn}{Filename for the stratified times.}

  \item{distsfn}{Filename for the changing distances.}

  \item{dispersal_multipliers_fn}{Filename for the changing
  hard-coded dispersal multipliers}

  \item{area_of_areas_fn}{Filename for the area of each
  area}

  \item{areas_allowed_fn}{Filename for the allowed
  connections between areas for single-species ranges.}

  \item{detects_fn}{Filename for the counts of detections
  of OTUs of interest. See \code{\link{calc_obs_like}}.}

  \item{controls_fn}{Filename for the counts of taphonomic
  controls (which INCLUDE the OTUs of interest). See
  \code{\link{calc_obs_like}}.}

  \item{max_range_size}{The maximum rangesize, in number of
  areas.  Having a smaller maximum range size means that
  you can have more areas (the size of the state space is
  greatly reduced; see
  \code{\link[cladoRcpp]{numstates_from_numareas}}.}

  \item{states_list}{A list of the possible
  states/geographic ranges, in 0-based index form.}

  \item{force_sparse}{Should sparse matrix exponentiation
  be used?  Default \code{FALSE}, which means dense matrix
  exponentiation is always used. If \code{NA}, the program
  will use sparse matrix exponentiation for transition
  matrices above rank 128 (size 128x128).  NOTE: Sparse
  matrix exponentiation seems to give correlated, but not
  exact, results, and these errors may accumulate.
  Presumably the problems become less with larger matrices,
  but I have not explored this in detail.}

  \item{use_detection_model}{If TRUE, use the detection
  model (with parameters mf, dp, and fdp) and counts of
  detections and counts of taphonomic controls to calculate
  the \code{tip_condlikes_of_data_on_each_state}.}

  \item{print_optim}{If TRUE (default), print the
  optimization steps as ML estimation progresses.}

  \item{tmpwd}{The working directory in which the input and
  output files will be placed. Default is
  \code{\link[base]{getwd}}. This is stored mostly for
  future reference; users are responsible for manually
  navigating to the appropriate directory ahead of time,
  using \code{\link[base]{setwd}}.}

  \item{num_cores_to_use}{If >1, parallel processing will
  be attempted. \bold{Note:} parallel processing via
  \code{library (parallel)} will work in Mac command-line
  R, but not in Mac GUI \code{R.app}.}

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

  \item{use_optimx}{If TRUE, use \code{optimx} rather that
  \code{optim}.}

  \item{return_condlikes_table}{If \code{TRUE}, return the
  table of ALL conditional likelihood results, including at
  branch subsections (only some should be used in
  calculating the final log-likelihood of the geography
  range data on the tree!)}

  \item{calc_TTL_loglike_from_condlikes_table}{If TRUE,
  force making of the condlikes table, and use it to
  calculate the log-likelihood (default=TRUE; matches
  LAGRANGE).}

  \item{calc_ancprobs}{If \code{TRUE} (default), calculate
  and return the necessary pieces (uppass and downpass
  probs) for ancestral states.}

  \item{fixnode}{If the state at a particular node is going
  to be fixed (e.g. for ML marginal ancestral states), give
  the node number.}

  \item{fixlikes}{The state likelihoods to be used at the
  fixed node.  I.e. 1 for the fixed state, and 0 for the
  others.}

  \item{speedup}{If \code{TRUE} (default), set the maximum
  number of iterations to \code{itnmax=50*(number of free
  parameters)}, instead of the \code{optimx} default, 250.
  Also set \code{optimx} \code{reltol} parameter to 0.001
  (instead of the default, ~1e-8).}
}
\value{
  \code{inputs} Inputs for ML search.
}
\description{
  Set up the inputs object for an ML search.  See parameter
  descriptions for defaults.
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
  \code{\link{readfiles_BioGeoBEARS_run}},
  \code{\link[BioGeoBEARS]{define_BioGeoBEARS_model_object}},
  \code{\link[base]{setwd}}, \code{\link[base]{getwd}}
}

