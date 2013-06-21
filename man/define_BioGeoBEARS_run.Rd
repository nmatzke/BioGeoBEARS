\name{define_BioGeoBEARS_run}
\alias{define_BioGeoBEARS_run}
\title{Define a MAXIMUM LIKELIHOOD search, perhaps stratified}
\usage{
  define_BioGeoBEARS_run(abbr = "default",
    description = "defaults",
    BioGeoBEARS_model_object = define_BioGeoBEARS_model_object(),
    trfn = "Psychotria_5.2.newick",
    geogfn = "Psychotria_geog.data", timesfn = NA,
    distsfn = NA, dispersal_multipliers_fn = NA,
    area_of_areas_fn = NA, areas_allowed_fn = NA,
    max_range_size = NA, states_list = NA,
    force_sparse = NA, print_optim = TRUE,
    num_cores_to_use = NA, cluster_already_open = FALSE,
    use_optimx = TRUE)
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

  \item{max_range_size}{The maximum rangesize, in number of
  areas.  Having a smaller maximum range size means that
  you can have more areas (the size of the state space is
  greatly reduced; see
  \code{\link[cladoRcpp]{numstates_from_numareas}}.}

  \item{states_list}{A list of the possible
  states/geographic ranges, in 0-based index form.}

  \item{force_sparse}{Should sparse matrix exponentiation
  be used?}

  \item{print_optim}{If TRUE (default), print the
  optimization steps as ML estimation progresses.}

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
}
\value{
  \code{inputs} Inputs for ML search.
}
\description{
  Set up the inputs object for an ML search.
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
  \code{\link[BioGeoBEARS]{define_BioGeoBEARS_model_object}}
}

