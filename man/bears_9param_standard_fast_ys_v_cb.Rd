\name{bears_9param_standard_fast_ys_v_cb}
\alias{bears_9param_standard_fast_ys_v_cb}
\title{6-parameter model, adding j (founder-event speciation), v (vicariance proportion), and both maxent_constraint_01 and maxent_constraint_01v}
\usage{
  bears_9param_standard_fast_ys_v_cb(trfn = "Psychotria_5.2.newick",
    geogfn = "Psychotria_geog.data", max_range_size = NULL,
    num_cores_to_use = NULL)
}
\arguments{
  \item{trfn}{The filename of the phylogenetic tree, in
  NEWICK format
  (\url{http://evolution.genetics.washington.edu/phylip/newicktree.html}).
  Tipnames should match the names in geogfn.  See
  \code{\link[ape]{read.tree}} in APE for reading in
  phylogenetic trees.}

  \item{geogfn}{A PHYLIP-style file with geographic range
  data (see \code{\link{getranges_from_LagrangePHYLIP}})
  for each tipname. This is the same format used by C++
  LAGRANGE (\cite{SmithRee2010_CPPversion}).}

  \item{max_range_size}{The maximum rangesize, in number of
  areas.  Having a smaller maximum range size means that
  you can have more areas (the size of the state space is
  greatly reduced; see
  \code{\link[cladoRcpp]{numstates_from_numareas}}.}

  \item{num_cores_to_use}{If >1, parallel processing will
  be attempted. \bold{Note:} parallel processing via
  \code{library (parallel)} will work in Mac command-line
  R, but not in Mac GUI \code{R.app}.}
}
\value{
  \code{bears_output} A list of outputs.
  bears_output$optim_result
}
\description{
  This implements a 6-parameter model, basically
  \code{LAGRANGE} or
  \code{\link{bears_2param_standard_fast}}, but with a
  parameter \emph{j} controlling the relative weight of
  "founder-event speciation" (\cite{Matzke (2012)}), and
  another parameter \emph{v} controlling the relative
  weight of vicariance.  The remainder of the weight
  (weights must sum to 1) is taken up by the range-copying
  and range-subset forms of sympatric speciation utilized
  in \code{LAGRANGE} and
  \code{\link{bears_2param_standard_fast}}.  A fifth
  parameter, \emph{maxent_constraint_01}, controls the
  relative probability of daughter lineages of different
  rangesizes.  A sixth parameter,
  \emph{maxent_constraint_01v}, controls the relative
  probability of daughter lineages of different rangesizes,
  but only for the vicariance events, which are rather
  different from other types of speciation events. If
  maxent_constraint_01v=0.0001, the smaller-ranged daughter
  lineage will have size 1 area, with probability 1.  If
  maxent_constraint_01v=0.5, all different rangesizes will
  have equal probability, and if
  maxent_constraint_01v=0.9999, the largest possible range
  will have probability 1 -- but note that in a vicariance
  context, this would mean at maximum rangesize of 50% of
  the areas.
}
\note{
  Go BEARS!
}
\examples{
test=1

# Get the example files directory
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
# tmp hard code: extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"

# Set the filenames (Hawaiian Psychotria from Ree & Smith 2008)
trfn = np(paste(extdata_dir, "/Psychotria_5.2.newick", sep=""))
tr = read.tree(file=trfn)

geogfn = np(paste(extdata_dir, "/Psychotria_geog.data", sep=""))

# Look at the tree and ranges, for kicks
getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tr

\dontrun{
# Run the ML search
bears_output = bears_6param_standard_fast_ys_v(trfn=trfn, geogfn=geogfn)
bears_output
}
# Add c (direct changes), b (branch length exponent), x (distance exponent)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  Felsenstein, Joe.  The Newick tree format.
  \url{http://evolution.genetics.washington.edu/phylip/newicktree.html}
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
  \url{https://code.google.com/p/lagrange/}

  Matzke_2012_IBS

  ReeSmith2008

  Ree2009configurator

  SmithRee2010_CPPversion

  Landis_Matzke_etal_2013_BayArea
}
\seealso{
  \code{\link{bears_2param_standard_fast}},
  \code{\link[cladoRcpp]{numstates_from_numareas}},
  \code{\link{getranges_from_LagrangePHYLIP}},
  \code{\link[ape]{read.tree}},
  \code{\link{calc_loglike_sp}}
}

