\name{bears_2param_standard_fast_symOnly}
\alias{bears_2param_standard_fast_symOnly}
\title{2-parameter model, no cladogenesis model (as in BayArea or other purely continuous-time model)}
\usage{
  bears_2param_standard_fast_symOnly(trfn = "Psychotria_5.2.newick",
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
  This implements a 2-parameter model, as in LAGRANGE or
  \code{\link{bears_2param_standard_fast}}, but omits the
  speciation/cladogenesis model.  This means that the model
  is purely continuous-time, as when biogeographic range is
  treated as a discrete character in software designed for
  inference on morphological () or molecular data ().  This
  model is that implemented in \code{BayArea}, if no
  distance-dependent effect on dispersal probability is
  assumed.  Such distance-dependence could easily be added
  with a third parameter, however.
}
\details{
  \code{BayArea} is a new program by Landis, Matzke, Moore,
  and Huelsenbeck; see \cite{Landis et al. (2013)}.
  However, BayArea does not currently implement
  cladogenesis models; it only has continuous-time model
  for evolutionary change along branches.  In effect, this
  means that the cladogenesis model is sympatric speciation
  with complete range copying with probability 1.
}
\note{
  Go BEARS!
}
\examples{
test=1

# Get the example files directory
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
# tmp hard code:
# extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"

# Set the filenames (Hawaiian Psychotria from Ree & Smith 2008)
trfn = np(paste(extdata_dir, "/Psychotria_5.2.newick", sep=""))
tr = read.tree(file=trfn)

geogfn = np(paste(extdata_dir, "/Psychotria_geog.data", sep=""))

# Look at the tree and ranges, for kicks
getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tr

\dontrun{
# Run the ML search
bears_output = bears_2param_standard_fast_symOnly(trfn=trfn, geogfn=geogfn)
bears_output
}
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  Felsenstein, Joe.  The Newick tree format.
  \url{http://evolution.genetics.washington.edu/phylip/newicktree.html}
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
  \url{https://code.google.com/p/lagrange/}

  Landis_Matzke_etal_2013_BayArea

  Matzke_2012_IBS

  ReeSmith2008

  Ree2009configurator

  SmithRee2010_CPPversion
}
\seealso{
  \code{\link[cladoRcpp]{numstates_from_numareas}},
  \code{\link{getranges_from_LagrangePHYLIP}},
  \code{\link[ape]{read.tree}},
  \code{\link{calc_loglike_sp}}
}

