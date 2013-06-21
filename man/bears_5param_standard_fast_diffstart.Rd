\name{bears_5param_standard_fast_diffstart}
\alias{bears_5param_standard_fast_diffstart}
\title{5-parameter model, with different starting points for optimization}
\usage{
  bears_5param_standard_fast_diffstart(trfn = "Psychotria_5.2.newick",
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
  This implements the same model as
  \code{\link{bears_5param_standard_fast}}, but uses
  different starting points and slightly different
  constraints.
}
\details{
  As the number of parameters increases, the importance of
  starting ML optimization runs from different places
  increases.  Several starting points should be tried,
  especially if the likelihood surface seems flat.
}
\note{
  Go BEARS!
}
\examples{
test=1

# Get the example files directory
extdata_dir = system.file("extdata/", package="BioGeoBEARS")
# tmp hard code: extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"

# Set the filenames (Hawaiian Psychotria from Ree & Smith 2008)
trfn = paste(extdata_dir, "Psychotria_5.2.newick", sep="")
tr = read.tree(file=trfn)

geogfn = paste(extdata_dir, "Psychotria_geog.data", sep="")

# Look at the tree and ranges, for kicks
getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tr

\dontrun{
# Run the ML search
bears_output = bears_5param_standard_fast_diffstart(trfn=trfn, geogfn=geogfn)
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

