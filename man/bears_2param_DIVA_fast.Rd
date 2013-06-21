\name{bears_2param_DIVA_fast}
\alias{bears_2param_DIVA_fast}
\title{2-parameter model, fixed cladogenesis model (as in LAGRANGE)}
\usage{
  bears_2param_DIVA_fast(trfn = "Psychotria_5.2.newick",
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
  This function implements a biogeographical model with 2
  free parameters (\emph{d}, rate of dispersal/range
  addition, and \emph{e}, rate of extinction/range
  contraction), and a fixed cladogenesis model copying the
  DIVA model (\cite{Ronquist (1997)}.  This has: equal
  probability of vicariance at all range sizes, but NO
  sympatric-subset speciation, no jump/founder-event
  speciation, and sympatric-range-copying events are
  limited to the smaller descendant always having a range
  size of 1 area(\cite{Ronquist et al. (2011)}).
}
\details{
  Once the model is set up, it is input into the
  optimization routine \code{\link[optimx]{optimx}} (the
  more common \code{\link[stats]{optim}} can also be used
  by editing the function), and
  \code{\link{calc_loglike_sp}} is used to calculate the
  log-likelihood of each set of parameters.  Once the
  parameter values that give the data the maximum
  likelihood are found, they are reported back to the
  function and returned to the user.

  This duplicates the model used in the standard DIVA
  implementation (\cite{Ree et al. (2008)}, \cite{Ree
  (2009)}, \cite{Smith et al. (2010)}, with no constraints
  on dispersal or range size.

  Here, all of the fastest processing options have been
  used.

  Model implementations are provided to show the user how a
  specific model can be set up and optimized.  This is
  preferable compared to the "black-box" nature of most
  other inference packages.  Users are encouraged to
  experiment.  Useful models can be added to later versions
  of \code{BioGeoBEARS}.
}
\note{
  Go BEARS!
}
\examples{
test=1

# Get the example files directory
extdata_dir = system.file("extdata/", package="BioGeoBEARS")
# tmp hard code:
# extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"

# Set the filenames (Hawaiian Psychotria from Ree & Smith 2008)
trfn = paste(extdata_dir, "Psychotria_5.2.newick", sep="")
tr = read.tree(file=trfn)

geogfn = paste(extdata_dir, "Psychotria_geog.data", sep="")

# Look at the tree and ranges, for kicks
getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tr

\dontrun{
# Run the ML search
bears_output = bears_2param_standard_fast(trfn=trfn, geogfn=geogfn)
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

  Ronquist_1997_DIVA

  Ronquist_Sanmartin_2011

  ReeSmith2008

  Ree2009configurator

  SmithRee2010_CPPversion

  Matzke_2012_IBS

  Ronquist1996_DIVA
}
\seealso{
  \code{\link[cladoRcpp]{numstates_from_numareas}},
  \code{\link{getranges_from_LagrangePHYLIP}},
  \code{\link[ape]{read.tree}},
  \code{\link{calc_loglike_sp}}
}

