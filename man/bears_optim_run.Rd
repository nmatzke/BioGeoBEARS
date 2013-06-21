\name{bears_optim_run}
\alias{bears_optim_run}
\title{Run ML search from \code{BioGeoBEARS_run} object}
\usage{
  bears_optim_run(BioGeoBEARS_run_object = define_BioGeoBEARS_run())
}
\arguments{
  \item{BioGeoBEARS_run_object}{Contains all inputs}
}
\value{
  \code{bears_output} A list of outputs.
  bears_output$optim_result
}
\description{
  Uses a BioGeoBEARS_run_object to simplify input.
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
bears_output = bears_optim_run(trfn=trfn, geogfn=geogfn)
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

