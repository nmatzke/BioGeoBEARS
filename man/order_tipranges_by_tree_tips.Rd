\name{order_tipranges_by_tree_tips}
\alias{order_tipranges_by_tree_tips}
\title{Reorder the rows in a \code{tipranges} object, to correspond to tree tips}
\usage{
  order_tipranges_by_tree_tips(tipranges, tr)
}
\arguments{
  \item{tipranges}{An object of class
  \code{\link{tipranges}}.}

  \item{tr}{A \code{\link[ape]{phylo}} tree object.}
}
\value{
  \code{tipranges} An object of class \code{tipranges}
}
\description{
  The tipranges object, as read from a
  LAGRANGE/PHYLIP-style geography file, may not have the
  species names as the same order as they are in the tips
  of the tree. This function allows the user to reorder
  them to match the tree
}
\note{
  Go BEARS!
}
\examples{
testval=1

# Get the example files directory
extdata_dir = system.file("extdata/", package="BioGeoBEARS")
# tmp hard code:
# extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"
# Set the filename (Hawaiian Psychotria from Ree & Smith 2008)

trfn = paste(extdata_dir, "Psychotria_5.2.newick", sep="")
tr = read.tree(trfn)

fn = paste(extdata_dir, "Psychotria_geog.data", sep="")
tipranges1 = getranges_from_LagrangePHYLIP(lgdata_fn=fn)
tipranges1

# Reorder the tipranges object
tipranges2 = order_tipranges_by_tree_tips(tipranges1, tr)
tipranges2
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
  \url{https://code.google.com/p/lagrange/}

  Matzke_2012_IBS

  ReeSmith2008

  SmithRee2010_CPPversion
}
\seealso{
  \code{\link{tipranges_to_area_strings}},
  \code{\link{define_tipranges_object}},
  \code{\link{save_tipranges_to_LagrangePHYLIP}}
}

