\name{tipranges_to_area_strings}
\alias{tipranges_to_area_strings}
\title{Convert tipranges binary coding to range strings}
\usage{
  tipranges_to_area_strings(tipranges, areaabbr = NULL)
}
\arguments{
  \item{tipranges}{An object of class
  \code{\link{tipranges}}.}

  \item{areaabbr}{A vector of the abbreviations (preferably
  1 character each).}
}
\value{
  \code{tiprange_names} A vector of strings.
}
\description{
  This function converts the 0110-type format of the
  tipranges object into a list of strings describing the
  geographic ranges.  E.g., 1100 becomes AB, 0111 become
  BCD (assuming the regions are abbreviated A, B, C...).
  Users can input their preferred abbreviations with
  \code{areaabbr}.
}
\details{
  Note that you will HAVE to use
  \code{\link{order_tipranges_by_tree_tips}} on the
  tipranges object first, to make sure the tipranges are in
  the correct order on the tree tips.
}
\note{
  Go BEARS!
}
\examples{
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
tipranges_to_area_strings(tipranges=tipranges1, areaabbr=NULL)
tipranges_to_area_strings(tipranges=tipranges1, areaabbr=c("K", "O", "M", "H"))

# Reorder the tipranges object
tipranges2 = order_tipranges_by_tree_tips(tipranges1, tr)
tipranges2
tipranges_to_area_strings(tipranges=tipranges2, areaabbr=NULL)
tipranges_to_area_strings(tipranges=tipranges2, areaabbr=c("K", "O", "M", "H"))
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
  \code{\link{order_tipranges_by_tree_tips}},
  \code{\link{define_tipranges_object}},
  \code{\link{save_tipranges_to_LagrangePHYLIP}}
}

