\name{getranges_from_LagrangePHYLIP}
\alias{getranges_from_LagrangePHYLIP}
\title{Read a LAGRANGE PHYLIP-style file containing geographic ranges into a \code{tipranges} object}
\usage{
  getranges_from_LagrangePHYLIP(lgdata_fn = "lagrange_area_data_file.data")
}
\arguments{
  \item{lgdata_fn}{The LAGRANGE geographic data file to be
  read.}
}
\value{
  \code{tipranges_object} An object of class
  \code{tipranges}
}
\description{
  Given some geographic range data for tips in the Lagrange
  C++/PHYLIP format (\cite{Smith et al. (2010)}), this
  function imports the range data into a
  \code{tipranges}-class \code{data.frame} structure.
}
\details{
  LAGRANGE C++ geographic range files are ASCII text files
  with the format:

  \code{19 4 (A B C D)}\cr \code{P_mariniana_Kokee2
  1000}\cr \code{P_mariniana_Oahu 0100}\cr
  \code{P_mariniana_MauiNui 0010}\cr
  \code{P_hawaiiensis_Makaopuhi 0001}\cr
  \code{P_wawraeDL7428 1000}\cr [...]\cr \cr

  The first row specifies the number of taxa (here, 19),
  the number of areas (here, 4), and finally, the
  names/abbreviations of the areas. The rest of the rows
  give the taxon names, followed by a tab and then the
  presence/absence in each range with 1s/0s.

  The file above is part of the geographic range data for
  the Hawaiian \emph{Psychotria} dataset used by \cite{Ree
  et al. (2008)}.
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
fn = paste(extdata_dir, "Psychotria_geog.data", sep="")
getranges_from_LagrangePHYLIP(lgdata_fn=fn)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
  \url{https://code.google.com/p/lagrange/}

  SmithRee2010_CPPversion

  ReeSmith2008

  Matzke_2012_IBS
}
\seealso{
  \code{\link{define_tipranges_object}},
  \code{\link{save_tipranges_to_LagrangePHYLIP}}
}

