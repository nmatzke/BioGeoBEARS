\name{save_tipranges_to_LagrangePHYLIP}
\alias{save_tipranges_to_LagrangePHYLIP}
\title{Save a tipranges object to a LAGRANGE PHYLIP-style file containing binary-encoded geographic ranges}
\usage{
  save_tipranges_to_LagrangePHYLIP(tipranges_object,
    lgdata_fn = "lagrange_area_data_file.data",
    areanames = colnames(tipranges_object@df))
}
\arguments{
  \item{tipranges_object}{An object of class
  \code{tipranges}.}

  \item{lgdata_fn}{The LAGRANGE geographic data file to be
  output.}

  \item{areanames}{A list of the names of the areas.}
}
\value{
  \code{tipranges_object} An object of class
  \code{tipranges}
}
\description{
  Given some geographic range data for tips in the
  \code{tipranges} object, this function exports them to an
  ASCII text file in the Lagrange C++/PHYLIP format
  (\cite{Smith et al. (2010)}).  This file can then be read
  by \code{\link{getranges_from_LagrangePHYLIP}}.
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

# Create an example tipranges object
tipranges = define_tipranges_object()

# See current directory
getwd()

\dontrun{
# Save the file
# Set the filename
fn = "example_tipranges.data"
save_tipranges_to_LagrangePHYLIP(tipranges_object=tipranges, lgdata_fn=fn)

# Show the file
tmplines = scan(file=fn, what="character", sep="\\n")
cat(tmplines, sep="\\n")

# Again, with areanames
save_tipranges_to_LagrangePHYLIP(tipranges_object=tipranges,
lgdata_fn=fn, areanames=c("area1","area2","area3"))

# Show the file
tmplines = scan(file=fn, what="character", sep="\\n")
cat(tmplines, sep="\\n")
} # End dontrun
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
  \code{\link{getranges_from_LagrangePHYLIP}}
}

