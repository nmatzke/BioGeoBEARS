\name{read_controls}
\alias{read_controls}
\title{Read a file with the total number of detections in a taphonomic control}
\usage{
  read_controls(controls_fn, OTUnames = NULL,
    areanames = NULL, tmpskip = 0, phy = NULL)
}
\arguments{
  \item{controls_fn}{The filename of the file containing
  the counts of taphonomic control detections.}

  \item{OTUnames}{Default \code{NULL}, in which case the
  first column of the text file is used as row names/OTU
  names.}

  \item{areanames}{Default \code{NULL}, in which case the
  text file column headings are used.}

  \item{tmpskip}{How many lines should be skipped before
  reading the text file?  Default \code{0}.}

  \item{phy}{An ape phylo object. If included, the rows
  will be sorted to match the order of tree tip labels.}
}
\value{
  \code{dtf}
}
\description{
  This function reads in a tab-delimited text file
  containing counts of detections of the taphonomic
  controls in each region.  These numbers should always be
  equal to or larger than the counts in the detections
  file.
}
\details{
  The idea of taphonomic controls dates back at least to
  work of Bottjer & Jablonski (1988).  The basic idea is
  that if you have taxa of roughly similar detectability,
  then detections of other taxa give some idea of overall
  detection effort.  Obviously this is a very simple model
  that can be criticized in any number of ways (different
  alpha diversity in each region, different detectability
  of individual taxa, etc.), but it is a useful starting
  point as there has been no implementation of any
  detection model in historical/phylogenetic biogeography
  to date.

  One could imagine (a) every OTU and area has a different
  count of detections and taphonomic control detections, or
  (b) the taphonomic control detections are specified by
  area, and shared across all OTUs. Situation (b) is likely
  more common, but this function implements (a).  Behavior
  (b) could be reproduced by summing each column, and/or
  copying this sum to all cells for a particular area.
}
\note{
  Go BEARS!
}
\examples{
testval=1
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS

  Bottjer_Jablonski_1988
}
\seealso{
  \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}}
}

