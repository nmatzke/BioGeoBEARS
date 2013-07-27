\name{readfiles_BioGeoBEARS_run}
\alias{readfiles_BioGeoBEARS_run}
\title{Read in the extra input files, if any}
\usage{
  readfiles_BioGeoBEARS_run(inputs)
}
\arguments{
  \item{inputs}{The inputs list}
}
\value{
  \code{inputs} The modified inputs list
}
\description{
  This function reads input files for stratification,
  constraints, and detection, i.e., everything except the
  tree and geography files. E.g., \code{areas_allowed_fn}
  file is just a list of distance matrices, separated by
  blank lines, from youngest to oldest.
}
\note{
  Go BEARS!
}
\examples{
test=1
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link{define_BioGeoBEARS_run}},
  \code{\link{read_times_fn}},
  \code{\link{read_distances_fn}},
  \code{\link{read_dispersal_multipliers_fn}},
  \code{\link{read_area_of_areas_fn}},
  \code{\link{read_areas_allowed_fn}},
  \code{\link{read_detections}},
  \code{\link{read_controls}}
}

