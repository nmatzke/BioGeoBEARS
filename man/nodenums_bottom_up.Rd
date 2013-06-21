\name{nodenums_bottom_up}
\alias{nodenums_bottom_up}
\title{Assign node labels in bottom-up, left-first format (as in e.g. r8s)}
\usage{
  nodenums_bottom_up(tr)
}
\arguments{
  \item{tr}{A tree object in \code{\link[ape]{phylo}}
  format.}
}
\value{
  \code{traverse_records}
}
\description{
  This function assigns node numbers by tracing up from the
  root.  This corresponds to the node numbers in e.g. r8s
  (\cite{Sanderson (2003)}).
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

  Sanderson_2003_r8s

  Matzke_2012_IBS

  r8s

  Marazzi_etal_Sanderson_2012_r8s_morph
}
\seealso{
  \code{\link[phylobase]{phylo4}},
}

