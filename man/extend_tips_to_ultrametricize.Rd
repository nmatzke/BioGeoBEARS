\name{extend_tips_to_ultrametricize}
\alias{extend_tips_to_ultrametricize}
\title{Take a tree, extend all tips (including fossils) up to 0.0 my before present}
\usage{
  extend_tips_to_ultrametricize(obj, age_of_root = 0,
    tips_end_at_this_date = NA)
}
\arguments{
  \item{obj}{An \code{\link[ape]{ape}}
  \code{\link[ape]{phylo}} object.}

  \item{age_of_root}{The length of the branch below the
  root. Default 0.}

  \item{tips_end_at_this_date}{The tips can be set to
  something other than 0, if desired.  (This could produce
  negative branclengths, however.)}
}
\value{
  \code{obj} The corrected phylogeny
}
\description{
  Makes tree precisely ultrametric by extending the
  terminal branches up to the highest tip (which is treated
  as 0 my before present).
}
\details{
  This function ADDS the time_before_present to everything,
  including fossils.  You have been warned.
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
  \code{\link[ape]{read.tree}}, \code{\link{prt}},
  \code{\link{average_tr_tips}}
}

