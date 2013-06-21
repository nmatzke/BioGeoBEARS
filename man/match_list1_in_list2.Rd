\name{match_list1_in_list2}
\alias{match_list1_in_list2}
\title{Return TRUE for list1 items when they occur in list2}
\usage{
  match_list1_in_list2(list1, list2)
}
\arguments{
  \item{list1}{The list of things you want to check}

  \item{list2}{The list of things you want to check
  against}
}
\value{
  \code{matchlist} The TRUE/FALSE list for list1
}
\description{
  Return matching TRUE/FALSE values.  E.g. list1 (e.g. a
  big list) TRUE if it is found in list2 (e.g. a smaller
  list)
}
\details{
  Utility function for %in%, when one's brain gets
  confused.
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
  \url{https://code.google.com/p/lagrange/}

  Matzke_2012_IBS
}
\seealso{
  \code{\link[base]{match}}
}

