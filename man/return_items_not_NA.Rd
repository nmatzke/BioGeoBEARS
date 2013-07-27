\name{return_items_not_NA}
\alias{return_items_not_NA}
\title{Remove NAs from a vector/list}
\usage{
  return_items_not_NA(x)
}
\arguments{
  \item{x}{The vector of items to check for being not NA.}
}
\value{
  \code{y} The surviving, non-NA cells of a vector.
}
\description{
  Utility function. This function returns the non-NA values
  from a vector.
}
\details{
  This is used by
  \code{\link{get_indices_where_list1_occurs_in_list2_noNA}},
  which is used by
  \code{\link{get_indices_of_branches_under_tips}}, which
  is used by \code{\link{extend_tips_to_ultrametricize}},
  which can be used by \code{\link{section_the_tree}}.
}
\note{
  Go BEARS!
}
\examples{
list1 = c("N", "I", NA, "C", "K")
return_items_not_NA(list1)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link{prt}}, \code{\link[base]{LETTERS}},
  \code{\link{get_indices_where_list1_occurs_in_list2_noNA}},
  \code{\link{get_indices_where_list1_occurs_in_list2}},
  \code{\link{extend_tips_to_ultrametricize}},
  \code{\link{section_the_tree}}
}

