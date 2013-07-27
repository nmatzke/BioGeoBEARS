\name{get_indices_where_list1_occurs_in_list2_noNA}
\alias{get_indices_where_list1_occurs_in_list2_noNA}
\title{Return (first!) indices in second list matching the first list, excluding NAs}
\usage{
  get_indices_where_list1_occurs_in_list2_noNA(list1,
    list2)
}
\arguments{
  \item{list1}{The first list.}

  \item{list2}{The second list list.}
}
\value{
  \code{match_indices} The match indices.
}
\description{
  This function will return one match (the first) for each
  item in the list; i.e. the second-list index for each
  item in the first list.  Only the first hit in the second
  list is returned.  Unlike
  \code{\link{get_indices_where_list1_occurs_in_list2}},
  non-hits (NAs) are excluded.
}
\details{
  This is used by get_indices_of_branches_under_tips, which
  is used by \code{\link{extend_tips_to_ultrametricize}},
  which can be used by section_the_tree.
}
\note{
  Go BEARS!
}
\examples{
list1 = c("N", "I", "C", "K")
list2 = LETTERS
get_indices_where_list1_occurs_in_list2_noNA(list1, list2)
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
  \code{\link{get_indices_where_list1_occurs_in_list2}},
  \code{\link{extend_tips_to_ultrametricize}},
  \code{\link{section_the_tree}},
  \code{\link{return_items_not_NA}}
}

