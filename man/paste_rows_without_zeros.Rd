\name{paste_rows_without_zeros}
\alias{paste_rows_without_zeros}
\title{Concatenate cells in each row of a text-based transition matrix, excluding zeros}
\usage{
  paste_rows_without_zeros(tmpmat)
}
\arguments{
  \item{tmpmat}{A cladogenesis/speciation probability
  matrix (text-based) to collapse each row of.}
}
\value{
  \code{tmpcol} A list containing each row, concatenated
}
\description{
  This is a utility function for
  \code{\link{make_relprob_txtmatrix_sp1}}.
}
\details{
  Convert e.g.:\cr \cr \code{ A|A A|B A|C A|A,B A|B,C A|A,C
  A|A,B,C}\cr \code{ A s j j 0 0 0 0}\cr \code{ B 0 j 0 0 0
  0 0}\cr \code{ C 0 0 j 0 0 0 0}\cr \code{ A,B 0 v 0 b1 0
  0 0}\cr \code{ B,C 0 0 0 0 j 0 0}\cr \code{ A,C 0 0 v 0 0
  b1 0}\cr \code{ A,B,C 0 0 0 0 v 0 b1}\cr \cr ...to...\cr
  \cr \code{ A B C A,B B,C A,C A,B,C }\cr \code{"s+j+j" "j"
  "j" "v+b1" "j" "v+b1" "v+b1" }\cr \cr
}
\note{
  Go BEARS!
}
\examples{
testval=1

spmat = make_relprob_matrix_bi(states_list=list("_", c("A"),
c("B"), c("C"), c("A","B"), c("B","C"), c("A","C"), c("A","B","C")),
split_ABC=FALSE, splitval="", code_for_overlapping_subsets=NA, printwarn=1)
spmat
tmpcol = paste_rows_without_zeros(tmpmat=spmat)
tmpcol
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link{make_relprob_txtmatrix_sp1}}
}

