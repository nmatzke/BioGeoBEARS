\name{binary_range_to_letter_code_txt}
\alias{binary_range_to_letter_code_txt}
\title{Convert binary presence/absence codes (1/0) to text area names}
\usage{
  binary_range_to_letter_code_txt(tipranges_row, areanames)
}
\arguments{
  \item{tipranges_row}{row of a \code{tipranges} object.}

  \item{areanames}{a list of the names of the areas}
}
\value{
  \code{statename} The corresponding name(s), collapsed
  into a string
}
\description{
  Given a row of a \code{tipranges} object, converts to the
  corresponding name(s), collapsed into a string.  E.g., if
  the areas were \code{(A,B,C,D)}, and the tipranges row
  had \code{(1 0 1 0)}, the output statename would be "AC".
}
\note{
  Go BEARS!
}
\examples{
testval=1
tipranges_row = c(1, 0, 1, 0)
areanames = c("A", "B", "C", "D")
statename = binary_range_to_letter_code_txt(tipranges_row, areanames)
statename
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
  \code{\link{binary_range_to_letter_code_list}},
  \code{\link{tipranges_to_tip_condlikes_of_data_on_each_state}}
}

