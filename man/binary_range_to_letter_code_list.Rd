\name{binary_range_to_letter_code_list}
\alias{binary_range_to_letter_code_list}
\title{Convert binary presence/absence codes (1/0) to a list of text area names}
\usage{
  binary_range_to_letter_code_list(tipranges_row,
    areanames)
}
\arguments{
  \item{tipranges_row}{row of a \code{tipranges} object.}

  \item{areanames}{a list of the names of the areas}
}
\value{
  \code{list_of_areas_in_the_state} A list of the name(s)
  of the areas corresponding to the presence/absence coding
  in the row
}
\description{
  Given a row of a \code{tipranges} object, converts to a
  list of the corresponding name(s). E.g., if the areas
  were \code{(A,B,C,D)}, and the tipranges row had \code{(1
  0 1 0)}, the output statename would be ("A","C").
}
\note{
  Go BEARS!
}
\examples{
testval=1
tipranges_row = c(1, 0, 1, 0)
areanames = c("A", "B", "C", "D")
list_of_areas_in_the_state = binary_range_to_letter_code_list(tipranges_row,
areanames)
list_of_areas_in_the_state
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link{binary_ranges_to_letter_codes}},
  \code{\link{binary_range_to_letter_code_txt}},
  \code{\link{tipranges_to_tip_condlikes_of_data_on_each_state}}
}

