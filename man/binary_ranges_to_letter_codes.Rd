\name{binary_ranges_to_letter_codes}
\alias{binary_ranges_to_letter_codes}
\title{Convert binary presence/absence codes (1/0) to a list of text area names}
\usage{
  binary_ranges_to_letter_codes(tipranges, areanames)
}
\arguments{
  \item{tipranges}{a \code{tipranges} object.}

  \item{areanames}{a list of the names of the areas}
}
\value{
  \code{letter_code_ranges} A list of the states -- there
  will be as many states as there are rows/tips in
  \code{tipranges}. Each state will be a list of area
  names.
}
\description{
  Given a row of a \code{tipranges} object, converts to a
  list of the corresponding statenames for each row.
}
\note{
  Go BEARS!
}
\examples{
testval=1
# Define a tipranges object
tipranges_object = define_tipranges_object()
tipranges_object

areanames = getareas_from_tipranges_object(tipranges_object)
areanames

letter_code_ranges = binary_ranges_to_letter_codes(tipranges=tipranges_object,
areanames)
letter_code_ranges
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
  \code{\link{letter_string_to_binary}},
  \code{\link{letter_strings_to_tipranges_df}},
  \code{\link{tipranges_to_tip_condlikes_of_data_on_each_state}}
}

