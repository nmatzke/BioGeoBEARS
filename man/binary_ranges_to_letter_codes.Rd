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

letter_code_ranges = binary_ranges_to_letter_codes(tipranges=tipranges_object, areanames)
letter_code_ranges
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
  \url{https://code.google.com/p/lagrange/}

  Matzke N (2012). "Founder-event speciation in BioGeoBEARS
  package dramatically improves likelihoods and alters
  parameter inference in Dispersal-Extinction-Cladogenesis
  (DEC) analyses." _Frontiers of Biogeography_, *4*(suppl.
  1), pp. 210. ISSN 1948-6596, Poster abstract published in
  the Conference Program and Abstracts of the International
  Biogeography Society 6th Biannual Meeting, Miami,
  Florida. Poster Session P10: Historical and
  Paleo-Biogeography. Poster 129B. January 11, 2013, <URL:
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}>.
}
\seealso{
  \code{\link{binary_range_to_letter_code_list}},
  \code{\link{tipranges_to_tip_condlikes_of_data_on_each_state}}
}

