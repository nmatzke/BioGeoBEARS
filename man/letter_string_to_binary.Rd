\name{letter_string_to_binary}
\alias{letter_string_to_binary}
\title{Convert ranges in the form of letters (A, AB, BFG, etc.) to binary state number codes}
\usage{
  letter_string_to_binary(letter_string,
    letter_codes_in_desired_order = "alphabet")
}
\arguments{
  \item{letter_string}{A string of letters (e.g. "ABD")}

  \item{letter_codes_in_desired_order}{The letter codes in
  the desired order. The default keyword, "alphabet", uses
  the standard 26 capital letters; the output binary codes
  will thus have 26 positions.  If the user inputs fewer
  letters here, or puts them in another order, those will
  be used.}
}
\value{
  \code{numcodes} A list with the binary codes.
}
\description{
  This function takes a letter string (e.g. ABD) and
  converts to binary encoding (e.g. 1101).
}
\note{
  Go BEARS!
}
\examples{
testval=1
letter_string = "ABD"
letter_string_to_binary(letter_string, letter_codes_in_desired_order="alphabet")

letter_string = "ABD"
letter_string_to_binary(letter_string,
letter_codes_in_desired_order=c("A","B","C","D","E","F"))

letter_string = "ABD"
letter_string_to_binary(letter_string,
letter_codes_in_desired_order=strsplit("ABCDEF", split="")[[1]])
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link{binary_ranges_to_letter_codes}}
}

