\name{letter_strings_to_tipranges_df}
\alias{letter_strings_to_tipranges_df}
\title{Convert ranges in the form of letters (A, AB, BFG, etc.) to a \code{tipranges} object}
\usage{
  letter_strings_to_tipranges_df(letter_strings,
    letter_codes_in_desired_order = "alphabet",
    tipnames_in_order = NULL)
}
\arguments{
  \item{letter_strings}{A list of ranges in concatenated
  letter form ("A", "AB", "BFG", etc.)}

  \item{letter_codes_in_desired_order}{The letter codes in
  the desired order. The default keyword, "alphabet", uses
  the standard 26 capital letters; the output binary codes
  will thus have 26 positions.  If the user inputs fewer
  letters here, or puts them in another order, those will
  be used.}

  \item{tipnames_in_order}{If given, the input tipnames
  will be applied as rownames in the tipranges object.
  Default is \code{NULL}, which results in numbering the
  rows.}
}
\value{
  \code{tipranges} An object of class \code{tipranges}.
}
\description{
  This function converts ranges in the form of concatenated
  letters (A, AB, BFG, etc.) to binary state number codes.
  Via \code{\link{apply}}, this is done to each member of
  the entire input vector of strings.  It outputs
  \code{\link{tipranges}} object.
}
\note{
  Go BEARS!
}
\examples{
testval=1
letter_strings = c("A", "B", "C", "AB", "AC", "BC", "ABC")
letter_strings_to_tipranges_df(letter_strings)

letter_strings = c("A", "B", "C", "AB", "AC", "BC", "ABC")
letter_strings_to_tipranges_df(letter_strings,
tipnames_in_order=paste("tip", seq(1,7), sep=""))
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link{letter_string_to_binary}},
  \code{\link{getranges_from_LagrangePHYLIP}}
}

