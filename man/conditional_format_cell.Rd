\name{conditional_format_cell}
\alias{conditional_format_cell}
\title{Conditionally format a number (mostly)}
\usage{
  conditional_format_cell(cellval,
    numbers_below_this_get_scientific = 1e-04,
    numdigits_for_superlow_scientific = 1,
    numbers_above_this_get_scientific = 1e+07,
    numdigits_for_superhigh_scientific = 2,
    numdigits_inbetween_have_fixed_digits = 4)
}
\arguments{
  \item{cellval}{The cell value to format.}

  \item{numbers_below_this_get_scientific}{When the
  absolute value of a number is below this value,
  scientific notation is used.}

  \item{numdigits_for_superlow_scientific}{Number of digits
  after the '.' for scientific notation of small numbers.}

  \item{numbers_above_this_get_scientific}{When the
  absolute value of a number is above this value,
  scientific notation is used.}

  \item{numdigits_for_superhigh_scientific}{Number of
  digits after the '.' for scientific notation of large
  numbers.}

  \item{numdigits_inbetween_have_fixed_digits}{Numbers of
  medium size have this many fixed digits.  Note that other
  cutoffs are specified in the code, and
  \code{\link[base]{signif}} is used to make e.g. integers
  appear as 0, 1, 2..}
}
\value{
  \code{cellval} The value, reformatted and of class
  \code{\link[base]{character}}.
}
\description{
  When a table has numbers that range over many orders of
  magnitude, it can be very distracting if the display
  program forces each column to the same format.  This
  function formats a cell much like Excel would.
}
\details{
  The defaults seem to work well, but could be modified.
  The current function also extracts just the filename, if
  a full path is given.
}
\note{
  Go BEARS!
}
\examples{
test=1

cellval = 143514514514532
conditional_format_cell(cellval)

cellval = -42.235235
conditional_format_cell(cellval)

cellval = -42.0000000
conditional_format_cell(cellval)

cellval = 0.0000
conditional_format_cell(cellval)

cellval = 0.0001
conditional_format_cell(cellval)

cellval = 0.00001
conditional_format_cell(cellval)

cellval = 0.0000111
conditional_format_cell(cellval)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link[base]{signif}}, \code{\link[base]{sprintf}}
}

