\name{divide_probs_by_number_of_options_txt}
\alias{divide_probs_by_number_of_options_txt}
\title{Divide each type of event by its frequency}
\usage{
  divide_probs_by_number_of_options_txt(probmat)
}
\arguments{
  \item{probmat}{A character matrix of probabilities in the
  form of formulas, not normalized by the sum of each row.}
}
\value{
  \code{probmat} A matrix of strings, where each cell
  contains the parameters describing the conditional
  probability of that ancestor-->(Left descendant,Right
  descendant) range inheritance scenario.
}
\description{
  In a speciation/cladogenesis matrix, the conditional
  probabilities of each row must sum to 1.  This function
  sums the number of events of each category and scales
  them accordingly.
}
\details{
  This function returns the strings, which can then be
  processed in other functions by e.g. find/replace or
  \code{\link[base]{eval}}.
}
\note{
  Go BEARS!
}
\examples{
testval=1
probmat = make_relprob_matrix_bi()
probmat

probmat2 = divide_probs_by_number_of_options_txt(probmat)
probmat2
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link{make_relprob_matrix_bi}},
  \code{\link{divide_probs_by_number_of_options_nums}}
}

