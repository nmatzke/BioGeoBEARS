\name{divide_probs_by_number_of_options_nums}
\alias{divide_probs_by_number_of_options_nums}
\title{Divide each type of event by its frequency, return calculated probabilities}
\usage{
  divide_probs_by_number_of_options_nums(spPmat, probmat)
}
\arguments{
  \item{spPmat}{A matrix of numbers, where each cell
  contains the conditional probability of that
  ancestor-->(Left descendant,Right descendant) range
  inheritance scenario.}

  \item{probmat}{A matrix of text, describing each of the
  allowed range-inheritance events.}
}
\value{
  \code{spPmat} A matrix of numbers, where each cell
  contains the conditional probability of that
  ancestor-->(Left descendant,Right descendant) range
  inheritance scenario.
}
\description{
  In a speciation/cladogenesis matrix, the conditional
  probabilities of each row must sum to 1.  This function
  sums the number of events of each category and scales
  them accordingly.
}
\details{
  This function returns the calculated conditional
  probabilities.
}
\note{
  Go BEARS!
}
\examples{
testval=1
spmat = make_relprob_matrix_bi()
spmat

spmat1 = divide_probs_by_number_of_options_txt(spmat)
spmat1


probmat = spmat
spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\\\+",
mergesym="*", ys=1, j=0, v=1, maxent_constraint_01=0.0001,
maxent_constraint_01v=0.0001, max_numareas=3)
spPmat
probmat2 = divide_probs_by_number_of_options_nums(spPmat, probmat)
probmat2

probmat = spmat1
spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\\\+",
mergesym="*", ys=1, j=0, v=1, maxent_constraint_01=0.0001,
maxent_constraint_01v=0.0001, max_numareas=3)
spPmat
probmat3 = divide_probs_by_number_of_options_nums(spPmat, probmat)
probmat3
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
  \code{\link{divide_probs_by_number_of_options_txt}}
}

