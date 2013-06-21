\name{get_ML_state_indices}
\alias{get_ML_state_indices}
\title{Extract the indices for the ML states at each node, given a row of relative probabilities}
\usage{
  get_ML_state_indices(relprobs_row, nums, maxprob)
}
\arguments{
  \item{relprobs_row}{A row from a \code{relprobs}, a
  numeric matrix of relative probabilities}

  \item{nums}{Numbers indexing the states from 1 to
  numstates}

  \item{maxprob}{The value of the maximum probability for
  the row.}
}
\value{
  \code{index_of_ML_state_s}
}
\description{
  Given a table with the rows representing nodes, and the
  columns representing the relative probabilities of each
  state, this function finds the ML (maximum likelihood)
  state(s) for each node;
  \code{\link{get_ML_state_indices}} does this for a row,
  \code{\link{get_ML_states}} iterates over all the rows.
}
\note{
  Go BEARS!
}
\examples{
testval=1
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
  \code{\link{get_ML_states}}
}

