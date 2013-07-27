\name{get_ML_states_from_relprobs}
\alias{get_ML_states_from_relprobs}
\title{Extract the ML states at each node, from a table of relative probabilities -- old version}
\usage{
  get_ML_states_from_relprobs(relprobs, statenames,
    returnwhat = "states", if_ties = "takefirst")
}
\arguments{
  \item{relprobs}{A numeric matrix of relative
  probabilities}

  \item{statenames}{The names of the states/geographic
  ranges (e.g., A, AB, CDE, ABD, etc...)}

  \item{returnwhat}{If "indices", return the 0-based
  indices of the states. If "states", return the name of
  the state, based on statenames.}

  \item{if_ties}{What to do with ties. Currently, the only
  option is to take the first (this will be shown in e.g. a
  pie chart, of course).}
}
\value{
  \code{ML_states} or \code{ML_states_indices}, depending
  on \code{returnwhat}.
}
\description{
  Given a table with the rows representing nodes, and the
  columns representing the relative probabilities of each
  state, this function finds the ML (maximum likelihood)
  state(s) for each node.
}
\details{
  If possible, the input matrix should be the actual ML
  estimate of the state probabilities at each node, rather
  than just the scaled conditional likelihoods at each
  node. The latter reflect only the tips-down information,
  whereas the former (the marginal ancestral state
  reconstruction) uses all of the information, and the
  probabilities of the states at the root and in the
  outgroup(s) can influence the estimates in the ingroups.
  This would not likely be particularly important in a pure
  continuous-time model, but in a model with cladogenesis
  it could matter quite a bit.

  See
  \url{http://blog.phytools.org/2013/03/marginal-ancestral-state-reconstruction.html}
  for more discussion of marginal ancestral state
  reconstructions, versus mere scaled conditional
  likelihoods.

  Revell and other sources (\cite{Felsenstein (2004)})
  advocate the "re-rooting" method for obtaining the
  marginal ancestral state reconstructions; however,
  re-rooting requires a time-reversible model and a tree
  with no root.  In biogeography we have a
  \emph{non}-reversible model, and typically a time-scaled
  chronogram.  However, the same result can be obtained by
  modifying the scaled conditional likelihoods obtained
  from a downpass from the tips, via an doing an up-pass
  from the root scaled conditional likelihoods, being
  careful to transfer probabilities via the time-forward
  version of the Q-matrix and cladogenesis/speciation
  matrix.

  \bold{Note:} further notes as this is implemented
  (required!)
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
  \url{http://blog.phytools.org/2013/03/marginal-ancestral-state-reconstruction.html}
  \url{https://code.google.com/p/lagrange/}

  Felsenstein2004

  Matzke_2012_IBS
}
\seealso{
  \code{\link{get_ML_state_indices}}
}

