\name{calc_prob_forward_onebranch_dense}
\alias{calc_prob_forward_onebranch_dense}
\title{Dense matrix exponentiation forward on a branch, with rexpokit}
\usage{
  calc_prob_forward_onebranch_dense(relprobs_branch_bottom,
    branch_length, Qmat)
}
\arguments{
  \item{relprobs_branch_bottom}{The relative probability of
  each state at the base of the branch (should sum to 1).}

  \item{branch_length}{The length of the branch.}

  \item{Qmat}{A Q transition matrix in square (dense)
  format}
}
\value{
  \code{actual_probs_after_forward_exponentiation} The
  probabilities of each state at the top of the branch.
}
\description{
  Take input probabilities, and get the probabilities at
  the end of a branch using matrix exponentiation.
}
\details{
  The \code{\link{calc_loglike_sp}} function calculates
  most transition probabilities internally via
  \code{\link[rexpokit]{rexpokit}}.  These are then stored
  and can be used again when an uppass is being done for
  ancestral state estimates.  However, if there is a root
  branch below the lowest fork, the uppass needs to
  calculate the forward probabilities.
}
\note{
  Go BEARS!
}
\examples{
# Make a square instantaneous rate matrix (Q matrix)
# This matrix is taken from Peter Foster's (2001) "The Idiot's Guide
# to the Zen of Likelihood in a Nutshell in Seven Days for Dummies,
# Unleashed" at:
# \\url{http://www.bioinf.org/molsys/data/idiots.pdf}
#
# The Q matrix includes the stationary base freqencies, which Pmat
# converges to as t becomes large.
Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168,
0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)

relprobs_branch_bottom = c(0.25, 0.25, 0.25, 0.25)

# Make a series of t values
branch_length = 0.1

calc_prob_forward_onebranch_dense(relprobs_branch_bottom, branch_length, Qmat)
calc_prob_forward_onebranch_dense(relprobs_branch_bottom, branch_length=0.5, Qmat)
calc_prob_forward_onebranch_dense(relprobs_branch_bottom, branch_length=1, Qmat)
calc_prob_forward_onebranch_dense(relprobs_branch_bottom, branch_length=2, Qmat)
calc_prob_forward_onebranch_dense(relprobs_branch_bottom, branch_length=10, Qmat)
calc_prob_forward_onebranch_dense(relprobs_branch_bottom, branch_length=20, Qmat)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
  \url{http://www.bioinf.org/molsys/data/idiots.pdf}

  Matzke_2012_IBS

  FosterIdiots
}
\seealso{
  \code{\link{expokit_dgpadm_Qmat2}},
  \code{\link[rexpokit]{expokit_dgpadm_Qmat}},
  \code{\link[rexpokit]{rexpokit}}
}

