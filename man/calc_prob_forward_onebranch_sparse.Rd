\name{calc_prob_forward_onebranch_sparse}
\alias{calc_prob_forward_onebranch_sparse}
\title{Sparse matrix exponentiation forward on a branch, with rexpokit}
\usage{
  calc_prob_forward_onebranch_sparse(relprobs_branch_bottom,
    branch_length, tmpQmat_in_REXPOKIT_coo_fmt,
    coo_n = coo_n, anorm = anorm, check_for_0_rows = TRUE,
    TRANSPOSE_because_forward = TRUE)
}
\arguments{
  \item{relprobs_branch_bottom}{The relative probability of
  each state at the base of the branch (should sum to 1).}

  \item{branch_length}{The length of the branch.}

  \item{tmpQmat_in_REXPOKIT_coo_fmt}{A Q transition matrix
  in sparse (COO) format.  See
  \code{\link[rexpokit]{mat2coo}}.}

  \item{coo_n}{If a COO matrix is input, \code{coo_n}
  specified the order (# rows, equals # columns) of the
  matrix.}

  \item{anorm}{\code{dgexpv} requires an initial guess at
  the norm of the matrix. Using the R function
  \code{\link{norm}} might get slow with large matrices. If
  so, the user can input a guess manually (\code{Lagrange}
  seems to just use 1 or 0, if I recall correctly).}

  \item{check_for_0_rows}{If TRUE or a numeric value, the
  input Qmat is checked for all-zero rows, since these will
  crash the FORTRAN \code{wrapalldmexpv} function. A small
  nonzero value set to check_for_0_rows or the default
  (0.0000000000001) is input to off-diagonal cells in the
  row (and the diagonal value is normalized), which should
  fix the problem.}

  \item{TRANSPOSE_because_forward}{For non-time-reversible
  models, the forward calculation is different than the
  backward one. Fortunately this just means switching the
  rows and columns of a transition matrix.}
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
require("rexpokit")

Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504,
0.168, 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
tmpQmat_in_REXPOKIT_coo_fmt = mat2coo(Qmat)

relprobs_branch_bottom = c(0.25, 0.25, 0.25, 0.25)

# Make a series of t values
branch_length = 0.1

calc_prob_forward_onebranch_sparse(relprobs_branch_bottom, branch_length,
tmpQmat_in_REXPOKIT_coo_fmt, coo_n=4, anorm=1, check_for_0_rows=TRUE,
TRANSPOSE_because_forward=TRUE)
calc_prob_forward_onebranch_sparse(relprobs_branch_bottom, branch_length=0.5,
tmpQmat_in_REXPOKIT_coo_fmt, coo_n=4, anorm=1, check_for_0_rows=TRUE,
TRANSPOSE_because_forward=TRUE)
calc_prob_forward_onebranch_sparse(relprobs_branch_bottom, branch_length=1,
tmpQmat_in_REXPOKIT_coo_fmt, coo_n=4, anorm=1, check_for_0_rows=TRUE,
TRANSPOSE_because_forward=TRUE)
calc_prob_forward_onebranch_sparse(relprobs_branch_bottom, branch_length=2,
tmpQmat_in_REXPOKIT_coo_fmt, coo_n=4, anorm=1, check_for_0_rows=TRUE,
TRANSPOSE_because_forward=TRUE)
calc_prob_forward_onebranch_sparse(relprobs_branch_bottom, branch_length=10,
tmpQmat_in_REXPOKIT_coo_fmt, coo_n=4, anorm=1, check_for_0_rows=TRUE,
TRANSPOSE_because_forward=TRUE)
calc_prob_forward_onebranch_sparse(relprobs_branch_bottom, branch_length=20,
tmpQmat_in_REXPOKIT_coo_fmt, coo_n=4, anorm=1, check_for_0_rows=TRUE,
TRANSPOSE_because_forward=TRUE)
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

