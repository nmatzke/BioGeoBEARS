#######################################################
# Calc loglike speciation -- separated to its own file, because it's such a huge function
#######################################################
require("ape")
require("rexpokit")
require("cladoRcpp")
#require("parallel")

#######################################################
# Helper functions to try to speed things up
#######################################################	


#######################################################
# expokit_dgpadm_Qmat2_prebyte
#######################################################
#' A version of expokit_dgpadm_Qmat to byte-compile
#' 
#' Byte-compiling is supposed to speed up functions; this is an 
#' attempt to do this on the \code{\link[rexpokit]{rexpokit}} function \code{\link[rexpokit]{expokit_dgpadm_Qmat}}.  It
#' is also possible to byte-compile everything during package installation (via \code{ByteCompile: true} in the
#' DESCRIPTION file), which is implemented in \code{BioGeoBEARS}, so this may be redundant.
#'
#' \code{\link{expokit_dgpadm_Qmat2_prebyte}} gets byte-compiled into \code{\link{expokit_dgpadm_Qmat2}}.
#' 
#' See \url{http://dirk.eddelbuettel.com/blog/2011/04/12/} for discussion of the \code{compiler} (\code{\link[compiler]{cmpfun}}) package.
#' 
#' @param times one or more time values to exponentiate by
#' @param Qmat an input Q transition matrix
#' @param transpose_needed If TRUE (default), matrix will be transposed (apparently
#' EXPOKIT needs the input matrix to be transposed compared to normal)
#' @return \code{tmpoutmat} the output matrix.
#' @export
#' @seealso \code{\link[rexpokit]{expokit_dgpadm_Qmat}}, \code{\link{expokit_dgpadm_Qmat2}}, \code{\link[compiler]{compile}}, \code{\link[compiler]{cmpfun}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' 
expokit_dgpadm_Qmat2_prebyte <- function(times, Qmat, transpose_needed=TRUE)
	{
	require("rexpokit")
	rexpokit::expokit_dgpadm_Qmat(Qmat=Qmat, t=times, transpose_needed)
	}


#######################################################
# expokit_dgpadm_Qmat2
#######################################################
#' A byte-compiled version of expokit_dgpadm_Qmat2_prebyte
#' 
#' Byte-compiling is supposed to speed up functions; this is an 
#' attempt to do this on the \code{\link[rexpokit]{rexpokit}} function \code{\link[rexpokit]{expokit_dgpadm_Qmat}}.  It
#' is also possible to byte-compile everything during package installation (via \code{ByteCompile: true} in the
#' DESCRIPTION file), which is implemented in \code{BioGeoBEARS}, so this may be redundant.
#'
#' \code{\link{expokit_dgpadm_Qmat2_prebyte}} gets byte-compiled into \code{\link{expokit_dgpadm_Qmat2}}.
#' 
#' See \url{http://dirk.eddelbuettel.com/blog/2011/04/12/} for discussion of the \code{\link[compiler]{compile}} package.
#' 
#' @param times one or more time values to exponentiate by
#' @param Qmat an input Q transition matrix
#' @param transpose_needed If TRUE (default), matrix will be transposed (apparently
#' EXPOKIT needs the input matrix to be transposed compared to normal)
#' @return \code{tmpoutmat} the output matrix.
#' @export
#' @seealso \code{\link[rexpokit]{expokit_dgpadm_Qmat}}, \code{\link{expokit_dgpadm_Qmat2}}, \code{\link[compiler]{compile}}, \code{\link[compiler]{cmpfun}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' 
expokit_dgpadm_Qmat2 = compiler::cmpfun(expokit_dgpadm_Qmat2_prebyte)





#######################################################
# mapply_likelihoods_prebyte
#######################################################
#' Use mapply on matrix exponentiations -- pre-byte-compiling
#' 
#' During the likelihood calculations from the tips to the root of a tree, the transition matrix Qmat needs to 
#' be exponentiated for each branch length in the tree.  This is the slowest step of the likelihood
#' calculation, especially for large matrices.  This function performs this with mapply.
#' 
#' Byte-compiling is supposed to speed up functions; this is an 
#' attempt to do this on the \code{\link[rexpokit]{rexpokit}} function \code{\link[rexpokit]{expokit_dgpadm_Qmat}}.  It
#' is also possible to byte-compile everything during package installation (via \code{ByteCompile: true} in the
#' DESCRIPTION file), which is implemented in \code{BioGeoBEARS}, so this may be redundant.
#'
#' \code{\link{mapply_likelihoods_prebyte}} gets byte-compiled into \code{\link{mapply_likelihoods}}.
#' 
#' See \url{http://dirk.eddelbuettel.com/blog/2011/04/12/} for discussion of the \code{\link[compiler]{compile}} package.
#' 
#' @param Qmat an input Q transition matrix.
#' @param phy2 A phylogenetic tree.
#' @param transpose_needed If TRUE (default), matrix will be transposed (apparently EXPOKIT needs the input matrix to be transposed compared to normal).
#' @return \code{independent_likelihoods_on_each_branch} The output matrix of the likelihoods for each state on each branch.
#' @export
#' @seealso \code{\link[base]{mapply}}, \code{\link[rexpokit]{expokit_dgpadm_Qmat}}, \code{\link{expokit_dgpadm_Qmat2}}, \code{\link[compiler]{compile}}, \code{\link[compiler]{cmpfun}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' 
mapply_likelihoods_prebyte <- function(Qmat, phy2, transpose_needed)
	{
	dims_Qmat = dim(Qmat)

	# Not parallel processing
	independent_likelihoods_on_each_branch = mapply(FUN=expokit_dgpadm_Qmat2, Qmat=list(Qmat), times=phy2$edge.length, transpose_needed=TRUE, SIMPLIFY="array")	

	# 2016-05-28_bug_fix: dims_Qmat
	# Fix this error, e.g. when DEC* model + areas_allowed means that
	# ranges_list = NULL, Kauai is just
	# ranges_list = Kauai
	# This means that:
	#   subtree_tip_relative_probs_of_each_state 
	# and thus
	#   tip_condlikes_of_data_on_each_state
	# ...are just a list of numbers, not a matrix, thus 
	# rowSums fails in calc_loglike_sp() in that time-stratum.
	# 
	#
	# This was the error message:
	# 
	# Error in rowSums(tip_condlikes_of_data_on_each_state) : 
	#  'x' must be an array of at least two dimensions
	# Calls: bears_optim_run ... calc_loglike_sp_stratified -> calc_loglike_sp -> rowSums
	#
	tmp_dims = dim(matrix(data=1,ncol=1,nrow=1))
	if ( identical(x=dims_Qmat, y=tmp_dims) )
		{
		independent_likelihoods_on_each_branch = array(data=independent_likelihoods_on_each_branch, dim=c(1,1,length(phy2$edge.length)))
		}

	return(independent_likelihoods_on_each_branch)
	}

#######################################################
# mapply_likelihoods
#######################################################
#' Use mapply on matrix exponentiations -- post-byte-compiling
#' 
#' During the likelihood calculations from the tips to the root of a tree, the transition matrix Qmat needs to 
#' be exponentiated for each branch length in the tree.  This is the slowest step of the likelihood
#' calculation, especially for large matrices.  This function performs this with mapply.
#' 
#' Byte-compiling is supposed to speed up functions; this is an 
#' attempt to do this on the \code{\link[rexpokit]{rexpokit}} function \code{\link[rexpokit]{expokit_dgpadm_Qmat}}.  It
#' is also possible to byte-compile everything during package installation (via \code{ByteCompile: true} in the
#' DESCRIPTION file), which is implemented in \code{BioGeoBEARS}, so this may be redundant.
#'
#' \code{\link{mapply_likelihoods_prebyte}} gets byte-compiled into \code{\link{mapply_likelihoods}}.
#' 
#' See \url{http://dirk.eddelbuettel.com/blog/2011/04/12/} for discussion of the \code{\link[compiler]{compile}} package.
#' 
#' @param Qmat an input Q transition matrix.
#' @param phy2 A phylogenetic tree.
#' @param transpose_needed If TRUE (default), matrix will be transposed (apparently EXPOKIT needs the input matrix to be transposed compared to normal).
#' @return \code{independent_likelihoods_on_each_branch} The output matrix of the likelihoods for each state on each branch.
#' @export
#' @seealso \code{\link[base]{mapply}}, \code{\link[rexpokit]{expokit_dgpadm_Qmat}}, \code{\link{expokit_dgpadm_Qmat2}}, \code{\link[compiler]{compile}}, \code{\link[compiler]{cmpfun}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#'
mapply_likelihoods = compiler::cmpfun(mapply_likelihoods_prebyte)










#######################################################
# calc_prob_forward_onebranch_dense
#######################################################
#' Dense matrix exponentiation forward on a branch, with rexpokit
#' 
#' Take input probabilities, and get the probabilities at the end of a branch using matrix exponentiation.
#' 
#' The \code{\link{calc_loglike_sp}} function calculates most transition probabilities internally
#' via \code{\link[rexpokit]{rexpokit}}.  These are then stored and can be used again
#' when an uppass is being done for ancestral state estimates.  However, if there is a root
#' branch below the lowest fork, the uppass needs to calculate the forward probabilities.
#' 
#' @param relprobs_branch_bottom The relative probability of each state at the base of the branch (should sum to 1).
#' @param branch_length The length of the branch.
#' @param Qmat A Q transition matrix in square (dense) format
#' @return \code{actual_probs_after_forward_exponentiation} The probabilities of each state at the top of the branch.
#' @export
#' @seealso \code{\link{expokit_dgpadm_Qmat2}}, \code{\link[rexpokit]{expokit_dgpadm_Qmat}}, \code{\link[rexpokit]{rexpokit}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.bioinf.org/molsys/data/idiots.pdf}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite FosterIdiots
#' @examples
#' # Make a square instantaneous rate matrix (Q matrix)
#' # This matrix is taken from Peter Foster's (2001) "The Idiot's Guide
#' # to the Zen of Likelihood in a Nutshell in Seven Days for Dummies,
#' # Unleashed" at:
#' # \url{http://www.bioinf.org/molsys/data/idiots.pdf}
#' #
#' # The Q matrix includes the stationary base freqencies, which Pmat 
#' # converges to as t becomes large.
#' Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 
#' 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
#' 
#' relprobs_branch_bottom = c(0.25, 0.25, 0.25, 0.25)
#' 
#' # Make a series of t values
#' branch_length = 0.1
#' 
#' calc_prob_forward_onebranch_dense(relprobs_branch_bottom, branch_length, Qmat)
#' calc_prob_forward_onebranch_dense(relprobs_branch_bottom, branch_length=0.5, Qmat)
#' calc_prob_forward_onebranch_dense(relprobs_branch_bottom, branch_length=1, Qmat)
#' calc_prob_forward_onebranch_dense(relprobs_branch_bottom, branch_length=2, Qmat)
#' calc_prob_forward_onebranch_dense(relprobs_branch_bottom, branch_length=10, Qmat)
#' calc_prob_forward_onebranch_dense(relprobs_branch_bottom, branch_length=20, Qmat)
#' 
calc_prob_forward_onebranch_dense <- function(relprobs_branch_bottom, branch_length, Qmat)
	{
	defaults='
	relprobs_branch_bottom = rep(1/nrow(Qmat), nrow(Qmat))
	branch_length = 1
	'
	
	# Dense matrix exponentiation, FORWARD
	# To get forward behavior, in this case, do transpose; what matters is the order in which you use %*%
	
	# cond_probs, aka independent_likelihoods
	# branch_length = 1
	cond_probs_after_forward_exponentiation = expokit_dgpadm_Qmat2(times=branch_length, Qmat=Qmat, transpose_needed=TRUE)
	
	# Probabilities of states at the top of the branch
	actual_probs_after_forward_exponentiation = relprobs_branch_bottom %*% cond_probs_after_forward_exponentiation
	actual_probs_after_forward_exponentiation
	
	return(actual_probs_after_forward_exponentiation)
	}




#######################################################
# calc_prob_forward_onebranch_sparse
#######################################################
#' Sparse matrix exponentiation forward on a branch, with rexpokit
#' 
#' Take input probabilities, and get the probabilities at the end of a branch using matrix exponentiation.
#' 
#' The \code{\link{calc_loglike_sp}} function calculates most transition probabilities internally
#' via \code{\link[rexpokit]{rexpokit}}.  These are then stored and can be used again
#' when an uppass is being done for ancestral state estimates.  However, if there is a root
#' branch below the lowest fork, the uppass needs to calculate the forward probabilities.
#' 
#' @param relprobs_branch_bottom The relative probability of each state at the base of the branch (should sum to 1).
#' @param branch_length The length of the branch.
#' @param tmpQmat_in_REXPOKIT_coo_fmt A Q transition matrix in sparse (COO) format.  See \code{\link[rexpokit]{mat2coo}}.
#' @param coo_n If a COO matrix is input, \code{coo_n} specified the order (# rows, equals # columns) of the matrix.
#' @param anorm \code{dgexpv} requires an initial guess at the norm of the matrix. Using the
#' R function \code{\link{norm}} might get slow with large matrices. If so, the user
#' can input a guess manually (\code{Lagrange} seems to just use 1 or 0, if I
#' recall correctly).
#' @param check_for_0_rows If TRUE or a numeric value, the input Qmat is checked for all-zero rows, since these will crash the FORTRAN 
#' \code{wrapalldmexpv} function. A small nonzero value set to check_for_0_rows or the default (0.0000000000001) is input to  off-diagonal
#' cells in the row (and the diagonal value is normalized), which should fix the problem.
#' @param TRANSPOSE_because_forward For non-time-reversible models, the forward calculation is different than the backward one.  
#' Fortunately this just means switching the rows and columns of a transition matrix.
#' @return \code{actual_probs_after_forward_exponentiation} The probabilities of each state at the top of the branch.
#' @export
#' @seealso \code{\link{expokit_dgpadm_Qmat2}}, \code{\link[rexpokit]{expokit_dgpadm_Qmat}}, \code{\link[rexpokit]{rexpokit}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.bioinf.org/molsys/data/idiots.pdf}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite FosterIdiots
#' @examples
#' # Make a square instantaneous rate matrix (Q matrix)
#' # This matrix is taken from Peter Foster's (2001) "The Idiot's Guide
#' # to the Zen of Likelihood in a Nutshell in Seven Days for Dummies,
#' # Unleashed" at:
#' # \url{http://www.bioinf.org/molsys/data/idiots.pdf}
#' #
#' # The Q matrix includes the stationary base freqencies, which Pmat 
#' # converges to as t becomes large.
#' require("rexpokit")
#' 
#' Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 
#' 0.168, 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
#' tmpQmat_in_REXPOKIT_coo_fmt = mat2coo(Qmat)
#' 
#' relprobs_branch_bottom = c(0.25, 0.25, 0.25, 0.25)
#' 
#' # Make a series of t values
#' branch_length = 0.1
#' 
#' calc_prob_forward_onebranch_sparse(relprobs_branch_bottom, branch_length, 
#' tmpQmat_in_REXPOKIT_coo_fmt, coo_n=4, anorm=1, check_for_0_rows=TRUE, 
#' TRANSPOSE_because_forward=TRUE)
#' calc_prob_forward_onebranch_sparse(relprobs_branch_bottom, branch_length=0.5, 
#' tmpQmat_in_REXPOKIT_coo_fmt, coo_n=4, anorm=1, check_for_0_rows=TRUE, 
#' TRANSPOSE_because_forward=TRUE)
#' calc_prob_forward_onebranch_sparse(relprobs_branch_bottom, branch_length=1, 
#' tmpQmat_in_REXPOKIT_coo_fmt, coo_n=4, anorm=1, check_for_0_rows=TRUE, 
#' TRANSPOSE_because_forward=TRUE)
#' calc_prob_forward_onebranch_sparse(relprobs_branch_bottom, branch_length=2, 
#' tmpQmat_in_REXPOKIT_coo_fmt, coo_n=4, anorm=1, check_for_0_rows=TRUE, 
#' TRANSPOSE_because_forward=TRUE)
#' calc_prob_forward_onebranch_sparse(relprobs_branch_bottom, branch_length=10, 
#' tmpQmat_in_REXPOKIT_coo_fmt, coo_n=4, anorm=1, check_for_0_rows=TRUE, 
#' TRANSPOSE_because_forward=TRUE)
#' calc_prob_forward_onebranch_sparse(relprobs_branch_bottom, branch_length=20, 
#' tmpQmat_in_REXPOKIT_coo_fmt, coo_n=4, anorm=1, check_for_0_rows=TRUE, 
#' TRANSPOSE_because_forward=TRUE)
#' 
calc_prob_forward_onebranch_sparse <- function(relprobs_branch_bottom, branch_length, tmpQmat_in_REXPOKIT_coo_fmt, coo_n, anorm, check_for_0_rows=TRUE, TRANSPOSE_because_forward=TRUE, printlevel=1)
	{
	defaults='

	relprobs_branch_bottom = rep(1/nrow(Qmat), nrow(Qmat))
	branch_length = 1

	# tmpQmat_in_REXPOKIT_coo_fmt = Qmat
	tmpQmat_in_REXPOKIT_coo_fmt = mat2coo(Qmat)
	tmpQmat_in_REXPOKIT_coo_fmt_original = tmpQmat_in_REXPOKIT_coo_fmt
	coo_n = 16
	anorm = 1
	'
	
	# Store the number of states/order of the matrix
	numstates = coo_n
	
	# Sparse matrix exponentiation, FORWARD
	
	#condlikes_Left = try (
	#			expokit_dgexpv_Qmat(Qmat=tmpQmat_in_REXPOKIT_coo_fmt, t=phy2$edge.length[i], inputprobs_for_fast=relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,], transpose_needed=FALSE, transform_to_coo_TF=FALSE, coo_n=coo_n, anorm=anorm, check_for_0_rows=TRUE)
	
	# Note 2013-04-16: transpose_needed only works on dense matrices!
	# For sparse matrices, switch columns manually
	if (TRANSPOSE_because_forward == TRUE)
		{
		tmpQmat_in_REXPOKIT_coo_fmt_transposed = tmpQmat_in_REXPOKIT_coo_fmt
		tmpQmat_in_REXPOKIT_coo_fmt_transposed[,1:2] = tmpQmat_in_REXPOKIT_coo_fmt[,2:1]
		}
	
	
# 	cond_probs_top = try(
# 		expokit_dgexpv_Qmat(Qmat=tmpQmat_in_REXPOKIT_coo_fmt_transposed, 
# 			t=branch_length, 
# 			inputprobs_for_fast=relprobs_branch_bottom, 
# 			transpose_needed=FALSE, 
# 			transform_to_coo_TF=FALSE, 
# 			coo_n=coo_n, 
# 			anorm=anorm, 
# 			check_for_0_rows=TRUE)
# 		)
# 	cond_probs_top
	
	# Trying with kexpmv
	tmpQmat_in_KEXPMV_crs_fmt_transposed = coo2crs(
			ia=tmpQmat_in_REXPOKIT_coo_fmt_transposed[,"ia"], 
			ja=tmpQmat_in_REXPOKIT_coo_fmt_transposed[,"ja"], 
			a =tmpQmat_in_REXPOKIT_coo_fmt_transposed[,"a"],
			n=coo_n, transpose_needed=FALSE)
	
	#print(tmpQmat_in_KEXPMV_crs_fmt_transposed)
	#print(branch_length)
		
	#print(relprobs_branch_bottom)
	#print(dim(relprobs_branch_bottom))
	#print(length(relprobs_branch_bottom))	
	
	if (printlevel >= 1)
		{
		cat("U")
		}
	cond_probs_top_try = try(
		kexpmv::expokit_dmexpv(mat=tmpQmat_in_KEXPMV_crs_fmt_transposed, 
			t=branch_length, 
			vector=relprobs_branch_bottom, 
			transpose_needed=FALSE, 
			transform_to_crs=FALSE, 
			crs_n=length(relprobs_branch_bottom), 
			anorm=NULL, 
			mxstep=10000,
			tol=as.numeric(1e-10))
		)
	
	# Error check
	if (class(cond_probs_top_try) == "try-error")
		{
		cat("\n\nIn calc_prob_forward_onebranch_sparse():\n\n")
		cat("\n\ntry-error on kexpmv::expokit_dmexpv():\n\n")
		cat("t=", branch_length, "\n")
		print("tmpQmat_in_KEXPMV_crs_fmt_transposed")
		print(tmpQmat_in_KEXPMV_crs_fmt_transposed)
		print(coo_n)
		print(anorm)
		print("Input vector: relprobs_branch_bottom")
		print(relprobs_branch_bottom)
		print("cond_probs_top_try:")
		print(cond_probs_top_try)
		stop("In calc_prob_forward_onebranch_sparse(), try-error on kexpmv::expokit_dmexpv(); see screen output above.")
		}
	
	# Extract the probabilities vector
	cond_probs_top = c(cond_probs_top_try$output_probs)
	cond_probs_top[cond_probs_top<0] = 0.0
	cond_probs_top[cond_probs_top>1] = 1.0
	
	if (any(is.nan(cond_probs_top)))
		{
		cat("\n\nIn calc_prob_forward_onebranch_sparse():\n\n")
		cat("\n\nnan error on kexpmv::expokit_dmexpv():\n\n")
		cat("\n\nBranch is the root branch:\n\n")
		cat("t=", branch_length, "\n")
		print(tmpQmat_in_KEXPMV_crs_fmt_transposed)
		print(coo_n)
		print(anorm)
		print("cond_probs_top:")
		print(cond_probs_top)
		}
	
	actual_probs_after_forward_exponentiation = cond_probs_top
	
	return(actual_probs_after_forward_exponentiation)
	} # END calc_prob_forward_onebranch_sparse <- function(relprobs_branch_bottom, branch_length, tmpQmat_in_REXPOKIT_coo_fmt, coo_n=coo_n, anorm=anorm, check_for_0_rows=TRUE, TRANSPOSE_because_forward=TRUE, printlevel=1)







#######################################################
# check_if_state_is_allowed
#######################################################
#' Check if a geographic range/state is allowed, given an areas-allowed matrix.
#' 
#' If the user has specified a matrix stating which areas are allowed to exist
#' (and thus have a species with a range including that area), this function checks if the 
#' input list of areas (as a 0-based vector of areas) in a single state/geographic range
#' is consistent with the areas-allowed matrix.
#' 
#' This function may be used by e.g. \code{\link[base]{apply}}.
#' 
#' @param state_0based_indexes The input state is a 0-based vector of area indices.
#' @param areas_allowed_mat A matrix (number of areas x number of areas) with 1s indicating allowed 
#' connections between areas, and 0s indicating disallowed connections.
#' @return \code{TRUE} or \code{FALSE}
#' @export
#' @seealso \code{\link[base]{apply}} \code{\link{check_if_state_is_allowed_by_adjacency}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
check_if_state_is_allowed <- function(state_0based_indexes, areas_allowed_mat)
	{
	# Check that areas-allowed matrix is square!!
	tmp_areas_allow_mat = areas_allowed_mat
	if (nrow(tmp_areas_allow_mat) != ncol(tmp_areas_allow_mat))
		{
		errortxt = paste0("\n\nSTOP ERROR in check_BioGeoBEARS_run():\n\nAreas-allowed matrices should be square, but yours ain't square:\n\nncol=", ncol(tmp_areas_allow_mat), "\nnrow=", nrow(tmp_areas_allow_mat), "\n\n")
		cat(errortxt)
		cat("Printing the offending areas-allowed matrix:\n\n")
		print(tmp_areas_allow_mat)
		cat("\n\n")
		stop(errortxt)
		} # END if (nrow(tmp_areas_allow_mat) != ncol(tmp_areas_allow_mat))

	# Check for NAs, let those through
	if ( (length(state_0based_indexes)==1) && (is.na(state_0based_indexes)) )
		{
		return(TRUE)
		} # END if (is.na(state_0based_indexes))
	
	state_1based_indexes = state_0based_indexes + 1
	tmp_numareas = length(state_1based_indexes)
	
	# Compare to areas_allowed_mat
	first_area_1based_index = state_1based_indexes[1]
	is_each_area_in_this_state_allowed = areas_allowed_mat[first_area_1based_index, state_1based_indexes]
	
	
# 	cat("\nstate_0based_indexes:")
# 	print(state_0based_indexes)
# 	
# 	cat("\ntmp_numareas =", tmp_numareas)
# 	
# 	cat("\nis_each_area_in_this_state_allowed =\n")
# 	print(is_each_area_in_this_state_allowed)
# 	
# 	cat("\nsum(is_each_area_in_this_state_allowed) =", sum(is_each_area_in_this_state_allowed))
# 	cat("\n")
	
	# Make sure is_each_area_in_this_state_allowed has positive length
	if ( (tmp_numareas == sum(is_each_area_in_this_state_allowed)) && (length(is_each_area_in_this_state_allowed) > 0) )
		{
		return(TRUE)	# This state/geographic range is allowed
		} else {
		return(FALSE)	# This state/geographic range is DISallowed
		} # END if (tmp_numareas == sum(is_each_area_in_this_state_allowed))
	
	return(stop("ERROR in check_if_state_is_allowed(): you shouldn't get here."))
	} # END check_if_state_is_allowed


#######################################################
# check_if_state_is_allowed_by_adjacency
#######################################################
#' Check if a geographic range/state is allowed, given an areas-adjacency matrix.
#' 
#' If the user has specified a matrix stating which areas are allowed to be connected
#' (and thus have a species with a range in both areas), this function checks if the 
#' input list of areas (as a 0-based vector of areas) in a single state/geographic range
#' is consistent with the areas-adjacency matrix.
#' 
#' This function may be used by e.g. \code{\link[base]{apply}}.
#' 
#' @param state_0based_indexes The input state is a 0-based vector of area indices.
#' @param areas_adjacency_mat A matrix (number of areas x number of areas) with 1s indicating allowed 
#' connections between areas, and 0s indicating disallowed connections.
#' @return \code{TRUE} or \code{FALSE}
#' @export
#' @seealso \code{\link[base]{apply}} \code{\link{check_if_state_is_allowed}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
check_if_state_is_allowed_by_adjacency <- function(state_0based_indexes, areas_adjacency_mat)
	{
	# Check that areas-allowed matrix is square!!
	tmp_areas_adjacency_mat = areas_adjacency_mat
	if (nrow(tmp_areas_adjacency_mat) != ncol(tmp_areas_adjacency_mat))
		{
		errortxt = paste0("\n\nSTOP ERROR in check_BioGeoBEARS_run():\n\nAreas-adjacency matrices should be square, but yours ain't square:\n\nncol=", ncol(tmp_areas_adjacency_mat), "\nnrow=", nrow(tmp_areas_adjacency_mat), "\n\n")
		cat(errortxt)
		cat("Printing the offending areas-adjacency matrix:\n\n")
		print(tmp_areas_adjacency_mat)
		cat("\n\n")
		stop(errortxt)
		} # END if (nrow(tmp_areas_adjacency_mat) != ncol(tmp_areas_adjacency_mat))

	# Check for NAs, let those through
	if ( (length(state_0based_indexes)==1) && (is.na(state_0based_indexes)) )
		{
		return(TRUE)
		} # END if (is.na(state_0based_indexes))
	
	
	# Compare to areas_adjacency_mat
	#state_1based_indexes = state_0based_indexes + 1
	#tmp_numareas = length(state_1based_indexes)
	#first_area_1based_index = state_1based_indexes[1]
	#is_each_area_in_this_state_allowed = areas_adjacency_mat[first_area_1based_index, state_1based_indexes]
	#if (tmp_numareas == sum(is_each_area_in_this_state_allowed))
	#	{
	#	return(TRUE)	# This state/geographic range is allowed
	#	} else {
	#	return(FALSE)	# This state/geographic range is DISallowed
	#	} # END if (tmp_numareas == sum(is_each_area_in_this_state_allowed))
	
	# Make sure this is a matrix, not a data.frame
	areas_adjacency_mat = as.matrix(areas_adjacency_mat)
	
	# Compare to areas_adjacency_mat
	areas_in_this_state_range_1based = state_0based_indexes + 1
	adjacent_if_1s = areas_adjacency_mat[areas_in_this_state_range_1based, areas_in_this_state_range_1based]

	if ( (length(adjacent_if_1s) == sum(adjacent_if_1s)) && (length(adjacent_if_1s) > 0) )
		{
		return(TRUE)	# This state/geographic range is allowed
		} else {
		return(FALSE)	# This state/geographic range is DISallowed
		} # END if ( length(adjacent_if_1s) == sum(adjacent_if_1s) )
	
	
# 	cat("\nstate_0based_indexes:")
# 	print(state_0based_indexes)
# 	
# 	cat("\ntmp_numareas =", tmp_numareas)
# 	
# 	cat("\nis_each_area_in_this_state_allowed =\n")
# 	print(is_each_area_in_this_state_allowed)
# 	
# 	cat("\nsum(is_each_area_in_this_state_allowed) =", sum(is_each_area_in_this_state_allowed))
# 	cat("\n")
	
	
	return(stop("ERROR in check_if_state_is_allowed_by_adjacency(): you shouldn't get here."))
	} # END check_if_state_is_allowed_by_adjacency








#######################################################
# Calc loglike speciation
#######################################################

#######################################################
# calc_loglike_sp_prebyte:
#######################################################
#' Calculate log-likelihood with a transition matrix and speciation events -- pre-byte-compiled
#'
#' This function is the pre-byte-compiled version of \code{\link{calc_loglike_sp}}.
#' 
#' Byte-compiling is supposed to speed up functions; this is an 
#' attempt to do this on the \code{\link[rexpokit]{rexpokit}} function \code{\link[rexpokit]{expokit_dgpadm_Qmat}}.  It
#' is also possible to byte-compile everything during package installation (via \code{ByteCompile: true} in the
#' DESCRIPTION file), which is implemented in \code{BioGeoBEARS}, so this may be redundant.
#'
#' \code{\link{calc_loglike_sp_prebyte}} gets byte-compiled into \code{\link{calc_loglike_sp}}.
#' 
#' See \url{http://dirk.eddelbuettel.com/blog/2011/04/12/} for discussion of the \code{\link[compiler]{compile}} package.
#' 
#' @param tip_condlikes_of_data_on_each_state A numeric matrix with rows representing tips, and columns representing states/geographic ranges.  The cells
#' give the likelihood of the observation data under the assumption that the tip has that state; typically this means that the known geographic range gets a 
#' '1' and all other states get a 0.
#' @param phy A phylogeny object.  The function converts it to pruningwise order.
#' @param Qmat A Q transition matrix representing the along-branch model for the evolution of geographic range, using parameters \emph{d} (dispersal/range expansion), 
#' \emph{e} (extinction/range contraction/local extirpation), and perhaps others (e.g. distance).  This matrix can be input in either dense or sparse (COO) format, 
#' as specified by \code{input_is_COO}.
#' @param spPmat Default is \code{NULL}; users should usually use \code{spPmat_inputs}.  \code{spPmat} is A numeric matrix representing the probability of each
#' ancestor range-->(Left range, Right range) transition at cladogenesis events.  There are 
#' different ways to represent this matrix.  In the simplest representation, this is just a rectangular matrix with numstates rows (representing the ancestral
#' states) and numstates^2 columns (representing all possible descendant pairs).  Use of this type of matrix is specified by \code{cppSpMethod=1}. It is calculated
#' from a textual speciation matrix (typically \code{spmat} in the code) via \code{\link{symbolic_to_relprob_matrix_sp}}. However, this matrix gets huge and
#' slow for large numbers of states/ranges.  \code{cppSpMethod=2} and \code{cppSpMethod=3} implement successively more efficient and faster 
#' representation and processing of this matrix in COO-like formats.  See \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOprobs}} for the \code{cppSpMethod=2} 
#' method, and \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}} for the \code{cppSpMethod=3} method (the fastest).
#' @param min_branchlength Nodes with branches below this branchlength will not be treated as cladogenesis events; instead, they will be treated as 
#' if an OTU had been sampled from an anagenetic lineage, i.e. as if you had a direct ancestor.  This is useful for putting fossils into the biogeography analysis,
#' when you have fossil species that range through time. (Note: the proper way to obtain such trees, given that most phylogenetic methods force all OTUs to be tips 
#' rather than direct ancestors, is another question subject to active research.  However, one method might be to just set a branch-length cutoff, and treat any
#' branches sufficiently small as direct ancestors.)
#' @param return_what What should be returned to the user? Options are "loglike" (the log-likelihood of the data under the tree, model, and model parameters), 
#' "nodelikes" (the scaled conditional likelihoods at the nodes), "rootprobs" (the relative probability of the geographic ranges/states at the root), or "all"
#' (all of the above in a list).  Typically the user will only want to return "loglike" while doing ML optimization, but then return "all" once the ML parameter
#' values have been found.
#' @param probs_of_states_at_root The prior probability of the states/geographic ranges at the root.  The default, \code{NULL}, effectively means an equal probability
#' for each state (this is also what \code{LAGRANGE} assumes; and running with NULL will reproduce exactly the \code{LAGRANGE} parameter inferences and
#' log-likelihood).
#' @param rootedge  Should the root edge be included in the calculation (i.e., calculate to the bottom of the root), if a root edge is present?  Default \code{FALSE}.
#' @param sparse Should sparse matrix exponentiation be performed?  This should be faster for very large matrices (> 100-200 states), however, the calculations 
#' appear to be less accurate.  The function will transform a dense matrix to COO format (see \code{\link[rexpokit]{mat2coo}}) if necessary according to 
#' the \code{input_is_COO} parameter.
#' @param printlevel If >= 1, various amounts of intermediate output will be printed to screen.  Note: Intermediate outputs from C++ and FORTRAN functions have been
#' commented out, to meet CRAN guidelines.
#' @param use_cpp Should the C++ routines from \code{\link[cladoRcpp]{cladoRcpp}} be used to speed up calculations?  Default \code{TRUE}.
#' @param input_is_COO Is the input Q matrix a sparse, COO-formatted matrix (\code{TRUE}) or a standard dense matrix (\code{FALSE}). Default \code{FALSE}.
#' @param spPmat_inputs A list of parameters so that \code{spPmat} (the speciation transition probability matrix) can be calculated on-the-fly, according
#' to the method in \code{cppSpMethod}.  See example.
#' @param cppSpMethod Three C++ methods from cladoRcpp for calculating and using the cladogenesis probability matrix.  1 is slowest but easiest to understand; 3 is fastest.
#' If \code{spPmat_inputs} is given, the program will generate the appropriate spPmat on-the-fly, and the user does not have to input the full \code{spPmat} manually.
#' @param cluster_already_open If the user wants to distribute the matrix exponentiation calculations from all the branches across a number of processors/nodes on 
#' a cluster, specify the cluster here.  E.g. \code{cluster_already_open = makeCluster(rep("localhost",num_cores_to_use), type = "SOCK")}.  Note: this will work on 
#' most platforms, including Macs running R from command line, but will NOT work on Macs running the R GUI \code{R.app}, because parallel processing functions like
#' \code{MakeCluster} from e.g. \code{library(parallel)} for some reason crash R.app.  The program runs a check for R.app and will just run on 1 node if found. 
#' @param calc_ancprobs Should ancestral state estimation be performed (adds an uppass at the end).
#' @param include_null_range Does the state space include the null range?
#' Default is \code{NULL} which means running on a single processor.
#' @param fixnode If the state at a particular node is going to be fixed (e.g. for ML marginal ancestral states), give the node number. Can now also be a vector of 
#' node numbers, to allow fixing of multiple nodes, but if this is done, \code{fixlikes} must be a matrix with a number of rows equal to \code{length(fixnode)}.
#' @param fixlikes The state likelihoods to be used at the fixed node.  I.e. 1 for the fixed state, and 0 for the others. If only one node is fixed, \code{fixlikes}
#' can be either a vector with length=numstates, or a matrix with 1 row and numstates columns.  For multiple fixed nodes, fixlikes must be a matrix with
#' \code{nrow(fixlikes)==length(fixnode)}.
#' @param stratified Default FALSE. If TRUE, you are running a stratified analysis, in which case uppass probs should be calculated elsewhere.
#' @param states_allowed_TF Default NULL. If user gives a vector of TRUE and FALSE values, these states will be set to 0 likelihood throughout the calculations.
#' @param m This is a vector of rate/weight multipliers for dispersal, conditional on the values of some
#' (non-biogeographical) trait. For example, one might hypothesize that flight/flightlessness effects
#' dispersal probability, and manually put a multiplier of 0.001 on the flightlessness state. Or, 
#' one might attempt to estimate this. The strategy used in cladoRcpp is to expand the default cladogenetic
#' rate matrix by length(m) times. I.e., if \emph{m} is not \code{NULL}, then loop through the values of \emph{m} and apply the 
#' multipliers to \emph{d} (and \emph{j}, and \emph{a}) events. Default is \code{NULL}.
#' @jts_matrix A numtraits x numtraits matrix containing the proportions for 
#' trait transitions during j events. E.g., for a sudden switch from 
#' trait 1 (flight) to trait 2 (flightlessness) during a jump event.
#' @BioGeoBEARS_model_object For printing the parameters in the event that the likelihood calculation produces NaNs. 
#' Ignored if NULL (default)
#' @on_NaN_error Default is NULL, which prints a long stop error if NaNs are hit during the likelihood calculation 
#' (probably due to extreme parameter values and underflow problems). If on_NaN_error is numeric, the value is the 
#' log-likelihood that is returned (e.g. a log-likelihood of -10e100 should be so incredibly low that it 
#' drives any optimizer away from that solution).
#' @return Return whatever is specified by \code{return_what}.
#' @export
#' @seealso \code{\link{calc_loglike_sp}}, \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp}}, \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOprobs}}, 
#' \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}}, \code{\link[rexpokit]{mat2coo}}, 
#' \code{\link{rcpp_calc_anclikes_sp_COOweights_faster}}
#' @bibliography /Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#'	 @cite Landis_Matzke_etal_2013_BayArea
#' @note Go BEARS!
#' @note (COO = Coordinate list format for a matrix, see \url{http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29}
#' @author Nicholas Matzke \email{matzke@@berkeley.edu}
#' @examples
#' testval=1
#'
calc_loglike_sp_prebyte <- function(tip_condlikes_of_data_on_each_state, phy, Qmat, spPmat=NULL, min_branchlength=0.000001, return_what="loglike", probs_of_states_at_root=NULL, rootedge=FALSE, sparse=FALSE, printlevel=1, use_cpp=TRUE, input_is_COO=FALSE, spPmat_inputs=NULL, cppSpMethod=3, cluster_already_open=NULL, calc_ancprobs=FALSE, include_null_range=TRUE, fixnode=NULL, fixlikes=NULL, stratified=FALSE, states_allowed_TF=NULL, m=NULL, jts_matrix=NULL, BioGeoBEARS_model_object=NULL, on_NaN_error=NULL)
	{
	defaults='
	# Phylogeny
	newick_str = "(((Humans:6.0, Chimps:6.0):1.0, Gorillas:7.0):5.0, Orangs:12.0):1.0;"
	phy = read.tree(file="", text=newick_str)

	# Areas: A=Africa, B=Not Africa
	areas = c("A","B")
	states_list = areas_list_to_states_list_old(areas, include_null_range=FALSE)
	dedf = make_relprob_matrix_de(states_list)
	Qmat = symbolic_to_Q_matrix(dedf, d=0.1, e=0.1)
	
	numnodes = length(phy$tip.label) + phy$Nnode
	tip_condlikes_of_data_on_each_state = matrix(data=0, nrow=length(phy$tip.label), ncol=length(states_list))
	# The columns are tip likelihoods of: A, B, AB
	tip_condlikes_of_data_on_each_state[1:3,1] = 1
	tip_condlikes_of_data_on_each_state[4,2] = 1
	tip_condlikes_of_data_on_each_state
	

	spPmat=NULL
	
	
	min_branchlength=0.000001
	return_what = "all"
	probs_of_states_at_root=NULL
	rootedge=FALSE
	cppSpMethod=3
	cluster_already_open=NULL
	printlevel=1
	input_is_COO=FALSE
	
	stratified=FALSE
	use_cpp=TRUE
	#sparse=TRUE
	sparse=FALSE
	
	#include_null_range=TRUE
	include_null_range=FALSE
	m=NULL
	
	states_allowed_TF=NULL
	probs_of_states_at_root = NULL
	' # end defaults
	
	
	# Store sparse in force_sparse
	force_sparse = sparse
	
	#print("printlevel")	
	#print(printlevel)
	
	#print("on_NaN_error:")
	#print(on_NaN_error)
	#print("Qmat#1")
	#print(Qmat)

	#######################################################
	# Error check on fixnode / fixlikes
	#######################################################
	if (!is.null(fixnode))
		{
		if (( is.null(dim(fixlikes)) == TRUE) && (length(fixnode)==1))
			{
			pass_fixlikes = TRUE
			} else {
			if ( (dim(fixlikes)[1]) == length(fixnode) )
				{
				pass_fixlikes = TRUE
				
				# Another error check: Multiple nodes in 'fixnode' MUST be sorted in increasing order
				if ( (all(c(order(fixnode) == 1:length(fixnode)))) == TRUE)
					{
					pass_fixlikes = TRUE
					} else {
					pass_fixlikes = FALSE
					error_msg = "ERROR in calc_loglike_sp(): \n             Multiple nodes in 'fixnode' MUST be sorted in increasing order.\n"
					cat(error_msg)
					stop(error_msg)
					} # end if order check

				} else {
				pass_fixlikes = FALSE
				error_msg = "ERROR: Either:\n             (1) fixnode must be a single node number, and fixlikes must be a vector, or\n             (2) fixlikes like must be a matrix with the # of rows equal to length(fixnode).\n"
				cat(error_msg)
				stop(error_msg)
				} # end 2nd if()
			} # end 1st if()
		}
	
	
	#######################################################
	# Calculating ancestral probs via an uppass requires the
	# downpass splitlikes be saved.
	#######################################################
	if (calc_ancprobs == TRUE)
		{
		# Also, cppSpMethod must be 3 (I'm not going to program the other ones!!)
		if (cppSpMethod != 3)
			{
			stop("ERROR: You must have cppSpMethod=3 if calc_ancprobs=TRUE")
			}
		}
	
#	print("m:")
#	print(m)
	
	#######################################################
	# ERROR CHECK
	#######################################################
	if (use_cpp == TRUE && !is.null(spPmat_inputs))
		{
		numstates_in_tip_condlikes_of_data_on_each_state = ncol(tip_condlikes_of_data_on_each_state)
		
		# Number of states in the spPmat, *including* the null range if desired
		if (include_null_range == TRUE)
			{
			numstates_in_spPmat_inputs_states_indices = 1 + length(spPmat_inputs$l)
			} else {
			numstates_in_spPmat_inputs_states_indices = 0 + length(spPmat_inputs$l)
			} # END if (include_null_range == TRUE)
		
		if (is.null(m) == FALSE)
			{
			numstates_in_spPmat_inputs_states_indices = numstates_in_spPmat_inputs_states_indices * length(m)
			}
		
		if (input_is_COO == TRUE)
			{
			#print("Qmat:")
			#print(Qmat)
			numstates_in_Qmat = max( max(Qmat[,1]), max(Qmat[,2]) )
			} else {
			numstates_in_Qmat = ncol(Qmat)
			} # END if (input_is_COO == TRUE)
		
		if ( all(numstates_in_tip_condlikes_of_data_on_each_state==numstates_in_spPmat_inputs_states_indices, numstates_in_tip_condlikes_of_data_on_each_state==numstates_in_Qmat ) == TRUE )
			{
			# Continue
			all_inputs_correct_size = TRUE
			} else {
			
			
			if (force_sparse == FALSE)
				{
				# Stop if not everything is equal
				stop_function <- function(numstates_in_tip_condlikes_of_data_on_each_state, numstates_in_spPmat_inputs_states_indices, numstates_in_Qmat)
					{
					cat("\n\nSTOP ERROR in calc_loglike_sp_prebyte(): Some inputs have incorrect size -- \n")
					cat("numstates_in_tip_condlikes_of_data_on_each_state:	", numstates_in_tip_condlikes_of_data_on_each_state, "\n")
					cat("numstates_in_spPmat_inputs_states_indices+1:	", numstates_in_spPmat_inputs_states_indices, "\n")
					cat("numstates_in_Qmat:								", numstates_in_Qmat, "\n")
					cat("\n")
					cat("This means either\n(a) you have a conflict between the length of states_list\nand the dimensions of the tip likelihoods, Qmat, etc.; or\n(b) it means 'd' and 'e' were so close to zero that rcpp_states_list_to_DEmat()'s\ndecided they were effectively 0; default min_precision=1e-16.\n(c) include_null_range setting is wrong, you have include_null_range=", include_null_range, "\n...Probably...(a) or (c).", sep="")
					cat("\n")
					} # END stop_function

				stop(stop_function(numstates_in_tip_condlikes_of_data_on_each_state, numstates_in_spPmat_inputs_states_indices, numstates_in_Qmat))
				} else {
				# Warning if not everything is equal
				warning_function <- function(numstates_in_tip_condlikes_of_data_on_each_state, numstates_in_spPmat_inputs_states_indices, numstates_in_Qmat)
					{
					cat("\n\nWarning in calc_loglike_sp_prebyte(): this would be a stop error under dense matrix exponentiation, but you have force_sparse=TRUE, so probably this message just means you have a zero rate in the bottom row of Q so it thinks the COO-formatted matrix is the wrong size. Probably ignore. Dense matrix message follows.\n\nSome inputs have incorrect size -- \n")
					cat("numstates_in_tip_condlikes_of_data_on_each_state:	", numstates_in_tip_condlikes_of_data_on_each_state, "\n")
					cat("numstates_in_spPmat_inputs_states_indices+1:	", numstates_in_spPmat_inputs_states_indices, "\n")
					cat("numstates_in_Qmat:								", numstates_in_Qmat, "\n")
					cat("\n")
					cat("This means either\n(a) you have a conflict between the length of states_list\nand the dimensions of the tip likelihoods, Qmat, etc.; or\n(b) it means 'd' and 'e' were so close to zero that rcpp_states_list_to_DEmat()'s\ndecided they were effectively 0; default min_precision=1e-16.\n(c) include_null_range setting is wrong, you have include_null_range=", include_null_range, "\n...Probably...(a) or (c).", sep="")
					cat("\n")
					} # END warning_function
				warning(warning_function(numstates_in_tip_condlikes_of_data_on_each_state, numstates_in_spPmat_inputs_states_indices, numstates_in_Qmat))
				} # END if (force_sparse == FALSE)
			}
		} else {
		numstates_in_tip_condlikes_of_data_on_each_state = ncol(tip_condlikes_of_data_on_each_state)
		
		# Number of states in the spPmat, *including* the null range if desired
		if (include_null_range == TRUE)
			{
			numstates_in_spPmat_ie_nrows = 1 + nrow(spPmat)
			} else {
			numstates_in_spPmat_ie_nrows = 0 + nrow(spPmat)
			} # END if (include_null_range == TRUE)

		if (is.null(m) == FALSE)
			{
			numstates_in_spPmat_ie_nrows = numstates_in_spPmat_ie_nrows * length(m)
			}

		
		if (input_is_COO == TRUE)
			{
			numstates_in_Qmat = max( max(Qmat[,1]), max(Qmat[,2]) )
			} else {
			numstates_in_Qmat = ncol(Qmat)
			} # END if (input_is_COO == TRUE)
		
		if ( all(numstates_in_tip_condlikes_of_data_on_each_state==numstates_in_spPmat_ie_nrows, numstates_in_tip_condlikes_of_data_on_each_state==numstates_in_Qmat ) == TRUE )
			{
			# Continue
			all_inputs_correct_size = TRUE
			} else {

			if (force_sparse == FALSE)
				{
				# Stop if not everything is equal
				stop_function2 <- function(numstates_in_tip_condlikes_of_data_on_each_state, numstates_in_spPmat_ie_nrows, numstates_in_Qmat)
					{
					cat("\n\nERROR: Some inputs have incorrect size -- \n")
					cat("numstates_in_tip_condlikes_of_data_on_each_state:	", numstates_in_tip_condlikes_of_data_on_each_state, "\n")
					cat("numstates_in_spPmat_ie_nrows+1:				", numstates_in_spPmat_ie_nrows, "\n")
					cat("numstates_in_Qmat:								", numstates_in_Qmat, "\n")
					cat("\n")
					cat("This means either\n(a) you have a conflict between the length of states_list\nand the dimensions of the tip likelihoods, Qmat, etc.; or\n(b) it means 'd' and 'e' were so close to zero that rcpp_states_list_to_DEmat()'s\ndecided they were effectively 0; default min_precision=1e-16.\n(c) include_null_range setting is wrong, you have include_null_range=", include_null_range, "\n...Probably...(a) or (c).", sep="")
					cat("\n")
					} # END stop_function2
				stop(stop_function2(numstates_in_tip_condlikes_of_data_on_each_state, numstates_in_spPmat_ie_nrows, numstates_in_Qmat))
				} else {
				warning_function2 <- function(numstates_in_tip_condlikes_of_data_on_each_state, numstates_in_spPmat_ie_nrows, numstates_in_Qmat)
					{
					cat("\n\nWarning only - this would be a stop error under dense matrix exponentiation, but you have force_sparse=TRUE, so probably this message just means you have a zero rate in the bottom row of Q so it thinks the COO-formatted matrix is the wrong size. Probably ignore. Dense matrix message follows.\n\nSome inputs have incorrect size -- \n")
					cat("numstates_in_tip_condlikes_of_data_on_each_state:	", numstates_in_tip_condlikes_of_data_on_each_state, "\n")
					cat("numstates_in_spPmat_ie_nrows+1:				", numstates_in_spPmat_ie_nrows, "\n")
					cat("numstates_in_Qmat:								", numstates_in_Qmat, "\n")
					cat("\n")
					cat("This means either\n(a) you have a conflict between the length of states_list\nand the dimensions of the tip likelihoods, Qmat, etc.; or\n(b) it means 'd' and 'e' were so close to zero that rcpp_states_list_to_DEmat()'s\ndecided they were effectively 0; default min_precision=1e-16.\n(c) include_null_range setting is wrong, you have include_null_range=", include_null_range, "\n...Probably...(a) or (c).", sep="")
					cat("\n")
					} # END stop_function2
				warning(warning_function2(numstates_in_tip_condlikes_of_data_on_each_state, numstates_in_spPmat_ie_nrows, numstates_in_Qmat))
				
				} # END if (force_sparse == FALSE)
			}
		} # END if (use_cpp == TRUE && !is.null(spPmat_inputs))
		# END error check
	
	
	
	
	
	
	# Fix "l" (which is states_indices, i.e. a list of lists of state_indices)
	# so that there is no "null geographic range", i.e. no "_", no "-", no "NA"
	if ( is.null(spPmat_inputs)==FALSE )
		{
		spPmat_inputs$l[spPmat_inputs$l == c("_")] = NULL
		spPmat_inputs$l[spPmat_inputs$l == c("-")] = NULL
		spPmat_inputs$l[spPmat_inputs$l == c("-1")] = NULL
		#spPmat_inputs$l[spPmat_inputs$l == c(-1)] = NULL
		}
	
	
	
	# Calculate likelihoods down tree
	#numstates = nrow(Qmat)
	numstates = ncol(tip_condlikes_of_data_on_each_state)

	
	#######################################################
	# Check if the phylogeny is actually just a number (i.e., a branch length)
	#######################################################
	# Hmm, this seems harder...


	#######################################################
	# The rest assumes that you've got a phylo3 object, in default order
	#######################################################
	
	
	
	
	edgelengths = phy$edge.length
	num_branches_below_min = sum(edgelengths < min_branchlength)

	if ( (printlevel >= 1) && (num_branches_below_min > 0) )
		{
		cat("\n")
		cat("Running calc_loglike_sp():\n")
		cat("This run of calc_loglike_sp() has a min_branchlength of: ", min_branchlength, "\n", sep="")
		cat("Branches shorter than this will be assumed to be connected to the tree with\n")
		cat("sympatric events (i.e., members of fossil lineages on ~0 length branches.)\n")
		cat("This tree has ", num_branches_below_min, " branches < ", min_branchlength, ".\n", sep="")
		cat("\n")
		}	

	
	num_internal_nodes = phy$Nnode
	numtips = length(phy$tip.label)
	
	# likelihoods are computed at all nodes
	# make a list to store the 
	numnodes = numtips + num_internal_nodes
	computed_likelihoods_at_each_node = numeric(length=numnodes)
	
	
	
	
	# Reorder the edge matrix into pruningwise order
	# This is CRUCIAL!!
	phy2 <- reorder(phy, "pruningwise")
	
	
	tipnums <- 1:numtips

	# Put in the sums of the probabilities of the states at each tip
	# (This only works if the tip data are 0000100000, etc...)
	#computed_likelihoods_at_each_node[tipnums] = rowSums(tip_condlikes_of_data_on_each_state)
	
	# 2014-05-07_NJM: change to divide by the sum
	computed_likelihoods_at_each_node[tipnums] = rowSums(tip_condlikes_of_data_on_each_state) / rowSums(tip_condlikes_of_data_on_each_state)
	
	#######################################################
	# Initialize matrices for downpass and uppass storage of state probabilities
	#######################################################
	
	# This is all you need for a standard likelihood calculation
	# relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS = rel probs AT A NODE
	# BE SURE TO DISTINGUISH STATE PROBABILITIES (VS LIKELIHOODS, WHICH ARE NOT NORMALIZED)
	condlikes_of_each_state <- matrix(data=0, nrow=numnodes, ncol=numstates)
	relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS <- matrix(data=0, nrow=numnodes, ncol=numstates)
	
	# In THEORY, this should be OK to divide, but BECAUSE I pass relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS to 
	# the matrix exponentiation calculations, these need to be likelihoods at the tips;
	# SO DON'T DIVIDE, if you do it normalizes the probabilities at the tips, and this screws up things that are being passed 
	# likelihoods that are not 1/0000, e.g. stratification analyses, and causes M0 and M0 strat to disagree on Psychotria (e.g.)
	#
	# At some point we should more rigorously separate relprobs and condlikes, but don't screw with it now!
	# 
	# For e.g. uppass calculations of tip probs, just do a hack later.
	# 
	relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[tipnums, ] = tip_condlikes_of_data_on_each_state #/ rowSums(tip_condlikes_of_data_on_each_state)
	
	# BUT, DO NOT DIVIDE THIS BY rowSums(tip_condlikes_of_data_on_each_state), it forces normalization which prevents e.g.
	# passing down likelihoods during stratification
	condlikes_of_each_state[tipnums, ] = tip_condlikes_of_data_on_each_state
	#relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[tipnums, ] <- 1
	
	
	# But, if you want to do ancestral states properly, and get marginal
	# reconstructions, you've gotta store:
	# 1. Probability of the states at the bottom of each branch on downpass
	# 2. Probability of the states at the bottom of each branch on UP-pass
	# 3. Probability of the states at the top of each branch on UP-pass
	# Combine 0-3 to get the ML marginal probability of states at
	# 4. The bottom of each branch
	# 5. The top of each branch
	# Combine 4 & 5 to get:
	# 6. MAYBE The ML marginal probability of each split scenario at each node - i.e. #5(top) * #4(left bot) * #4(right bot)
	# 
	if (calc_ancprobs == TRUE)
		{
		# Every node (except maybe the root) has a branch below it, and there is also a 
		# relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS at the bottom of this branch
		relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS <- matrix(data=NA, nrow=numnodes, ncol=numstates)
		relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS <- matrix(data=NA, nrow=numnodes, ncol=numstates)
		relative_probs_of_each_state_at_branch_top_AT_node_UPPASS <- matrix(data=NA, nrow=numnodes, ncol=numstates)

		ML_marginal_prob_each_state_at_branch_bottom_below_node <- matrix(data=NA, nrow=numnodes, ncol=numstates)
		ML_marginal_prob_each_state_at_branch_top_AT_node <- matrix(data=NA, nrow=numnodes, ncol=numstates)
		ML_marginal_prob_each_split_at_branch_top_AT_node = list()
		}



	
	# If you want to use standard dense matrix exponentiation, rexpokit's dgpadm is best,
	# and you can run it through mapply on the branch lengths of all branches at once
	if (sparse==FALSE)
		{
		# Get probmats for each branch, put into a big array
		# Creates empty array to store results
		#print("Qmat")
		#print(dput(Qmat))

		
		independent_likelihoods_on_each_branch = calc_independent_likelihoods_on_each_branch(phy2, Qmat, cluster_already_open, Qmat_is_sparse=FALSE)
		} else {
		# Sparse matrices
		# Here, if you want to use sparse matrix exponentiation, you are NOT
		# going to store the entire Pmat for each branch; you are going to directly
		# calculate the output probabilities w (w), via w=exp(Qt)v via myDMEXPV.
		#
		# This is most efficiently done if you transpose and COO-ify the matrix here,
		# ahead of time.
		if (input_is_COO==FALSE)
			{
			original_Qmat = Qmat
			
			# number of states in the original matrix
			coo_n = numstates
			anorm = as.numeric(norm(original_Qmat, type="O"))
			matvec = original_Qmat
			
			# DO NOT TRANSPOSE; we want to go BACKWARDS in time, NOT FORWARDS!
			#tmatvec = base::t(matvec)
			tmatvec = matvec
			tmpQmat_in_REXPOKIT_coo_fmt = mat2coo(tmatvec)
			} else {
			# The input Qmat is already in COO format (untransposed, as we are doing 
			# likelihoods backwards in time)
			
			# CHECK THAT ITS IN COO FORMAT
			if (class(Qmat) != "data.frame")
				{
				txt = paste0("ERROR: calc_loglike_sp is attempting to use a sparse COO-formated Q matrix, but you provided a Qmat not in data.frame form")
				cat("\n\n")
				cat(txt)
				cat("\n\n")
				print("class(Qmat)")
				print(class(Qmat) )
				print("Qmat")
				print(Qmat)
				stop(txt)
				}
			
			
			if ( (ncol(Qmat) != 3) )
				{
				txt = paste0("ERROR: calc_loglike_sp is attempting to use a sparse COO-formated Q matrix, but you provided a Qmat that does't have 3 columns")
				cat("\n\n")
				cat(txt)
				cat("\n\n")
				print("class(Qmat)")
				print(class(Qmat) )
				print("Qmat")
				print(Qmat)
				stop(txt)
				}
			coo_n = numstates
			anorm = 1
			tmpQmat_in_REXPOKIT_coo_fmt = Qmat
			#tmpQmat_in_REXPOKIT_coo_fmt = cooQmat
			} # END if (input_is_COO==FALSE)
		
		# Make a CRS-formatted matrix, for kexpmv
		# DO THE TRANSPOSE HERE, trait+geog matrices assembled transposed
		tmpQmat_in_kexpmv_crs_fmt = coo2crs(
			ia=tmpQmat_in_REXPOKIT_coo_fmt[,"ia"], 
			ja=tmpQmat_in_REXPOKIT_coo_fmt[,"ja"], 
			a =tmpQmat_in_REXPOKIT_coo_fmt[,"a"],
			n=numstates, transpose_needed=FALSE)

# 		print("numstates")
# 		print(numstates)
		} # END if (sparse==FALSE)


	# DEFINE DOWNPASS THROUGH THE BRANCHES	
	i = 1
	edges_to_visit = seq(from=1, by=2, length.out=num_internal_nodes)
	
	
	# Calculate the rowsums of the speciation matrix ONCE
	if ( use_cpp == TRUE )
		{
		if ( is.null(spPmat_inputs)==FALSE )
			{
			# Calculate the rowsums (for input into rcpp_calc_anclikes_sp()
			# Actually, just do this ONCE

			# (above) Fix "l" (which is states_indices, i.e. a list of lists of state_indices)
			# so that there is no "null geographic range", i.e. no "_", no "-", no "NA"
			l = spPmat_inputs$l		# states_indices
			s = spPmat_inputs$s
			v = spPmat_inputs$v
			j = spPmat_inputs$j
			y = spPmat_inputs$y
			
			# dmat means combined dispersal multipliers matrix, NOT dmat_times_d
			dmat = spPmat_inputs$dmat
			
			# Take the max of the indices of the possible areas, and add 1
			# numareas = max(unlist(spPmat_inputs$l), na.rm=TRUE) + 1 # old, bogus
			numareas = max(sapply(X=spPmat_inputs$l, FUN=length), na.rm=TRUE) + 0
			
			maxent01s_param = spPmat_inputs$maxent01s_param
			maxent01v_param = spPmat_inputs$maxent01v_param
			maxent01j_param = spPmat_inputs$maxent01j_param
			maxent01y_param = spPmat_inputs$maxent01y_param
			
			maxent01s = relative_probabilities_of_subsets(max_numareas=numareas, maxent_constraint_01=maxent01s_param, NA_val=0)
			maxent01v = relative_probabilities_of_vicariants(max_numareas=numareas, maxent_constraint_01v=maxent01v_param, NA_val=0)
			maxent01j = relative_probabilities_of_subsets(max_numareas=numareas, maxent_constraint_01=maxent01j_param, NA_val=0)
			maxent01y = relative_probabilities_of_subsets(max_numareas=numareas, maxent_constraint_01=maxent01y_param, NA_val=0)

			# You really need a list of sizes here:
			
			# Matrix of probs for each ancsize
			maxprob_as_function_of_ancsize_and_decsize = mapply(FUN=max, maxent01s, maxent01v, maxent01j, maxent01y, MoreArgs=list(na.rm=TRUE))
			maxprob_as_function_of_ancsize_and_decsize = matrix(data=maxprob_as_function_of_ancsize_and_decsize, nrow=nrow(maxent01s), ncol=ncol(maxent01s))
			maxprob_as_function_of_ancsize_and_decsize[maxprob_as_function_of_ancsize_and_decsize > 0] = 1
			maxprob_as_function_of_ancsize_and_decsize[maxprob_as_function_of_ancsize_and_decsize <= 0] = 0
			
			# Now, go through, and make a list of the max minsize for each decsize
			max_minsize_as_function_of_ancsize = apply(X=maxprob_as_function_of_ancsize_and_decsize, MARGIN=1, FUN=maxsize)

			# -1 for null range
			if (include_null_range == TRUE)
				{
				state_space_size_Qmat_to_cladoMat = -1
				} else {
				state_space_size_Qmat_to_cladoMat = 0
				}
			
			tmpca_1 = rep(1, (ncol(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS) + state_space_size_Qmat_to_cladoMat))
			tmpcb_1 = rep(1, (ncol(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS) + state_space_size_Qmat_to_cladoMat))

			# Print the matrix to screen from C++
			printmat = FALSE
			if (printlevel >= 2)
				{
				printmat = TRUE
				}
			
			# Calculate the rowSums of the speciation matrix, for future reference (i.e., setting all input likelihoods to 1)
			# Only need be done once.
			#
			# But, actually, what would make this all REALLY efficient would be to just calculate the conditional
			# probability matrix ONCE, storing it in a COO-like format.  The Rsp_rowsums would be easily derived
			# from that, and we wouldn't have to calculate the speciation model Nnodes times independently.
			#Rsp_rowsums = rcpp_calc_anclikes_sp_rowsums(Rcpp_leftprobs=tmpca_1, Rcpp_rightprobs=tmpcb_1, l=l, s=s, v=v, j=j, y=y, dmat=dmat, maxent01s=maxent01s, maxent01v=maxent01v, maxent01j=maxent01j, maxent01y=maxent01y, max_minsize_as_function_of_ancsize=max_minsize_as_function_of_ancsize, printmat=printmat)
			
			# Get the speciation matrix conditional probabilities in a COO-like format
			# [[1]] = inums = indexes of left descendant state in speciation matrix, by ancestral rowsnums 1-15
			# [[2]] = jnums = indexes of right descendant state in speciation matrix, by ancestral rowsnums 1-15
			# [[3]] = probs = probvals of this left|right combination in speciation matrix, by ancestral rowsnums 1-15
			if (cppSpMethod == 2)
				{
				# Error check
				if (is.null(m) == FALSE)
					{
					txt = "STOP ERROR in calc_loglike_sp_prebyte(): if m is not NULL, only option cppSpMethod==3 can be used. You have cppSpMethod==2."
					cat("\n\n")
					cat(txt)
					cat("\n\n")
					stop(txt)
					} # END if (is.null(m) == FALSE)
				
				COO_probs_list = rcpp_calc_anclikes_sp_COOprobs(Rcpp_leftprobs=tmpca_1, Rcpp_rightprobs=tmpcb_1, l=l, s=s, v=v, j=j, y=y, dmat=dmat, maxent01s=maxent01s, maxent01v=maxent01v, maxent01j=maxent01j, maxent01y=maxent01y, max_minsize_as_function_of_ancsize=max_minsize_as_function_of_ancsize, printmat=printmat)

				# Sum the probabilities (list [[3]]) in each row ([[3]][[1]] list of probs
				# through [[3]][[15]] list of probs)
				Rsp_rowsums = sapply(X=COO_probs_list[[3]], FUN=sum)
				} # END if (cppSpMethod == 2)
	
			if (cppSpMethod == 3)
				{
				if (printlevel >= 2)
					{
					params_to_print = c("tmpca_1", "tmpcb_1", "l", "s", "v", "j", "y", "dmat", "maxent01s", "maxent01v", "maxent01j", "maxent01y", "max_minsize_as_function_of_ancsize")

					for (tmppval in params_to_print)
						{
						cmdstr = paste(	"cat('", tmppval, "', ':\n', sep='')", sep="")
						eval(parse(text=cmdstr))
						
						# Get the value
						cmdstr = paste("tmppval = ", tmppval, sep="")
						eval(parse(text=cmdstr))
						
						# Print it
						print(tmppval)
						}
					} # END if (printlevel >= 2)
				
				#print(m)
				
				COO_weights_columnar = rcpp_calc_anclikes_sp_COOweights_faster(Rcpp_leftprobs=tmpca_1, Rcpp_rightprobs=tmpcb_1, l=l, s=s, v=v, j=j, y=y, dmat=dmat, maxent01s=maxent01s, maxent01v=maxent01v, maxent01j=maxent01j, maxent01y=maxent01y, max_minsize_as_function_of_ancsize=max_minsize_as_function_of_ancsize, printmat=printmat, m=m, m_null_range=include_null_range, jts_matrix=jts_matrix)
				
				# combine with C++ function
				# This causes an error with spPmat=NULL; spPmat_inputs=NULL; use_cpp=TRUE; sparse=FALSE
				# i.e. gives 16 states with a 0 on the end, rather than 15 states
				#Rsp_rowsums = rcpp_calc_rowsums_for_COOweights_columnar(COO_weights_columnar=COO_weights_columnar, numstates=numstates)
				
				# This gives 15 states
				Rsp_rowsums = rcpp_calc_rowsums_for_COOweights_columnar(COO_weights_columnar=COO_weights_columnar)
				} # END if (cppSpMethod == 3)
			} # END if ( is.null(spPmat_inputs)==FALSE )
		} # END if ( use_cpp == TRUE )


	#######################################################
	#######################################################
	# THIS IS A DOWNPASS FROM THE TIPS TO THE ROOT
	#######################################################
	#######################################################
	for (i in edges_to_visit)
		{
		# First edge visited is i
		#print(i)
		
		# Its sister is j 
		j <- i + 1
		#print(j)

		# Get the node numbers at the tips of these two edges		
		left_desc_nodenum <- phy2$edge[i, 2]
		right_desc_nodenum <- phy2$edge[j, 2]
		
# 		if (left_desc_nodenum == 53)
# 			{
# 			print("left_desc_nodenum")
# 			print("right_desc_nodenum")
# 			print(left_desc_nodenum)
# 			print(right_desc_nodenum)
# 			}
		
		# And for the ancestor edge (i or j shouldn't matter, should produce the same result!!!)
		anc <- phy2$edge[i, 1]
		
		txt = paste("anc:", anc, " left:", left_desc_nodenum, " right:", right_desc_nodenum, sep="")
		#print(txt)
		
		# Is sparse is FALSE, input the pre-calculated likelihoods;
		# If sparse is TRUE, dynamically calculate using expokit_dgexpv_Qmat
		if (sparse==FALSE)
			{
			
			#print("Checking Q matrix")
			#cat("p=", p, ", rate=", rate, "\n", sep=" ")
			#print(Q)
			#print(phy$edge.length[i])
			#condlikes_Left <- matexpo(Qmat * phy2$edge.length[i]) %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,]
			#condlikes_Right <- matexpo(Qmat * phy2$edge.length[j]) %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum,]


			if (printlevel >= 2)
				{
				print("dense matrix exponentiation")
				}
			
			
			if (is.null(cluster_already_open))
				{
				# Conditional likelihoods of states at the bottom of left branch
				condlikes_Left = independent_likelihoods_on_each_branch[,,i] %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,]
							
				# Conditional likelihoods of states at the bottom of right branch
				condlikes_Right = independent_likelihoods_on_each_branch[,,j] %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum,]
				} else {
				
				#cat("dim(independent_likelihoods_on_each_branch[[i]]):\n", sep="")
				#cat(dim(independent_likelihoods_on_each_branch))
				#cat("\n\n")
				
				#cat("length(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,]):\n", sep="")
				#cat(length(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,]))
				#cat("\n\n")
				
				condlikes_Left = independent_likelihoods_on_each_branch[[i]] %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,]
				
				# Conditional likelihoods of states at the bottom of right branch
				condlikes_Right = independent_likelihoods_on_each_branch[[j]] %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum,]

				}
			
# 			if (left_desc_nodenum == 53)
# 				{
# 				print("saving relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,]")
# 				#save(Qmat, file="Qmat.Rdata")
# 				# 
# 				print("dput(Qmat):")
# 				#print(dput(Qmat))
# 				save(Qmat, file="Qmat.Rdata")
# 				print(Qmat[1:15,1:15])
# 
# 				
# 				print("condlikes_Left")
# 				print(condlikes_Left)
# 				print("condlikes_Right")
# 				print(condlikes_Right)
# 				
# 				
# 				tip_probs_left = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,]
# 				tip_probs_right = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum,]
# 				save(tip_probs_left, file="left_tipprobs.Rdata")
# 				save(tip_probs_right, file="right_tipprobs.Rdata")
# 				condlikes_Left_row = condlikes_Left
# 				condlikes_Right_row = condlikes_Right
# 				save(condlikes_Left_row, file="condlikes_Left.Rdata")
# 				save(condlikes_Right_row, file="condlikes_Right.Rdata")
# 
# 				save(tip_probs_left, file="tip_probs_left.Rdata")
# 				save(tip_probs_right, file="tip_probs_right.Rdata")
# 
# 				}
# 			
			
			
			# Error check
			if (any(is.nan(condlikes_Left)))
				{
				if (printlevel > 0)
					{
					cat("\n\nindependent likelihood %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS returned NaNs:\n\n")

					cat("\n\nPrinting parameter values being attempted when this occurred:\n\n")
					if (is.null(BioGeoBEARS_model_object))
						{
						print("BioGeoBEARS_model_object was not passed to calc_loglike_sp() for printing.")
						} else {
						print(BioGeoBEARS_model_object)
						} # END if (is.null(BioGeoBEARS_model_object))
				
					cat("i=", i, "\n")
					cat("phy2$edge.length[i]=", phy2$edge.length[i], "\n")
					cat("\n\nprint(condlikes_Left):\n\n")
					print(condlikes_Left)

					cat("\n\nprint(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS):\n\n")
					print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS)

					cat("\n\nprint(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,]):\n\n")
					print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,])
					#print(coo_n)
					#print(anorm)
					}
				
				if (is.null(on_NaN_error))
					{
					stop("\n\nStopping on error in dense exponentiation downpass (left branch): NaNs produced in likelihood calculation. This may mean your transition matrix disallows necessary transitions.  E.g., if your ranges are 'A' and 'B', and your model is DEC, then allowing range 'AB' as a possible state is required, so that you can get from 'A' to 'B' via 'AB' as the intermediate. Alternatively, NaNs can be produced sometimes if your Maximum Likelihood (ML) search proposes weird parameter values (such as a negative rate or weight) or a parameter so small that required transitions have a probability that machine precision rounds to zero or negative.  Sometimes this seems to occur because optim, optimx, etc. propose parameters slightly outside the user-specified upper and lower (min/max) boundaries for some reason. One solution is often to narrow the min/max limits. \n\nAnother solution: To have this error report an extremely low log-likelihood,, set BioGeoBEARS_run_object$on_NaN_error to something like -1e50.\n\n")
					}
				
				if ( (is.numeric(on_NaN_error)) && (return_what == "loglike") )
					{
					if (printlevel > 0)
						{
						warning(paste0("\n\nWarning on error in dense exponentiation downpass (left branch): NaNs produced in likelihood calculation. This may mean your transition matrix disallows necessary transitions.  E.g., if your ranges are 'A' and 'B', and your model is DEC, then allowing range 'AB' as a possible state is required, so that you can get from 'A' to 'B' via 'AB' as the intermediate. Alternatively, NaNs can be produced sometimes if your Maximum Likelihood (ML) search proposes weird parameter values (such as a negative rate or weight) or a parameter so small that required transitions have a probability that machine precision rounds to zero or negative.  Sometimes this seems to occur because optim, optimx, etc. propose parameters slightly outside the user-specified upper and lower (min/max) boundaries for some reason. One solution is often to narrow the min/max limits. \n\nYou are using another solution: Normally, this would be a stop error, but you specified that BioGeoBEARS_run_object$on_NaN_error=", on_NaN_error, "\n\n"))
						}
					
					return(on_NaN_error)
					}
				
				}

			if (any(is.nan(condlikes_Right)))
				{
				
				if (printlevel > 0)
					{
					cat("\n\nindependent likelihood %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS returned NaNs:\n\n")
				
					cat("\n\nPrinting parameter values being attempted when this occurred:\n\n")
					if (is.null(BioGeoBEARS_model_object))
						{
						print("BioGeoBEARS_model_object was not passed to calc_loglike_sp() for printing.")
						} else {
						print(BioGeoBEARS_model_object)
						} # END if (is.null(BioGeoBEARS_model_object))
				
					cat("j=", j, "\n")
					cat("phy2$edge.length[j]=", phy2$edge.length[j], "\n")
					cat("\n\nprint(condlikes_Right):\n\n")
					print(condlikes_Right)

					cat("\n\nprint(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS):\n\n")
					print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS)

					cat("\n\nprint(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum,]):\n\n")
					print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum,])
					}
				#print(coo_n)
				#print(anorm)

				if (is.null(on_NaN_error))
					{
					stop("\n\nStopping on error in dense exponentiation downpass (right branch): NaNs produced in likelihood calculation. This may mean your transition matrix disallows necessary transitions.  E.g., if your ranges are 'A' and 'B', and your model is DEC, then allowing range 'AB' as a possible state is required, so that you can get from 'A' to 'B' via 'AB' as the intermediate. Alternatively, NaNs can be produced sometimes if your Maximum Likelihood (ML) search proposes weird parameter values (such as a negative rate or weight) or a parameter so small that required transitions have a probability that machine precision rounds to zero or negative.  Sometimes this seems to occur because optim, optimx, etc. propose parameters slightly outside the user-specified upper and lower (min/max) boundaries for some reason. One solution is often to narrow the min/max limits. \n\nAnother solution: To have this error report an extremely low log-likelihood,, set BioGeoBEARS_run_object$on_NaN_error to something like -1e50.\n\n")
					}
				
				if (printlevel > 0)
					{
					print("print(on_NaN_error):")
					print(on_NaN_error)
					}

				if ( (is.numeric(on_NaN_error)) && (return_what == "loglike") )
					{
					if (printlevel > 0)
						{
						warning(paste0("\n\nWarning on error in dense exponentiation downpass (right branch): NaNs produced in likelihood calculation. This may mean your transition matrix disallows necessary transitions.  E.g., if your ranges are 'A' and 'B', and your model is DEC, then allowing range 'AB' as a possible state is required, so that you can get from 'A' to 'B' via 'AB' as the intermediate. Alternatively, NaNs can be produced sometimes if your Maximum Likelihood (ML) search proposes weird parameter values (such as a negative rate or weight) or a parameter so small that required transitions have a probability that machine precision rounds to zero or negative.  Sometimes this seems to occur because optim, optimx, etc. propose parameters slightly outside the user-specified upper and lower (min/max) boundaries for some reason. One solution is often to narrow the min/max limits. \n\nYou are using another solution: Normally, this would be a stop error, but you specified that BioGeoBEARS_run_object$on_NaN_error=", on_NaN_error, "\n\n"))
						}

					return(on_NaN_error)
					}
					
				} # END if (any(is.nan(condlikes_Right)))
			
			if (printlevel >= 2) {
			txt = paste("condlikes at bottom of L: ", paste(round(condlikes_Left, 4), collapse=" ", sep=""), sep="")
			print(txt)
			}
			
			if (printlevel >= 2) {
			txt = paste("condlikes at bottom of R: ", paste(round(condlikes_Right, 4), collapse=" ", sep=""), sep="")
			print(txt)
			}
			#condlikes_Left <- expm(Qmat * phy$edge.length[i], method="Ward77") %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum, 
			#  ]
			#condlikes_Right <- expm(Qmat * phy$edge.length[j], method="Ward77") %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum, 
			#  ]
			
			
			#condlikes_Left <- exp(Qmat * phy$edge.length[i]) %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum, 
			#  ]
			#condlikes_Right <- exp(Qmat * phy$edge.length[j]) %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum, 
			#  ]
			} else {
			#######################################################
			# Sparse matrix exponentiation
			#######################################################
			# sparse == TRUE

			# For rapid exponentiation of sparse matrices, we use myDMEXPV, and input
			# starting probabilities/likelihoods, and output ending probabilities/likelihoods
			#
			# When there ARE inputprobs_for_fast, myDMEXPV is invoked, which appears to do a 
			# forward-probability matrix exponentiation calculation, not a backward probability
			# matrix exponentiation calculation.  This produces a positive value for a state which is 
			# impossible in the ancestor (e.g. the null range is impossible in an ancestor)
			#
			# To fix this, multiple the output probabilities by 0 if check_for_0_rows[i[ == TRUE

			if (printlevel >= 2)
				{
				print("sparse matrix exponentiation")
				}
			
			# Conditional likelihoods of data given the states at the bottom of left branch
			#condlikes_Left = try (
			#expokit_dgexpv_Qmat(Qmat=tmpQmat_in_REXPOKIT_coo_fmt, t=phy2$edge.length[i], inputprobs_for_fast=relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,], transpose_needed=TRUE, transform_to_coo_TF=FALSE, coo_n=numstates, anorm=anorm, check_for_0_rows=TRUE)
			#)
			#txt = paste0("Note: sparse matrix exponentiation on left branch, length t=", round(phy2$edge.length[i], 3), ". Transition matrix: ", numstates, " states.")
			#cat("\n")
			#print("relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,]")
			tip_probs_left = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,]
			
			try_result_Left = try (
			kexpmv::expokit_dgexpv(mat=tmpQmat_in_kexpmv_crs_fmt, t=phy2$edge.length[i], vector=relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,], transpose_needed=TRUE, transform_to_crs=FALSE, crs_n=length(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,]), anorm=NULL, mxstep=10000, tol=as.numeric(1e-10))
			)
			
			if (printlevel >= 1)
				{
				txt = "L"
				cat(txt)
				}
			
			
			# Error check
			if (class(try_result_Left) == "try-error")
				{
				if (printlevel > 0)
					{
					cat("\n\ntry-error on kexpmv::expokit_dgexpv():\n\n")
					cat("i=", i, "\n")
					cat("phy2$edge.length[i]=", phy2$edge.length[i], "\n")
					print(tmpQmat_in_kexpmv_crs_fmt)
					print(phy2)
					print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS)
					print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,])
					print(coo_n)
					print(anorm)
					print("BioGeoBEARS_model_object")
					print(BioGeoBEARS_model_object)
					}
				
				stoptxt = "\n\nStopping on error in sparse exponentiation downpass (treepiece, aka a branch segment): NaNs produced in likelihood calculation. This may mean your transition matrix disallows necessary transitions.  E.g., if your ranges are 'A' and 'B', and your model is DEC, then allowing range 'AB' as a possible state is required, so that you can get from 'A' to 'B' via 'AB' as the intermediate. Alternatively, NaNs can be produced sometimes if your Maximum Likelihood (ML) search proposes weird parameter values (such as a negative rate or weight) or a parameter so small that required transitions have a probability that machine precision rounds to zero or negative.  Sometimes this seems to occur because optim, optimx, etc. propose parameters slightly outside the user-specified upper and lower (min/max) boundaries for some reason. One solution is often to narrow the min/max limits. \n\nAnother solution: To have this error report an extremely low log-likelihood, set BioGeoBEARS_run_object$on_NaN_error to something like -1e50.\n\n"
				
				if (is.null(on_NaN_error))
					{
					cat(stoptxt)
					stop(stoptxt)
					}
				
				if (printlevel > 0)
					{
					print("print(on_NaN_error):")
					print(on_NaN_error)
					}
					
				if ( (is.numeric(on_NaN_error)) && (return_what == "loglike") )
					{
					if (printlevel > 0)
						{
						warning(paste0("\n\nWarning  on error in sparse exponentiation downpass (left branch): NaNs produced in likelihood calculation. This may mean your transition matrix disallows necessary transitions.  E.g., if your ranges are 'A' and 'B', and your model is DEC, then allowing range 'AB' as a possible state is required, so that you can get from 'A' to 'B' via 'AB' as the intermediate. Alternatively, NaNs can be produced sometimes if your Maximum Likelihood (ML) search proposes weird parameter values (such as a negative rate or weight) or a parameter so small that required transitions have a probability that machine precision rounds to zero or negative.  Sometimes this seems to occur because optim, optimx, etc. propose parameters slightly outside the user-specified upper and lower (min/max) boundaries for some reason. One solution is often to narrow the min/max limits. \n\nYou are using another solution: Normally, this would be a stop error, but you specified that BioGeoBEARS_run_object$on_NaN_error=", on_NaN_error, "\n\n"))
						}
					return(on_NaN_error)
					} else {
					stop(stoptxt)
					}
				} # END if (any(is.nan(condlikes_Left)))

			# Load up the conditional likelihoods
			condlikes_Left = c(try_result_Left$output_probs)
			condlikes_Left[condlikes_Left<0] = 0.0
			
			# Conditional likelihoods of data given the states at the bottom of right branch
			# condlikes_Right = try(
	 		#	expokit_dgexpv_Qmat(Qmat=tmpQmat_in_REXPOKIT_coo_fmt, t=phy2$edge.length[j], inputprobs_for_fast=relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum,], transpose_needed=TRUE, transform_to_coo_TF=FALSE, coo_n=coo_n, anorm=anorm, check_for_0_rows=TRUE)
 			# )

			#print("relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum,]")
			#tip_probs_right = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum,]

			try_result_Right = try (
			kexpmv::expokit_dgexpv(mat=tmpQmat_in_kexpmv_crs_fmt, t=phy2$edge.length[j], vector=relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum,], transpose_needed=TRUE, transform_to_crs=FALSE, crs_n=length(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum,]), anorm=NULL, mxstep=10000, tol=as.numeric(1e-10))
			)

			if (printlevel >= 1)
				{
				txt = "R"
				cat(txt)
				}
			

			# Error check
			if (class(try_result_Right) == "try-error")
				{
				if (printlevel > 0)
					{
					cat("\n\ntry-error on kexpmv::expokit_dgexpv():\n\n")
					cat("j=", j, "\n")
					cat("phy2$edge.length[j]=", phy2$edge.length[j], "\n")
					print(tmpQmat_in_kexpmv_crs_fmt)
					print(phy2)
					print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS)
					print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum,])
					print(coo_n)
					print(anorm)
					print("BioGeoBEARS_model_object")
					print(BioGeoBEARS_model_object)
					}

				if (is.null(on_NaN_error))
					{
					stop("\n\nStopping on error in sparse exponentiation downpass (right branch): NaNs produced in likelihood calculation. This may mean your transition matrix disallows necessary transitions.  E.g., if your ranges are 'A' and 'B', and your model is DEC, then allowing range 'AB' as a possible state is required, so that you can get from 'A' to 'B' via 'AB' as the intermediate. Alternatively, NaNs can be produced sometimes if your Maximum Likelihood (ML) search proposes weird parameter values (such as a negative rate or weight) or a parameter so small that required transitions have a probability that machine precision rounds to zero or negative.  Sometimes this seems to occur because optim, optimx, etc. propose parameters slightly outside the user-specified upper and lower (min/max) boundaries for some reason. One solution is often to narrow the min/max limits. \n\nAnother solution: To have this error report an extremely low log-likelihood, set BioGeoBEARS_run_object$on_NaN_error to something like -1e50.\n\n")
					}
				
				if (printlevel > 0)
					{
					print("print(on_NaN_error):")
					print(on_NaN_error)
					}
				
				if ( (is.numeric(on_NaN_error)) && (return_what == "loglike") )
					{
					if (printlevel > 0)
						{
						warning(paste0("\n\nWarning  on error in sparse exponentiation downpass (right branch): NaNs produced in likelihood calculation. This may mean your transition matrix disallows necessary transitions.  E.g., if your ranges are 'A' and 'B', and your model is DEC, then allowing range 'AB' as a possible state is required, so that you can get from 'A' to 'B' via 'AB' as the intermediate. Alternatively, NaNs can be produced sometimes if your Maximum Likelihood (ML) search proposes weird parameter values (such as a negative rate or weight) or a parameter so small that required transitions have a probability that machine precision rounds to zero or negative.  Sometimes this seems to occur because optim, optimx, etc. propose parameters slightly outside the user-specified upper and lower (min/max) boundaries for some reason. One solution is often to narrow the min/max limits. \n\nYou are using another solution: Normally, this would be a stop error, but you specified that BioGeoBEARS_run_object$on_NaN_error=", on_NaN_error, "\n\n"))
						}
					return(on_NaN_error)
					}
				} # END if (any(is.nan(condlikes_Right)))

			# Load up the conditional likelihoods
			condlikes_Right = c(try_result_Right$output_probs)
			condlikes_Right[condlikes_Right<0] = 0.0

			#print(c(condlikes_Left))
			#print(c(condlikes_Right))

# 			if (left_desc_nodenum == 53)
# 				{
# 				print("saving relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,]")
# 				#print(dput(tmpQmat_in_REXPOKIT_coo_fmt))
# 				save(try_result_Left, file="try_result_Left.Rdata")
# 				save(try_result_Right, file="try_result_Right.Rdata")
# 
# 				save(tip_probs_left, file="tip_probs_left.Rdata")
# 				save(tip_probs_right, file="tip_probs_right.Rdata")
# 				save(condlikes_Left, file="condlikes_Left.Rdata")
# 				save(condlikes_Right, file="condlikes_Right.Rdata")
# 				save(Qmat, file="Qmat.Rdata")
# 				save(tmpQmat_in_REXPOKIT_coo_fmt, file="tmpQmat_in_REXPOKIT_coo_fmt.Rdata")
# 				}

			} # END if (sparse==FALSE) (ENDING MATRIX EXPONENTIATION)	
		
		# Zero out impossible states
		if (!is.null(states_allowed_TF))
			{
			condlikes_Left[states_allowed_TF==FALSE] = 0
			condlikes_Right[states_allowed_TF==FALSE] = 0
			}
	

	

		
		
		# Save the conditional likelihoods of the data at the bottoms of each branch 
		# (needed for marginal ancestral state probabilities)
		if (calc_ancprobs == TRUE)
			{
			# Every node (except maybe the root) has a branch below it, and there is also a 
			# relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS at the bottom of this branch
			relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[left_desc_nodenum,] = condlikes_Left / sum(condlikes_Left)
			relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[right_desc_nodenum,] = condlikes_Right / sum(condlikes_Right)
			}			
			
		
		
		
		# If there is no speciational model, you are assuming 100% sympatry (range duplication)
		# at each speciation event
		#
		# In this case, you can just multiply the two conditional likelihood matrices together
		#
		# Also, if a branch is extremely short (a "hook"), this is essentially a zero-length
		# branch, we are assuming that this represents the range of a lineage at that 
		# point.  There is no speciation event here -- both "lineages" inherit
		# the same range.  This allows fossils to closely influence ancestral states.
		#
		# This was developed with Kaitlin Maguire over several years of screwing around.
		
		# Check for a short "hook" branch; if found, use just allopatric speciational model

		# get the correct edge
		left_edge_TF = phy2$edge[,2] == left_desc_nodenum
		right_edge_TF = phy2$edge[,2] == right_desc_nodenum
		
		# Check the branchlength of each edge
		is_leftbranch_hook_TF = phy2$edge.length[left_edge_TF] < min_branchlength
		is_rightbranch_hook_TF = phy2$edge.length[right_edge_TF] < min_branchlength
		hooknode_TF = (is_leftbranch_hook_TF + is_rightbranch_hook_TF) > 0
		
# 		print(left_desc_nodenum)
# 		print(right_desc_nodenum)
# 		
# 		print(phy2$edge.length[left_desc_nodenum])
# 		print(phy2$edge.length[right_desc_nodenum])
# 		
# 		print(is_leftbranch_hook_TF)
# 		print(is_rightbranch_hook_TF)
# 		
# 		print(hooknode_TF)
		if (use_cpp == FALSE)
			{
			if ( (is.null(spPmat) || hooknode_TF==TRUE) )
				{
				##################################
				# no speciational model of range change
				##################################
				node_likelihood <- condlikes_Left * condlikes_Right
				if (printlevel >= 2)
					{
					print("use_cpp=FALSE, direct multiplication")
					}
				} else {
				##################################
				# combine the likelihoods from each branch bottom with this speciational model
				# of range change
				##################################

				if (printlevel >= 2)
					{
					print("use_cpp=FALSE, speciation model")
					}
				
				# for each ancestral state, get prob of branch pairs 
				outmat = matrix(0, nrow=nrow(spPmat), ncol=ncol(spPmat))
				
				# Get the probs of each pair, list by row
				# (This might be a slightly slow step for large matrices)
				
				# Exclude "_" ranges from this calculation, as if there is speciation,
				# they are not going extinct
				
				if (include_null_range == TRUE)
					{
					probs_of_each_desc_pair = c(base::t(outer(c(condlikes_Left[-1]), c(condlikes_Right[-1]), "*")))
					} else {
					probs_of_each_desc_pair = c(base::t(outer(c(condlikes_Left), c(condlikes_Right), "*")))
					}
				#sum(probs_of_each_desc_pair)
	
				# Make a matrix with the probability of each pair, in each cell
				# (duplicated for each ancestral state)
				for (spi in 1:nrow(spPmat))
					{
					outmat[spi,] = probs_of_each_desc_pair
					}
				# Multiply the probabilities of each pair, by the probability of each
				# type of speciation event, to get the probabilities at the 
				# ancestral node
				outmat2 = outmat * spPmat
				
				node_likelihood_with_speciation = rowSums(outmat2)
				
				
				# THIS ZERO IS ALREADY OVER-WRITING THE NULL STATE LIKELIHOOD!!
				
				# Add the 0 back in, representing 0 probability of "_"
				# range just below speciation event
				if (include_null_range == TRUE)
					{
					node_likelihood = c(0, node_likelihood_with_speciation)
					} else {
					node_likelihood = node_likelihood_with_speciation
					}
				
				print("node_likelihood:")
				print(node_likelihood)
				
				#######################################################
				# If the states/likelihoods have been fixed at a particular node
				#######################################################
				if (!is.null(fixnode))
					{
					# For multiple fixnodes
					if (length(fixnode) > 1)
						{
						# Get the matching node
						TF = (anc == fixnode)
						temporary_fixnode = fixnode[TF]
						temporary_fixlikes = c(fixlikes[TF,])
						} else {
						temporary_fixnode = fixnode
						temporary_fixlikes = c(fixlikes)
						}

					
					if ((length(temporary_fixnode) > 0) && (anc == temporary_fixnode))
						{
						# If the node is fixed, ignore the calculation for this node, and
						# instead use the fixed likelihoods (i.e., the "known" state) for
						# this node.
						# fix the likelihoods of the (NON-NULL) states
						node_likelihood = node_likelihood * temporary_fixlikes
						}
					} # end if (!is.null(fixnode))
				}
			} else {
			##################################
			# use_cpp == TRUE
			##################################
			if ( hooknode_TF==TRUE )
				{
				##################################
				# no speciational model of range change
				##################################
				if (printlevel >= 2) 
					{
					print("use_cpp=TRUE, direct multiplication")
					}
# 				print("use_cpp=TRUE, direct multiplication")
# 				print("condlikes_Left:")
# 				print(condlikes_Left)
# 				print("condlikes_Right:")
# 				print(condlikes_Right)
				node_likelihood <- condlikes_Left * condlikes_Right			
				} else {
				if ( is.null(spPmat_inputs)==TRUE )
					{
					if ( is.null(spPmat) == TRUE )
						{
						##################################
						# no speciational model of range change
						##################################
						node_likelihood <- condlikes_Left * condlikes_Right
						print ("WARNING #2: (use_cpp==TRUE && is.null(spPmat_inputs)==TRUE)")
						print("use_cpp=TRUE, no spPmat_inputs, direct multiplication")

						} else {
						# otherwise ...
						if (printlevel >= 2)
							{
							print("use_cpp=TRUE, no spPmat_inputs, speciation model")
							}
						##################################
						# combine the likelihoods from each branch bottom with this speciational model
						# of range change
						##################################

						# for each ancestral state, get prob of branch pairs 
						outmat = matrix(0, nrow=nrow(spPmat), ncol=ncol(spPmat))
						
						# Get the probs of each pair, list by row
						# (This might be a slightly slow step for large matrices)
						
						# Exclude "_" ranges from this calculation, as if there is speciation,
						# they are not going extinct
						
						if (include_null_range == TRUE)
							{
							probs_of_each_desc_pair = c(base::t(outer(c(condlikes_Left[-1]), c(condlikes_Right[-1]), "*")))
							} else {
							probs_of_each_desc_pair = c(base::t(outer(c(condlikes_Left), c(condlikes_Right), "*")))
							}
						#sum(probs_of_each_desc_pair)
			
						# Make a matrix with the probability of each pair, in each cell
						# (duplicated for each ancestral state)
						for (spi in 1:nrow(spPmat))
							{
							outmat[spi,] = probs_of_each_desc_pair
							}
						# Multiply the probabilities of each pair, by the probability of each
						# type of speciation event, to get the probabilities at the 
						# ancestral node
						outmat2 = outmat * spPmat
						
						node_likelihood_with_speciation = rowSums(outmat2)
						
						# THIS ZERO IS ALREADY OVER-WRITING THE NULL STATE LIKELIHOOD!!
						
						# Add the 0 back in, representing 0 probability of "_"
						# range just below speciation event
						if (include_null_range == TRUE)
							{
							node_likelihood = c(0, node_likelihood_with_speciation)
							} else {
							node_likelihood = node_likelihood_with_speciation
							}

						print("node_likelihood:")
						print(node_likelihood)


						#######################################################
						# If the states/likelihood have been fixed at a particular node
						#######################################################
						if (!is.null(fixnode))
							{
							# For multiple fixnodes
							if (length(fixnode) > 1)
								{
								# Get the matching node
								TF = (anc == fixnode)
								temporary_fixnode = fixnode[TF]
								temporary_fixlikes = c(fixlikes[TF,])
								} else {
								temporary_fixnode = fixnode
								temporary_fixlikes = c(fixlikes)
								}
					
							if ((length(temporary_fixnode) > 0) && (anc == temporary_fixnode))
								{
								# If the node is fixed, ignore the calculation for this node, and
								# instead use the fixed likelihoods (i.e., the "known" state) for
								# this node.
								# fix the likelihoods of the (NON-NULL) states
								node_likelihood = node_likelihood * temporary_fixlikes
								}
							}
						}
					} else {
					## if ( is.null(spPmat_inputs)==FALSE )
					if (printlevel >= 2)
						{
						print("use_cpp=TRUE, yes spPmat_inputs, speciation model")
						}
					################################
					# Use the C++ function!
					################################
					##################################
					# combine the likelihoods from each branch bottom with this speciational model
					# of range change
					##################################
					
					# Calculate the likelihoods at each node,
					# give the speciation model, and input probs for 
					# each branch
					if (include_null_range == TRUE)
						{
						ca = condlikes_Left[-1]
						cb = condlikes_Right[-1]
						} else {
						ca = condlikes_Left
						cb = condlikes_Right						
						}
					
					# Rcpp does weird alterations to the input variables, so use each only once!
					#tmpca_1 = ca
					#tmpcb_1 = cb
					tmpca_2 = ca
					tmpcb_2 = cb
					
					# (above) Fix "l" (which is states_indices, i.e. a list of lists of state_indices)
					# so that there is no "null geographic range", i.e. no "_", no "-", no "NA"
					l = spPmat_inputs$l		# states_indices
					s = spPmat_inputs$s
					v = spPmat_inputs$v
					j = spPmat_inputs$j
					y = spPmat_inputs$y
					
					# dmat means combined dispersal multipliers matrix, NOT d
					dmat = spPmat_inputs$dmat
	
					# Print the matrix to screen from C++
					printmat = FALSE
					if (printlevel >= 2) {
					printmat = TRUE
					}
					
					# Calculate the rowsums (for input into rcpp_calc_anclikes_sp()
					# Actually, just do this ONCE (above)
					#Rsp_rowsums = rcpp_calc_anclikes_sp_rowsums(Rcpp_leftprobs=tmpca_1, Rcpp_rightprobs=tmpcb_1, l=l, s=s, v=v, j=j, y=y, dmat=dmat, maxent01s, maxent01v, maxent01j, maxent01y, printmat=printmat)
					
					if (printlevel >= 2)
						{
						print("Rsp_rowsums:")
						print(Rsp_rowsums)
						}
					
					#node_likelihood_with_speciation = rep(0.0, length(tmpca_2))
					# This calculates the speciation probabilities again at each node; this is inefficient
					if (cppSpMethod == 1)
						{
						node_likelihood_with_speciation = rcpp_calc_anclikes_sp(Rcpp_leftprobs=tmpca_2, Rcpp_rightprobs=tmpcb_2, l=l, s=s, v=v, j=j, y=y, dmat=dmat, maxent01s=maxent01s, maxent01v=maxent01v, maxent01j=maxent01j, maxent01y=maxent01y, max_minsize_as_function_of_ancsize=max_minsize_as_function_of_ancsize, Rsp_rowsums=Rsp_rowsums, printmat=printmat)						
						#print(node_likelihood_with_speciation[1:5])
						} # END if (cppSpMethod == 1)
					
					# Really, we should just iterate through the COO_probs_list using a C++ function
					# 2012-12-04 NJM says: this KICKS ASS in C++!!
					
					if (cppSpMethod == 2)
						{
						node_likelihood_with_speciation2 = rep(0.0, length(tmpca_2))
						node_likelihood_with_speciation2 = rcpp_calc_anclikes_sp_using_COOprobs(Rcpp_leftprobs=tmpca_2, Rcpp_rightprobs=tmpcb_2, RCOO_left_i_list=COO_probs_list[[1]], RCOO_right_j_list=COO_probs_list[[2]], RCOO_probs_list=COO_probs_list[[3]], Rsp_rowsums=Rsp_rowsums, printmat=printmat)
						#print(node_likelihood_with_speciation2[1:5])
						node_likelihood_with_speciation = node_likelihood_with_speciation2
						} # END if (cppSpMethod == 2)

					if (cppSpMethod == 3)
						{
						node_likelihood_with_speciation3 = rep(0.0, length(tmpca_2))
						node_likelihood_with_speciation3 = rcpp_calc_splitlikes_using_COOweights_columnar(Rcpp_leftprobs=tmpca_2, Rcpp_rightprobs=tmpcb_2, COO_weights_columnar=COO_weights_columnar, Rsp_rowsums=Rsp_rowsums, printmat=printmat)
						#print(node_likelihood_with_speciation2[1:5])
						node_likelihood_with_speciation = node_likelihood_with_speciation3
						
						# Null ranges can leave NaNs
						# when the downpass is calculated on 
						# ranges x trait m's
						if (length(m) > 0)
							{
							nanTF = is.nan(node_likelihood_with_speciation)
							
							# IF this is the issue, the number of NaNs will be length(m)-1
							# (don't want to zero out any weird NaNs)
							if (sum(nanTF) == (length(m)-1))
								{
								node_likelihood_with_speciation[nanTF] = 0.0
								}
							}
						} # END if (cppSpMethod == 3)


					

					
					if (printlevel >= 2)
						{
						print("node_likelihood_with_speciation:")
						print(node_likelihood_with_speciation)
						}
	
					if (include_null_range == TRUE)
						{
						node_likelihood = c(0, node_likelihood_with_speciation)
						} else {
						node_likelihood = node_likelihood_with_speciation
						}



					#######################################################
					# If the states/likelihood have been fixed at a particular node
					#######################################################
					if (!is.null(fixnode))
						{
						# For multiple fixnodes
						if (length(fixnode) > 1)
							{
							# Get the matching node
							TF = (anc == fixnode)
							temporary_fixnode = fixnode[TF]
							temporary_fixlikes = c(fixlikes[TF,])
							} else {
							temporary_fixnode = fixnode
							temporary_fixlikes = c(fixlikes)
							}
					
						if ((length(temporary_fixnode) > 0) && (anc == temporary_fixnode))
							{
							# If the node is fixed, ignore the calculation for this node, and
							# instead use the fixed likelihoods (i.e., the "known" state) for
							# this node.
							# fix the likelihoods of the (NON-NULL) states
							node_likelihood = node_likelihood * temporary_fixlikes
							}
						} # END if (!is.null(fixnode))

					#formatted_table = conditional_format_table(t(node_likelihood))
					#print(formatted_table)

					#print("!is.null(states_allowed_TF)")
					#print(!is.null(states_allowed_TF))
					# Zero out impossible states					
					if (!is.null(states_allowed_TF))
						{
						
						node_likelihood[states_allowed_TF==FALSE] = 0
						} # END if (!is.null(states_allowed_TF))
					} # END if ( is.null(spPmat_inputs)==TRUE )
				} # END if ( hooknode_TF==TRUE )
			} # END if (use_cpp == FALSE)


		########################################################################
		# If you are at the root node, you also multiply by the probabilities
		# of the starting states
		########################################################################
		if (i == max(edges_to_visit))
			{
			# If not specified, you are assuming even probabilities of each state
			if (is.null(probs_of_states_at_root) == TRUE)
				{
				node_likelihood = node_likelihood
				} else {
				# Otherwise, use user-specified probs_of_states_at_root
				# These could be:
				# 1. User-specified base frequencies
				# 2. Observed frequencies at tips
				#    (WARNING: likely have states with 0 tip observations!)
				# 3. Observed frequencies in some reference set of species
				# 4. Some distribution on range sizes, even within range sizes
				# 5. Equilibrium probabilities assuming Qmat rate matrix is
				#    run to infinite time
				node_likelihood = probs_of_states_at_root * node_likelihood 
				}
			}

		
		#txt = paste("left node: ", left_desc_nodenum, ", right node:", right_desc_nodenum, sep="")
		#print(txt)
		#print(node_likelihood)
		
		nanTF = is.nan(node_likelihood)
		
		if (sum(nanTF) > 0)
			{
			anc <- phy2$edge[i, 1]

			
			# Default printlevel is 1
			if (printlevel >= 2)
				{
				errortxt = paste("STOP ERROR in calc_loglike_sp -- LONG error message (because printlevel >= 2).\n\nWith 'if (sum(nanTF) > 0)': you have ", sum(nanTF), "NaNs in the 'node_likelihood' at edges #", i, " & ", j, ", node #", anc, "!", sep="")

				cat("\n\n")
				cat(errortxt)
				cat("\n\n")
			
				cat("\n\nnanTF:\n\n")
				print(c(nanTF))
				cat("\n\nprint(node_likelihood):\n\n")
				print(c(node_likelihood))
				
				cat("\n","i=", i, ", j=", j, "\n", sep="")
				cat("\n","left_desc_nodenum=", left_desc_nodenum, ", right_desc_nodenum=", right_desc_nodenum, "\n", sep="")

				cat("\nFree parameters from BioGeoBEARS_model_object:\n")
				free_params = BioGeoBEARS_model_object@params_table$est[BioGeoBEARS_model_object@params_table$type=="free"]
				free_params = as.data.frame(matrix(free_params, nrow=1), stringsAsFactors=FALSE)
				names(free_params) = row.names(BioGeoBEARS_model_object@params_table)[BioGeoBEARS_model_object@params_table$type=="free"]
				print(free_params)
				} else {
				errortxt = "STOP ERROR in calc_loglike_sp: NaN error(s) in likelihood calculation. No 'on_NaN_error' given, so stopping. Set printlevel > 1 for longer error message."
				} # END if (printlevel >= 2)

			if (is.numeric(on_NaN_error) == FALSE)
				{
				stop(errortxt)
				} # END if (is.numeric(on_NaN_error) == FALSE)
		
			# And for the ancestor edge (i or j shouldn't matter, should produce the same result!!!)
			#anc <- phy2$edge[i, 1]
			#plot(phy2)
			#nodelabels()
			
			#printall(prt(phy2))
			
			#print(tip_condlikes_of_data_on_each_state[left_desc_nodenum,])
			#print(tip_condlikes_of_data_on_each_state[right_desc_nodenum,])
			
			if ( (is.numeric(on_NaN_error)) && (return_what == "loglike") )
					{
					# Default printlevel is 1
					if (printlevel >= 2)
						{
						warning(paste0("\n\nWarning  on error in combining likelihoods at node, at 'if (sum(nanTF) > 0)': NaNs produced in likelihood calculation. This may mean your transition matrix disallows necessary transitions.  E.g., if your ranges are 'A' and 'B', and your model is DEC, then allowing range 'AB' as a possible state is required, so that you can get from 'A' to 'B' via 'AB' as the intermediate. Alternatively, NaNs can be produced sometimes if your Maximum Likelihood (ML) search proposes weird parameter values (such as a negative rate or weight) or a parameter so small that required transitions have a probability that machine precision rounds to zero or negative.  Sometimes this seems to occur because optim, optimx, etc. propose parameters slightly outside the user-specified upper and lower (min/max) boundaries for some reason. One solution is often to narrow the min/max limits. \n\nYou are using another solution: Normally, this would be a stop error, but you specified that BioGeoBEARS_run_object$on_NaN_error=", on_NaN_error, "\n\n"))
						} else {
						if (printlevel > 0)
							{
							txt = paste0("WARNING in calc_loglike_sp: NaN error(s). Returning $on_NaN_error=", on_NaN_error, ". $printlevel > 1 for more info.")
							warning(txt)
							}
						}
					return(on_NaN_error)
					}
			} # END if (sum(nanTF) > 0)

		
		total_likelihood_for_node = sum(node_likelihood)
		
		#print(total_likelihood_for_node)
		
		computed_likelihoods_at_each_node[anc] = total_likelihood_for_node
		
		relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[anc, ] = node_likelihood / total_likelihood_for_node
		condlikes_of_each_state[anc, ] = node_likelihood
		
		#print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS)
		} # end downpass
	
		if (force_sparse == TRUE)
			{
			if (printlevel >= 1)
				{
				txt = paste0("...segments done\n")
				cat(txt)
				}
			} # END if (force_sparse == TRUE)
	
	
	#######################################################
	# END DOWNPASS FROM THE TIPS TO THE ROOT
	#######################################################
	rootnode = anc
	
	relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS
	# relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS2 = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS

	# You are now at the anc (ancestor node; in Psychotria, node 20)
	# If you want, you could calculate the likelihood down to the bottom of a root edge
	# below the root node
	# check if the rootedge option is desired, and the phylogeny HAS a root edge
	#phy2$root.edge = 0.1
	if ( (rootedge == TRUE) && (!is.null(phy2$root.edge)) && (phy2$root.edge > 0) )
		{
		#phy2$root.edge = 1.1
		root_edge_length = phy2$root.edge
		
		if (force_sparse == FALSE)
			{
			# Dense matrix exponentiation on root branch segment
			independent_likelihoods_on_root_edge = expokit_dgpadm_Qmat(Qmat=Qmat, t=root_edge_length, transpose_needed=TRUE)
		
			# Get the likelihoods at the bottom of the branch, condition on the relative likelihoods at the top 
			# (here, the node at the "top" is the Last Common Ancestor node on the tree)
			conditional_likelihoods_at_bottom_of_root_branch = c(independent_likelihoods_on_root_edge %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[anc, ])
			} else {
			# Sparse matrix exponentiation on root branch segment
			try_result_root_segment = try (
				kexpmv::expokit_dgexpv(mat=tmpQmat_in_kexpmv_crs_fmt, t=root_edge_length, vector=relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[anc, ], transpose_needed=TRUE, transform_to_crs=FALSE, crs_n=length(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[anc, ]), anorm=NULL, mxstep=10000, tol=as.numeric(1e-10))
				)
			
			conditional_likelihoods_at_bottom_of_root_branch = c(try_result_root_segment$output_probs)
			}

		total_likelihood_for_bottom_of_root_branch = sum(conditional_likelihoods_at_bottom_of_root_branch)
		total_likelihood_for_bottom_of_root_branch
		
		
		# Get the relative probability of each state at the bottom of the root branch
		relative_probs_of_each_state_at_bottom_of_root_branch = c(conditional_likelihoods_at_bottom_of_root_branch) / total_likelihood_for_bottom_of_root_branch
		
		
		# Add the likelihood of the bottom of the root branch to the list of node likelihoods
		computed_likelihoods_at_each_node = c(computed_likelihoods_at_each_node, total_likelihood_for_bottom_of_root_branch)
		
		# Add a bottom row to the relative_probs
		relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS = rbind(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS, relative_probs_of_each_state_at_bottom_of_root_branch)
		condlikes_of_each_state = rbind(condlikes_of_each_state, conditional_likelihoods_at_bottom_of_root_branch)

		} else {
		# Otherwise, the root relative probabilities are just the last relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[anc, ]
		relative_probs_of_each_state_at_bottom_of_root_branch = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[anc, ]
		}
	
	# Zero out impossible states
	if (!is.null(states_allowed_TF))
		{
		relative_probs_of_each_state_at_bottom_of_root_branch[states_allowed_TF==FALSE] = 0
		relative_probs_of_each_state_at_bottom_of_root_branch = relative_probs_of_each_state_at_bottom_of_root_branch / sum(relative_probs_of_each_state_at_bottom_of_root_branch)
		}

	
	# why times 2?? -- probably this is some AIC thing from the original APE function
	#output_loglike = -2 * sum(log(comp[-TIPS]))
	#output_loglike = 2 * sum(log(comp[-TIPS]))
	
	total_loglikelihood = sum(log(computed_likelihoods_at_each_node))
	#print(total_loglikelihood)
		
	#######################################################
	# Calculating the probabilities of ancestral states at
	# nodes (and immediately after nodes) is not just a 
	# matter of normalizing the downpass conditional
	# likelihoods so they sum to 1.
	#
	# With traditional reversible models, one would re-root 
	# the phylogeny at each node (or just do an uppass),
	# but here we need to do an uppass with the FORWARD
	# model to get the ancprobs.
	#######################################################
	if (calc_ancprobs == TRUE)
		{
		if ((printlevel >= 1) && (stratified == FALSE))
			{
			cat("\nUppass starting for marginal ancestral states estimation!\n", sep="")
			}


		# First, figure out the probs AT THE LOWEST FORK (different from bottom of the root branch,
		# which may or may not exist)
		starting_probs = NULL
		
		# If there was a root edge, start with that:
		if ( (rootedge == TRUE) && (!is.null(phy2$root.edge)) && (phy2$root.edge > 0) )
			{
			#######################################################
			# THERE IS A ROOT EDGE
			#######################################################

			# Do sparse or dense matrix exponentiation
			if (sparse==FALSE)
				{
				# Dense matrix exponentiation
				# Need to do a forward matrix exponentiation
				actual_probs_after_forward_exponentiation = calc_prob_forward_onebranch_dense(relprobs_branch_bottom=relative_probs_of_each_state_at_bottom_of_root_branch, branch_length=root_edge_length, Qmat)
				# *IF* the null range is excluded, ensure that it has
				# 0 probability on the UPPASS
				# (The end of a fossil range could be a tip with null range, but
				#  this seems like a dubious idea as it elevates all rates and leads 
				#  to large uncertainty -- see Matzke & Maguire, SVP 2011)
				if (include_null_range == TRUE)
					{
					actual_probs_after_forward_exponentiation[1] = 0 	# NULL range is impossible
					} # END if (include_null_range == TRUE)
				
				# Normalize the uppass probabilities by the sum
				actual_probs_after_forward_exponentiation = actual_probs_after_forward_exponentiation / sum(actual_probs_after_forward_exponentiation)
				} else {
				# Sparse matrix exponentiation
				actual_probs_after_forward_exponentiation = calc_prob_forward_onebranch_sparse(relprobs_branch_bottom=relative_probs_of_each_state_at_bottom_of_root_branch, branch_length=root_edge_length, tmpQmat_in_REXPOKIT_coo_fmt=tmpQmat_in_REXPOKIT_coo_fmt, coo_n=coo_n, anorm=anorm, check_for_0_rows=TRUE)

				# *IF* the null range is excluded, ensure that it has
				# 0 probability on the UPPASS
				# (The end of a fossil range could be a tip with null range, but
				#  this seems like a dubious idea as it elevates all rates and leads 
				#  to large uncertainty -- see Matzke & Maguire, SVP 2011)
				if (include_null_range == TRUE)
					{
					actual_probs_after_forward_exponentiation[1] = 0 	# NULL range is impossible
					} # END if (include_null_range == TRUE)
				
				# Normalize the uppass probabilities by the sum
				actual_probs_after_forward_exponentiation = actual_probs_after_forward_exponentiation / sum(actual_probs_after_forward_exponentiation)
				}
			
			# Combine this uppass probability with the downpass condlikes for the states just below the root node.
			# anc is the ancestral node
			
			# relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS = results_object$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS
			# anc=20
			
			# Multiply the probs, then divide by the sum
			multiplied_probs = actual_probs_after_forward_exponentiation * relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[anc, ]
			starting_probs = multiplied_probs / sum(multiplied_probs)
			} else {
			#######################################################
			# NO ROOT EDGE
			#######################################################
			# Otherwise, just start with the bottom fork
			starting_probs = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[anc, ]
			}
		
		# Zero out impossible states
		if (!is.null(states_allowed_TF))
			{
			starting_probs[states_allowed_TF==FALSE] = 0
			starting_probs = starting_probs / sum(starting_probs)
			}
		
		
		#######################################################
		#######################################################
		# THIS IS AN UPPASS FROM THE ROOT TO THE TIPS
		#######################################################
		#######################################################

		# Check to make sure you have the necessary inputs
		if (exists("COO_weights_columnar") == FALSE)
			{
			stop("\nERROR_A: calc_loglike_sp requires 'COO_weights_columnar', 'Rsp_rowsums', and cppSpMethod==3 for marginal ancestral state estimations.\n")
			}
		if (exists("Rsp_rowsums") == FALSE)
			{
			stop("\nERROR_B: calc_loglike_sp requires 'COO_weights_columnar', 'Rsp_rowsums', and cppSpMethod==3 for marginal ancestral state estimations.\n")
			}
		if (cppSpMethod != 3)
			{
			stop("\nERROR_C: calc_loglike_sp requires 'COO_weights_columnar', 'Rsp_rowsums', and cppSpMethod==3 for marginal ancestral state estimations.\n")
			}

		# OK, from the starting probs at the bottom fork, we just need to do an uppass from anc, reversing the downpass
		# Get the descendant node nums
		
		# Keep track of the uppass probabilities
		# NJM 2015-04-07
		#relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc, ] = starting_probs
		relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc, ] = 1/length(starting_probs)

		# Visit edges in reverse order from the downpass
		edges_to_visit_uppass = seq(from=(num_internal_nodes*2), by=-2, length.out=num_internal_nodes)
		# Since we are going backwards
		#print(edges_to_visit_uppass)
		#print(i)
		#print(j)
		#cat("\n")
		
		#for (i in edges_to_visit_uppass)
		#j=edges_to_visit_uppass[1]
		for (j in edges_to_visit_uppass)		# Since we are going backwards
			{
			# First edge visited is i
			#print(i)
			
			# Its sister is j 
			#j <- i - 1
			i <- j - 1		# Since we are going backwards
			
			# Get the node numbers at the tips of these two edges		
			left_desc_nodenum <- phy2$edge[i, 2]
			right_desc_nodenum <- phy2$edge[j, 2]

			# And for the ancestor edge (i or j shouldn't matter, should produce the same result!!!)
			if (phy2$edge[i, 1] != phy2$edge[j, 1])
				{
				errortxt = "STOP ERROR in UPPASS: (phy2$edge[i, 1] != phy2$edge[i, 1]). Probably you didn't reorder() the tree in pruningwise fashion"
				cat("\n\n")
				cat(errortxt)
				cat("\n\n")
				stop(errortxt)
				}
			# anc is the focal node
			anc <- phy2$edge[i, 1]
			# ancedge
			anc_edgenum_TF = phy2$edge[,2] == anc
			anc_edgenum = (1:length(anc_edgenum_TF))[anc_edgenum_TF]
			
			# For the marginal state probability, uppass calculations
			# Get the mother and sister of "anc" (which is the focal node)
			mother_of_anc_TF = phy2$edge[,2] == anc
			mother_of_anc = phy2$edge[mother_of_anc_TF,1]

			sister_of_anc_TF = phy2$edge[,1] == mother_of_anc
			sister_of_anc_TF2 = (sister_of_anc_TF + mother_of_anc_TF) == 1
			sister_of_anc = phy2$edge[sister_of_anc_TF2,2]
			mother_of_anc
			sister_of_anc
			
			# Is the sister left or right?
			# (note: these are reversed from what you would get with:
			#  plot(tr); nodelabels()
			sister_is_LR = "rootnode"
			if (anc != rootnode)
				{
				if (sister_of_anc > anc)
					{
					sister_is_LR = "right"
					} else {
					sister_is_LR = "left"
					} # END if (sister_of_anc > anc)
				} # END if (anc != rootnode)
			
			# get the correct edges
			left_edge_TF = phy2$edge[,2] == left_desc_nodenum
			right_edge_TF = phy2$edge[,2] == right_desc_nodenum
			left_edgenum = (1:length(left_edge_TF))[left_edge_TF]
			right_edgenum = (1:length(right_edge_TF))[right_edge_TF]
			
			
			# Check the branchlength of each edge
			# It's a hook if either branch is super-short
			is_leftbranch_hook_TF = phy2$edge.length[left_edge_TF] < min_branchlength
			is_rightbranch_hook_TF = phy2$edge.length[right_edge_TF] < min_branchlength
			hooknode_TF = (is_leftbranch_hook_TF + is_rightbranch_hook_TF) > 0


			#cat(i, j, left_desc_nodenum, right_desc_nodenum, hooknode_TF, "\n", sep="	")
			
			
			# You start with these uppass probs
			relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc, ]
			
			############################################################
			# Apply speciation model to get the uppass probs at the 
			# bottom of the branch below the anc node
			############################################################
			if (hooknode_TF == TRUE)
				{
				# If you have a "hooknode" (short branch = direct ancestor), for
				# the uppass, it is simpler to convert the cladogenesis model
				# to an all-1s model
				temp_COO_weights_columnar = COO_weights_columnar

				# NJM 2016-02-24 -- see:
				# /drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/AAAB_M3_ancestor_check
				# ...for an example that traces the -2 issue in detail. 
				# Or:
				# http://phylo.wikidot.com/fossil-data-in-biogeographical-analysis-in-biogeobears#toc3
				#
				# Basically, this converts 1-based state numbers (e.g., 1-16, with
				# 1: null range
				# 2: A
				# 3: B
				# ...etc..
				#
				# ...to the 0-based state names in a cladogenesis matrix, where the
				# null range is automatically excluded (if used in the first place).
				#
				# Then the states in the cladogenesis matrix are numbered, starting from 0.
				# E.g.:
				# 0: A
				# 1: B
				# 2: C
				# ...etc...
				
				if (include_null_range == TRUE)
					{
					# NJM 2016-02-05 bug fix: -1 to -2
					# Uppass starting for marginal ancestral states estimation!
					# Error in if (sum(relprob_each_split_scenario) > 0) { : 
					#  missing value where TRUE/FALSE needed
					highest_clado_state_0based_considering_null_range = numstates - 2
					} else {
					# NJM 2016-02-05 bug fix: -0 to -1
					# Uppass starting for marginal ancestral states estimation!
					# Error in if (sum(relprob_each_split_scenario) > 0) { : 
					#  missing value where TRUE/FALSE needed
					highest_clado_state_0based_considering_null_range = numstates - 1
					} # if (include_null_range == TRUE)
				
				# Ancestral, left, and right states all the same
				temp_COO_weights_columnar[[1]] = 0:highest_clado_state_0based_considering_null_range
				temp_COO_weights_columnar[[2]] = 0:highest_clado_state_0based_considering_null_range
				temp_COO_weights_columnar[[3]] = 0:highest_clado_state_0based_considering_null_range
				temp_COO_weights_columnar[[4]] = rep(1, times=highest_clado_state_0based_considering_null_range+1)
				} else {
				temp_COO_weights_columnar = COO_weights_columnar
				} # END if (hooknode_TF == TRUE)

			#print("temp_COO_weights_columnar:")
			#print(temp_COO_weights_columnar)
							
			# Apply regular speciation model, with the weights 
			# given in COO_weights_columnar, and the 
			# normalization factor (sum of the weights across 
			# each row/ancestral state) in Rsp_rowsums.
			num_nonzero_split_scenarios = length(temp_COO_weights_columnar[[1]])
			
			relprobs_just_after_speciation_UPPASS_Left = rep(0, numstates)
			relprobs_just_after_speciation_UPPASS_Right = rep(0, numstates)

			# Go through each ancestral state, add up the probability contributed to the
			# left descendent corner and right descendent corner for each state
			
			# ORIGINAL UPPASS function (2014-03-ish)
			#relprobs_just_after_speciation = calc_uppass_probs_orig(anc, relative_probs_of_each_state_at_branch_top_AT_node_UPPASS, COO_weights_columnar, numstates, include_null_range)
			
			# NEW UPPASS function (2015-01-10) -- slight differences
			#########################################################################
			# Uppass probabilities EXCLUDING downpass likelihoods
			# Normal way, non-JOINT
			#########################################################################
			probs_ancstate = relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc,]
			#probs_ancstate = 1
			#prob_each_split_scenario_df2 = calc_uppass_scenario_probs_new(probs_ancstate, COO_weights_columnar, numstates, include_null_range, left_branch_downpass_likes=NULL, right_branch_downpass_likes=NULL, Rsp_rowsums=NULL)
			
			# Probs at the mother have been predetermined, in the uppass
			# 1. Get uppass probabilities at the base of the branch below the 
			#    focal (anc) node, including the probabilities coming down
			#    from the sister, and up from the mother.

			if (anc == rootnode)
				{
				# You are at the root ancestor node
				probs_at_mother = 1/length(starting_probs)
				likes_at_sister = 1/length(starting_probs)
				left_branch_downpass_likes = NULL
				right_branch_downpass_likes = NULL
				#relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[anc,] = NA		
				probs_of_mother_and_sister_uppass_to_anc = 1/length(starting_probs)
				relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc,] = probs_of_mother_and_sister_uppass_to_anc
				} else {
				# You are NOT at the root ancestor node
				probs_at_mother = relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[mother_of_anc,]
				likes_at_sister_branch_bottom = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[sister_of_anc,]

				if (sister_is_LR == "left")
					{	
					left_branch_downpass_likes = likes_at_sister_branch_bottom
					right_branch_downpass_likes = NULL
					}
				if (sister_is_LR == "right")
					{	
					left_branch_downpass_likes = NULL
					right_branch_downpass_likes = likes_at_sister_branch_bottom
					}
				
				# Calculate the uppass probs at the branch 	
				#print("calculation of uppass at split")
				#print(probs_at_mother)
				#print(left_branch_downpass_likes)
				#print(right_branch_downpass_likes)
				
				uppass_probs_at_bottom_below_anc_results = calc_uppass_probs_new2(probs_ancstate=probs_at_mother, COO_weights_columnar=temp_COO_weights_columnar, numstates=numstates, include_null_range=include_null_range, left_branch_downpass_likes=left_branch_downpass_likes, right_branch_downpass_likes=right_branch_downpass_likes, Rsp_rowsums=NULL)
				#print("...finished uppass at split")

				# Store
				if (sister_is_LR == "left")
					{	
					Rprobs_brbot_below_anc = uppass_probs_at_bottom_below_anc_results$relprobs_just_after_speciation_UPPASS_Right
					relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[anc,] = Rprobs_brbot_below_anc
					}
				if (sister_is_LR == "right")
					{
					Lprobs_brbot_below_anc = uppass_probs_at_bottom_below_anc_results$relprobs_just_after_speciation_UPPASS_Left
					relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[anc,] = Lprobs_brbot_below_anc
					}
				
				# 2. Exponentiate up from the mother to the focal/anc node
				probs_at_branch_bot = relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[anc,]
				
				if (force_sparse == FALSE)
					{
					# Dense matrix forward exponentiation
					probs_of_mother_and_sister_uppass_to_anc = probs_at_branch_bot %*% expokit_dgpadm_Qmat2(times=phy2$edge.length[anc_edgenum], Qmat=Qmat, transpose_needed=TRUE)
					} else {
					# Sparse matrix forward exponentiation
					probs_of_mother_and_sister_uppass_to_anc = calc_prob_forward_onebranch_sparse(relprobs_branch_bottom=probs_at_branch_bot, branch_length=phy2$edge.length[anc_edgenum], tmpQmat_in_REXPOKIT_coo_fmt=tmpQmat_in_REXPOKIT_coo_fmt, coo_n=coo_n, anorm=anorm, check_for_0_rows=FALSE, TRANSPOSE_because_forward=TRUE)
					}				
				
				# Store in uppass
				relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc,] = probs_of_mother_and_sister_uppass_to_anc
				} # END if (anc == rootnode)

			##################################################################
			# Finish uppass to tips
			##################################################################
			# Check if either the left or right descendant nodes are tips;
			# if so, do the exponentiation here, so as to completely fill
			# in the UPPASS table
			##################################################################
			# If Left descendant is a tip
			if (left_desc_nodenum <= length(phy2$tip.label))
				{
				probs_at_anc = relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc,]

				left_branch_downpass_likes = NULL
				right_branch_downpass_likes = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[right_desc_nodenum,]

				uppass_probs_at_bottom_below_tip_results = calc_uppass_probs_new2(probs_ancstate=probs_at_anc, COO_weights_columnar=temp_COO_weights_columnar, numstates=numstates, include_null_range=include_null_range, left_branch_downpass_likes=left_branch_downpass_likes, right_branch_downpass_likes=right_branch_downpass_likes, Rsp_rowsums=NULL)
				
				# The UPPASS probabilities below the tip:
				Lprobs_brbot_below_tip = uppass_probs_at_bottom_below_tip_results$relprobs_just_after_speciation_UPPASS_Left

				if (force_sparse == FALSE)
					{
					# Dense matrix forward exponentiation
					# The UPPASS probabilities AT the tip:
					Lprobs_brtop_AT_tip = Lprobs_brbot_below_tip %*% expokit_dgpadm_Qmat2(times=phy2$edge.length[left_edgenum], Qmat=Qmat, transpose_needed=TRUE)
					} else {
					# Sparse matrix forward exponentiation
					# The UPPASS probabilities AT the tip:
					Lprobs_brtop_AT_tip = calc_prob_forward_onebranch_sparse(relprobs_branch_bottom=Lprobs_brbot_below_tip, branch_length=phy2$edge.length[left_edgenum], tmpQmat_in_REXPOKIT_coo_fmt=tmpQmat_in_REXPOKIT_coo_fmt, coo_n=coo_n, anorm=anorm, check_for_0_rows=FALSE, TRANSPOSE_because_forward=TRUE)
					}				
				
				# Store: branch bottoms
				relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[left_desc_nodenum,] = Lprobs_brbot_below_tip
				# Store: branch tops
				relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[left_desc_nodenum,] = Lprobs_brtop_AT_tip
				} # END if (left_desc_nodenum <= length(phy2$tip.label))
			
			# If Right descendant is a tip
			if (right_desc_nodenum <= length(phy2$tip.label))
				{
				probs_at_anc = relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc,]

				left_branch_downpass_likes = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[left_desc_nodenum,]
				right_branch_downpass_likes = NULL

				uppass_probs_at_bottom_below_tip_results = calc_uppass_probs_new2(probs_ancstate=probs_at_anc, COO_weights_columnar=temp_COO_weights_columnar, numstates=numstates, include_null_range=include_null_range, right_branch_downpass_likes=right_branch_downpass_likes, left_branch_downpass_likes=left_branch_downpass_likes, Rsp_rowsums=NULL)
				
				# The UPPASS probabilities below the tip:
				Rprobs_brbot_below_tip = uppass_probs_at_bottom_below_tip_results$relprobs_just_after_speciation_UPPASS_Right

				if (force_sparse == FALSE)
					{
					# Dense matrix forward exponentiation
					# The UPPASS probabilities AT the tip:
					Rprobs_brtop_AT_tip = Rprobs_brbot_below_tip %*% expokit_dgpadm_Qmat2(times=phy2$edge.length[right_edgenum], Qmat=Qmat, transpose_needed=TRUE)
					} else {
					# Sparse matrix forward exponentiation
					# The UPPASS probabilities AT the tip:
					Rprobs_brtop_AT_tip = calc_prob_forward_onebranch_sparse(relprobs_branch_bottom=Rprobs_brbot_below_tip, branch_length=phy2$edge.length[right_edgenum], tmpQmat_in_REXPOKIT_coo_fmt=tmpQmat_in_REXPOKIT_coo_fmt, coo_n=coo_n, anorm=anorm, check_for_0_rows=FALSE, TRANSPOSE_because_forward=TRUE)
					}				
				
				# Store: branch bottoms
				relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[right_desc_nodenum,] = Rprobs_brbot_below_tip
				# Store: branch tops
				relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[right_desc_nodenum,] = Rprobs_brtop_AT_tip
				} # END if (right_desc_nodenum <= length(phy2$tip.label))				

			##################################################################
			# END finish uppass to tips
			##################################################################
			
			# Zero out impossible states
			if (!is.null(states_allowed_TF))
				{
relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc,][states_allowed_TF==FALSE] = 0
relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[anc,][states_allowed_TF==FALSE] = 0

relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[left_desc_nodenum,][states_allowed_TF==FALSE] = 0
relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[right_desc_nodenum,][states_allowed_TF==FALSE] = 0

relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[left_desc_nodenum,][states_allowed_TF==FALSE] = 0				
relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[right_desc_nodenum,][states_allowed_TF==FALSE] = 0

				} # END if (!is.null(states_allowed_TF))	
			
			
			# 2014-03-20_NJM
			# FIXNODES FOR UPPASS
			# NON-BUG BUG BUG (I don't think anyone has ever addressed UPPASS 
			# calculations with fixed ancestral states)
			
			use_fixnodes_on_uppass = TRUE	# turn off if desired/buggy
			if (use_fixnodes_on_uppass)
				{
				#######################################################
				# If the states/likelihoods have been fixed at a particular node
				# (check top of anc branch)
				#######################################################
				if (!is.null(fixnode))
					{
					# For multiple fixnodes
					if (length(fixnode) > 1)
						{
						# Get the matching node
						TF = (anc == fixnode)
						temporary_fixnode = fixnode[TF]
						temporary_fixlikes = c(fixlikes[TF,])
						} else {
						temporary_fixnode = fixnode
						temporary_fixlikes = c(fixlikes)
						}

				
					if ((length(temporary_fixnode) > 0) && (anc == temporary_fixnode))
						{
						# If the node is fixed, ignore the calculation for this node, and
						# instead use the fixed likelihoods (i.e., the "known" state) for
						# this node.
						# fix the likelihoods of the (NON-NULL) states
						relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc,] = relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc,] * temporary_fixlikes
						}
					} # end if (!is.null(fixnode))
				} # end if (use_fixnodes_on_uppass)

			#######################################################
			# End of UPPASS loop for this ancnode.  Move to next ancnode.
			#######################################################
			} # End uppass loop

		#######################################################
		# End of calculating ancestral states
		#######################################################
		if ((printlevel >= 1) && (stratified == FALSE))
			{
			cat("\nUppass completed for marginal ancestral states estimation!\n", sep="")
			}
		} # End IF calc_ancprobs==TRUE


	if (printlevel >= 2)
		{
		cat("total_loglikelihood:	", total_loglikelihood, "\n", sep="")
		}
	
	if (return_what == "loglike")
		{
		#print(total_loglikelihood)
		return(total_loglikelihood)
		} 

	if (return_what == "nodelikes")
		{
		return(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS)
		}

	if (return_what == "rootprobs")
		{
		return(total_loglikelihood)
		}
	if (return_what == "all")
		{
		calc_loglike_sp_results = list()
		calc_loglike_sp_results$computed_likelihoods_at_each_node = computed_likelihoods_at_each_node
		calc_loglike_sp_results$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS
		calc_loglike_sp_results$condlikes_of_each_state = condlikes_of_each_state

		if (calc_ancprobs == TRUE)
			{
			calc_loglike_sp_results$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS
			calc_loglike_sp_results$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS = relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS
			calc_loglike_sp_results$relative_probs_of_each_state_at_branch_top_AT_node_UPPASS = relative_probs_of_each_state_at_branch_top_AT_node_UPPASS


			
			# Calculate marginal estimates of ancestral states
			
			#######################################################
			# For branch bottoms
			#######################################################
			ML_marginal_prob_each_state_at_branch_bottom_below_node = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS * relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS

			# print
			ML_marginal_prob_each_state_at_branch_bottom_below_node
			rowSums(ML_marginal_prob_each_state_at_branch_bottom_below_node)
			
			ML_marginal_prob_each_state_at_branch_bottom_below_node = ML_marginal_prob_each_state_at_branch_bottom_below_node / rowSums(ML_marginal_prob_each_state_at_branch_bottom_below_node)

			# print
			ML_marginal_prob_each_state_at_branch_bottom_below_node
			rowSums(ML_marginal_prob_each_state_at_branch_bottom_below_node)
			
			#######################################################
			# Probabilities at branch tops
			#######################################################
			if (stratified == FALSE)
				{
				ML_marginal_prob_each_state_at_branch_top_AT_node = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS * relative_probs_of_each_state_at_branch_top_AT_node_UPPASS
				
				# NJM - 2015-01-06: But, at the root, the root probabilities 
				# should NOT be multiplied, that would
				# be downpass * downpass, which 
				# results in focusing probability on the 
				# most-probable downpass state
				# Instead, the probabilities are just the downpass probabilities
				ML_marginal_prob_each_state_at_branch_top_AT_node[rootnode, ] = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[rootnode, ]
				
				} else {
				# In stratified analyses, we DON'T multiply through the
				# uppass probabilities, because these are not even 
				# calculated until a separate step in the 
				# stratified function
				ML_marginal_prob_each_state_at_branch_top_AT_node = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS * 1
				}
			
			# print
			ML_marginal_prob_each_state_at_branch_top_AT_node
			rowSums(ML_marginal_prob_each_state_at_branch_top_AT_node)
			
			ML_marginal_prob_each_state_at_branch_top_AT_node = ML_marginal_prob_each_state_at_branch_top_AT_node / rowSums(ML_marginal_prob_each_state_at_branch_top_AT_node)
			
			# print
			ML_marginal_prob_each_state_at_branch_top_AT_node
			rowSums(ML_marginal_prob_each_state_at_branch_top_AT_node)
			
			# Save them
			calc_loglike_sp_results$ML_marginal_prob_each_state_at_branch_bottom_below_node = ML_marginal_prob_each_state_at_branch_bottom_below_node
			calc_loglike_sp_results$ML_marginal_prob_each_state_at_branch_top_AT_node = ML_marginal_prob_each_state_at_branch_top_AT_node
			
			}

		calc_loglike_sp_results$relative_probs_of_each_state_at_bottom_of_root_branch = relative_probs_of_each_state_at_bottom_of_root_branch
		calc_loglike_sp_results$total_loglikelihood = total_loglikelihood

		class(calc_loglike_sp_results) = "calc_loglike_sp_results"
		return(calc_loglike_sp_results)
		}
	} # END calc_loglike_sp_prebyte

##########################################################
# END calculation of log-likelihood (non-stratified)
# (pre-byte version)
##########################################################































# Byte-compiling functions can speed them up by several times
#
# e.g.: http://blog.revolutionanalytics.com/2011/08/with-byte-compiler-r-214-will-be-even-faster.html
# 
#######################################################
# calc_loglike_sp:
#######################################################
#' Calculate log-likelihood with a transition matrix and speciation events -- byte-compiled
#'
#' This is the workhorse function of \code{\link[BioGeoBEARS]{BioGeoBEARS}}.  It calculates the likelihood of the tip data (the geographic ranges
#' observed at the tips) given a phylogenetic tree, a Q transition matrix specifying the model of range evolution along branches, and a speciation probability
#' matrix specifying the probability of the various possible ancestor-->(Left descendant, Right descendant) range evolution events at phylogenetic nodes/speciation
#' events.
#'
#' This likelihood calculation will be repeated many hundreds or thousands of times in any ML (maximum likelihood) or Bayesian estimation procedure.  Thus, if the 
#' calculation of the log-likelihood of the data under one set of parameter values is too slow, inference takes days or becomes impossible.  However, by using fast
#' matrix exponentiation (package \code{\link[rexpokit]{rexpokit}}) and fast C++ routines for calculating the probabilities of range inheritance scenarios
#' at cladogenesis (package \code{\link[cladoRcpp]{cladoRcpp}}), major speed gains can be achieved. Most of the complexity in the input parameters and the code serves
#' these more rapid alternatives.
#' 
#' However, note that due to the explosion of the geographic range state space with more geographic areas (see \code{\link[cladoRcpp]{numstates_from_numareas}}), 
#' any computational method that explicitly calculates the likelihood of all states will eventually become unusable between 8-20 areas, depending on details.  An
#' alternative method, which is fast for large numbers of areas, is BayArea, by Landis, Matzke, Moore, and Huelsenbeck; see \cite{Landis_Matzke_etal_2013_BayArea}. 
#' However, BayArea does not currently implement cladogenesis models; it only has continuous-time model for evolutionary change along branches.  In effect,
#' this means that the cladogenesis model is sympatric speciation with complete range copying with probability 1.
#' 
#' @param tip_condlikes_of_data_on_each_state A numeric matrix with rows representing tips, and columns representing states/geographic ranges.  The cells
#' give the likelihood of the observation data under the assumption that the tip has that state; typically this means that the known geographic range gets a 
#' '1' and all other states get a 0.
#' @param phy A phylogeny object.  The function converts it to pruningwise order.
#' @param Qmat A Q transition matrix representing the along-branch model for the evolution of geographic range, using parameters \emph{d} (dispersal/range expansion), 
#' \emph{e} (extinction/range contraction/local extirpation), and perhaps others (e.g. distance).  This matrix can be input in either dense or sparse (COO) format, 
#' as specified by \code{input_is_COO}.
#' @param spPmat Default is \code{NULL}; users should usually use \code{spPmat_inputs}.  \code{spPmat} is A numeric matrix representing the probability of each
#' ancestor range-->(Left range, Right range) transition at cladogenesis events.  There are 
#' different ways to represent this matrix.  In the simplest representation, this is just a rectangular matrix with numstates rows (representing the ancestral
#' states) and numstates^2 columns (representing all possible descendant pairs).  Use of this type of matrix is specified by \code{cppSpMethod=1}. It is calculated
#' from a textual speciation matrix (typically \code{spmat} in the code) via \code{\link{symbolic_to_relprob_matrix_sp}}. However, this matrix gets huge and
#' slow for large numbers of states/ranges.  \code{cppSpMethod=2} and \code{cppSpMethod=3} implement successively more efficient and faster 
#' representation and processing of this matrix in COO-like formats.  See \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOprobs}} for the \code{cppSpMethod=2} 
#' method, and \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}} for the \code{cppSpMethod=3} method (the fastest).
#' @param min_branchlength Nodes with branches below this branchlength will not be treated as cladogenesis events; instead, they will be treated as 
#' if an OTU had been sampled from an anagenetic lineage, i.e. as if you had a direct ancestor.  This is useful for putting fossils into the biogeography analysis,
#' when you have fossil species that range through time. (Note: the proper way to obtain such trees, given that most phylogenetic methods force all OTUs to be tips 
#' rather than direct ancestors, is another question subject to active research.  However, one method might be to just set a branch-length cutoff, and treat any
#' branches sufficiently small as direct ancestors.)
#' @param return_what What should be returned to the user? Options are "loglike" (the log-likelihood of the data under the tree, model, and model parameters), 
#' "nodelikes" (the scaled conditional likelihoods at the nodes), "rootprobs" (the relative probability of the geographic ranges/states at the root), or "all"
#' (all of the above in a list).  Typically the user will only want to return "loglike" while doing ML optimization, but then return "all" once the ML parameter
#' values have been found.
#' @param probs_of_states_at_root The prior probability of the states/geographic ranges at the root.  The default, \code{NULL}, effectively means an equal probability
#' for each state (this is also what \code{LAGRANGE} assumes; and running with NULL will reproduce exactly the \code{LAGRANGE} parameter inferences and
#' log-likelihood).
#' @param rootedge  Should the root edge be included in the calculation (i.e., calculate to the bottom of the root), if a root edge is present?  Default \code{FALSE}.
#' @param sparse Should sparse matrix exponentiation be performed?  This should be faster for very large matrices (> 100-200 states), however, the calculations 
#' appear to be less accurate.  The function will transform a dense matrix to COO format (see \code{\link[rexpokit]{mat2coo}}) if necessary according to 
#' the \code{input_is_COO} parameter.
#' @param printlevel If >= 1, various amounts of intermediate output will be printed to screen.  Note: Intermediate outputs from C++ and FORTRAN functions have been
#' commented out, to meet CRAN guidelines.
#' @param use_cpp Should the C++ routines from \code{\link[cladoRcpp]{cladoRcpp}} be used to speed up calculations?  Default \code{TRUE}.
#' @param input_is_COO Is the input Q matrix a sparse, COO-formatted matrix (\code{TRUE}) or a standard dense matrix (\code{FALSE}). Default \code{FALSE}.
#' @param spPmat_inputs A list of parameters so that \code{spPmat} (the speciation transition probability matrix) can be calculated on-the-fly, according
#' to the method in \code{cppSpMethod}.  See example.
#' @param cppSpMethod Three C++ methods from cladoRcpp for calculating and using the cladogenesis probability matrix.  1 is slowest but easiest to understand; 3 is fastest.
#' If \code{spPmat_inputs} is given, the program will generate the appropriate spPmat on-the-fly, and the user does not have to input the full \code{spPmat} manually.
#' @param cluster_already_open If the user wants to distribute the matrix exponentiation calculations from all the branches across a number of processors/nodes on 
#' a cluster, specify the cluster here.  E.g. \code{cluster_already_open = makeCluster(rep("localhost",num_cores_to_use), type = "SOCK")}.  Note: this will work on 
#' most platforms, including Macs running R from command line, but will NOT work on Macs running the R GUI \code{R.app}, because parallel processing functions like
#' \code{MakeCluster} from e.g. \code{library(parallel)} for some reason crash R.app.  The program runs a check for R.app and will just run on 1 node if found. 
#' @param calc_ancprobs Should ancestral state estimation be performed (adds an uppass at the end).
#' @param include_null_range Does the state space include the null range?
#' @return Return whatever is specified by \code{return_what}.
#' @param fixnode If the state at a particular node is going to be fixed (e.g. for ML marginal ancestral states), give the node number. Can now also be a vector of 
#' node numbers, to allow fixing of multiple nodes, but if this is done, \code{fixlikes} must be a matrix with a number of rows equal to \code{length(fixnode)}.
#' @param fixlikes The state likelihoods to be used at the fixed node.  I.e. 1 for the fixed state, and 0 for the others. If only one node is fixed, \code{fixlikes}
#' can be either a vector with length=numstates, or a matrix with 1 row and numstates columns.  For multiple fixed nodes, fixlikes must be a matrix with
#' \code{nrow(fixlikes)==length(fixnode)}.
#' @param stratified Default FALSE. If TRUE, you are running a stratified analysis, in which case uppass probs should be calculated elsewhere.
#' @param states_allowed_TF Default NULL. If user gives a vector of TRUE and FALSE values, these states will be set to 0 likelihood throughout the calculations.
#' @param m This is a vector of rate/weight multipliers for dispersal, conditional on the values of some
#' (non-biogeographical) trait. For example, one might hypothesize that flight/flightlessness effects
#' dispersal probability, and manually put a multiplier of 0.001 on the flightlessness state. Or, 
#' one might attempt to estimate this. The strategy used in cladoRcpp is to expand the default cladogenetic
#' rate matrix by length(m) times. I.e., if \emph{m} is not \code{NULL}, then loop through the values of \emph{m} and apply the 
#' multipliers to \emph{d} (and \emph{j}, and \emph{a}) events. Default is \code{NULL}.
#' @jts_matrix A numtraits x numtraits matrix containing the proportions for 
#' trait transitions during j events. E.g., for a sudden switch from 
#' trait 1 (flight) to trait 2 (flightlessness) during a jump event.
#' @BioGeoBEARS_model_object For printing the parameters in the event that the likelihood calculation produces NaNs. 
#' Ignored if NULL (default)
#' @on_NaN_error Default is NULL, which prints a long stop error if NaNs are hit during the likelihood calculation 
#' (probably due to extreme parameter values and underflow problems). If on_NaN_error is numeric, the value is the 
#' log-likelihood that is returned (e.g. a log-likelihood of -10e100 should be so incredibly low that it 
#' drives any optimizer away from that solution).
#' @export
#' @seealso \code{\link{calc_loglike_sp}}, \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp}}, \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOprobs}}, 
#' \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}}, \code{\link[rexpokit]{mat2coo}}, 
#' \code{\link{rcpp_calc_anclikes_sp_COOweights_faster}}
#' @bibliography /Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#'	 @cite Landis_Matzke_etal_2013_BayArea
#' @note Go BEARS!
#' @note (COO = Coordinate list format for a matrix, see \url{http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29}
#' @author Nicholas Matzke \email{matzke@@berkeley.edu}
#' @examples
#' testval=1
#'
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)





# calc_independent_likelihoods_on_each_branch 
calc_independent_likelihoods_on_each_branch_prebyte <- function(phy2, Qmat, cluster_already_open=NULL, Qmat_is_sparse=FALSE)
	{
	# phy2 must have been reordered, e.g.
	phy2 = reorder(phy2, "pruningwise")
	
	if (Qmat_is_sparse==FALSE)
		{
		# Get probmats for each branch, put into a big array
		# Create empty array to store results
		#independent_likelihoods_on_each_branch = array(0, dim=c(nrow(Qmat), ncol(Qmat), length(phy2$edge.length)))
		
		independent_likelihoods_on_each_branch = vector("list", length(phy2$edge.length))
		tmpmatrix = matrix(data=0, nrow=nrow(Qmat), ncol=ncol(Qmat))
		for (m in 1:length(phy2$edge.length))
			{
			independent_likelihoods_on_each_branch[[m]] = tmpmatrix
			}
		# Calculate the conditional likelihoods for each branch
		# dgexpv NOT ALLOWED when you have a null range state
		# (maybe try very very small values here)

		# clusterApply and other multicore stuff (e.g. doMC) are apparently dangerous on R.app
		if (!is.null(cluster_already_open))
			{
			# 
			if (.Platform$GUI == "AQUA")
				{
				cat("In calc_loglike_sp(), cluster_already_open=", cluster_already_open, " which means you want to calculate likelihoods on branches using a multicore option.\n", sep="")
				cat("But .Platform$GUI='AQUA', which means you are running the Mac GUI R.app version of R.  Parallel multicore functions, e.g. as accessed via \n", sep="")
				cat("library(parallel), are apparently dangerous/will crash R.app (google multicore 'R.app').  So, changing to cluster_already_open=NULL.\n", sep="")
				cluster_already_open=NULL
				} # END if (.Platform$GUI == "AQUA")
			} # END if (!is.null(cluster_already_open))

		
		# Run on the cluster of nodes, if one is open
		# clusterApply etc. appear to NOT work on R.app
		if (!is.null(cluster_already_open))
			{
			# mcmapply
			#library(parallel)
			#independent_likelihoods_on_each_branch = mcmapply(FUN=expokit_dgpadm_Qmat, Qmat=list(Qmat), t=phy2$edge.length, transpose_needed=TRUE, SIMPLIFY="array", mc.cores=Ncores)
			independent_likelihoods_on_each_branch = clusterApply(cl=cluster_already_open, x=phy2$edge.length, fun=expokit_dgpadm_Qmat2, Qmat=Qmat, transpose_needed=TRUE)
			} else {
			# Not parallel processing
			#independent_likelihoods_on_each_branch = mapply(FUN=expokit_dgpadm_Qmat, Qmat=list(Qmat), t=phy2$edge.length, transpose_needed=TRUE, SIMPLIFY="array")
			independent_likelihoods_on_each_branch = mapply_likelihoods(Qmat, phy2, transpose_needed=TRUE)
			#independent_likelihoods_on_each_branch
			}
		} else {
		errortxt = paste("\n\nError in calc_independent_likelihoods_on_each_branch():\n Cannot be used with sparse matrix exponentiation! (Qmat cannot be sparse)\n\n", sep="")
		cat(errortxt)
		stop("\nStopping on error.\n")
		} # END if (Qmat_is_sparse==FALSE)
		
	return(independent_likelihoods_on_each_branch)
	} # END calc_independent_likelihoods_on_each_branch()





# Byte-compile the function
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)

