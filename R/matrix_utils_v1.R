
# Note: experiments with large matrices
# are here:
# 


#sourcedir = "/Dropbox/_njm/"
#source8 = '_matrix_utils_v1.R'
#source(paste(sourcedir, source8, sep=""))



# for e.g. calc_loglike
# sourcedir = '/Dropbox/_njm/'
# source3 = '_R_tree_functions_v1.R'
# source(paste(sourcedir, source3, sep=""))

# Probability matrix fun

# For background, see:
# The Complete Idiot's Guide to the Zen of Likelihood in a Nutshell in Seven Days for Dummies, Unleashed.
# http://www.bioinf.org/molsys/data/idiots.pdf

# Rmat to Qmat
# This is an Rmat where the offdiag AND diag are appropriately scaled, e.g.:
# > Rmat
#      [,1]   [,2]  [,3]       [,4]
# [1,] -2.9  0.300  0.40  0.3000000
# [2,]  0.3 -0.525  0.30  0.4000000
# [3,]  0.4  0.300 -1.25  0.3000000
# [4,]  0.3  0.400  0.30 -0.8333333
#
# More usually, only the offdiag will matter
# > Rmat
#      [,1]   [,2]  [,3]       [,4]
# [1,]  --   0.300  0.40  0.3000000
# [2,]  0.3    --   0.30  0.4000000
# [3,]  0.4  0.300   --   0.3000000
# [4,]  0.3  0.400  0.30     --
#
# basefreqs=rep(1/ncol(Rmat), ncol(Rmat))
# 
Rmat_perfect_to_Qmat <- function(Rmat, basefreqs=rep(1/ncol(Rmat), nrow(Rmat)), check=FALSE)
	{
	# Repeat basefreqs each row
	basefreqs_mat = basefreqs_to_basefreqs_mat(basefreqs)
	
	# element-by-element multiplication
	tmpQmat = basefreqs_mat * Rmat
	
	# Get the sum of the diagonal
	tmp_sumdiag = sumdiag(tmpQmat)
	
	# check if the diag equals the off-diag
	if (check == TRUE)
		{
		tmp_sumoffdiag = sumoffdiag(tmpQmat)
		
		if ( round(abs(tmp_sumdiag),2) != round(tmp_sumoffdiag,2) )
			{
			cat("\n")
			cat("ERROR: Rmat * basefreqs_mat did not produce matrix where\n")
			cat("abs(sum(diag)) equals sum(offdiag).\n")
			cat("\n")
			cat("Rmat:\n")
			print(Rmat)
			cat("\nbasefreqs:\n")
			print(basefreqs)
			cat("\ntmpQmat:\n")
			print(tmpQmat)
			cat("\n")
			cat("sumdiag=", tmp_sumdiag, ", sumoffdiag=", tmp_sumoffdiag, ", diff=", abs(tmp_sumdiag)-tmp_sumoffdiag, "\n")
			cat("\n")
			} else {
			cat("Check passed, with Rmat*basefreqs_mat, abs(sum(diag)) equals sum(offdiag).\n")
			}
		}
	
	# This is sort of a bogus scaling factor (any Rmat matrix would have it)
	scaling_factor = abs(tmp_sumdiag)
	
	Qmat = tmpQmat / scaling_factor
	
	return(Qmat)
	}


# Here, we ignore the diagonal
# Equal base frequencies = default
Rmat_to_Qmat <- function(Rmat, basefreqs=rep(1/ncol(Rmat)) )
	{
	# Repeat basefreqs each row
	basefreqs_mat = basefreqs_to_basefreqs_mat(basefreqs)
	
	# element-by-element multiplication
	tmpQmat = basefreqs_mat * Rmat
	
	# Get the sum of the diagonal
	# This is sort of a bogus scaling factor (any Rmat matrix would have it)
	scaling_factor = sumoffdiag(tmpQmat)
	
	Qmat = tmpQmat / scaling_factor
	
	# Normalize the diagonal to be -rowSum(row_without_diagonal)
	Qmat = normat(relative_matrix=Qmat)
	
	return(Qmat)
	}


# Here, we ignore the diagonal
# Equal base frequencies = default
Rmat_to_Qmat_noScaling <- function(Rmat, basefreqs=rep(1/ncol(Rmat)), use_scaling=FALSE, use_basefreqs=FALSE)
	{
	if (use_basefreqs)
		{
		# Repeat basefreqs each row
		basefreqs_mat = basefreqs_to_basefreqs_mat(basefreqs)
		# element-by-element multiplication
		# (ie, multiply the first COLUMN by the first basefreq, etc.)
		tmpQmat = basefreqs_mat * Rmat
		} else {
		tmpQmat = Rmat
		} # END if (use_basefreqs)
	
	if (use_scaling)
		{
		# Get the sum of the diagonal
		# This is sort of a bogus scaling factor (any Rmat matrix would have it)
		scaling_factor = sumoffdiag(tmpQmat)
	
		Qmat = tmpQmat / scaling_factor
		} else {
		Qmat = tmpQmat
		} # END if (use_scaling)
	
	# Normalize the diagonal to be -rowSum(row_without_diagonal)
	Qmat = normat(relative_matrix=Qmat)
	
	return(Qmat)
	}



Qmat_to_Rmat <- function(Qmat, basefreqs, return_scalefactor=FALSE)
	{
	# Repeat basefreqs each row
	basefreqs_mat = basefreqs_to_basefreqs_mat(basefreqs)
	
	# Get a matrix with 1s on off-diagonal, 0s on diagonal
	offdiags_1 = zeros_diag_ones_offdiag(size=ncol(Qmat))
	
	# Make sure the off-diagonals IN EACH ROW sum to 1
	Rmat = Qmat / basefreqs_mat

	if (return_scalefactor == FALSE)
		{
		Rmat = Rmat / rowSums(Rmat * offdiags_1)
		return(Rmat)
		} else {
		return(rowSums(Rmat * offdiags_1))
		}

	return(Rmat)
	}






# The Pmat should have each row add to 1
Pmat_to_Qmat <- function(Pmat)
	{
	require(expm) 	# for logm, the *matrix* log
	
	# Take the *matrix* log of the Pmat
	logPmat = logm::logm(Pmat)

	# Get a matrix with 1s on off-diagonal, 0s on diagonal
	offdiags_1 = zeros_diag_ones_offdiag(size=ncol(Qmat))
	
	# Scale so the sum of all off-diagonals IN THE WHOLE MATRIX sums to 1
	tmpQmat = logPmat / sum(logPmat * offdiags_1)
	
	Qmat = normat(tmpQmat)
		
	return(Qmat)
	}


# This repeats the basefreqs each row;
# useful for multiplying Rmat
basefreqs_to_basefreqs_mat <- function(basefreqs)
	{
	basefreqs_mat = matrix(data=rep(basefreqs, length(basefreqs)), nrow=length(basefreqs), byrow=TRUE)
	return(basefreqs_mat)	
	}



sumdiag <- function(mat)
	{
	mat_sumdiag = sum(diag(mat))
	return(mat_sumdiag)
	}

sumoffdiag <- function(mat)
	{
	mat_sumoffdiag = sum( mat * zeros_diag_ones_offdiag(size=nrow(mat)), na.rm=TRUE )
	return(mat_sumoffdiag)
	}

zeros_diag_ones_offdiag <- function(size)
	{
	mat = matrix(data=1, nrow=size, ncol=size)
	
	# Put zeros into the diagonal
	diag(mat) = 0
	
	return(mat)
	}


normat <- function(relative_matrix)
	{
	# make the rows sum to zero (not 1), as in here:
	# http://en.wikipedia.org/wiki/Substitution_model
	m = as.matrix(relative_matrix)
	
	diag(m) = 0
	rowsums = rowSums(m)
	diag(m) = -rowsums
	
	return(m)
	}

#' Return a matrix as a list item
#'
#' Puts a matrix into a 1-item list.  This can be useful
#' for adding a matrix to a list of other matrices.
#'
#' @param m, a matrix
#' @return a list with the matrix as the only list item
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples # Example code:
#' 
mat_to_Qmat_as_list <- function(m)
	{
	tmpmat3 = list(m)
	return(tmpmat3)	
	}

#' Put zeros into the diagonal of a matrix
#'
#' The diagonal is ignored here, e.g. in ACE
#'
#' @param m, a matrix
#' @return the matrix with 0 on the diagonal
#' @export
#' @seealso \code{\link{fermat.test}}
#' @references
#' \url{http://www.bioinf.org/molsys/data/idiots.pdf}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
zerosmat <- function(m)
	{
	diag(m) = 0
	return(m)
	}




ones <- function(size)
	{
	# The diag() function creates a 1s matrix (an identity matrix)
	ones_matrix = diag(size)
	
	return(ones_matrix)
	}



############################################################################
#         R FUNCTIONS FOR MATRIX POWERS
#
#         Copied by Jeffrey S. Rosenthal, probability.ca, 2009, from:
#         https://stat.ethz.ch/pipermail/r-help/2007-May/131330.html
#
# The base package of the statistical software R does not seem to include
# built-in functions for computing powers of matrices.  So, I provide
# them here, copied from the above-mentioned web page.
############################################################################

# a simple version, by Ron E. VanNimwegen:
matpow <- function(mat,n){
  ans <- mat
  for ( i in 1:(n-1)){
    ans <- mat %*% ans
  }
  return(ans)
}



# a faster version, by Alberto Monteiro
matpowfast <- function(mat, n)
{
  if (n == 1) return(mat)
  result <- diag(1, ncol(mat))
  while (n > 0) {
    if (n %% 2 != 0) {
      result <- result %*% mat
      n <- n - 1
    }
    mat <- mat %*% mat
    n <- n / 2
  }
  return(result)
}



# Here, each result is normalized so the rows sum to 
# 1, making a relative probability matrix
# Modified for probability matrices (rows sum to 1)
# by Nick Matzke
matpow2 <- function(mat,n){
  ans <- mat
  for ( i in 1:(n-1)){
    ans <- mat %*% ans
    ans = ans/rowSums(ans)
  }
  return(ans)
}

# a faster version, by Alberto Monteiro
# Modified for probability matrices (rows sum to 1)
# by Nick Matzke
matpowfast2 <- function(mat, n)
{
  if (n == 1) return(mat)
  result <- diag(1, ncol(mat))
  while (n > 0) {
    if (n %% 2 != 0) {
      result <- result %*% mat
      result = result/rowSums(result)
      n <- n - 1
    }
    mat <- mat %*% mat
    mat = mat/rowSums(mat)
    n <- n / 2
  }
  return(result)
}



#######################################################
# MATRIX MANIPULATION FUNCTIONS
#######################################################

#' Convert a square matrix to a Q matrix
#'
#' Converts any square matrix containing relative rates
#' to an instantaneous rate matrix (Q matrix) where all
#' off-diagonal elements are positive, and where the
#' diagonal elements are the negative of the sum of the
#' off-diagonal elements in that row.
#'
#' This means that each row sums to 0, which is how an
#' instantaneous rate matrix should be structured.
#' 
#' @param relative_matrix a matrix with the relative rates of transitions
#' @return Q, the instantaneous transition matrix
#' @export
#' @seealso \code{\link{fermat.test}}
#' @references
#' \url{http://www.bioinf.org/molsys/data/idiots.pdf}
#' @bibliography /Dropbox/_njm/__packages/rexpokit_setup/rexpokit_refs.bib
#'   @cite FosterIdiots
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples # Example code:
#' foo(c(1,2,3,4)) # 2.5
mat2q <- function(relative_matrix)
	{
	# make the rows sum to zero (not 1), as in here:
	# http://en.wikipedia.org/wiki/Substitution_model
	Q = relative_matrix
	
	diag(Q) = 0
	rowsums = rowSums(Q)
	diag(Q) = -rowsums
	
	return(Q)
	}


#' Assign a rate to matrix cells symmetrically
#'
#' Assign a rate to matrix cells symmetrically (mirrored across diagonal)
#'
#' @param m, a matrix
#' @return m, a modified matrix
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
assign_rates_sym <- function(m, rownums, rate)
	{
	for (i in rownums)
		{
		for (j in rownums)
			{
			cat(i, j, rate, "\n", sep=" ")
			m[i,j] = rate
			m[j,i] = rate
			}
		}
	return(m)
	}


#' Assign a rate to matrix cells symmetrically
#'
#' Assign a rate to matrix cells asymmetrically (not mirrored across diagonal)
#'
#' @param m, a matrix
#' @return m, a modified matrix
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
assign_rates_asym <- function(m, rownums_from, rownums_to, rate)
	{
	for (i in rownums_from)
		{
		for (j in rownums_to)
			{
			cat(i, j, rate, "\n", sep=" ")
			m[i,j] = rate
			m[j,i] = rate
			}
		}
	return(m)
	}




#' Get the unique categories out of a rate matrix
#'
#' descrip
#'
#' @param Qmat the Q transition matrix
#' @return \code{Quniqs} the unique rate categories
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @seealso \code{\link{counts_of_each_rate}}, \code{\link{counts_of_each_rate_minus_rowcol}}
#' 
uniq_rates <- function(Qmat)
	{
	Qtmp = Qmat
	diag(Qtmp) = NA
	Qlist = c(Qtmp[ !is.na(Qtmp) ])
	
	# return just the unique values
	Quniqs = rev(sort(unique(Qlist)))
	return(Quniqs)
	}


#' Assign rate categories by CRP
#'
#' Assigns rate categories to an input matrix bye CRP (Chinese Restaurant
#' Process, a form of clustering where the number of clusters is estimated
#' and can range from 1 to infinity (or really, the maximum number of clusters
#' is the number of things being clustered).
#'
#' @param nrows number of rows in the matrix (order of the square matrix)
#' @param alpha the clustering (or concentration) parameter of the Chinese Restaurant Process (CRP)
#' @param priorlambda a hyperparameter specifying the mean of the exponential for
#' drawing a rate for each category.  This could itself be drawn from a hyperprior
#' in another function.
#' @return Qmat, the instantaneous transition matrix (Q matrix)
#' @references
#' \url{http://en.wikipedia.org/wiki/Chinese_restaurant_process}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
assign_rate_categories_rchinese <- function(nrows, alpha, priorlambda = 0.1)
	{
	# empty list of categories
	Qcategories = ones(nrows)
	Qmat = Qcategories
	diag(Qcategories) = 0
	numrates = nrows * nrows - nrows

	categories = rchinese(numrates, alpha)
	
	# randomly assign the categories to cells
	Qcategories[Qcategories == 1] = categories
	
	# randomly draw rates for each category
	uniq_categories = sort(unique(categories))
	
	# using a uniform as the prior
	numcats = length(uniq_categories)
	cat("# rate categories = ", numcats, "\n", sep="")
	#uniq_rates = runif(numcats, priormin, priormax)
	uniq_rates = rexp(numcats, priorlambda)
	
	for (i in 1:numcats)
		{
		tmpcat = uniq_categories[i]
		Qmat[Qcategories == tmpcat] = uniq_rates[i]
		}
	
	# normalize the diagonals so that the rows sum to 1:
	Qmat = mat2q(Qmat)
	return(Qmat)
	}

#' Count the number of each category in a Q matrix
#'
#' For each unique category, count the number of cells in that category.
#' 
#' In CRP (Chinese Restaurant Process) terms, this is the number of "patrons"
#' (items or observations) at each "table" (each cluster).
#'
#' @param Qmat, the Q matrix
#' @return Quniqs_counts, counts in the same order as the uniques returned by \code{\link{uniq_rates}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @seealso \code{\link{uniq_rates}}
counts_of_each_rate <- function(Qmat)
	{
	Quniqs = uniq_rates(Qmat)
	Qtmp = Qmat
	diag(Qtmp) = -1
	Qlist = c(Qtmp[Qtmp != -1])
	
	Quniqs_counts = Quniqs
	
	for (i in 1:length(Quniqs))
		{
		rateval = Quniqs[i]
		
		Quniqs_counts[i] = sum(Qlist == rateval)
		}
	return(Quniqs_counts)
	}



#' Count the number of each category, minus the rate of interest
#'
#' count the number of each category, minus the rate of interest (typically the one being changed) which is specified by rowcol, a list of 2 numbers: c(rownum, colnum) 
#'
#' @param Qmat, the Q transition matrix
#' @param rowcol, a list of 2 numbers: c(rownum, colnum) 
#' @return Quniqs_counts, counts in the same order as the uniques returned by \link{uniq_rates}
#' @export
#' @seealso \code{\link{counts_of_each_rate}}
#' @references
#' \url{http://www.bioinf.org/molsys/data/idiots.pdf}
#' @bibliography /Dropbox/_njm/__packages/rexpokit_setup/rexpokit_refs.bib
#'   @cite FosterIdiots
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples # Example code:
#' # Example code
counts_of_each_rate_minus_rowcol <- function(Qmat, rowcol)
	{
	Qtmp = Qmat
	
	# exclude diagonals
	diag(Qtmp) = -1

	# also exclude the rate that is being changed
	Qtmp[rowcol[1], rowcol[2]] = -1
	
	# get the rest and get/count uniques
	Qlist = c(Qtmp[Qtmp != -1])
	
	Quniqs = unique(Qlist)
	Quniqs_counts = Quniqs
	
	for (i in 1:length(Quniqs))
		{
		rateval = Quniqs[i]
		
		Quniqs_counts[i] = sum(Qlist == rateval)
		}
	return(Quniqs_counts)
	}

#' Return the unique rate categories
#'
#' Return the unique rate categories.
#'
#' @param Qmat the Q transition matrix
#' @param rowcol a list of 2 numbers: c(rownum, colnum) 
#' @return \code{Quniqs} the list of unique rate categories
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
uniqs_of_each_rate_minus_rowcol <- function(Qmat, rowcol)
	{
	Qtmp = Qmat
	
	# exclude diagonals
	diag(Qtmp) = -1

	# also exclude the rate that is being changed
	Qtmp[rowcol[1], rowcol[2]] = -1
	
	# get the rest and get/count uniques
	Qlist = c(Qtmp[Qtmp != -1])
	
	Quniqs = unique(Qlist)
	return(Quniqs)
	}


#' Change a randomly-selected rate
#'
#' Randomly selects one of the unique rate categories
#' in a rate matrix and changes it by a draw from a 
#' normal distribution of mean 0 and standard deviation
#' specified by stddev (stddev = standard deviation as a 
#' fraction of the value being modified; a check step disallows
#' moving a parameter below 0)
#'
#' @param Qmat, the Q transition matrix
#' @param stddev, the standard deviation for the Normal(0, stddev*value) draw to 
#' modify the 
#' @return newQmat, the modified Q transition matrix
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
modify_a_rate <- function(Qmat, stddev=0.3)
	{
	newQmat = Qmat
	
	uniq_vals = uniq_rates(Qmat)
	
	val_to_change = sample(uniq_vals, 1)
	
	new_val = val_to_change + rnorm(1, val_to_change, sd=stddev*val_to_change)
	
	# don't let new_val go below 0!!
	while(new_val <= 0)
		{
		new_val = val_to_change + rnorm(1, val_to_change, sd=stddev*val_to_change)
		}
	
	newQmat[Qmat == val_to_change] = new_val
	newQmat = mat2q(newQmat)
	return(newQmat)
	}



#' Put one cell in a new category and assess probability
#'
#' Put one cell in a new category return Q and the "conditional
#' prior probability" of this new category assignment compared
#' to the old assignment
#'
#' @param Qmat, the Q transition matrix
#' @param alpha the clustering (or concentration) parameter of the
#' Chinese Restaurant Process (CRP)
#' @param priorlambda, a hyperparameter specifying the mean of the exponential for
#' drawing a rate for each category.  This could itself be drawn from a hyperprior
#' in another function.
#' @return modQ_vals, a list containing the modified new matrix (modQ_vals$newQmat)
#' and the ratio of the conditional probabilities of the old vs. new 
#' assignment of categories, modQ_vals$ratio_new_to_old_cond_probs
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
modify_rate_matrix_categories <- function(Qmat, alpha, priorlambda = 0.1)
	{
	# Put one cell in a new category
	# return Q and the 
	# "conditional prior probability" of this 
	# new category assignment
	# compared to the old assignment
	printflag=FALSE

	modQ_vals = c()
	modQ_vals$newQmat = Qmat
	
	nrows = nrow(Qmat)
	items = 1:nrows
	
	# Get the current categories
	current_categories = uniq_rates(Qmat)
	
	# randomly sample a (non-identical) row/column pair:
	rowcol = sample(items, 2, replace=FALSE)

	# identify the value in that cell
	cellval_to_change = Qmat[rowcol[1], rowcol[2]]
	
	prflag("modify_rate_matrix_categories #1a", printflag)
	prflag(cellval_to_change, printflag)
	prflag(rowcol, printflag)
	
	
	########################################################
	# calculate the conditional prior of the starting state
	########################################################
	Quniqs = uniqs_of_each_rate_minus_rowcol(Qmat, rowcol)
	Quniqs_counts = counts_of_each_rate_minus_rowcol(Qmat, rowcol)
	
	conditional_probs = diri_probs(Quniqs_counts, alpha)
	
	# Calculate the probability of assignment to the present category 
	matchval = match(cellval_to_change, Quniqs)
	
	# this could return NA, if this rate is unique;
	# In this case, assign the conditional prob
	# of a new category
	if (is.na(matchval))
		{
		matchval = length(conditional_probs)
		}
	
	conditional_prob_oldval = conditional_probs[matchval]

	prflag("modify_rate_matrix_categories #1", printflag)
	
	# Now generate a new value (category) from the prior
	#outcdf = empcdf(conditional_probs)

	# sample one of the options
	numcats = length(conditional_probs)
	
	# don't force a change
	#new_cat = sample((1:numcats)[-matchval], 1, replace=TRUE, prob=conditional_probs[-matchval])
	
	prflag("modify_rate_matrix_categories #2", printflag)
	
	
	# (force a change away from matchval)
	prflag(conditional_probs, printflag)
	prflag(numcats, printflag)	
	prflag(matchval, printflag)
	
	if (length(conditional_probs) <= 2)
		{
		new_cat = 2
		}
	else
		{
		new_cat = sample((1:numcats)[-matchval], 1, replace=TRUE, prob=conditional_probs[-matchval])
		}

	prflag("modify_rate_matrix_categories #3", printflag)


	if (new_cat < numcats)
		{
		newval = Quniqs[new_cat]
		# This could cause an issue if the categories
		# are very focused on 1 group
		# (e.g. you just get the new value back)
		}
	else
		{
		# pick a new value
		newval = rexp(1, priorlambda)
		}
	
	prflag("modify_rate_matrix_categories #3", printflag)

	ratio_new_to_old_cond_probs = ( conditional_probs[new_cat] / conditional_probs[matchval] )
	
	# update the Qmatrix
	modQ_vals$newQmat[rowcol[1], rowcol[2]] = newval
	modQ_vals$newQmat = mat2q(modQ_vals$newQmat)
	
	modQ_vals$ratio_new_to_old_cond_probs = ratio_new_to_old_cond_probs
	
	return(modQ_vals)
	}


#' Calculate conditional probability of clustering under a Dirichlet Process
#'
#' Calculate the conditional probabilities of assignment
#' of a value to the cluster given the current distribution
#' of counts.
#'
#' @param Quniqs_counts, counts in the same order as the uniques returned by \link{uniq_rates}.  Obtained from \code{\link{counts_of_each_rate_minus_rowcol}} or \code{\link{counts_of_each_rate}}
#' @param alpha the clustering (or concentration) parameter of the
#' Chinese Restaurant Process (CRP)
#' @return conditional_probs, the ratio of the conditional probabilities
#' of the old vs. new assignment of categories
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
diri_probs <- function(Quniqs_counts, alpha)
	{
	ntotal = sum(Quniqs_counts)
	n_in_cat = Quniqs_counts
	
	#conditional_probs_oldval = n_in_cat / (ntotal - 1 + alpha)	
	#conditional_probs_newval = alpha / (ntotal - 1 + alpha)

	conditional_probs_oldval = n_in_cat / (ntotal - 0 + alpha)	
	conditional_probs_newval = alpha / (ntotal - 0 + alpha)
	
	conditional_probs = c(conditional_probs_oldval, conditional_probs_newval)
	return(conditional_probs)
	}


#' Empirical CDF
#'
#' Returns an actual empirical CDF (cumulative distribution function) of data \code{x}.
#'
#' @param x, the data
#' @return outcdf, the empirical cumulative distribution function
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
empcdf <- function(x)
	{
	outcdf = x
	for (i in 2:length(x))
		{
		outcdf[i] = x[i] + outcdf[i-1]
		}
	
	return(outcdf)
	}





#' Summary of a series of Q matrices
#'
#' descrip
#'
#' @param Qmat3D a 3D matrix of Q matrices
#' @param burnin the number of MCMC samples to skip as burnin
#' @return \code{summary_mat} summary (mean) of matrix
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
summarize_Qmats <- function(Qmat3D, burnin)
	{
	print("put something in function 'summarize_Qmats'")
	summary_mat = Qmat3D[,,1]
	return(summary_mat)
	}










#######################################################
# TEST SOME MATRIX EXPONENTIALS
#######################################################



#' Run speed tests on various matrix exponentiation algorithms
#'
#' Runs speed tests on various matrix exponentiation algorithms, which are 
#' hard-coded as strings specifying the specific command.
#'
#' The user will have to install each used package with \code{\link{install.packages}}
#' and \code{\link{library}} or \code{\link{require}}, OR comment out the functions they
#' don't want to bother with.
#'
#' @param Qmat the Q transition matrix
#' @param t the time to exponentiate by
#' @return \code{exponentiation_times} a \code{data.frame} of the exponentiation times
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
run_exp_timetests <- function(Qmat, t)
	{
	exp_cmds = c("base::exp(Qmat_test*t)", 
	"msm::MatrixExp(Qmat*t, method='pade')",
	"msm::MatrixExp(Qmat*t, method='series')",
	"expm::expm(Qmat*t)",
	"expm::expm(Qmat*t, method='Higham08.b')",
	"expm::expm(Qmat*t, method='Ward77')",
	"expm::expm(Qmat*t, method='R_Eigen')",
	"expm::expm(Qmat*t, method='Taylor')",
	"Matrix::expm(Qmat*t)",
	"ape::matexpo(Qmat*t)",
	"ctarma::zexpm(Qmat*t)",
	"rexpokit::expokit_dgpadm_Qmat(Qmat_test, t)",
	"rexpokit::expokit_dmexpv_Qmat(Qmat_test, t)")
	
	exponentiation_times = NULL
	
	for (i in 1:length(exp_cmds))
		{
		# these should be right-ish
		expr = paste("Pmat = try(", exp_cmds[i], ")", sep="")
		tmptime = system.time(eval(parse(text=expr)))

		# Check for an error running it
		if (class(Pmat) == "try-error")
			{
			tmptime = rep(NA, times=5)
			}
			
		# Create the row for this command
		tmprow = c(nrow(Qmat), t, exp_cmds[i], tmptime)
		exponentiation_times = rbind(exponentiation_times, tmprow)
		#Pmat = msm::MatrixExp(Qmat*t)
		#(Pmat)
		(apply(Pmat, 2, sum))
		}
	exponentiation_times = adf2(exponentiation_times)
	exponentiation_times = dfnums_to_numeric(exponentiation_times)
	names(exponentiation_times) = c("order", "t", "exp_method", "user.self", "sys.self", "elapsed", "user.child", "sys.child")
	
	exponentiation_times[order(exponentiation_times$elapsed),]
	return(exponentiation_times)
	}





###################################################
# Testing speed of calculation of the exponential
# of matrices
###################################################

#' Run an exponentiation command
#'
#' Exponentiate a rate matrix according to a specific command string
#' that specifies the package, the function, and the options.
#'
#' @param Qmat the Q transition matrix
#' @param t the time to exponentiate by
#' @param cmdstr the exponentiation command to run.  This must have 
#' \code{Pmat = } on the left side, and ideally will specify the 
#' package on the right side, e.g. \code{Pmat = expm::expm(Qmat*t, method='Higham08.b')}
#' @param printflag should details be printed?
#' @return \code{Pmat} the transition probability matrix for time \code{t}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
exp_cmdstr <- function(Qmat, t=1, cmdstr="Pmat = expm::expm(Qmat*t, method='Higham08.b')", printflag=FALSE)
	{
	# print cmdstr if desired
	if (printflag == TRUE)
		{
		cat("exp_cmdstr(Qmat, t, cmdstr) is running: ", cmdstr, "\n", sep="")
		}

	# Evaluate the string
	#print(1)
	# R: evaluate a string without printing
	eval(parse(text = cmdstr))
	#print(2)	
	return(Pmat)
	}



#' A version of exp_cmdstr that accepts lists (for lapply)
#'
#' \code{lapply} takes lists of inputs, so this function 
#' is written to take lists
#'
#' @param i index value
#' @param Qmats_list a list of Q matrices
#' @param t_list a list of times
#' @param cmdstr the exponentiation command to run.  This must have 
#' \code{Pmat = } on the left side, and ideally will specify the 
#' package on the right side, e.g. \code{Pmat = expm::expm(Qmat*t, method='Higham08.b')}
#' @param printflag should details be printed?
#' @return Pmat, the transition probability matrix
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
exp_cmdstr_list <- function(i, Qmats_list, t_list, cmdstr="Pmat = expm::expm(Qmat*t, method='Higham08.b')", printflag=FALSE)
	{
	defaults = '
	Qmats_list = A_list
	i = 1
	'

	#print(3)
	Pmat = exp_cmdstr(Qmats_list[[i]], t_list[i], cmdstr)
	#print(4)
	return(Pmat)
	}





#' Run lots of exponentiations with lapply
#'
#' Takes lists of indices, matrices, and times, and runs an
#' exponentiation on each.
#'
#' @param Qmats_list a list of Q matrices
#' @param t_list a list of times
#' @param cmdstr the exponentiation command to run.  This must have 
#' \code{Pmat = } on the left side, and ideally will specify the 
#' package on the right side, e.g. \code{Pmat = expm::expm(Qmat*t, method='Higham08.b')}
#' @param printflag should details be printed?
#' @return something
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
run_lotsa_exponentiations <- function(Qmats_list, t_list, cmdstr="Pmat = expm::expm(Qmat*t, method='Higham08.b')", printflag=FALSE)
	{
	defaults='
	Qmats_list = NULL
	t_list = NULL
	printflag=FALSE
	Qmat = Qmats_list[[1]]
	'
	#cmdstr="Pmat = expm::expm(Qmat*t, method='Higham08.b')"

	indices = 1:length(Qmats_list)
	Pmat_list = lapply(indices, exp_cmdstr_list, Qmats_list, t_list, cmdstr)
	return(Pmat_list)
	}




#' Run lots of exponentiations and get timing
#'
#' For a series of exponentiation commands, run each through a 
#' series of matrices and times.
#'
#' @param A_list a list of Q matrices
#' @param t_list a list of times to exponentiate by
#' @param tmp_cmdstr_list a list of several exponentiation commands
#' @return results_table commands and the respective times
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 

run_lotsa_exp_functions <- function(A_list, t_list, tmp_cmdstr_list, inputprobs_for_fast=NULL)
	{
	
	# Initialize results table
	results_table = matrix(NA, nrow=length(tmp_cmdstr_list), ncol=4)
	
	
	for (i in 1:length(tmp_cmdstr_list))
		{
		tmp_cmdstr = tmp_cmdstr_list[i]
		cat("run_lotsa_exp_functions() running & timing: '", tmp_cmdstr, "'\n", sep="")
		
		
		# Run the matrix exponentiation, with error trapping
		# set time limit to 1 second
		timelimit = 1
		setTimeLimit(elapsed = timelimit)
		z = try (( time_result = system.time((Pmat_list = run_lotsa_exponentiations(A_list, t_list, cmdstr=tmp_cmdstr))) ))
		setTimeLimit(elapsed = Inf)
				
		# If exponentiation failed		
		if (class(z) == "try-error")
			{
			error_condition = attr(z, "condition")
			
			if (grepl(pattern="elapsed time limit", x=error_condition) == TRUE)
				{
				time_result = c(Inf, Inf, Inf, Inf, Inf)
				cat("TIME-LIMIT ERROR in run_lotsa_exponentiations():\n\n'", z[1], "'\n\n...exceeded ", timelimit, " second limit...\n\n\n", sep="")
				} else {
				time_result = c(NA, NA, NA, NA, NA)
				cat("CAPTURED ERROR in run_lotsa_exponentiations():\n\n'", z[1], "'\n\n...returning NAs to time_result...\n\n\n", sep="")
				}
			}
		
		# print the result
		print(time_result)
		cat("\n\n")
		
		tmp_row = c(tmp_cmdstr, time_result[1], time_result[2], time_result[3])
		
		results_table[i, ] = tmp_row
		}
	
	results_table = adf2(results_table)
	names(results_table) = c("cmd_str", "user", "system", "elapsed")
	return(results_table)
	}


#' Element-wise exponentiation of a SparseM matrix
#'
#' Taking code from \code{\link{http://emdbolker.wikidot.com/seeddisp}}, this
#' was alleged to be a fast way to exponentiate a matrix.  However, this is 
#' exponentiation of each individual element, not matrix exponentiation.
#' 
#' Assumes no true 0s.
#'
#' @param Qmat the Q transition matrix
#' @param t the time to exponentiate by
#' @return x the element-wise exponentiated matrix
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
exp_SparseM <- function(x, t=1)
	{
	defaults = '
	x = Qmat
	t = 1
	'

	# Exponentiate a SparseM matrix
	# code from: http://emdbolker.wikidot.com/seeddisp
	if (inherits(x,"matrix.csr"))
		{
		x@ra <- exp(t*x@ra) ## Sparse matrix HACK
		}
	else if (inherits(x,"matrix.csc"))
		{
		x@ra <- exp(t*x@ra) ## Sparse matrix HACK
		}
	else if (inherits(x,"matrix.ssr"))
		{
		x@ra <- exp(t*x@ra) ## Sparse matrix HACK
		}
	else if (inherits(x,"matrix.ssc"))
		{
		x@ra <- exp(t*x@ra) ## Sparse matrix HACK
		}
	else if (inherits(x,"matrix.coo"))
		{
		x@ra <- exp(t*x@ra) ## Sparse matrix HACK
		}
	else if (is.matrix(x))
		{
		x <- exp(t*x)
		}
	else {
		stop("x should be a matrix")
		}
	return(x)
	}


#' Run lots of matrix exponentiation methods, including SparseM which we shouldn't
#'
#' SparseM and some other methods do element-wise exponentiation, not 
#' matrix exponentiation, so we don't need them.
#'
#' Spam matrices also don't work.
#'
#' @param matrix_dim the order (nrows and ncols) of the Q matrix to simulate
#' @param num_to_make number of matrices to make
#' @param tmp_newvals_str user-input string specifying default value of newvals
#' @param tmp_numvals_to_change number of values to change; default equals the order of the matrix (\code{matrix_dim})
#' @return results_all the time results for each method of exponentiation
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
run_all_exponentiation_methods <- function(matrix_dim=10, num_to_make=1000, tmp_newvals_str = "newvals=0.5", tmp_numvals_to_change=matrix_dim, inputprobs_for_fast=NULL)
	{
	defaults='
	matrix_dim=10
	num_to_make=1000
	tmp_newvals_str = "newvals=0.5"
	tmp_numvals_to_change=matrix_dim
	inputprobs_for_fast=NULL
	'
	
	t_list = rep(1, length(A_list))
	
	A_list = make_lotsa_sparse_matrices(num_to_make, matrix_dim, newvals_str=tmp_newvals_str,  numvals_to_change=matrix_dim)
	
	inputprobs_for_fast = runif(n=matrix_dim, min=0, max=1)
	
	# Standard matrix format
	tmp_cmdstr_list = c("msm::MatrixExp(Qmat*t, method='pade')",
	"msm::MatrixExp(Qmat*t, method='series')",
	"ape::matexpo(Qmat*t)",
	"expm::expm(Qmat*t)",
	"expm::expm(Qmat*t, method='Higham08.b')",
	"expm::expm(Qmat*t, method='Ward77')",
	"expm::expm(Qmat*t, method='R_Eigen')",
	"expm::expm(Qmat*t, method='Taylor')",
	"Matrix::expm(Qmat*t)",
	"ctarma::zexpm(Qmat*t)",
	"rexpokit::expokit_dgpadm_Qmat(Qmat_test, t)",
	"rexpokit::expokit_dmexpv_Qmat(Qmat_test, t)",
	"rexpokit::expokit_dmexpv_Qmat(Qmat_test, t, inputprobs_for_fast=inputprobs_for_fast)")
	
	# Old
	# Standard matrix format
	# tmp_cmdstr_list = c("Pmat = msm::MatrixExp(Qmat*t)", 
	# 	"Pmat = ape::matexpo(Qmat*t)", 
	# 	"Pmat = expm::expm(Qmat*t, method='Higham08.b')", 
	# 	"Pmat = expm::expm(Qmat*t, method='Ward77')", 
	# 	"Pmat = expm::expm(Qmat*t, method='R_Eigen')", 
	# 	"Pmat = expm::expm(Qmat*t, method='Taylor')", 
	# 	"Pmat = Matrix::expm(Qmat*t)")



	results_std = run_lotsa_exp_functions(A_list, t_list, tmp_cmdstr_list, inputprobs_for_fast=inputprobs_for_fast)

	
	# SPAM matrix format
	# make a list of SPAM matrices
# 	Aspam_list = lapply(A_list, matrix_to_spam)
# 	
# 	tmp_cmdstr_list = c("Pmat = expm::expm(Qmat*t, method='Higham08.b')") 
# 	results_spam = run_lotsa_exp_functions(A_list, t_list, tmp_cmdstr_list)
# 	
# 	
# 	
# 	
# 	# SparseM matrix format
# 	
# 	# exponentiate a SparseM matrix
# 	# code from: http://emdbolker.wikidot.com/seeddisp
# 	z = as.matrix.csr(A_list[[1]])
# 	exp_SparseM(z)
# 	
# 	
# 	# make a list of SPAM matrices
# 	A_SparseM_csr_list = lapply(A_list, as.matrix.csr)
# 	A_SparseM_csc_list = lapply(A_list, as.matrix.csc)
# 	
# 	# for Symmetric matrices only
# 	# A_SparseM_ssr_list = lapply(A_list, as.matrix.ssr)
# 	# A_SparseM_ssc_list = lapply(A_list, as.matrix.ssc)
# 	A_SparseM_coo_list = lapply(A_list, as.matrix.coo)
# 	
# 	tmp_cmdstr_list = c("Pmat = exp_SparseM(Qmat, t)") 
# 	results_SparseM_csr = run_lotsa_exp_functions(A_SparseM_csr_list, t_list, tmp_cmdstr_list)
# 	
# 	tmp_cmdstr_list = c("Pmat = exp_SparseM(Qmat, t)") 
# 	results_SparseM_csc = run_lotsa_exp_functions(A_SparseM_csc_list, t_list, tmp_cmdstr_list)
# 	
# 	#tmp_cmdstr_list = c("Pmat = exp_SparseM(Qmat, t)") 
# 	#results_SparseM_ssr = run_lotsa_exp_functions(A_SparseM_ssr_list, t_list, tmp_cmdstr_list)
# 	
# 	#tmp_cmdstr_list = c("Pmat = exp_SparseM(Qmat, t)") 
# 	#results_SparseM_ssc = run_lotsa_exp_functions(A_SparseM_ssc_list, t_list, tmp_cmdstr_list)
# 	
# 	tmp_cmdstr_list = c("Pmat = exp_SparseM(Qmat, t)") 
# 	results_SparseM_coo = run_lotsa_exp_functions(A_SparseM_coo_list, t_list, tmp_cmdstr_list)
# 
# 	
# 	
# 	results_all = rbind(results_std, results_spam, results_SparseM_csr, results_SparseM_csc, results_SparseM_coo)

	results_all = results_std
	
	results_all$matrix_dim = matrix_dim
	results_all$num_to_make = num_to_make
	results_all$newvals_str = tmp_newvals_str
	results_all$numvals_nonzero = tmp_numvals_to_change
	
	return(results_all)
	}






###########################################
# Make matrices
###########################################


#' Convert a matrix to spam format
#'
#' spam is one of the sparse matrix formats, however the spam
#' package does not include matrix exponentiation specifically
#'
#' @param A a matrix
#' @return \code{Aspam}, a matrix in spam format
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
matrix_to_spam <- function(A)
	{
	Aspam = spam(as.numeric(A))
	return(Aspam)
	}

#' Return the number of cells in the triangle, excluding the diagonal
#'
#' returns matrix_dim * (matrix_dim + 1)/2 - matrix_dim
#'
#' @param matrix_dim the order (nrows or ncols) of the matrix
#' @return num_tri_cells the number of cells in the triangle (excluding the diagonal)
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
size_half_tri_nodiag <- function(matrix_dim)
	{
	num_tri_cells = matrix_dim * (matrix_dim + 1)/2 - matrix_dim
	return(num_tri_cells)
	}


#' Make a sparse matrix
#'
#' Makes a sparse matrix with all zeros, except for a user-specified number
#' of cells (numvals_to_change, default is set to the order of the matrix)
#' which are set to eval(parse(text=newvals_str)).
#'
#' @param matrix_dim the order of the matrix (order = number of rows/columns)
#' @param newvals_str the string which will input new values via eval(parse(newvals_str))
#' @param numvals_to_change the number of values to change; default is the order of the matrix
#' @return \code{A} a sparse matrix
#' @export
#' @seealso \code{\link{make_sparse_matrix_by_i}}, \code{\link{make_lotsa_sparse_matrices}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
make_sparse_matrix <- function(matrix_dim, newvals_str="newvals=0.5", numvals_to_change=matrix_dim)
	{
	# Generate "newvals" from the string
	# newvals = 
	eval(parse(text = newvals_str))
	
	
	# make one or zero non-diagonals to be 0.5
	# numvals_to_change = matrix_dim
	num_tri_cells = size_half_tri_nodiag(matrix_dim)
	
	k = sample(1:num_tri_cells, numvals_to_change, replace=FALSE)
	i = ceiling(k/matrix_dim) + 0
	j = (k - (i-1)*matrix_dim) + 0
	
	
	# blank matrix
	A <- matrix(0, nrow=matrix_dim, ncol=matrix_dim)
	# (A <- spam(0, matrix_dim, matrix_dim))
	
	# make the i,j cells 0.5
	(A[cbind(i, j)] <- newvals)
	
	return(A)
	}

#' A lapply-usable version of make_sparse_matrix
#'
#' This function allows the use of make_sparse_matrix by \code{lapply}
#'
#' @param i an index value
#' @param matrix_dim the order of the matrix (order = number of rows/columns)
#' @param newvals_str the string which will input new values via eval(parse(newvals_str))
#' @param numvals_to_change the number of values to change; default is the order of the matrix
#' @return \code{A} a sparse matrix
#' @export
#' @seealso \code{\link{make_sparse_matrix}}, \code{\link{make_lotsa_sparse_matrices}}#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
make_sparse_matrix_by_i <- function(i, matrix_dim, newvals_str="newvals=0.5", numvals_to_change=matrix_dim)
	{
	A = make_sparse_matrix(matrix_dim, newvals_str, numvals_to_change)
	return(A)
	}




#' Runs make_sparse_matrix a lot of times
#'
#' Combines \code{lapply} and \code{make_sparse_matrix_by_i} to make a list
#' of sparse matrices for speed testing purposes.
#'
#' @param num_to_make number of sparse matrices to make
#' @param matrix_dim the order of the matrix (order = number of rows/columns)
#' @param newvals_str the string which will input new values via eval(parse(newvals_str))
#' @param numvals_to_change the number of values to change; default is the order of the matrix
#' @return \code{A_list} a list of sparse matrices
#' @export
#' @seealso \code{\link{make_sparse_matrix}}, \code{\link{make_sparse_matrix_by_i}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
make_lotsa_sparse_matrices <- function(num_to_make, matrix_dim, newvals_str="newvals=0.5", numvals_to_change=matrix_dim)
	{
	defaults = '
	num_to_make
	matrix_dim
	newvals_str="newvals=0.5"
	numvals_to_change=matrix_dim
	'

	# Make the appropriate number of indices
	indices = 1:num_to_make
	
	# Initialize list of matrices
	A_list = lapply(rep.int(0, time=num_to_make), matrix, nrow=matrix_dim, ncol=matrix_dim)
	
	# Make a list of matrices
	A_list = lapply(indices, make_sparse_matrix_by_i, matrix_dim, newvals_str)
	return(A_list)
	}




#' Get the order of a big matrix produced by combining several small matrices
#'
#' When combining several transition matrices (e.g. geography, elevation, and climate),
#' it is useful to know the size of the output matrix, which will inflate very rapidly.
#'
#' The new matrix order is the product of the order of each input matrix. E.g., a 
#' 2x2 matrix for elevation (high/low), a 2x2 matrix for precipitation (wet/dry), and 
#' a 4x4 matrix (for 4 islands) will produce a matrix of order=2x2x4=16, i.e. 16 rows
#' and columns.  This would soon become impractical for exponentiation, but certain
#' cases may produce sparse matrices.
#'
#' @param list_of_matrices a list of input matrices
#' @return \code{matrix_dim} the order of the big matrix
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
get_bigmatrix_dim <- function(list_of_matrices)
	{
	matrix_dim = 1
	for (i in 1:length(list_of_matrices))
		{
		matrix_dim = matrix_dim * nrow(list_of_matrices[[i]])
		}
	return(matrix_dim)
	}


#' Combine small matrices into one large compound matrix
#'
#' For certain purposes it is interesting to combine several transition matrices
#' (e.g. geography, elevation, and climate) into one large matrix.  However, 
#' the size of this matrix will inflate very rapidly.
#'
#' The new matrix order is the product of the order of each input matrix. E.g., a 
#' 2x2 matrix for elevation (high/low), a 2x2 matrix for precipitation (wet/dry), and 
#' a 4x4 matrix (for 4 islands) will produce a matrix of order=2x2x4=16, i.e. 16 rows
#' and columns.  This would soon become impractical for exponentiation, but certain
#' cases may produce sparse matrices.
#'
#' @param list_of_matrices the small input matrices to combine
#' @return \code{tmpmat} the large output matrix
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
combine_matrices_into_supermatrix <- function(list_of_matrices)
	{
	# Get the matrix dimensions
	matrix_dim = 1
	for (i in 1:length(list_of_matrices))
		{
		matrix_dim = matrix_dim * nrow(list_of_matrices[[i]])
		}
	
	# Set up the empty matrix
	tmpmat = matrix(data=1, nrow=matrix_dim, ncol=matrix_dim)
	tmpmat_names = matrix(data=NA, nrow=matrix_dim, ncol=matrix_dim)
	
	# recursive algorithm
	tmpmat = nesting_mats(tmpmat, list_of_matrices, level=0, maxlevel=length(list_of_matrices))
	
	return(tmpmat)
	}


#' recursively nest matrices to construct a large compound matrix
#'
#' recursively fill in the cells of a large output matrix produced
#' by combining an input matrix.
#'
#' This is a dynamic function that calls itself until \code{maxlevel} is reached.
#'
#' @param tmpmat a blank matrix with the correct size as the output matrix
#' @param list_of_matrices the small input matrices to combine
#' @param level the current level of nesting (this is kept track of with each level of recursion
#' @param maxlevel recursion stops when \code{maxlevel} is reached; typically \code{maxlevel}
#' is the number of input matrices
#' @return \code{tmpmat} the large output matrix
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
nesting_mats <- function(tmpmat, list_of_matrices, level=0, maxlevel=3)
	{
	if (level >= maxlevel)
		{
		return(tmpmat)
		}

	level = level + 1
	
	# Figure out the steps for each level
	nrows_of_matrices_left = figure_out_step_val(level, list_of_matrices)
	
	byval = nrows_of_matrices_left

	# If you are not at the last small matrix
	starts = seq(1, nrow(tmpmat), by=byval)
	ends = starts-1 + byval
	#cat("level=", level, ", step value=", byval, "\n", sep="")
	#cat("starts", starts, "\n", sep=", ")
	#cat("ends", ends, "\n", sep=", ")
	
	# edit the matrix
	tmpsmall_mat = list_of_matrices[[level]]
	nstates_small = nrow(tmpsmall_mat)

	# Oscillate through the number of states in the small matrix
	statenum = 0
	for (m in 1:nstates_small)
		{
		for (n in 1:nstates_small)
			{
			indices_to_use_rows = seq(m, length(starts), by=nstates_small)
			indices_to_use_cols = seq(n, length(starts), by=nstates_small)

			#print(indices_to_use_rows)			
			#print(indices_to_use_cols)
			
			tmpmat_starts_to_use_rows = starts[indices_to_use_rows]
			tmpmat_ends_to_use_rows = ends[indices_to_use_rows]

			tmpmat_starts_to_use_cols = starts[indices_to_use_cols]
			tmpmat_ends_to_use_cols = ends[indices_to_use_cols]
			
			# update the big matrix by multiplying the small matrix at appropriate cells
			row_indices = make_list_of_indices_from_2_lists(tmpmat_starts_to_use_rows, tmpmat_ends_to_use_rows)

			col_indices = make_list_of_indices_from_2_lists(tmpmat_starts_to_use_cols, tmpmat_ends_to_use_cols)
			
			#cat("m=", m, ", rows=", row_indices, "\n", sep="")
			#cat("n=", n, ", rows=", col_indices, "\n", sep="")
			#print(tmpsmall_mat)
			
			tmpmat[row_indices, col_indices] = tmpsmall_mat[m,n] * tmpmat[row_indices, col_indices]
			
			#print(tmpmat)
			}		
		}
	
	# Go one level deeper
	tmpmat = nesting_mats(tmpmat, list_of_matrices, level)
	
	return(tmpmat)
	}
	











#' Combine matrices using the variable names
#'
#' Fill in a supermatrix with variable names in the cells, by 
#' combining input matrices with names.
#'
#' @param list_of_matrices_w_variable_names_w_variable_names A list of input matrices,
#' where the cells have been filled in with variable names which will be
#' multiplied together in the supermatrix.
#' @return \code{tmpmat} A supermatrix with variable names in the cells
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
combine_matrices_into_supermatrix_w_names <- function(list_of_matrices_w_variable_names_w_variable_names)
	{
	# Get the matrix dimensions
	matrix_dim = 1
	for (i in 1:length(list_of_matrices_w_variable_names))
		{
		matrix_dim = matrix_dim * nrow(list_of_matrices_w_variable_names[[i]])
		}
	
	# Set up the empty matrix
	tmpmat = matrix(data="1", nrow=matrix_dim, ncol=matrix_dim)
	tmpmat_names = matrix(data=NA, nrow=matrix_dim, ncol=matrix_dim)
	
	# recursive algorithm
	tmpmat = nesting_mats_w_names(tmpmat, list_of_matrices_w_variable_names, level=0, maxlevel=length(list_of_matrices_w_variable_names))
	
	return(tmpmat)
	}

#' Recursive (dynamic) function for filling in the variable names of a supermatrix
#'
#' This is a dynamic function which calls itself in order to recursively fill in
#' variable names in a combined (or compound) supermatrix made by combining
#' several small matrices.
#'
#' @param tmpmat a blank matrix with the correct size as the output matrix
#' @param list_of_matrices_w_variable_names the small input matrices to combine (which
#' have variable names, not numbers
#' @param level the current level of nesting (this is kept track of with each level of recursion
#' @param maxlevel recursion stops when \code{maxlevel} is reached; typically \code{maxlevel}
#' is the number of input matrices
#' @return \code{tmpmat} the large output matrix
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
nesting_mats_w_names <- function(tmpmat, list_of_matrices_w_variable_names, level=0, maxlevel=3)
	{
	if (level >= maxlevel)
		{
		return(tmpmat)
		}

	level = level + 1
	
	# Figure out the steps for each level
	nrows_of_matrices_left = figure_out_step_val(level, list_of_matrices_w_variable_names)
	
	byval = nrows_of_matrices_left

	# If you are not at the last small matrix
	starts = seq(1, nrow(tmpmat), by=byval)
	ends = starts-1 + byval
	#cat("level=", level, ", step value=", byval, "\n", sep="")
	#cat("starts", starts, "\n", sep=", ")
	#cat("ends", ends, "\n", sep=", ")
	
	# edit the matrix
	tmpsmall_mat = list_of_matrices_w_variable_names[[level]]
	nstates_small = nrow(tmpsmall_mat)

	# Oscillate through the number of states in the small matrix
	statenum = 0
	for (m in 1:nstates_small)
		{
		for (n in 1:nstates_small)
			{
			indices_to_use_rows = seq(m, length(starts), by=nstates_small)
			indices_to_use_cols = seq(n, length(starts), by=nstates_small)

			#print(indices_to_use_rows)			
			#print(indices_to_use_cols)
			
			tmpmat_starts_to_use_rows = starts[indices_to_use_rows]
			tmpmat_ends_to_use_rows = ends[indices_to_use_rows]

			tmpmat_starts_to_use_cols = starts[indices_to_use_cols]
			tmpmat_ends_to_use_cols = ends[indices_to_use_cols]
			
			# update the big matrix by multiplying the small matrix at appropriate cells
			row_indices = make_list_of_indices_from_2_lists(tmpmat_starts_to_use_rows, tmpmat_ends_to_use_rows)

			col_indices = make_list_of_indices_from_2_lists(tmpmat_starts_to_use_cols, tmpmat_ends_to_use_cols)
			
			#cat("m=", m, ", rows=", row_indices, "\n", sep="")
			#cat("n=", n, ", rows=", col_indices, "\n", sep="")
			#print(tmpsmall_mat)
			
			tmpmat[row_indices, col_indices] = sapply(tmpmat[row_indices, col_indices], paste, "*", tmpsmall_mat[m,n], sep="")
			
			#print(tmpmat)
			}		
		}
	
	# Go one level deeper
	tmpmat = nesting_mats_w_names(tmpmat, list_of_matrices_w_variable_names, level, maxlevel=length(list_of_matrices_w_variable_names))

	
	return(tmpmat)
	}
	














#' Take 2 lists of indices, make a single list of indices
#'
#' Given e.g. a list of start indices and a list of end indices
#' this function produces one big list, e.g.:
#'
#' c(start1:end1, start2:end2, start3:end3...etc...)
#'
#' @param list1 list of starting indices
#' @param list2 list of ending indices
#' @return \code{newlist} an output list of indices
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
make_list_of_indices_from_2_lists <- function(list1, list2)
	{
	defaults = '
	list1=c(1,3,5,7)
	list2 = list1+1
	'
	
	# initialize a list
	newlist = NULL
	
	# error check
	if (length(list1) != length(list2))
		{
		stop("make_list_of_indices_from_2_lists() says: ERROR, lists must be of the same length...")
		}
	
	# Go through the lists and make the new list
	for (i in 1:length(list1))
		{
		# get the new indices
		new_indices = list1[i] : list2[i]
		
		# add them to the list
		newlist = c(newlist, new_indices)
		}
	
	return(newlist)
	}



#' Get the step size for recursively filling in a supermatrix
#'
#' When recursively filling in a matrix, each input matrix 
#' will be distributed differently in the output matrix. 
#' This function figures out the step size based on the level of recursion.
#'
#' E.g., when combining:
#'
#' 0 A and 0 B
#' A 0     B 0
#'
#' These basically get expanded to:
#' 
#' 0 0 A A 
#' 0 0 A A
#' A A 0 0
#' A A 0 0
#'
#' ...and...
#'
#' 0 B 0 B
#' B 0 B 0
#' 0 B 0 B
#' B 0 B 0
#'
#' ...resulting in...
#'
#' 00 0B A0 AB
#' 0B 00 AB A0
#' A0 AB 00 0B
#' AB A0 0B 00
#' 
#'
#' @param list_of_matrices the small input matrices to combine
#' @param level the current level of nesting (this is kept track of with each level of recursion
#' @return \code{nrows_of_matrices_left} the number of rows left to subdivide
#' with more matrices; when this equals 1, the recursion stops
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
figure_out_step_val <- function(level, list_of_matrices)
	{
	# Figure out the step value (byval)
	nrows_of_matrices_left = 1

	# if you're on the last small matrix
	if (level >= length(list_of_matrices))
		{
		nrows_of_matrices_left = 1
		return(nrows_of_matrices_left)
		}
	
	# Otherwise
	for (j in ((level+1):length(list_of_matrices)))
		{
		nrows_of_matrices_left = nrows_of_matrices_left * nrow(list_of_matrices[[j]])
		}
	return(nrows_of_matrices_left)
	
	}



#' Combine 3 small matrices.
#'
#' Unnecessarily elaborate function to combine 3 relative-probability
#' transition matricies into one big one.
#'
#' @param dmat a small input matrix, e.g. transition probability based on distance
#' @param emat a small input matrix, e.g. transition probability based on elevation
#' @param cmat a small input matrix, e.g. transition probability based on climate
#' @return \code{result_list} a list containing \code{result_list$tmpmat}, the output supermatrix with numbers, and \code{result_list$tmpmat_names}, the output supermatrix with variable names
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
combine_3matrices <- function(dmat, emat, cmat)
	{
	dnum = nrow(dmat)
	enum = nrow(emat)
	cnum = nrow(cmat)

	matrix_dim = dnum * enum * cnum
	tmpmat = matrix(data=1, nrow=matrix_dim, ncol=matrix_dim)
	tmpmat_names = matrix(data=NA, nrow=matrix_dim, ncol=matrix_dim)
	
	# populate with ds
	size_non_d = enum*cnum
	for (i in 1:dnum)
		{
		start1 = 1+(i-1)*size_non_d 
		end1 = start1 + size_non_d - 1
		indices_to_change1 = start1:end1
		
		for (j in 1:dnum)
			{
			start2 = 1+(j-1)*size_non_d 
			end2 = start2 + size_non_d - 1
			indices_to_change2 = start2:end2
			
			cat(indices_to_change1, " -- ", indices_to_change2, "\n", sep="")
			
			tmpmat[indices_to_change1, indices_to_change2] = tmpmat[indices_to_change1, indices_to_change2] * dmat[i, j]
			
			tmpmat_names[indices_to_change1, indices_to_change2] = paste("d[", i, ",", j, "]*", sep="")
			

			size_non_e = cnum
			for (ei in 1:enum)
				{
				estart1 = start1+(ei-1)*size_non_e 
				eend1 = estart1 + size_non_e - 1
				eindices_to_change1 = estart1:eend1
								
				for (ej in 1:enum)
					{
					estart2 = start2+(ej-1)*size_non_e 
					eend2 = estart2 + size_non_e - 1
					eindices_to_change2 = estart2:eend2
					
					cat(eindices_to_change1, " -- ", eindices_to_change2, "\n", sep="")
					
					tmpmat[eindices_to_change1, eindices_to_change2] = tmpmat[eindices_to_change1, eindices_to_change2] * emat[ei, ej]
					
					tmpmat_names[eindices_to_change1, eindices_to_change2] = paste(tmpmat_names[eindices_to_change1, eindices_to_change2], "e[", ei, ",", ej, "]*", sep="")
					
					
					
					size_non_c = 1
					for (ci in 1:cnum)
						{
						cstart1 = estart1+(ci-1)*size_non_c 
						cend1 = cstart1 + size_non_c - 1
						cindices_to_change1 = cstart1:cend1
										
						for (cj in 1:cnum)
							{
							cstart2 = estart2+(cj-1)*size_non_c 
							cend2 = cstart2 + size_non_c - 1
							cindices_to_change2 = cstart2:cend2
							
							cat(cindices_to_change1, " -- ", cindices_to_change2, "\n", sep="")
							
							tmpmat[cindices_to_change1, cindices_to_change2] = tmpmat[cindices_to_change1, cindices_to_change2] * cmat[ci, cj]
							
							tmpmat_names[cindices_to_change1, cindices_to_change2] = paste(tmpmat_names[cindices_to_change1, cindices_to_change2], "c[", ci, ",", cj, "]", sep="")
							}
						}	
					}
				}		
			}
		}
	result_list = list(tmpmat, tmpmat_names)
	return(result_list)
	}



#' Combine 3 small matrices, second version
#'
#' Unnecessarily elaborate function to combine 3 relative-probability
#' transition matricies into one big one.
#'
#' @param dmat a small input matrix, e.g. transition probability based on distance
#' @param emat a small input matrix, e.g. transition probability based on elevation
#' @param cmat a small input matrix, e.g. transition probability based on climate
#' @return \code{result_list} a list containing \code{result_list$tmpmat}, the output supermatrix with numbers, and \code{result_list$tmpmat_names}, the output supermatrix with variable names
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
combine_3matrices2 <- function(dmat, emat, cmat)
	{
	# Unnecessarily elaborate function to combine 3 relative-prob transition 
	# matricies into one big one.
	
	dnum = nrow(dmat)
	enum = nrow(emat)
	cnum = nrow(cmat)

	matrix_dim = dnum * enum * cnum
	tmpmat = matrix(data=1, nrow=matrix_dim, ncol=matrix_dim)
	tmpmat_names = matrix(data=NA, nrow=matrix_dim, ncol=matrix_dim)
	
	# populate with ds
	size_non_d = enum*cnum
	for (i in 1:dnum)
		{
		start1 = 1+(i-1)*size_non_d 
		end1 = start1 + size_non_d - 1
		indices_to_change1 = start1:end1
		
		for (j in 1:dnum)
			{
			start2 = 1+(j-1)*size_non_d 
			end2 = start2 + size_non_d - 1
			indices_to_change2 = start2:end2
			
			cat(indices_to_change1, " -- ", indices_to_change2, "\n", sep="")
			
			tmpmat[indices_to_change1, indices_to_change2] = tmpmat[indices_to_change1, indices_to_change2] * dmat[i, j]
			
			tmpmat_names[indices_to_change1, indices_to_change2] = paste("d*", sep="")
			

			size_non_e = cnum
			for (ei in 1:enum)
				{
				estart1 = start1+(ei-1)*size_non_e 
				eend1 = estart1 + size_non_e - 1
				eindices_to_change1 = estart1:eend1
								
				for (ej in 1:enum)
					{
					estart2 = start2+(ej-1)*size_non_e 
					eend2 = estart2 + size_non_e - 1
					eindices_to_change2 = estart2:eend2
					
					cat(eindices_to_change1, " -- ", eindices_to_change2, "\n", sep="")
					
					tmpmat[eindices_to_change1, eindices_to_change2] = tmpmat[eindices_to_change1, eindices_to_change2] * emat[ei, ej]
					
					tmpmat_names[eindices_to_change1, eindices_to_change2] = paste(tmpmat_names[eindices_to_change1, eindices_to_change2], "e*", sep="")
					
					
					
					size_non_c = 1
					for (ci in 1:cnum)
						{
						cstart1 = estart1+(ci-1)*size_non_c 
						cend1 = cstart1 + size_non_c - 1
						cindices_to_change1 = cstart1:cend1
										
						for (cj in 1:cnum)
							{
							cstart2 = estart2+(cj-1)*size_non_c 
							cend2 = cstart2 + size_non_c - 1
							cindices_to_change2 = cstart2:cend2
							
							cat(cindices_to_change1, " -- ", cindices_to_change2, "\n", sep="")
							
							tmpmat[cindices_to_change1, cindices_to_change2] = tmpmat[cindices_to_change1, cindices_to_change2] * cmat[ci, cj]
							
							tmpmat_names[cindices_to_change1, cindices_to_change2] = paste(tmpmat_names[cindices_to_change1, cindices_to_change2], "c", sep="")
							}
						}	
					}
				}		
			}
		}
	result_list = list(tmpmat, tmpmat_names)
	return(result_list)
	}













#' JUNK - Test speed of exponentiation - spam
#'
#' Given a number of areas (\code{n}), make a large spam matrix.
#'
#' Pointless probably.
#'
#' @param n number of areas
#' @param numtimes number of exponentiations to run
#' @param expwhich option "A" (spam matrix, not real exponentiation really) or "B" (standard matrix)
#' @return \code{expAB}, the exponentiated matrix
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 

testexp <- function(n, numtimes, expwhich="A")
	{
	# test timing of matrix exponentiation
	#n = 30
	number_of_areas = (2^n) - 1
	cat("number_of_areas = ", number_of_areas, "\n", sep="")
	
	matrix_dim = number_of_areas
	
	for (l in 1:numtimes)
		{
		# make one or zero non-diagonals to be 0.5
		k = sample(1:matrix_dim, matrix_dim, replace=FALSE)
		i = ceiling(k/matrix_dim) + 0
		j = (k - (i-1)*matrix_dim) + 0
		
		
		# blank matrix
		(A <- spam(0, matrix_dim, matrix_dim))
		
		# make the i,j cells 0.5
		(A[cbind(i, j)] <- rep(.5, length(i)))
		
		# t(A) reflects across the matrix
		# A is original
		# diag.spam makes the diagonals = 1
		(A <- t(A) + A + diag.spam(matrix_dim))
		
		
		
		#stA = system.time( (expA = exp(A)) )
		if (expwhich == "A")
			{
			expA = exp(A)
			#return(expA)
			expAB = expA
			} else {
			#stB = system.time( (expB = exp(B)) )
			B <- as.matrix(A)
			expB = exp(B)
			#return(expB)
			expAB = expB
			}
		}
	return(expAB)
	#cat("number_of_areas = ", number_of_areas, "\n", sep="")
	# cat("time to exponentiate B, user: ", stB[1], ", system: ", stB[2], ", elapsed: ", stB[3], "\n", sep="")
	}
	




#######################################################
# THIS STUFF MAY BE JUNK -- PERHAPS EVERYTHING INVOLVING 
# SPARSE MATRIX EXPONENTIATION
#######################################################

#' JUNK - testexp enabled for \code{lapply}
#'
#' descrip
#'
#' @param index The index of the input (just for \code{\link{lapply}} purposes.
#' @param matrix_dim The order of the matrix (number of rows or columns)
#' @param expwhich option "A" (spam matrix, not real exponentiation really) or "B" (standard matrix) or "C" (exponentiate_sparse_matrix_w_lots_of_zeros)
#' @return \code{index} Just returns the input index value, since we are just
#' burning cycles for timing purposes here.
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
testexp_sub <- function(index, matrix_dim, expwhich="A")
		{
		# make one or zero non-diagonals to be 0.5
		k = sample(1:matrix_dim, matrix_dim, replace=FALSE)
		i = ceiling(k/matrix_dim) + 0
		j = (k - (i-1)*matrix_dim) + 0
		
		
		# blank matrix
		(A <- spam(0, matrix_dim, matrix_dim))
		
		# make the i,j cells 0.5
		(A[cbind(i, j)] <- rep(.5, length(i)))
		
		# t(A) reflects across the matrix
		# A is original
		# diag.spam makes the diagonals = 1
		(A <- t(A) + A + diag.spam(matrix_dim))
		
		# A is a SPAM matrix
		
		# B is a standard matrix
		B <- as.matrix(A)
		
		# make a CSR matrix for SparseM
		C = as.matrix.csr(B)
		
		#stA = system.time( (expA = exp(A)) )
		if (expwhich == "A")
			{
			expA = exp(A)
			#return(expA)
			expABC = expA
			}
		else if (expwhich == "B")
			{
			#stB = system.time( (expB = exp(B)) )
			# B <- as.matrix(A)
			expB = exp(B)
			#return(expB)
			expABC = expB
			}
		else if (expwhich == "C")
			{
			# C = as.matrix.csr(as.matrix(A))
			expABC = exponentiate_sparse_matrix_w_lots_of_zeros(C)
			}
		return(index)
		}



#' JUNK - Run \code{testexp_sub} with \code{lapply}
#'
#' Runs \code{\link{testexp_sub}} with \code{\link{lapply}}
#'
#' @param n number of areas
#' @param numtimes number of times to do the exponentiation
#' @param expwhich which type of matrix exponentiation to run (see \code{\link{testexp_sub}})
#' @return something
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
testexp2 <- function(n, numtimes, expwhich="A")
	{
	# test timing of matrix exponentiation
	#n = 30
	number_of_areas = (2^n) - 1
	cat("number_of_areas = ", number_of_areas, "\n", sep="")
	
	matrix_dim = number_of_areas
	
	indices = 1:numtimes
	lapply(indices, testexp_sub, matrix_dim, expwhich)
	
	return(indices)
	#cat("number_of_areas = ", number_of_areas, "\n", sep="")
	# cat("time to exponentiate B, user: ", stB[1], ", system: ", stB[2], ", elapsed: ", stB[3], "\n", sep="")
	}


#' JUNK? - Exponentiate a sparse matrix with lots of zeros
#'
#' This text is taken from: http://emdbolker.wikidot.com/seeddisp
#'
#' I'm not at all sure if it is legit - NJM.
#'
#' "The clever way is to use an R package that supports sparse matrices, and hence
#' automatically ignores all the out-of-plot interactions. However, in order to 
#' do this, we have to hack things a bit. First, we have to tell R how to
#' exponentiate just the non-zero elements of a sparse matrix: this requires 
#' hacking into the internal structure of the matrix a bit. Second, this approach
#' is mathematically wrong (but works) if we set the out-of-plot element distances
#' to zero - they should really be infinite (and then reduced to zero by the
#' seed-shadow function), but in this case "two wrongs make a right". As long as
#' there are no true zero distances in the data set, and as long as the out-of-plot
#' seed shadow contributions are always zero, this approach gives us the right
#' answer. The expected-seed function:"
#'
#' (i.e. exponentiate only the nonzero cells)
#' 
#' @param x the matrix (in SparseM format I think)
#' @param alpha the time or distance power to exponentiate by
#' @return something
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
exponentiate_sparse_matrix_w_lots_of_zeros <- function(x, alpha=1)
	{
	if (inherits(x,"matrix.csr"))
		{
		print("Sparse matrix HACK")
		x@ra <- exp(alpha * x@ra) ## Sparse matrix HACK
		}
	else if (is.matrix(x))
		{
		print("Regular matrix exponentiation")
		x <- exp(alpha * x)
		}
	else
		{
		stop("x should be a matrix")
		}
	return(x)
	}













make_list_of_indices_from_2_lists <- function(list1, list2)
	{
	defaults = '
	list1=c(1,3,5,7)
	list2 = list1+1
	'
	# initialize a list
	newlist = NULL
	
	# error check
	if (length(list1) != length(list2))
		{
		stop("make_list_of_indices_from_2_lists() says: ERROR, lists must be of the same length...")
		}
	
	# Go through the lists and make the new list
	for (i in 1:length(list1))
		{
		# get the new indices
		new_indices = list1[i] : list2[i]
		
		# add them to the list
		newlist = c(newlist, new_indices)
		}
	
	return(newlist)
	}


figure_out_step_val <- function(level, list_of_matrices)
	{
	# Figure out the step value (byval)
	nrows_of_matrices_left = 1

	# if you're on the last small matrix
	if (level >= length(list_of_matrices))
		{
		nrows_of_matrices_left = 1
		return(nrows_of_matrices_left)
		}
	
	# Otherwise
	for (j in ((level+1):length(list_of_matrices)))
		{
		nrows_of_matrices_left = nrows_of_matrices_left * nrow(list_of_matrices[[j]])
		}
	return(nrows_of_matrices_left)
	
	}




combine_3matrices <- function(dmat, emat, cmat)
	{
	# Unnecessarily elaborate function to combine 3 relative-prob transition 
	# matricies into one big one.
	
	dnum = nrow(dmat)
	enum = nrow(emat)
	cnum = nrow(cmat)

	matrix_dim = dnum * enum * cnum
	tmpmat = matrix(data=1, nrow=matrix_dim, ncol=matrix_dim)
	tmpmat_names = matrix(data=NA, nrow=matrix_dim, ncol=matrix_dim)
	
	# populate with ds
	size_non_d = enum*cnum
	for (i in 1:dnum)
		{
		start1 = 1+(i-1)*size_non_d 
		end1 = start1 + size_non_d - 1
		indices_to_change1 = start1:end1
		
		for (j in 1:dnum)
			{
			start2 = 1+(j-1)*size_non_d 
			end2 = start2 + size_non_d - 1
			indices_to_change2 = start2:end2
			
			cat(indices_to_change1, " -- ", indices_to_change2, "\n", sep="")
			
			tmpmat[indices_to_change1, indices_to_change2] = tmpmat[indices_to_change1, indices_to_change2] * dmat[i, j]
			
			tmpmat_names[indices_to_change1, indices_to_change2] = paste("d[", i, ",", j, "]*", sep="")
			

			size_non_e = cnum
			for (ei in 1:enum)
				{
				estart1 = start1+(ei-1)*size_non_e 
				eend1 = estart1 + size_non_e - 1
				eindices_to_change1 = estart1:eend1
								
				for (ej in 1:enum)
					{
					estart2 = start2+(ej-1)*size_non_e 
					eend2 = estart2 + size_non_e - 1
					eindices_to_change2 = estart2:eend2
					
					cat(eindices_to_change1, " -- ", eindices_to_change2, "\n", sep="")
					
					tmpmat[eindices_to_change1, eindices_to_change2] = tmpmat[eindices_to_change1, eindices_to_change2] * emat[ei, ej]
					
					tmpmat_names[eindices_to_change1, eindices_to_change2] = paste(tmpmat_names[eindices_to_change1, eindices_to_change2], "e[", ei, ",", ej, "]*", sep="")
					
					
					
					size_non_c = 1
					for (ci in 1:cnum)
						{
						cstart1 = estart1+(ci-1)*size_non_c 
						cend1 = cstart1 + size_non_c - 1
						cindices_to_change1 = cstart1:cend1
										
						for (cj in 1:cnum)
							{
							cstart2 = estart2+(cj-1)*size_non_c 
							cend2 = cstart2 + size_non_c - 1
							cindices_to_change2 = cstart2:cend2
							
							cat(cindices_to_change1, " -- ", cindices_to_change2, "\n", sep="")
							
							tmpmat[cindices_to_change1, cindices_to_change2] = tmpmat[cindices_to_change1, cindices_to_change2] * cmat[ci, cj]
							
							tmpmat_names[cindices_to_change1, cindices_to_change2] = paste(tmpmat_names[cindices_to_change1, cindices_to_change2], "c[", ci, ",", cj, "]", sep="")
							}
						}	
					}
				}		
			}
		}
	result_list = list(tmpmat, tmpmat_names)
	return(result_list)
	}




combine_3matrices2 <- function(dmat, emat, cmat)
	{
	# Unnecessarily elaborate function to combine 3 relative-prob transition 
	# matricies into one big one.
	
	dnum = nrow(dmat)
	enum = nrow(emat)
	cnum = nrow(cmat)

	matrix_dim = dnum * enum * cnum
	tmpmat = matrix(data=1, nrow=matrix_dim, ncol=matrix_dim)
	tmpmat_names = matrix(data=NA, nrow=matrix_dim, ncol=matrix_dim)
	
	# populate with ds
	size_non_d = enum*cnum
	for (i in 1:dnum)
		{
		start1 = 1+(i-1)*size_non_d 
		end1 = start1 + size_non_d - 1
		indices_to_change1 = start1:end1
		
		for (j in 1:dnum)
			{
			start2 = 1+(j-1)*size_non_d 
			end2 = start2 + size_non_d - 1
			indices_to_change2 = start2:end2
			
			cat(indices_to_change1, " -- ", indices_to_change2, "\n", sep="")
			
			tmpmat[indices_to_change1, indices_to_change2] = tmpmat[indices_to_change1, indices_to_change2] * dmat[i, j]
			
			tmpmat_names[indices_to_change1, indices_to_change2] = paste("d*", sep="")
			

			size_non_e = cnum
			for (ei in 1:enum)
				{
				estart1 = start1+(ei-1)*size_non_e 
				eend1 = estart1 + size_non_e - 1
				eindices_to_change1 = estart1:eend1
								
				for (ej in 1:enum)
					{
					estart2 = start2+(ej-1)*size_non_e 
					eend2 = estart2 + size_non_e - 1
					eindices_to_change2 = estart2:eend2
					
					cat(eindices_to_change1, " -- ", eindices_to_change2, "\n", sep="")
					
					tmpmat[eindices_to_change1, eindices_to_change2] = tmpmat[eindices_to_change1, eindices_to_change2] * emat[ei, ej]
					
					tmpmat_names[eindices_to_change1, eindices_to_change2] = paste(tmpmat_names[eindices_to_change1, eindices_to_change2], "e*", sep="")
					
					
					
					size_non_c = 1
					for (ci in 1:cnum)
						{
						cstart1 = estart1+(ci-1)*size_non_c 
						cend1 = cstart1 + size_non_c - 1
						cindices_to_change1 = cstart1:cend1
										
						for (cj in 1:cnum)
							{
							cstart2 = estart2+(cj-1)*size_non_c 
							cend2 = cstart2 + size_non_c - 1
							cindices_to_change2 = cstart2:cend2
							
							cat(cindices_to_change1, " -- ", cindices_to_change2, "\n", sep="")
							
							tmpmat[cindices_to_change1, cindices_to_change2] = tmpmat[cindices_to_change1, cindices_to_change2] * cmat[ci, cj]
							
							tmpmat_names[cindices_to_change1, cindices_to_change2] = paste(tmpmat_names[cindices_to_change1, cindices_to_change2], "c", sep="")
							}
						}	
					}
				}		
			}
		}
	result_list = list(tmpmat, tmpmat_names)
	return(result_list)
	}









# Convert a Qmat to an Rmat, when the Qmat is COO-formatted

Qmat_to_Rmat_COO <- function(cooQmat, n=NA, basefreqs=NA, return_scalefactor=FALSE)
	{
	defaults='
	vals = c(0,0.01,0.01,0.01,0.01,0,0,0,0,0,0,0,0,0,0,0,0,-0.04,0,0,0,0.01,0.01,0.01,0,0,0,0,0,0,0,0,0,0,-0.04,0,0,0.01,0,0,0.01,0.01,0,0,0,0,0,0,0,0,0,-0.04,0,0,0.01,0,0.01,0,0.01,0,0,0,0,0,0,0,0,0,-0.04,0,0,0.01,0,0.01,0.01,0,0,0,0,0,0,0.01,0.01,0,0,-0.06,0,0,0,0,0,0.01,0.01,0,0,0,0,0.01,0,0.01,0,0,-0.06,0,0,0,0,0.01,0,0.01,0,0,0,0.01,0,0,0.01,0,0,-0.06,0,0,0,0,0.01,0.01,0,0,0,0,0.01,0.01,0,0,0,0,-0.06,0,0,0.01,0,0,0.01,0,0,0,0.01,0,0.01,0,0,0,0,-0.06,0,0,0.01,0,0.01,0,0,0,0,0.01,0.01,0,0,0,0,0,-0.06,0,0,0.01,0.01,0,0,0,0,0,0,0.02,0.02,0,0.02,0,0,-0.06,0,0,0,0.01,0,0,0,0,0,0.02,0,0.02,0,0.02,0,0,-0.06,0,0,0.01,0,0,0,0,0,0,0.02,0.02,0,0,0.02,0,0,-0.06,0,0.01,0,0,0,0,0,0,0,0,0.02,0.02,0.02,0,0,0,-0.06,0.01,0,0,0,0,0,0,0,0,0,0,0,0.03,0.03,0.03,0.03,-0.04)
	
	Qmat = matrix(vals, nrow=16, ncol=16, byrow=FALSE)
	
	library(BioGeoBEARS)
	cooQmat = mat2coo(Qmat)
	cooQmat
	
	n=16
	basefreqs = rep(1, times=n)
	' # END defaults
	
	if (nrow(cooQmat) < 1)
		{
		stoptxt = paste0("STOP ERROR in Qmat_to_Rmat_COO(): the input 'cooQmat' has 0 rows! This may happen e.g. if you have a geography rate matrix with 0 transitions, i.e. a 1x1 rate matrix, i.e. a state space of size 1. The BioGeoBEARS default dense matrix routine can handle this, but the force_sparse=TRUE option cannot. And you don't need a sparse matrix anyway, with such a small state space!.")
		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
		stop(stoptxt)
		}
	
	if (is.na(n))
		{
		n = max(max(ia), max(ja))
		txt = paste0("WARNING in Qmat_to_Rmat_COO(). The input 'n', the assumed dimension of the transition matrix, was guessed from max(max(ia), max(ja)). This guess was n=", n, ". If this is wrong, input n manually.")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		warning(txt)
		} # END if (is.na(n))
	
	if ((length(basefreqs)==1) && is.na(basefreqs))
		{
		basefreqs = rep(1, times=n)
		} # END if ((length(basefreqs)==1) && is.na(basefreqs))
	
	if (length(basefreqs) != n)
		{
		stoptxt = paste0("STOP ERROR in Qmat_to_Rmat_COO(). The input 'n', the assumed dimension of the transition matrix, has to be the same as the length of basefreqs. Instead, you have n=", n, ", length(basefreqs)=", length(basefreqs), ".")
		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
		stop(stoptxt)
		} # END if (length(basefreqs) != n)


	
	
	# Sum the rows, removing diagonal values
	cooQmat_df = as.data.frame(cooQmat)

	# Before summing the rows, first, divide each column by the base frequencies
	cooQmat_df$a = cooQmat_df$a / basefreqs[cooQmat_df$ja]

	# Remove diagonals for COO-Rmat (since the diagonals are zero/NA in an Rmat
	notdiag_TF = (cooQmat_df$ia == cooQmat_df$ja) == FALSE
	cooQmat_df2 = cooQmat_df[notdiag_TF, ]
	cooQmat_df2_summed = aggregate(a~ia, data=cooQmat_df2, FUN=sum)
	
	
	# Check for rows that sum to 0, insert
	row_summed_to_zero_TF = ((1:n) %in% cooQmat_df2_summed$ia) == FALSE
	rownums_to_add = (1:n)[row_summed_to_zero_TF]
	rows_to_add = matrix(data=0, nrow=length(rownums_to_add), ncol=2)
	rows_to_add[,1] = rownums_to_add
	rows_to_add_df = as.data.frame(rows_to_add)
	names(rows_to_add_df) = names(cooQmat_df2_summed)
	cooQmat_df2_summed2 = rbind(rows_to_add_df, cooQmat_df2_summed)
	cooQmat_df2_summed2 = as.data.frame(cooQmat_df2_summed2, stringsAsFactors=FALSE)
	
	# Sort sums table by ia, just in case
# 	print(cooQmat_df2_summed2$ia)
# 	print(class(cooQmat_df2_summed2$ia))
# 	print("here2")
	cooQmat_df2_summed3 = cooQmat_df2_summed2[order(cooQmat_df2_summed2$ia), ]
	cooQmat_df2_summed3
	
	# scalefactor is a series of rowSums
	scalefactor = cooQmat_df2_summed3
	
	
	# Error check -- does the table of rowsums have the same number of rows as the state space length?
	if (nrow(cooQmat_df2_summed3) != n)
		{
		stoptxt = paste0("STOP ERROR in Qmat_to_Rmat_COO(). The cooQmat_df2_summed3 object, the sum of each row of the Q matrix (which is in COO format) must have a number of rows equal to 'n', the assumed dimension of the transition matrix. Instead, you have n=", n, ", nrow(cooQmat_df2_summed3)=", nrow(cooQmat_df2_summed3), ". Printing cooQmat_df2_summed3:")
		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
		print("cooQmat_df2_summed3")
		print(cooQmat_df2_summed3)
		stop(stoptxt)
		}
	

	# Divide the non-diagonals by the row-sums (ignoring diagonals)
	# Also, make the COO-formated Rmat (zeros on diagonal, then rows sum to 1)
	cooRmat_df_nDiags = cooQmat_df2
	cooRmat_df_nDiags$a = cooQmat_df2$a / cooQmat_df2_summed3$a[cooQmat_df2$ia]
	if (return_scalefactor == FALSE)
		{
		return(cooRmat_df_nDiags)
		} else {
#				print(cooRmat_df_nDiags)
#				print(class(cooRmat_df_nDiags))

		scalefactor = scalefactor

		return(scalefactor)
		} # END if (return_scalefactor == FALSE)
	
	return(cooRmat_df_nDiags)
	
	# Dense matrix method
# 	Repeat basefreqs each row
# 	basefreqs_mat = basefreqs_to_basefreqs_mat(basefreqs)
#   (ie, so multiply the first COLUMN by the first basefreq, etc.)
# 
# 	Get a matrix with 1s on off-diagonal, 0s on diagonal
# 	offdiags_1 = zeros_diag_ones_offdiag(size=ncol(Qmat))
# 	
# 	Make sure the off-diagonals IN EACH ROW sum to 1
# 	Rmat = Qmat / basefreqs_mat
# 
# 	if (return_scalefactor == FALSE)
# 		{
# 		Rmat = Rmat / rowSums(Rmat * offdiags_1)
# 		return(Rmat)
# 		} else {
# 		return(rowSums(Rmat * offdiags_1))
# 		}
# 
# 	return(Rmat)
	}


# scale_by_ja_instead_of_ia = TRUE is what works for trait-based rate matrices
Rmat_to_Qmat_noScaling_COO <- function(cooRmat_df, n=NA, basefreqs=NA, use_scaling=FALSE, use_basefreqs=FALSE, scale_by_ja_instead_of_ia)
	{
	defaults='
	vals = c(0,0.01,0.01,0.01,0.01,0,0,0,0,0,0,0,0,0,0,0,0,-0.04,0,0,0,0.01,0.01,0.01,0,0,0,0,0,0,0,0,0,0,-0.04,0,0,0.01,0,0,0.01,0.01,0,0,0,0,0,0,0,0,0,-0.04,0,0,0.01,0,0.01,0,0.01,0,0,0,0,0,0,0,0,0,-0.04,0,0,0.01,0,0.01,0.01,0,0,0,0,0,0,0.01,0.01,0,0,-0.06,0,0,0,0,0,0.01,0.01,0,0,0,0,0.01,0,0.01,0,0,-0.06,0,0,0,0,0.01,0,0.01,0,0,0,0.01,0,0,0.01,0,0,-0.06,0,0,0,0,0.01,0.01,0,0,0,0,0.01,0.01,0,0,0,0,-0.06,0,0,0.01,0,0,0.01,0,0,0,0.01,0,0.01,0,0,0,0,-0.06,0,0,0.01,0,0.01,0,0,0,0,0.01,0.01,0,0,0,0,0,-0.06,0,0,0.01,0.01,0,0,0,0,0,0,0.02,0.02,0,0.02,0,0,-0.06,0,0,0,0.01,0,0,0,0,0,0.02,0,0.02,0,0.02,0,0,-0.06,0,0,0.01,0,0,0,0,0,0,0.02,0.02,0,0,0.02,0,0,-0.06,0,0.01,0,0,0,0,0,0,0,0,0.02,0.02,0.02,0,0,0,-0.06,0.01,0,0,0,0,0,0,0,0,0,0,0,0.03,0.03,0.03,0.03,-0.04)
	
	Qmat = matrix(vals, nrow=16, ncol=16, byrow=FALSE)
	
	library(BioGeoBEARS)
	cooQmat = mat2coo(Qmat)
	cooRmat_df = as.data.frame(cooQmat)
	names(cooRmat_df) = c("ia", "ja", "a")
	
	n=16
	basefreqs = rep(1, times=n)
	use_scaling=FALSE
	use_basefreqs=FALSE
	scale_by_ja_instead_of_ia = TRUE
	' # END defaults


	if (is.na(n))
		{
		n = max(max(ia), max(ja))
		txt = paste0("WARNING in Rmat_to_Qmat_noScaling_COO(). The input 'n', the assumed dimension of the transition matrix, was guessed from max(max(ia), max(ja)). This guess was n=", n, ". If this is wrong, input n manually.")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		warning(txt)
		} # END if (is.na(n))
	
	if ((length(basefreqs)==1) && is.na(basefreqs))
		{
		basefreqs = rep(1, times=n)
		} # END if ((length(basefreqs)==1) && is.na(basefreqs))
	
	if (length(basefreqs) != n)
		{
		stoptxt = paste0("STOP ERROR in Rmat_to_Qmat_noScaling_COO(). The input 'n', the assumed dimension of the transition matrix, has to be the same as the length of basefreqs. Instead, you have n=", n, ", length(basefreqs)=", length(basefreqs), ".")
		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
		warning(stoptxt)
		} # END if (length(basefreqs) != n)

	
	# Starter
	tmp_cooQmat_df = cooRmat_df
	
	# Just in case, exclude diagonals
	nodiag_TF = (tmp_cooQmat_df$ia == tmp_cooQmat_df$ja) == FALSE
	
	# Sum the "rows" (ia) for the diagonal
	if (scale_by_ja_instead_of_ia == FALSE)
		{
		if (use_basefreqs)
			{
			# Repeat basefreqs each row
			#basefreqs_mat = basefreqs_to_basefreqs_mat(basefreqs)
			# element-by-element multiplication
			# (ie, multiply the first COLUMN by the first basefreq, etc.)
			tmp_cooQmat_df$a = tmp_cooQmat_df$a * basefreqs[tmp_cooQmat_df$ja]
			} else {
			junk=1
			} # END if (use_basefreqs)
	
		if (use_scaling)
			{
			# Get the sum of the diagonal
			# This is sort of a bogus scaling factor (any Rmat matrix would have it)
			scaling_factor = sum(tmp_cooQmat_df$a[nodiag_TF])
	
			tmp_cooQmat_df$a = tmp_cooQmat_df$a / scaling_factor
			} else {
			junk=1
			} # END if (use_scaling)
	
		# Normalize the diagonal to be -rowSum(row_without_diagonal)
		#Qmat = normat(relative_matrix=Qmat)
	
		# Sum each row of the COOmat
		tmp_cooQmat_df2 = tmp_cooQmat_df[nodiag_TF,]
		cooQmat_df2_summed = aggregate(a~ia, data=tmp_cooQmat_df2, FUN=sum)

		# Check for rows that sum to 0, insert
		row_summed_to_zero_TF = ((1:n) %in% cooQmat_df2_summed$ia) == FALSE
		rownums_to_add = (1:n)[row_summed_to_zero_TF]
		rows_to_add = matrix(data=0, nrow=length(rownums_to_add), ncol=2)
		rows_to_add[,1] = rownums_to_add
		rows_to_add_df = as.data.frame(rows_to_add)
		names(rows_to_add_df) = names(cooQmat_df2_summed)
		cooQmat_df2_summed2 = rbind(rows_to_add_df, cooQmat_df2_summed)
		} # END if (scale_by_ja_instead_of_ia == FALSE)


	# Sum the "columns" (ja) for the diagonal
	if (scale_by_ja_instead_of_ia == TRUE)
		{
		if (use_basefreqs)
			{
			# Repeat basefreqs each row
			#basefreqs_mat = basefreqs_to_basefreqs_mat(basefreqs)
			# element-by-element multiplication
			# (ie, multiply the first COLUMN by the first basefreq, etc.)
			tmp_cooQmat_df$a = tmp_cooQmat_df$a * basefreqs[tmp_cooQmat_df$ia]
			} else {
			junk=1
			} # END if (use_basefreqs)
	
		if (use_scaling)
			{
			# Get the sum of the djagonal
			# This is sort of a bogus scaling factor (any Rmat matrix would have it)
			scaling_factor = sum(tmp_cooQmat_df$a[nodiag_TF])
	
			tmp_cooQmat_df$a = tmp_cooQmat_df$a / scaling_factor
			} else {
			junk=1
			} # END if (use_scaling)
	
		# Normalize the djagonal to be -rowSum(row_without_djagonal)
		#Qmat = normat(relative_matrix=Qmat)
	
		# Sum each row of the COOmat
		tmp_cooQmat_df2 = tmp_cooQmat_df[nodiag_TF,]
		cooQmat_df2_summed = aggregate(a~ja, data=tmp_cooQmat_df2, FUN=sum)
		
		
		# Check for rows that sum to 0, insert
		row_summed_to_zero_TF = ((1:n) %in% cooQmat_df2_summed$ja) == FALSE
		rownums_to_add = (1:n)[row_summed_to_zero_TF]
		rows_to_add = matrix(data=0, nrow=length(rownums_to_add), ncol=2)
		rows_to_add[,1] = rownums_to_add
		rows_to_add_df = as.data.frame(rows_to_add)

		# change the colname "ja" to "ia" for consistency
		names(cooQmat_df2_summed) = c("ia", "a")
		names(rows_to_add_df) = names(cooQmat_df2_summed)
		cooQmat_df2_summed2 = rbind(rows_to_add_df, cooQmat_df2_summed)
		} # END if (scale_by_ja_instead_of_ia == TRUE)

	# Sort sums table by ia, just in case
	cooQmat_df2_summed2 = as.data.frame(cooQmat_df2_summed2, stringsAsFactors=FALSE)
# 	print("here1")
# 	print(cooQmat_df2_summed2$ia)
# 	print(class(cooQmat_df2_summed2$ia))
	cooQmat_df2_summed3 = cooQmat_df2_summed2[order(cooQmat_df2_summed2$ia), ]
	cooQmat_df2_summed3

	# Error check -- does the table of rowsums have the same number of rows as the state space length?
	if (nrow(cooQmat_df2_summed3) != n)
		{
		stoptxt = paste0("STOP ERROR in Rmat_to_Qmat_noScaling_COO(). The cooQmat_df2_summed3 object, the sum of each row of the Q matrix (which is in COO format) must have a number of rows equal to 'n', the assumed dimension of the transition matrix. Instead, you have n=", n, ", nrow(cooQmat_df2_summed3)=", nrow(cooQmat_df2_summed3), ". Printing cooQmat_df2_summed3:")
		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
		print("cooQmat_df2_summed3")
		print(cooQmat_df2_summed3)
		stop(stoptxt)
		}
	
	# Add these sums as diagonals
	tmp_ia = 1:n
	tmp_ja = 1:n
	tmp_a = -1 * cooQmat_df2_summed3$a
	diag_vals_to_add_to_cooQmat = cbind(tmp_ia, tmp_ja, tmp_a)
	diag_vals_to_add_to_cooQmat_df = as.data.frame(diag_vals_to_add_to_cooQmat)
	names(diag_vals_to_add_to_cooQmat_df) = names(tmp_cooQmat_df2)
	
	# Combine the diagonals and off-diagonals
	tmp_cooQmat_df3 = rbind(diag_vals_to_add_to_cooQmat_df, tmp_cooQmat_df2)
	# Sort by columns, then sort by rows for final product
	tmp_cooQmat_df3 = tmp_cooQmat_df3[order(tmp_cooQmat_df3$ja), ]
	tmp_cooQmat_df3 = tmp_cooQmat_df3[order(tmp_cooQmat_df3$ia), ]
	
	return(tmp_cooQmat_df3)
	}




# Convert a COO-formated matrix to CRS (Compressed sparse row (CSR, CRS or Yale format))
# https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format)
# 
# See also: kexpmv by Meabh, e.g. crs2mat, mat2crs
# See also: SparseM library (faster, but I am avoiding the dependencies and class structure)
# 
# Inputs:
# ia = row indexes (1-based) for non-zero values in the transition matrix
# ja = column indexes (1-based) for non-zero values in the transition matrix
# a = the values at those coordinates
coo2crs <- function(ia, ja, a, n=NA, transpose_needed=FALSE)
	{
	defaults='
	# Qmat from Psychotria (16x16)
	vals = c(0,0.01,0.01,0.01,0.01,0,0,0,0,0,0,0,0,0,0,0,0,-0.04,0,0,0,0.01,0.01,0.01,0,0,0,0,0,0,0,0,0,0,-0.04,0,0,0.01,0,0,0.01,0.01,0,0,0,0,0,0,0,0,0,-0.04,0,0,0.01,0,0.01,0,0.01,0,0,0,0,0,0,0,0,0,-0.04,0,0,0.01,0,0.01,0.01,0,0,0,0,0,0,0.01,0.01,0,0,-0.06,0,0,0,0,0,0.01,0.01,0,0,0,0,0.01,0,0.01,0,0,-0.06,0,0,0,0,0.01,0,0.01,0,0,0,0.01,0,0,0.01,0,0,-0.06,0,0,0,0,0.01,0.01,0,0,0,0,0.01,0.01,0,0,0,0,-0.06,0,0,0.01,0,0,0.01,0,0,0,0.01,0,0.01,0,0,0,0,-0.06,0,0,0.01,0,0.01,0,0,0,0,0.01,0.01,0,0,0,0,0,-0.06,0,0,0.01,0.01,0,0,0,0,0,0,0.02,0.02,0,0.02,0,0,-0.06,0,0,0,0.01,0,0,0,0,0,0.02,0,0.02,0,0.02,0,0,-0.06,0,0,0.01,0,0,0,0,0,0,0.02,0.02,0,0,0.02,0,0,-0.06,0,0.01,0,0,0,0,0,0,0,0,0.02,0.02,0.02,0,0,0,-0.06,0.01,0,0,0,0,0,0,0,0,0,0,0,0.03,0.03,0.03,0.03,-0.04)
	
	Qmat = matrix(vals, nrow=16, ncol=16, byrow=FALSE)
	library(SparseM)
	library(kexpmv)
	crsQmat = as.matrix.csr(Qmat)
	crsQmat
	
	library(BioGeoBEARS)
	cooQmat = mat2coo(Qmat)
	cooQmat
	
	# Test coo2crs
	ia=cooQmat[,"ia"]
	ja=cooQmat[,"ja"]
	a=cooQmat[,"a"]
	transpose_needed = FALSE
	n=16
	crsR_Qmat = coo2crs(ia=ia, ja=ja, a=a, n=16, transpose_needed=FALSE)
	crsR_Qmat
	
	# Check:
	crsQmat@ia
	crsR_Qmat$ia
	all.equal(crsR_Qmat$ia, crsQmat@ia)
	'
	
	
	# If the input is a zero-length vector (COO format with 0 rows)
	if (any(length(ia)==0, length(ja)==0, length(a)==0) == TRUE)
		{
		txt = paste0("WARNING from coo2crs: you have input a COO-formated matrix with 0 rows. This suggests you have a 1x1 transition matrix. If so, you really don't need to be doing sparse matrix exponentiation!  In BioGeoBEARS, change force_sparse to FALSE.")
		cat("\n")
		cat(txt)
		warning(txt)
		}
	
	
	# Calculate the dimension of the transition matrix, if needed
	if (is.na(n))
		{
		n = max(max(ia), max(ja))
		txt = paste0("WARNING in coo2crs(). The input 'n', the assumed dimension of the transition matrix, was guessed from max(max(ia), max(ja)). This guess was n=", n, ". If this is wrong, input n manually.")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		warning(txt)
		} # END if (is.na(n))
	
	
	# Check that the COO inputs are legit
	if (all.equal(length(ia), length(ja), length(a)) == FALSE)
		{
		txt = paste0("STOP ERROR in coo2crs(). The inputs 'ia', 'ja', and 'a' must all have the same length. Instead, you have length(ia)=", length(ia), ", length(ja)=", length(ja), ", length(a)=", length(a), ".")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		} # END if (all(equal(length(ia), length(ja), length(a))) == FALSE)
	
	# Tranpose if needed
	if (transpose_needed == TRUE)
		{
		tmp = ia 
		ia = ja
		ja = tmp
		} # END if (transpose_needed == TRUE)
	
	nrows = n
	ncols = n

	# Sort the values by row, then column
	tdf = data.frame(ia, ja, a)
	tdf2 = tdf[order(tdf[,2]),] # order by columns, then by rows
	tdf2 = tdf[order(tdf[,1]),] # order by rows last
	
	# Compressed Row Storage (CRS)
	# http://netlib.org/linalg/html_templates/node91.html	

	# Get the rows of the ordered-COO table where the row increments
	# (and the first row)
	if (nrow(tdf2) > 1)
		{
		all_ia_minus_last = tdf2$ia[1:(length(tdf2$ia)-1)]
		all_ia_minus_first = tdf2$ia[2:length(tdf2$ia)]
		TF = (all_ia_minus_last+1) == all_ia_minus_first 
		crs_ia = ((1:length(ia))[TF])+1
		} else {
		# If there is only 1 row
		# (e.g., for a 1-area problem)
		# (then all you need for ia is the first row)
		crs_ia = NULL
		}

	# Add the row for the first entry
	crs_ia = c(1, crs_ia, length(tdf2$a)+1)
	
# 	if ( (crs_ia[length(crs_ia)]) == (length(tdf2$a)+1) )
# 		{
# 		# Don't add another index after the last ia referencing the last a
# 		} else {
# 		# Add the row for the first entry and the last (dummy) row
# 		crs_ia = c(crs_ia, length(tdf2$a)+1)
# 		}
	
	# If the length of ia is not the number of rows+1, you have blank rows. These should be repeated once
	if (length(crs_ia) < (nrows+1))
		{
		TF = ((1:nrows) %in% ia) == FALSE
		missing_rows = (1:nrows)[TF]
		crs_ia = sort(c(missing_rows, crs_ia))
		
		if (length(crs_ia) != (nrows+1))
			{
			stoptxt = paste0("STOP ERROR in coo2crs: the output crs_ia (indexes to which 'a' values where the rows change) should have nrows+1 entries. However, the function is producing length(crs_ia)=", length(crs_ia), " and (nrows+1)=", (nrows+1), ". This indicates an error in the function dealing with missing rows. Check your inputs, n, and the function.")
			cat("\n\n")
			cat(stoptxt)
			cat("\n\n")
			
			print("crs_ia:")
			print(crs_ia)
			stop(stoptxt)
			} # END if (length(crs_ia) < (nrows+1))
		} # if (length(crs_ia) != (nrows+1))


# 	if (length(crs_ia) > (nrows+1))
# 		{
# 		stoptxt = paste0("STOP ERROR in coo2crs: the output crs_ia (indexes to which 'a' values where the rows change) should have nrows+1 entries. However, the function is producing length(crs_ia)=", length(crs_ia), " and (nrows+1)=", (nrows+1), ". This indicates somehow you've added extra rows. Check your inputs, n, and the function.")
# 		cat("\n\n")
# 		cat(stoptxt)
# 		cat("\n\n")
# 		
# 		print("crs_ia:")
# 		print(crs_ia)
# 		stop(stoptxt)
# 		} # END if (length(crs_ia) > (nrows+1))
	
	crsmat = list()
	crsmat$ia = crs_ia
	crsmat$ja = tdf2$ja # Cols are identical to COO (once sorted)
	crsmat$a = tdf2$a		# Values are identical to COO (once sorted)
	crsmat$dimension = c(nrows, ncols)
	return(crsmat)
	}




