#######################################################
# Calc loglike speciation -- separated to its own file, because it's such a huge function
#######################################################
require(ape)
require(rexpokit)
require(cladoRcpp)
require(parallel)

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
	# Not parallel processing
	independent_likelihoods_on_each_branch = mapply(FUN=expokit_dgpadm_Qmat2, Qmat=list(Qmat), t=phy2$edge.length, transpose_needed=TRUE, SIMPLIFY="array")	
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
#' Default is \code{NULL} which means running on a single processor.
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
calc_loglike_sp_prebyte <- function(tip_condlikes_of_data_on_each_state, phy, Qmat, spPmat=NULL, min_branchlength=1e-21, return_what="loglike", probs_of_states_at_root=NULL, rootedge=FALSE, sparse=FALSE, printlevel=1, use_cpp=TRUE, input_is_COO=FALSE, spPmat_inputs=NULL, cppSpMethod=3, cluster_already_open=NULL)
	{
	defaults='
	# Phylogeny
	newick_str = "(((Humans:6.0, Chimps:6.0):1.0, Gorillas:7.0):5.0, Orangs:12.0):1.0;"
	phy = read.tree(file="", text=newick_str)

	# Areas: A=Africa, B=Not Africa
	areas = c("A","B")
	states_list = areas_list_to_states_list(areas, include_null_range=FALSE)
	dedf = make_relprob_matrix_de(states_list)
	Qmat = symbolic_to_Q_matrix(dedf, d=0.1, e=0.1)
	
	numnodes = length(phy$tip.label) + phy$Nnode
	tip_condlikes_of_data_on_each_state = matrix(data=0, nrow=length(phy$tip.label), ncol=length(states_list))
	# The columns are tip likelihoods of: A, B, AB
	tip_condlikes_of_data_on_each_state[1:3,1] = 1
	tip_condlikes_of_data_on_each_state[4,2] = 1
	tip_condlikes_of_data_on_each_state
	

	spPmat=NULL
	
	
	min_branchlength=1e-21
	return_what = "all"
	probs_of_states_at_root=NULL
	rootedge=FALSE
	cppSpMethod=3
	cluster_already_open=NULL
	printlevel=1
	sparse=FALSE
	#sparse=TRUE
	
	' # end defaults
	
	#print("Qmat#1")
	#print(Qmat)
	
	

	
	
	#######################################################
	# ERROR CHECK
	#######################################################
	if (use_cpp == TRUE && !is.null(spPmat_inputs))
		{
		numstates_in_tip_condlikes_of_data_on_each_state = ncol(tip_condlikes_of_data_on_each_state)
		numstates_in_spPmat_inputs_states_indices = 1 + length(spPmat_inputs$l)
		
		if (input_is_COO == TRUE)
			{
			numstates_in_Qmat = max( max(Qmat[,1]), max(Qmat[,2]) )
			} else {
			numstates_in_Qmat = ncol(Qmat)
			}
		
		if ( all(numstates_in_tip_condlikes_of_data_on_each_state==numstates_in_spPmat_inputs_states_indices, numstates_in_tip_condlikes_of_data_on_each_state==numstates_in_Qmat ) == TRUE )
			{
			# Continue
			all_inputs_correct_size = TRUE
			} else {
			
			# Stop if not everything is equal
			stop_function <- function(numstates_in_tip_condlikes_of_data_on_each_state, numstates_in_spPmat_inputs_states_indices, numstates_in_Qmat)
				{
				cat("\n\nERROR: Some inputs have incorrect size -- \n")
				cat("numstates_in_tip_condlikes_of_data_on_each_state:	", numstates_in_tip_condlikes_of_data_on_each_state, "\n")
				cat("numstates_in_spPmat_inputs_states_indices+1:	", numstates_in_spPmat_inputs_states_indices, "\n")
				cat("numstates_in_Qmat:								", numstates_in_Qmat, "\n")
				cat("\n")
				cat("This probably means 'd' and 'e' were so close to zero that rcpp_states_list_to_DEmat()'s\n")
				cat("decided they were effectively 0; default min_precision=1e-16.\n")
				cat("\n")
				}
			stop(stop_function(numstates_in_tip_condlikes_of_data_on_each_state, numstates_in_spPmat_inputs_states_indices, numstates_in_Qmat))
			}
		} else {
		numstates_in_tip_condlikes_of_data_on_each_state = ncol(tip_condlikes_of_data_on_each_state)
		numstates_in_spPmat_ie_nrows = 1 + nrow(spPmat)
		
		if (input_is_COO == TRUE)
			{
			numstates_in_Qmat = max( max(Qmat[,1]), max(Qmat[,2]) )
			} else {
			numstates_in_Qmat = ncol(Qmat)
			}
		
		if ( all(numstates_in_tip_condlikes_of_data_on_each_state==numstates_in_spPmat_ie_nrows, numstates_in_tip_condlikes_of_data_on_each_state==numstates_in_Qmat ) == TRUE )
			{
			# Continue
			all_inputs_correct_size = TRUE
			} else {
			
			# Stop if not everything is equal
			stop_function2 <- function(numstates_in_tip_condlikes_of_data_on_each_state, numstates_in_spPmat_ie_nrows, numstates_in_Qmat)
				{
				cat("\n\nERROR: Some inputs have incorrect size -- \n")
				cat("numstates_in_tip_condlikes_of_data_on_each_state:	", numstates_in_tip_condlikes_of_data_on_each_state, "\n")
				cat("numstates_in_spPmat_ie_nrows+1:				", numstates_in_spPmat_ie_nrows, "\n")
				cat("numstates_in_Qmat:								", numstates_in_Qmat, "\n")
				cat("\n")
				}
			stop(stop_function2(numstates_in_tip_condlikes_of_data_on_each_state, numstates_in_spPmat_ie_nrows, numstates_in_Qmat))
			}
		
		}
	
	
	
	
	
	
	
	edgelengths = phy$edge.length
	num_branches_below_min = sum(edgelengths < min_branchlength)

	if ( (printlevel >= 2) || (num_branches_below_min > 0) )
		{
		cat("\n")
		cat("Running calc_loglike_sp():\n")
		cat("This run of calc_loglike_sp() has a min_branchlength of: ", min_branchlength, "\n", sep="")
		cat("Branches shorter than this will be assumed to be connected to the tree with\n")
		cat("sympatric events (i.e., members of fossil lineages on ~0 length branches.)\n")
		cat("This tree has ", num_branches_below_min, " branches < ", min_branchlength, ".\n", sep="")
		cat("\n")
		}	
	
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
	computed_likelihoods_at_each_node[tipnums] = rowSums(tip_condlikes_of_data_on_each_state)
	

	relative_probs_of_each_state <- matrix(data=0, nrow=numnodes, ncol=numstates)
	relative_probs_of_each_state[tipnums, ] = tip_condlikes_of_data_on_each_state
	condlikes_of_each_state = relative_probs_of_each_state
	#relative_probs_of_each_state[tipnums, ] <- 1
	
	
	# If you want to use standard dense matrix exponentiation, rexpokit's dgpadm is best,
	# and you can run it through mapply on the branch lengths of all branches at once
	if (sparse==FALSE)
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
					}
				}

			
			# Run on the cluster of nodes, if one is open
			# clusterApply etc. appear to NOT work on R.app
			if (!is.null(cluster_already_open))
				{
				# mcmapply
				#library(parallel)
				#independent_likelihoods_on_each_branch = mcmapply(FUN=expokit_dgpadm_Qmat, Qmat=list(Qmat), t=phy2$edge.length, transpose_needed=TRUE, SIMPLIFY="array", mc.cores=Ncores)
				independent_likelihoods_on_each_branch = clusterApply(cl=cluster_already_open, x=phy2$edge.length, fun=expokit_dgpadm_Qmat2,Qmat=Qmat, transpose_needed=TRUE)
				} else {
				# Not parallel processing
				#independent_likelihoods_on_each_branch = mapply(FUN=expokit_dgpadm_Qmat, Qmat=list(Qmat), t=phy2$edge.length, transpose_needed=TRUE, SIMPLIFY="array")
				independent_likelihoods_on_each_branch = mapply_likelihoods(Qmat, phy2, transpose_needed=TRUE)
				#independent_likelihoods_on_each_branch
				}
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
				if ( (class(Qmat) != "data.frame") || (ncol(Qmat) != 3) )
					{
					stop("ERROR: calc_loglike_sp is attempting to use a sparse COO-formated Q matrix, but you provided a regular square dense matrix")
					}
				
				
				coo_n = numstates
				anorm = 1
				tmpQmat_in_REXPOKIT_coo_fmt = Qmat
				#tmpQmat_in_REXPOKIT_coo_fmt = cooQmat
				}
			}
	
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
			
			dmat = spPmat_inputs$dmat
			
			# Take the max of the indices of the possible areas, and add 1
			numareas = max(unlist(spPmat_inputs$l), na.rm=TRUE) + 1
			
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


			tmpca_1 = rep(1, (ncol(relative_probs_of_each_state)-1))
			tmpcb_1 = rep(1, (ncol(relative_probs_of_each_state)-1))

			# Print the matrix to screen from C++
			printmat = FALSE
			if (printlevel >= 1)
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
				COO_probs_list = rcpp_calc_anclikes_sp_COOprobs(Rcpp_leftprobs=tmpca_1, Rcpp_rightprobs=tmpcb_1, l=l, s=s, v=v, j=j, y=y, dmat=dmat, maxent01s=maxent01s, maxent01v=maxent01v, maxent01j=maxent01j, maxent01y=maxent01y, max_minsize_as_function_of_ancsize=max_minsize_as_function_of_ancsize, printmat=printmat)

				# Sum the probabilities (list [[3]]) in each row ([[3]][[1]] list of probs
				# through [[3]][[15]] list of probs)
				Rsp_rowsums = sapply(X=COO_probs_list[[3]], FUN=sum)
				}
	
			if (cppSpMethod == 3)
				{
				if (printlevel >= 1)
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

					}
				
				COO_weights_columnar = rcpp_calc_anclikes_sp_COOweights_faster(Rcpp_leftprobs=tmpca_1, Rcpp_rightprobs=tmpcb_1, l=l, s=s, v=v, j=j, y=y, dmat=dmat, maxent01s=maxent01s, maxent01v=maxent01v, maxent01j=maxent01j, maxent01y=maxent01y, max_minsize_as_function_of_ancsize=max_minsize_as_function_of_ancsize, printmat=printmat)
				
				# combine with C++ function
				# This causes an error with spPmat=NULL; spPmat_inputs=NULL; use_cpp=TRUE; sparse=FALSE
				# i.e. gives 16 states with a 0 on the end, rather than 15 states
				#Rsp_rowsums = rcpp_calc_rowsums_for_COOweights_columnar(COO_weights_columnar=COO_weights_columnar, numstates=numstates)
				
				# This gives 15 states
				Rsp_rowsums = rcpp_calc_rowsums_for_COOweights_columnar(COO_weights_columnar=COO_weights_columnar)
				}

	
	
			}
		}
	
	
	
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
		
		# And for the ancestor edge
		anc <- phy2$edge[i, 1]
		
		txt = paste("anc:", anc, " left:", left_desc_nodenum, " right:", right_desc_nodenum, sep="")
		#print(txt)
		
		# Is sparse is FALSE, input the pre-calculated likelihoods;
		# If sparse is TRUE, dynamically calculate using expokit_dmexpv_Qmat
		if (sparse==FALSE)
			{
			
			#print("Checking Q matrix")
			#cat("p=", p, ", rate=", rate, "\n", sep=" ")
			#print(Q)
			#print(phy$edge.length[i])
			#condlikes_Left <- matexpo(Qmat * phy2$edge.length[i]) %*% relative_probs_of_each_state[left_desc_nodenum,]
			#condlikes_Right <- matexpo(Qmat * phy2$edge.length[j]) %*% relative_probs_of_each_state[right_desc_nodenum,]


			if (printlevel >= 1)
				{
				print("dense matrix exponentiation")
				}
			
			
			if (is.null(cluster_already_open))
				{
				# Conditional likelihoods of states at the bottom of left branch
				condlikes_Left = independent_likelihoods_on_each_branch[,,i] %*% relative_probs_of_each_state[left_desc_nodenum,]
							
				# Conditional likelihoods of states at the bottom of right branch
				condlikes_Right = independent_likelihoods_on_each_branch[,,j] %*% relative_probs_of_each_state[right_desc_nodenum,]
				} else {
				
				#cat("dim(independent_likelihoods_on_each_branch[[i]]):\n", sep="")
				#cat(dim(independent_likelihoods_on_each_branch))
				#cat("\n\n")
				
				#cat("length(relative_probs_of_each_state[left_desc_nodenum,]):\n", sep="")
				#cat(length(relative_probs_of_each_state[left_desc_nodenum,]))
				#cat("\n\n")
				
				condlikes_Left = independent_likelihoods_on_each_branch[[i]] %*% relative_probs_of_each_state[left_desc_nodenum,]
							
				# Conditional likelihoods of states at the bottom of right branch
				condlikes_Right = independent_likelihoods_on_each_branch[[j]] %*% relative_probs_of_each_state[right_desc_nodenum,]
				}


			
			if (printlevel >= 2) {
			txt = paste("condlikes at bottom of L: ", paste(round(condlikes_Left, 4), collapse=" ", sep=""), sep="")
			print(txt)
			}
			
			if (printlevel >= 2) {
			txt = paste("condlikes at bottom of R: ", paste(round(condlikes_Right, 4), collapse=" ", sep=""), sep="")
			print(txt)
			}
			#condlikes_Left <- expm(Qmat * phy$edge.length[i], method="Ward77") %*% relative_probs_of_each_state[left_desc_nodenum, 
			#  ]
			#condlikes_Right <- expm(Qmat * phy$edge.length[j], method="Ward77") %*% relative_probs_of_each_state[right_desc_nodenum, 
			#  ]
			
			
			#condlikes_Left <- exp(Qmat * phy$edge.length[i]) %*% relative_probs_of_each_state[left_desc_nodenum, 
			#  ]
			#condlikes_Right <- exp(Qmat * phy$edge.length[j]) %*% relative_probs_of_each_state[right_desc_nodenum, 
			#  ]
			} else {
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

			if (printlevel >= 1)
				{
				print("sparse matrix exponentiation")
				}
			
			# Conditional likelihoods of states at the bottom of left branch
			condlikes_Left = try (
			expokit_dmexpv_Qmat(Qmat=tmpQmat_in_REXPOKIT_coo_fmt, t=phy2$edge.length[i], inputprobs_for_fast=relative_probs_of_each_state[left_desc_nodenum,], transpose_needed=FALSE, transform_to_coo_TF=FALSE, coo_n=coo_n, anorm=anorm, check_for_0_rows=TRUE)
			)
			
			# Error check
			if (class(condlikes_Left) == "try-error")
				{
				cat("\n\ntry-error on expokit_dmexpv_Qmat():\n\n")
				cat("i=", i, "\n")
				cat("phy2$edge.length[i]=", phy2$edge.length[i], "\n")
				print(tmpQmat_in_REXPOKIT_coo_fmt)
				print(phy2)
				print(relative_probs_of_each_state)
				print(relative_probs_of_each_state[left_desc_nodenum,])
				print(coo_n)
				print(anorm)
				}
			
			if (any(is.nan(condlikes_Left)))
				{
				cat("\n\nexpokit_dmexpv_Qmat() returned NaNs:\n\n")
				cat("i=", i, "\n")
				cat("phy2$edge.length[i]=", phy2$edge.length[i], "\n")
				print(tmpQmat_in_REXPOKIT_coo_fmt)
				print(phy2)
				print(relative_probs_of_each_state)
				print(relative_probs_of_each_state[left_desc_nodenum,])
				print(coo_n)
				print(anorm)
				}
	
			# Conditional likelihoods of states at the bottom of right branch
			condlikes_Right = try(
			expokit_dmexpv_Qmat(Qmat=tmpQmat_in_REXPOKIT_coo_fmt, t=phy2$edge.length[j], inputprobs_for_fast=relative_probs_of_each_state[right_desc_nodenum,], transpose_needed=FALSE, transform_to_coo_TF=FALSE, coo_n=coo_n, anorm=anorm, check_for_0_rows=TRUE)
			)

			# Error check
			if (class(condlikes_Right) == "try-error")
				{
				cat("\n\ntry-error on expokit_dmexpv_Qmat():\n\n")
				cat("j=", j, "\n")
				cat("phy2$edge.length[j]=", phy2$edge.length[j], "\n")
				print(tmpQmat_in_REXPOKIT_coo_fmt)
				print(phy2)
				print(relative_probs_of_each_state)
				print(relative_probs_of_each_state[left_desc_nodenum,])
				print(coo_n)
				print(anorm)
				}


			if (any(is.nan(condlikes_Right)))
				{
				cat("\n\ntry-error on expokit_dmexpv_Qmat():\n\n")
				cat("j=", j, "\n")
				cat("phy2$edge.length[j]=", phy2$edge.length[j], "\n")
				print(tmpQmat_in_REXPOKIT_coo_fmt)
				print(phy2)
				print(relative_probs_of_each_state)
				print(relative_probs_of_each_state[left_desc_nodenum,])
				print(coo_n)
				print(anorm)
				}

			#print(c(condlikes_Left))
			#print(c(condlikes_Right))
			
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
				if (printlevel >= 1) {
				print("use_cpp=FALSE, direct multiplication")
				}
				} else {
				##################################
				# combine the likelihoods from each branch bottom with this speciational model
				# of range change
				##################################

				if (printlevel >= 1) {
				print("use_cpp=FALSE, speciation model")
				}
				
				# for each ancestral state, get prob of branch pairs 
				outmat = matrix(0, nrow=nrow(spPmat), ncol=ncol(spPmat))
				
				# Get the probs of each pair, list by row
				# (This might be a slightly slow step for large matrices)
				
				# Exclude "_" ranges from this calculation, as if there is speciation,
				# they are not going extinct
				probs_of_each_desc_pair = c(base::t(outer(c(condlikes_Left[-1]), c(condlikes_Right[-1]), "*")))
				#sum(probs_of_each_desc_pair)
	
				# Make a matrix with the probability of each pair, in each cell
				# (duplicated for each ancestral state)
				for (i in 1:nrow(spPmat))
					{
					outmat[i,] = probs_of_each_desc_pair
					}
				# Multiply the probabilities of each pair, by the probability of each
				# type of speciation event, to get the probabilities at the 
				# ancestral node
				outmat2 = outmat * spPmat
				
				node_likelihood_with_speciation = rowSums(outmat2)
				
				
				# THIS ZERO IS ALREADY OVER-WRITING THE NULL STATE LIKELIHOOD!!
				
				# Add the 0 back in, representing 0 probability of "_"
				# range just below speciation event
				node_likelihood = c(0, node_likelihood_with_speciation)
				}
			} else {
			# use_cpp == TRUE
			if ( hooknode_TF==TRUE )
				{
				##################################
				# no speciational model of range change
				##################################
				if (printlevel >= 1) 
					{
					print("use_cpp=TRUE, direct multiplication")
					}
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
						if (printlevel >= 1)
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
						probs_of_each_desc_pair = c(base::t(outer(c(condlikes_Left[-1]), c(condlikes_Right[-1]), "*")))
						#sum(probs_of_each_desc_pair)
			
						# Make a matrix with the probability of each pair, in each cell
						# (duplicated for each ancestral state)
						for (i in 1:nrow(spPmat))
							{
							outmat[i,] = probs_of_each_desc_pair
							}
						# Multiply the probabilities of each pair, by the probability of each
						# type of speciation event, to get the probabilities at the 
						# ancestral node
						outmat2 = outmat * spPmat
						
						node_likelihood_with_speciation = rowSums(outmat2)
						
						
						# THIS ZERO IS ALREADY OVER-WRITING THE NULL STATE LIKELIHOOD!!
						
						# Add the 0 back in, representing 0 probability of "_"
						# range just below speciation event
						node_likelihood = c(0, node_likelihood_with_speciation)
						}
					} else {
					if (printlevel >= 1)
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
					ca = condlikes_Left[-1]
					cb = condlikes_Right[-1]
					
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
					
					dmat = spPmat_inputs$dmat
	
					# Print the matrix to screen from C++
					printmat = FALSE
					if (printlevel >= 1) {
					printmat = TRUE
					}
					
					# Calculate the rowsums (for input into rcpp_calc_anclikes_sp()
					# Actually, just do this ONCE (above)
					#Rsp_rowsums = rcpp_calc_anclikes_sp_rowsums(Rcpp_leftprobs=tmpca_1, Rcpp_rightprobs=tmpcb_1, l=l, s=s, v=v, j=j, y=y, dmat=dmat, maxent01s, maxent01v, maxent01j, maxent01y, printmat=printmat)
					
					if (printlevel >= 2) {
					print("Rsp_rowsums:")
					print(Rsp_rowsums)
					}
					
					#node_likelihood_with_speciation = rep(0.0, length(tmpca_2))
					# This calculates the speciation probabilities again at each node; this is inefficient
					if (cppSpMethod == 1)
						{
						node_likelihood_with_speciation = rcpp_calc_anclikes_sp(Rcpp_leftprobs=tmpca_2, Rcpp_rightprobs=tmpcb_2, l=l, s=s, v=v, j=j, y=y, dmat=dmat, maxent01s=maxent01s, maxent01v=maxent01v, maxent01j=maxent01j, maxent01y=maxent01y, max_minsize_as_function_of_ancsize=max_minsize_as_function_of_ancsize, Rsp_rowsums=Rsp_rowsums, printmat=printmat)						
						#print(node_likelihood_with_speciation[1:5])
						}
					
					# Really, we should just iterate through the COO_probs_list using a C++ function
					# 2012-12-04 NJM says: this KICKS ASS in C++!!
					
					if (cppSpMethod == 2)
						{
						node_likelihood_with_speciation2 = rep(0.0, length(tmpca_2))
						node_likelihood_with_speciation2 = rcpp_calc_anclikes_sp_using_COOprobs(Rcpp_leftprobs=tmpca_2, Rcpp_rightprobs=tmpcb_2, RCOO_left_i_list=COO_probs_list[[1]], RCOO_right_j_list=COO_probs_list[[2]], RCOO_probs_list=COO_probs_list[[3]], Rsp_rowsums=Rsp_rowsums, printmat=printmat)
						#print(node_likelihood_with_speciation2[1:5])
						node_likelihood_with_speciation = node_likelihood_with_speciation2
						}

					if (cppSpMethod == 3)
						{
						node_likelihood_with_speciation3 = rep(0.0, length(tmpca_2))
						node_likelihood_with_speciation3 = rcpp_calc_splitlikes_using_COOweights_columnar(Rcpp_leftprobs=tmpca_2, Rcpp_rightprobs=tmpcb_2, COO_weights_columnar=COO_weights_columnar, Rsp_rowsums=Rsp_rowsums, printmat=printmat)
						#print(node_likelihood_with_speciation2[1:5])
						node_likelihood_with_speciation = node_likelihood_with_speciation3
						}

					
					if (printlevel >= 2) {
					print("node_likelihood_with_speciation:")
					print(node_likelihood_with_speciation)
					}
					node_likelihood = c(0, node_likelihood_with_speciation)				
					}
				}
			}


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
		
		total_likelihood_for_node = sum(node_likelihood)
		
		computed_likelihoods_at_each_node[anc] = total_likelihood_for_node
		
		relative_probs_of_each_state[anc, ] = node_likelihood / total_likelihood_for_node
		condlikes_of_each_state[anc, ] = node_likelihood
		
		#print(relative_probs_of_each_state)
		}
	relative_probs_of_each_state
	# relative_probs_of_each_state2 = relative_probs_of_each_state

	# You are now at the anc (ancestor node; in Psychotria, node 20)
	# If you want, you could calculate the likelihood down to the bottom of a root edge
	# below the root node
	# check if the rootedge option is desired, and the phylogeny HAS a root edge
	#phy2$root.edge = 0.1
	if ( (rootedge == TRUE) && (!is.null(phy2$root.edge)) && (phy2$root.edge > 0) )
		{
		phy2$root.edge = 1.1
		root_edge_length = phy2$root.edge
		independent_likelihoods_on_root_edge = expokit_dgpadm_Qmat(Qmat=Qmat, t=root_edge_length, transpose_needed=TRUE)
		
		# Get the likelihoods at the bottom of the branch, condition on the relative likelihoods at the top 
		# (here, the node at the "top" is the Last Common Ancestor node on the tree)
		conditional_likelihoods_at_bottom_of_root_branch = c(independent_likelihoods_on_root_edge %*% relative_probs_of_each_state[anc, ])
		total_likelihood_for_bottom_of_root_branch = sum(conditional_likelihoods_at_bottom_of_root_branch)
		total_likelihood_for_bottom_of_root_branch
		
		
		# Get the relative probability of each state at the bottom of the root branch
		relative_probs_of_each_state_at_bottom_of_root_branch = c(conditional_likelihoods_at_bottom_of_root_branch) / total_likelihood_for_bottom_of_root_branch
		
		
		# Add the likelihood of the bottom of the root branch to the list of node likelihoods
		computed_likelihoods_at_each_node = c(computed_likelihoods_at_each_node, total_likelihood_for_bottom_of_root_branch)
		
		# Add a bottom row to the relative_probs
		relative_probs_of_each_state = rbind(relative_probs_of_each_state, relative_probs_of_each_state_at_bottom_of_root_branch)
		condlikes_of_each_state = rbind(condlikes_of_each_state, conditional_likelihoods_at_bottom_of_root_branch)

		} else {
		# Otherwise, the root relative probabilities are just the last relative_probs_of_each_state[anc, ]
		relative_probs_of_each_state_at_bottom_of_root_branch = relative_probs_of_each_state[anc, ]
		}
	
	
	# why times 2?? -- probably this is some AIC thing from the original APE function
	#output_loglike = -2 * sum(log(comp[-TIPS]))
	#output_loglike = 2 * sum(log(comp[-TIPS]))
	
	total_loglikelihood = sum(log(computed_likelihoods_at_each_node))




	if (printlevel >= 1)
		{
		cat("total_loglikelihood:	", total_loglikelihood, "\n", sep="")
		}
	
	if (return_what == "loglike")
		{
		return(total_loglikelihood)
		} 

	if (return_what == "nodelikes")
		{
		return(relative_probs_of_each_state)
		}

	if (return_what == "rootprobs")
		{
		return(total_loglikelihood)
		}
	if (return_what == "all")
		{
		calc_loglike_sp_results = list()
		calc_loglike_sp_results$computed_likelihoods_at_each_node = computed_likelihoods_at_each_node
		calc_loglike_sp_results$total_loglikelihood = total_loglikelihood
		calc_loglike_sp_results$relative_probs_of_each_state = relative_probs_of_each_state
		calc_loglike_sp_results$condlikes_of_each_state = condlikes_of_each_state
		calc_loglike_sp_results$relative_probs_of_each_state_at_bottom_of_root_branch = relative_probs_of_each_state_at_bottom_of_root_branch
		class(calc_loglike_sp_results) = "calc_loglike_sp_results"
		return(calc_loglike_sp_results)
		}
	}























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
#' Default is \code{NULL} which means running on a single processor.
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
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)






