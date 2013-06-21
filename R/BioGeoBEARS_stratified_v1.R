require("ape")
require("rexpokit")
require("cladoRcpp")


#######################################################
# Do a stratified or other constrained analysis
#######################################################

# Based on:
# /_examples/changing_geog_v1.R




#######################################################
# section_the_tree
#######################################################
#' Section a tree for stratified analysis
#'
#' A utility function for stratified analysis.
#' 
#' @param inputs The list of inputs for stratified analysis
#' @return \code{inputs} with inputs$tree_sections_list added.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
section_the_tree <- function(inputs)
	{
	timeperiods = inputs$timeperiods
	phy_as_it_is_chopped_down = read.tree(inputs$trfn)
	
	plot(phy_as_it_is_chopped_down)
	#abline(v=timeperiods)
	axisPhylo()
	# CHECK THIS FUNCTION
	#phy_as_it_is_chopped_down$edge.length = phy_as_it_is_chopped_down$edge.length + 0.0001
	
	
	tree_sections_list = NULL
	tnum = 0
	
	if (length(timeperiods) <= 1)
		{
		chainsaw_result = list()
		chainsaw_result$tree_to_chainsaw = phy_as_it_is_chopped_down
		chainsaw_result$return_pieces_list[[1]] = phy_as_it_is_chopped_down
		chainsaw_result$return_pieces_basenames[[1]] = paste(sort(phy_as_it_is_chopped_down$tip.label), collapse="|", sep="")
		attr(chainsaw_result, "class") = "chainsaw_result"
		tree_sections_list[[1]] = chainsaw_result
		} else {
		
		for (i in 1:(length(timeperiods)))
		#for (i in 1:3))
			{
			cat("\n", i, ": ", timeperiods[i], "\n", "\n", sep="")
			# Chainsaw the top off the tree
			if (i == 1)
				{
				timepoint = timeperiods[i] - 0
				} else {
				timepoint = timeperiods[i]# - timeperiods[i-1]
				}
			# Update timepoints so you are subtracting the right amount!!!!!!!!
			timeperiods = timeperiods - timepoint
		  timeperiods
		  
		  # If it's the last piece, just use the remaining tree
		  if (i < length(timeperiods))
		  	{
			chainsaw_result = chainsaw2(phy_as_it_is_chopped_down, timepoint=timepoint, return_pieces=TRUE)
			#print(chainsaw_result)
			} else {
			chainsaw_result = list()
			chainsaw_result$tree_to_chainsaw = phy_as_it_is_chopped_down
			chainsaw_result$return_pieces_list[[1]] = phy_as_it_is_chopped_down
			chainsaw_result$return_pieces_basenames[[1]] = paste(sort(phy_as_it_is_chopped_down$tip.label), collapse="|", sep="")
			attr(chainsaw_result, "class") = "chainsaw_result"
			}
			
			# Store the chainsaw result
			tree_sections_list[[(tnum=tnum+1)]] = chainsaw_result
	
			# Convey the tree to the next round of chopping
			phy_as_it_is_chopped_down = chainsaw_result$tree_to_chainsaw
	
			plot(phy_as_it_is_chopped_down)
			#axisPhylo2(side = 1, roundlabels=TRUE, minage=timeperiods[i] 
			axisPhylo()
			
			}
		}

		
	# Append to inputs and return
	inputs$tree_sections_list = tree_sections_list
	
	return(inputs)
	}




#######################################################
# calc_loglike_sp_stratified:
#######################################################
#' Calculate log-likelihood with a transition matrix and speciation events, and with stratification
#'
#' This function is the stratified version of \code{\link{calc_loglike_sp}}.
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
#' @param null_range_allowed Does the state space include the null range?
#' Default is \code{NULL} which means running on a single processor.
#' @param fixnode If the state at a particular node is going to be fixed (e.g. for ML marginal ancestral states), give the node number.
#' @param fixlikes The state likelihoods to be used at the fixed node.  I.e. 1 for the fixed state, and 0 for the others.
#' @param inputs A list of inputs containing the dispersal matrix for each time period, etc.
#' @param allareas A list of all the areas in the total analysis
#' @param all_states_list A list of all the stats in the total analysis (0-based coding - ?)
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
calc_loglike_sp_stratified <- function(tip_condlikes_of_data_on_each_state, phy, Qmat=NULL, spPmat=NULL, min_branchlength=1e-21, return_what="loglike", probs_of_states_at_root=NULL, rootedge=TRUE, sparse=FALSE, printlevel=0, use_cpp=TRUE, input_is_COO=FALSE, spPmat_inputs=NULL, cppSpMethod=3, cluster_already_open=NULL, calc_ancprobs=FALSE, null_range_allowed=TRUE, fixnode=NULL, fixlikes=NULL, inputs=inputs, allareas=allareas, all_states_list=all_states_list)
	{
	defaults='
	Qmat=NULL; spPmat=NULL; min_branchlength=1e-21; return_what="loglike"; probs_of_states_at_root=NULL; rootedge=FALSE; sparse=FALSE; printlevel=1; use_cpp=TRUE; input_is_COO=FALSE; spPmat_inputs=NULL; cppSpMethod=3; cluster_already_open=NULL; calc_ancprobs=FALSE; null_range_allowed=TRUE; fixnode=NULL; fixlikes=NULL; inputs=inputs; allareas=allareas; all_states_list=all_states_list
	
	
	maxareas = 4
	phy = read.tree(inputs$trfn)
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=inputs$geogfn)
	tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges, phy, maxareas=maxareas)
	
	allareas = getareas_from_tipranges_object(tipranges)
	all_states_list = rcpp_areas_list_to_states_list(areas=allareas, include_null_range=TRUE, maxareas=maxareas)
	
	tmpres = calc_loglike_sp_stratified(tip_condlikes_of_data_on_each_state, phy, Qmat=NULL, spPmat=NULL, min_branchlength=1e-21, return_what="all", probs_of_states_at_root=NULL, rootedge=TRUE, sparse=FALSE, printlevel=0, use_cpp=TRUE, input_is_COO=FALSE, spPmat_inputs=NULL, cppSpMethod=3, cluster_already_open=NULL, calc_ancprobs=FALSE, null_range_allowed=TRUE, fixnode=NULL, fixlikes=NULL, inputs=inputs, allareas=allareas, all_states_list=all_states_list)
	tmpres
	'

	# Get the timeperiods; if 1 time period, run once; if multiple, run 
	if (is.null(inputs$timeperiods) || length(inputs$timeperiods) == 1)
		{
		num_iterations = 1
		} else {
		# Multiple timeperiods
		timeperiods = inputs$timeperiods
		num_iterations = length(timeperiods)
		}
	
	# All areas in the total analysis
	allareas=allareas
	allareas_list = seq(0, length(allareas)-1, 1)		# 0-base indexes

	# All states in the total analysis (after e.g. limitation on total # of areas)
	all_states_list=all_states_list	
	
	# Other variables
	BioGeoBEARS_model_object = inputs$BioGeoBEARS_model_object
	force_sparse = sparse
	
	#######################################################
	# Set up the starting probabilities etc.
	#######################################################
	# Starting tip_relative_probs_of_each_state
	current_condlikes_row = 0
	tip_relative_probs_of_each_state = tip_condlikes_of_data_on_each_state
	tip_relative_probs_of_each_state
	current_tip_relative_probs_of_each_state = tip_relative_probs_of_each_state
	current_condlikes_row = nrow(current_tip_relative_probs_of_each_state)
	
	
	# matrix to hold all of the relative probabilities; Making this purposely too big
	numnodes = phy$Nnode + length(phy$tip.label)
	all_relative_probs_of_each_state = matrix(0, ncol=length(all_states_list), nrow=(numnodes*length(timeperiods)))
	all_condlikes_of_each_state = matrix(0, ncol=length(all_states_list), nrow=(numnodes*length(timeperiods)))
	
	all_relative_probs_of_each_state[1:current_condlikes_row, ] = current_tip_relative_probs_of_each_state
	all_condlikes_of_each_state[1:current_condlikes_row, ] = current_tip_relative_probs_of_each_state
	
	
	#######################################################
	# Take the original tree and scale the branchlengths by b (branch-length exponent)
	# b=0, all branches=1; b=1, all branches normal
	#######################################################
	previous_timepoint = 0
	original_phy = phy
	b_branch_length_exponent = inputs$BioGeoBEARS_model_object@params_table["b", "est"]
	#original_phy$edge.length = original_phy$edge.length ^ b_branch_length_exponent
	phy_as_it_is_chopped_down = original_phy

	tiplikes_to_delete = list()

	
	for (i in 1:num_iterations)
		{
		#cat("\ni=",i, sep="")

		# Cut down the number of areas, by what is allowed
		areas_allowed_mat = inputs$list_of_areas_allowed_mats[[i]]

		states_allowed_TF = sapply(X=all_states_list, FUN=check_if_state_is_allowed, areas_allowed_mat)
		#states_to_use_TF = all_states_list %in% tmp_states_list
		states_to_use_TF = states_allowed_TF
		
		
		# States allowed in this timeperiod
		states_list = all_states_list[states_allowed_TF]
		
		# Add back NULL range, if needed
		if (null_range_allowed == TRUE)
			{
			states_list = c(NA, states_list)
			states_to_use_TF[1] = TRUE
			}
		
		# Make the dedf matrix for this time period
		dispersal_multipliers_matrix = inputs$list_of_dispersal_multipliers_mats[[i]]
		distances_mat = inputs$list_of_distances_mats[[i]]
		x_exponent = inputs$BioGeoBEARS_model_object@params_table["x", "est"]
		dispersal_multipliers_matrix = make_dispersal_multiplier_matrix(areas=allareas, states_list=NULL, dispersal_multipliers_matrix=dispersal_multipliers_matrix, distances_mat=distances_mat, x_exponent=x_exponent)
		d_current = inputs$BioGeoBEARS_model_object@params_table["d", "est"]
		dmat = d_current * dispersal_multipliers_matrix
		#print(dmat)
		
		# Calculate the extinction (local extipration) probability multipliers
		area_of_areas = inputs$list_of_area_of_areas[[i]]
		# Check elist for 0s
		if (any(elist <= 0))
			{
			stop("ERROR: Minimum distance between regions must be >= 1; correct this.\n\n", sep="")
			}
		u_extirpation_exponent = inputs$BioGeoBEARS_model_object@params_table["u", "est"]
		e_current = inputs$BioGeoBEARS_model_object@params_table["e", "est"]
		elist = e_current * (area_of_areas ^ (-1 * u_extirpation_exponent))
		#print(elist)
		
		# Calculate the Q matrix
		# someday we'll have to put "a" (anagenic range-switching) in here...
		if (is.null(Qmat))
			{
			Qmat_tmp = rcpp_states_list_to_DEmat(areas_list=allareas_list, states_list=states_list, dmat=dmat, elist=elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)
			} else {
			# If Qmat is pre-specified
			Qmat_tmp = Qmat
			}
		
		
		# Now. IF you have a subtree structure, you need to run this with a cladogenesis matrix, 
		# through calc_loglike_sp(), like normal.
		
		# If there's just one tree, store it in the object
		if (is.null(inputs$timeperiods) || length(inputs$timeperiods) == 1)
			{
			tr = read.tree(inputs$trfn)
			
			tree_to_chainsaw = NULL
			tree_to_chainsaw[[1]] = tr

			return_pieces_list = NULL
			return_pieces_list[[1]] = tr
			
			return_pieces_basenames = NULL
			return_pieces_basenames[[1]] = paste(sort(tr$tip.label), collapse="|", sep="")
			
			chainsaw_object = list()
			chainsaw_object$tree_to_chainsaw = tree_to_chainsaw
			chainsaw_object$return_pieces_list = return_pieces_list
			chainsaw_object$return_pieces_basenames = return_pieces_basenames
			attr(chainsaw_object, "class") = "chainsaw_result"
			
			inputs$tree_sections_list[[1]] = chainsaw_object
			}
		
		
		# OK, if you have a tree here, do that
		# if not, exp the branch
		
		#######################################################
		# Cladogenic model 
		#######################################################
		j = BioGeoBEARS_model_object@params_table["j","est"]
		ysv = BioGeoBEARS_model_object@params_table["ys","est"]
		v = BioGeoBEARS_model_object@params_table["v","est"]
		ys = BioGeoBEARS_model_object@params_table["ys","est"]
		sum_SPweights = ys + j + v
	
		maxent_constraint_01 = BioGeoBEARS_model_object@params_table["mx01","est"]
		
		# Text version of speciation matrix	
		maxent_constraint_01v = BioGeoBEARS_model_object@params_table["mx01v","est"]
		#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
			
		# Set the parameter controlling the size distribution of 
		# the smaller descendant species
		maxent01s_param = BioGeoBEARS_model_object@params_table["mx01s","est"]
		maxent01v_param = BioGeoBEARS_model_object@params_table["mx01v","est"]
		maxent01j_param = BioGeoBEARS_model_object@params_table["mx01j","est"]
		maxent01y_param = BioGeoBEARS_model_object@params_table["mx01y","est"]
	
	
		# Cladogenesis model inputs
		spPmat_inputs = NULL
		states_indices = states_list
		states_indices[1] = NULL	# shorten the states_indices by 1 (cutting the null range state from the speciation matrix)
		spPmat_inputs$l = states_indices
		spPmat_inputs$s = ys
		spPmat_inputs$v = v
		spPmat_inputs$j = j
		spPmat_inputs$y = ys
		spPmat_inputs$dmat = dmat
		spPmat_inputs$maxent01s_param = maxent01s_param
		spPmat_inputs$maxent01v_param = maxent01v_param
		spPmat_inputs$maxent01j_param = maxent01j_param
		spPmat_inputs$maxent01y_param = maxent01y_param
	
			
		#######################################################
		# Go through the tree pieces
		#######################################################
		chainsaw_result = inputs$tree_sections_list[[i]]
		
		# You will need the new tip likelihoods of the new tree:
		current_tip_relative_probs_of_each_state
		new_tip_likelihoods = matrix(0, nrow=length(chainsaw_result$return_pieces_list), ncol=length(all_states_list))

		for (jj in 1:length(chainsaw_result$return_pieces_list))
			{
			#cat("\njj=",jj, sep="")
			treepiece = chainsaw_result$return_pieces_list[[jj]]
		
			
			# If it's jjust a branch section
			if (is.numeric(treepiece))
				{
				tipname = chainsaw_result$return_pieces_basenames[[jj]]
				tip_TF = phy_as_it_is_chopped_down$tip.label == tipname
				relative_probs_of_each_state_at_the_tip_of_this_branch = current_tip_relative_probs_of_each_state[tip_TF, states_to_use_TF]
	
				# t = treepiece
				independent_likelihoods_at_branch_section_bottom = expokit_dgpadm_Qmat2(times=treepiece, Qmat=Qmat_tmp,  transpose_needed=TRUE)
				#independent_likelihoods_at_branch_section_bottom = expokit_dgpadm_Qmat(Qmat=Qmat_tmp,  t=treepiece, transpose_needed=FALSE)
				
				
				conditional_likelihoods_at_branch_section_bottom = matrix(independent_likelihoods_at_branch_section_bottom %*% relative_probs_of_each_state_at_the_tip_of_this_branch, nrow=1)
				
				
				# Test forward exponentiation instead...NO
				# independent_likelihoods_at_branch_section_bottom = expokit_dgpadm_Qmat2(times=treepiece, Qmat=Qmat_tmp,  transpose_needed=FALSE)
# 				#independent_likelihoods_at_branch_section_bottom = expokit_dgpadm_Qmat(Qmat=Qmat_tmp,  t=treepiece, transpose_needed=FALSE)
# 				
# 				
# 				conditional_likelihoods_at_branch_section_bottom = matrix(independent_likelihoods_at_branch_section_bottom %*% relative_probs_of_each_state_at_the_tip_of_this_branch, nrow=1)
# 				conditional_likelihoods_at_branch_section_bottom[1] = 0
# 				
				
				
				
				# Also, store the conditional likelihoods for all nodes in this subtree
				chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]] = conditional_likelihoods_at_branch_section_bottom
	
				# Relative probabilities -- jjust the new tip
				chainsaw_result$relative_probs_of_each_state_at_bottom_of_root_branch[[jj]] = conditional_likelihoods_at_branch_section_bottom / sum(conditional_likelihoods_at_branch_section_bottom)
	
				# Relative probabilities -- all nodes plus branch bottom (jjust branch bottom, here)
				chainsaw_result$relative_probabilities_for_nodes_plus_bottom_in_this_section[[jj]] = conditional_likelihoods_at_branch_section_bottom / sum(conditional_likelihoods_at_branch_section_bottom)
	
				} else {
				# Otherwise, treepiece is a subtree
				tmp_subtree = treepiece
				
				# Get the names of the tips in this subtree
				tipnames = tmp_subtree$tip.label
				
				# Use the tipnames to get the conditional likelihoods at these tips
				tips_for_subtree_TF = phy_as_it_is_chopped_down$tip.label %in% tipnames
				subtree_tip_relative_probs_of_each_state = current_tip_relative_probs_of_each_state[tips_for_subtree_TF,states_to_use_TF]
	
				# Calculate the likelihoods for this subtree
				calc_loglike_sp_results = calc_loglike_sp(
					tip_condlikes_of_data_on_each_state=subtree_tip_relative_probs_of_each_state, 
					phy=tmp_subtree, 
					Qmat=Qmat_tmp, 
					spPmat=NULL,
					return_what="all",
					probs_of_states_at_root=NULL,
					rootedge=TRUE,
					sparse=FALSE, 
					printlevel=printlevel,
					use_cpp=TRUE,
					input_is_COO=FALSE,
					spPmat_inputs=spPmat_inputs,
					cppSpMethod=cppSpMethod,
					cluster_already_open=cluster_already_open,
					calc_ancprobs=FALSE,
					null_range_allowed=null_range_allowed,
					fixnode=NULL,
					fixlikes=NULL
					)
	
				#chainsaw_result$conditional_likelihoods_at_branch_section_bottom[[jj]] = 
				
				# Also, store the conditional likelihoods for all nodes in this subtree
				# MINUS THE GODDAMN TIPS OF THE SUBTREE, THESE ARE ALREADY IN THERE
				tmp_tipnums = 1:length(tipnames)
				chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]] = matrix(data=calc_loglike_sp_results$condlikes_of_each_state[-tmp_tipnums, ], ncol=ncol(calc_loglike_sp_results$condlikes_of_each_state))
				
				# Matrix of tip likelihoods to delete so you don't repeat using them in the total
				# loglike
				tiplikes_to_delete[[jj]] = calc_loglike_sp_results$condlikes_of_each_state[tmp_tipnums, ]
				
				# Relative probabilities -- all nodes plus branch bottom (jjust branch bottom, here)
				chainsaw_result$relative_probabilities_for_nodes_plus_bottom_in_this_section[[jj]] = calc_loglike_sp_results$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[-tmp_tipnums, ]
	
				# Relative probabilities -- jjust the new tip
				chainsaw_result$relative_probs_of_each_state_at_bottom_of_root_branch[[jj]] = calc_loglike_sp_results$relative_probs_of_each_state_at_bottom_of_root_branch


				# ONLY for the nodes in original tree, store the condlikes
				# Add these to the overall list of conditional likelihoods
				numrows_to_add = nrow(chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]])
				
				
				
				startrow = current_condlikes_row + 1
				endrow = current_condlikes_row + numrows_to_add
				all_relative_probs_of_each_state[startrow:endrow, states_to_use_TF] = chainsaw_result$relative_probabilities_for_nodes_plus_bottom_in_this_section[[jj]]
				
				all_condlikes_of_each_state[startrow:endrow, states_to_use_TF] = chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]]

				current_condlikes_row = current_condlikes_row + numrows_to_add
				
				} # End if/then on branch vs. subtree
	
			# Also, store the relative probabilities for the new tip
			#new_tip_likelihoods[jj, states_to_use_TF] = chainsaw_result$relative_probs_of_each_state_at_bottom_of_root_branch[[jj]]
			new_tip_likelihoods[jj, states_to_use_TF] = chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]][nrow(chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]]), ]
	
			
# Add these to the overall list of conditional likelihoods
# 			numrows_to_add = nrow(chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]])
# 			
# 			startrow = current_condlikes_row + 1
# 			endrow = current_condlikes_row + numrows_to_add
# 			all_relative_probs_of_each_state[startrow:endrow, states_to_use_TF] = chainsaw_result$relative_probabilities_for_nodes_plus_bottom_in_this_section[[jj]]
# 			
# 			all_condlikes_of_each_state[startrow:endrow, states_to_use_TF] = chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]]
	
			#print(chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]])
			#print(log(sum(chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]])))
			
			#rowSums(all_condlikes_of_each_state) != 0
			#tmp_all_condlikes_of_each_state = all_condlikes_of_each_state[rowSums(all_condlikes_of_each_state) != 0,]
			#currLnL = sum(log(rowSums(tmp_all_condlikes_of_each_state)))
			#cat("i=", i, "; jj=", jj, "; currLnL=", currLnL, "\n")
	
# 			current_condlikes_row = current_condlikes_row + numrows_to_add
	
	
			} # End loop through jj tree pieces WITHIN a stratum



		# Update for the next loop
		# Tip likelihoods
		current_tip_relative_probs_of_each_state = new_tip_likelihoods
	
		# Store previous round
		#old_phy_as_it_is_chopped_down = phy_as_it_is_chopped_down
		#old_chainsaw_result = chainsaw_result
		#old_new_tip_likelihoods = new_tip_likelihoods
		
		# Convey the tree to the next iteration
		phy_as_it_is_chopped_down = chainsaw_result$tree_to_chainsaw

		} # End loop through i strata

	# Remove rows that have not been filled (till zero)
	rowSums(all_condlikes_of_each_state) != 0
	final_all_condlikes_of_each_state = all_condlikes_of_each_state[rowSums(all_condlikes_of_each_state) != 0,]
	
	rowSums(all_relative_probs_of_each_state) != 0
	all_relative_probs_of_each_state = all_relative_probs_of_each_state[rowSums(all_relative_probs_of_each_state) != 0,]
	
	#all_relative_probs_of_each_state
	
	if (rootedge == TRUE)
		{
		grand_total_likelihood = sum(log(rowSums(final_all_condlikes_of_each_state)))
		grand_total_likelihood
		} else {
		# Skip the last row
		grand_total_likelihood = sum(log(rowSums(final_all_condlikes_of_each_state[-nrow(final_all_condlikes_of_each_state),])))
		grand_total_likelihood
		}
	
	return(grand_total_likelihood)
	}








#######################################################
# calc_loglike_for_optim_stratified
#######################################################
#' Take model parameters and the data and calculate the log-likelihood -- stratified version
#' 
#' This is the stratified version of \code{\link{calc_loglike_for_optim}}. This function is an input to optim or optimx, the ML
#' estimation routines.
#' 
#' @param tip_condlikes_of_data_on_each_state A numeric matrix with rows representing tips, and columns representing states/geographic ranges.  The cells
#' give the likelihood of the observation data under the assumption that the tip has that state; typically this means that the known geographic range gets a 
#' '1' and all other states get a 0.
#' @param params A vector of parameters for optimization.
#' @param BioGeoBEARS_run_object Object containing the run parameters, and the model.
#' @param phy An ape tree object
#' @param force_sparse Should sparse matrix exponentiation be used?
#' @param print_optim If TRUE (default), print the optimization steps as ML estimation progresses.
#' @param areas_list A list of the desired area names/abbreviations/letters (?).
#' @param states_list A list of the possible states/geographic ranges, in 0-based index form.
#' @param cluster_already_open The cluster object, if it has already been started.
#' @return \code{ttl_loglike} The log-likelihood of the data under the input model and parameters.
#' @export
#' @seealso \code{\link[stats]{convolve}} chainsaw_result
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
calc_loglike_for_optim_stratified <- function(params, BioGeoBEARS_run_object, phy, tip_condlikes_of_data_on_each_state, print_optim=TRUE, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open)
	{
	
	BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
	
	# Put the parameters into the BioGeoBEARS_model_object, so that they can be universally read out
	# into any function
	BioGeoBEARS_model_object = params_into_BioGeoBEARS_model_object(BioGeoBEARS_model_object=BioGeoBEARS_model_object, params=params)
	
	# Update linked parameters
	BioGeoBEARS_model_object = calc_linked_params_BioGeoBEARS_model_object(BioGeoBEARS_model_object)
	
	# Set the dispersal and extinction rate
	d = BioGeoBEARS_model_object@params_table["d","est"]
	e = BioGeoBEARS_model_object@params_table["e","est"]

	# Equal dispersal in all directions (unconstrained)
	# Equal extinction probability for all areas
	areas = areas_list
# 	distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

	#dmat = matrix(d, nrow=length(areas), ncol=length(areas))
	#elist = rep(e, length(areas))
	
	# Set up the instantaneous rate matrix (Q matrix)
	#Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)

	# Cladogenic model
	j = BioGeoBEARS_model_object@params_table["j","est"]
	ysv = BioGeoBEARS_model_object@params_table["ys","est"]
	v = BioGeoBEARS_model_object@params_table["v","est"]
	ys = BioGeoBEARS_model_object@params_table["ys","est"]
	sum_SPweights = ys + j + v
	sum_SPweights
	
	# Store back in the run object
	BioGeoBEARS_run_object$BioGeoBEARS_model_object = BioGeoBEARS_model_object
	
# 	maxent_constraint_01 = BioGeoBEARS_model_object@params_table["mx01","est"]
# 	
# 	# Text version of speciation matrix	
# 	maxent_constraint_01v = BioGeoBEARS_model_object@params_table["mx01v","est"]
# 	#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
# 		
# 	# Set the parameter controlling the size distribution of 
# 	# the smaller descendant species
# 	maxent01s_param = BioGeoBEARS_model_object@params_table["mx01s","est"]
# 	maxent01v_param = BioGeoBEARS_model_object@params_table["mx01v","est"]
# 	maxent01j_param = BioGeoBEARS_model_object@params_table["mx01j","est"]
# 	maxent01y_param = BioGeoBEARS_model_object@params_table["mx01y","est"]
# 
# 
# 	# Cladogenesis model inputs
# 	spPmat_inputs = NULL
# 	states_indices = states_list
# 	states_indices[1] = NULL	# shorten the states_indices by 1 (cutting the null range state from the speciation matrix)
# 	spPmat_inputs$l = states_indices
# 	spPmat_inputs$s = ys
# 	spPmat_inputs$v = v
# 	spPmat_inputs$j = j
# 	spPmat_inputs$y = ys
# 	spPmat_inputs$dmat = distances_mat
# 	spPmat_inputs$maxent01s_param = maxent01s_param
# 	spPmat_inputs$maxent01v_param = maxent01v_param
# 	spPmat_inputs$maxent01j_param = maxent01j_param
# 	spPmat_inputs$maxent01y_param = maxent01y_param
# 

	if (print_optim == TRUE)
		{
		#outvars = as.data.frame(t(BioGeoBEARS_model_object@params_table$est))
		#names(outvars) = rownames(BioGeoBEARS_model_object@params_table)
		#outvars = c(BioGeoBEARS_model_object@params_table$est)
		
		#cat("\n")
		#cat(outvars, sep="	")
		
		# Before calculating the log likelihood, print it, in case there is e.g. a bug
		#cat("d=", d, "; e=", e, "; j=", j, "; ys=", ys, "; v=", v, "; maxent01=", maxent_constraint_01, "; maxent01v=", maxent_constraint_01v, "; sum=", sum_SPweights, "; LnL=", sep="")
		}

	# Calculate the log-likelihood of the data, given the model parameters during this iteration	
# 	ttl_loglike = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="loglike", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, cluster_already_open=cluster_already_open)
# 	ttl_loglike

	ttl_loglike = calc_loglike_sp_stratified(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=NULL, spPmat=NULL, min_branchlength=1e-21, return_what="loglike", probs_of_states_at_root=NULL, rootedge=TRUE, sparse=FALSE, printlevel=0, use_cpp=TRUE, input_is_COO=FALSE, spPmat_inputs=NULL, cppSpMethod=3, cluster_already_open=NULL, calc_ancprobs=FALSE, null_range_allowed=TRUE, fixnode=NULL, fixlikes=NULL, inputs=BioGeoBEARS_run_object, allareas=areas, all_states_list=states_list)

	if (print_optim == TRUE)
		{
		LnL = ttl_loglike
		# If the log likelihood is successful, print it
		outvars = adf(t(c(BioGeoBEARS_model_object@params_table$est, LnL)))
		#outvars = cbind(outvars, LnL)
		
		names(outvars) = c(rownames(BioGeoBEARS_model_object@params_table), "LnL")
		print(round(outvars,3))

		#cat(ttl_loglike, "\n", sep="")
		}
	
	return(ttl_loglike)
	}












#######################################################
# chainsaw2
#######################################################
#' Saw a tree off at a particular time before present
#' 
#' This function chops a tree like a hedge-trimmer, cutting straight across at a particular timepoint. 
#' The pieces are returned, as is the leftover tree, with branches shortened appropriately.  Pieces
#' that are mini-trees are returned as ape objects, whereas single branches are just lengths.
#'
#' This function is used during stratification, but could have other uses as well.
#' 
#' @param tr An ape phylo object.
#' @param timepoint The time at which the tree should be "chopped".
#' @param return_pieces Default TRUE, which means pieces should be returned
#' @return \code{chainsaw_result} (a list object with the pieces) or \code{tree_to_chainsaw}, just the leftover tree
#' @export
#' @seealso \code{\link{section_the_tree}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
chainsaw2 <- function(tr, timepoint=10, return_pieces=TRUE)
	{
	# Take a tree and saw it off evenly across a certain timepoint.
	# This removes any tips above the timepoint, and replaces them 
	# with a single tip representing the lineage crossing
	# the timepoint (with a new tip name).

	# Get the tree in a table
	tr_table = prt(tr, printflag=FALSE)
	
	# Find the tips that are less than 10 my old and drop them
	TF_exists_more_recently_than_10mya = tr_table$time_bp < timepoint
	
	# Get the corresponding labels
	labels_for_tips_existing_more_recently_than_10mya = tr_table$label[ TF_exists_more_recently_than_10mya == TRUE ]
	
	###########################################
	# Draft chainsaw function
	###########################################
	# loop through the branches that cross 10 mya
	
	# get a list of the edge start/stops in the phylogeny's edges
	edge_times_bp = get_edge_times_before_present(tr)
	
	# which of these branches cross 10 mya?
	edges_start_earlier_than_10mya = edge_times_bp[, 1] > timepoint
	edges_end_later_than_10mya = edge_times_bp[, 2] <= timepoint
	edges_to_chainsaw = edges_start_earlier_than_10mya + edges_end_later_than_10mya == 2
	
	# then, for each of these edges, figure out how many tips exist descending from it
	nodes_to_chainsaw = tr$edge[, 2][edges_to_chainsaw]
	
	# Take only internal nodes (? why ?)
	numtips = length(tr$tip.label)
	#nodes_to_chainsaw = nodes_to_chainsaw[nodes_to_chainsaw > numtips]
	
	# create a copy of the tree to chainsaw
	tree_to_chainsaw = tr
	
	if (return_pieces == TRUE)
		{
		return_pieces_list = as.list(rep(NA, length(nodes_to_chainsaw)))
		return_pieces_basenames = as.list(rep(NA, length(nodes_to_chainsaw)))
		}
	
	for (i in 1:length(nodes_to_chainsaw))
		{
		# If this is a tip node on the current tree, shorten the branch rather than cut it off
		if (nodes_to_chainsaw[i] <= numtips)
			{
			# Here, chainsaw is cutting an internal node, so extract the sectioned branch before you cut it down
			# (the cutting happens after the forloop)
			# (This is easy, it is just the length of the timeslab;
			#  which you should UPDATE as you move down the tree
			if (return_pieces == TRUE)
				{
				# Record the length of the branch section, and the name of that tip
				# (which is also the name of that base)
				return_pieces_list[[i]] = timepoint
				return_pieces_basenames[[i]] = tr$tip.label[nodes_to_chainsaw[i]]
				}
			# You don't have to do anything else, the chopping of single branches is 
			# covered after the forloop
			#cat("\ni=", i, "	ntips=", length(tree_to_chainsaw$tip.label), sep="")

			} else {
			# Here, it's an internal node, so extract the subtree before you drop it
			tmp_subtree = extract.clade(tr, nodes_to_chainsaw[i])
			#plot(tmp_subtree, root.edge=TRUE)
			# Also, record the branchlength below this node
			branchlength_below_subtree_LCA_node = timepoint - get_max_height_tree(tmp_subtree)
			# Add this to the bottom of the subtree
			tmp_subtree$root.edge = branchlength_below_subtree_LCA_node
			#plot(tmp_subtree, root.edge=TRUE)
			
			# Record the piece, if desired
			if (return_pieces == TRUE)
				{
				# Record the length of the branch section, and the name of that tip
				# (which is also the name of that base)
				return_pieces_list[[i]] = tmp_subtree
				new_labels = sort(tmp_subtree$tip.label)
				basename_after_cutting = paste(new_labels, collapse="|", sep="")
				return_pieces_basenames[[i]] = basename_after_cutting
				}

			#print(tmp_subtree$tip.label)

			tmp_number_of_tips = length(tmp_subtree$tip.label)
			#print(tmp_number_of_tips)
			
			# number of tips to drop = (numtips -1)
			numtips_to_drop = tmp_number_of_tips - 1 
			
			# tips_to_drop
			tmp_labels = tmp_subtree$tip.label
			
			labels_to_drop = tmp_labels[1:numtips_to_drop]
			ordered_labels_to_make_into_new_name = sort(tmp_labels)
			name_new_tip = paste(ordered_labels_to_make_into_new_name, collapse="|", sep="")
			
			# new label
			label_kept_num = length(tmp_labels)
			label_kept = tmp_labels[label_kept_num]
			#new_label = paste("CA_", label_kept, "+", numtips_to_drop, "_tips", sep="")
			new_label = name_new_tip
			tree_to_chainsaw$tip.label[tree_to_chainsaw$tip.label == label_kept] = new_label

			# chop off e.g. 2 of the 3 tips
			tree_to_chainsaw = drop.tip(tree_to_chainsaw, labels_to_drop)
			#cat("\ni=", i, "	ntips=", length(tree_to_chainsaw$tip.label), sep="")
			} # end else
		} # end for loop
	#plot(tree_to_chainsaw)
	#axisPhylo()
	
	tree_to_chainsaw_table = prt(tree_to_chainsaw, printflag=FALSE)
	
	tree_to_chainsaw_table_tips_TF_time_bp_LT_10my = tree_to_chainsaw_table$time_bp < timepoint
	
	
	tmp_edge_lengths =  tree_to_chainsaw_table$edge.length[tree_to_chainsaw_table_tips_TF_time_bp_LT_10my]
	
	times_bp_for_edges_to_chainsaw = tree_to_chainsaw_table$time_bp[tree_to_chainsaw_table_tips_TF_time_bp_LT_10my]
	
	adjustment = times_bp_for_edges_to_chainsaw - timepoint
	
	revised_tmp_edge_lengths = tmp_edge_lengths + adjustment
	
	tree_to_chainsaw_table$edge.length[tree_to_chainsaw_table_tips_TF_time_bp_LT_10my] = revised_tmp_edge_lengths
	
	# revised
	ordered_nodenames = get_nodenums(tree_to_chainsaw)
	parent_branches = get_indices_where_list1_occurs_in_list2(ordered_nodenames, tree_to_chainsaw$edge[,2])
	
	NA_false = is.not.na(tree_to_chainsaw_table$edge.length)
	
	tree_to_chainsaw$edge.length[parent_branches[NA_false]] = tree_to_chainsaw_table$edge.length[NA_false]

	if (return_pieces == TRUE)
		{
		chainsaw_result = NULL
		chainsaw_result$tree_to_chainsaw = tree_to_chainsaw
		chainsaw_result$return_pieces_list = return_pieces_list
		chainsaw_result$return_pieces_basenames = return_pieces_basenames
		class(chainsaw_result) = "chainsaw_result"
		return(chainsaw_result)
		} else {
		return(tree_to_chainsaw)
		}
	}


#######################################################
# get_daughters
#######################################################
#' Get all the direct daughters nodes of a node
#' 
#' @param nodenum The node number to get the daughters of
#' @param t An ape phylo object
#' @return \code{daughter_nodenums} List of the daughter node numbers
#' @export
#' @seealso \code{\link{findall}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_daughters <- function(nodenum, t)
	{
	daughter_edgenums = findall(nodenum, t$edge[,1])
	daughter_nodenums = t$edge[,2][daughter_edgenums]
	return(daughter_nodenums)
	}




# Get indices of all matches to a list
#######################################################
# findall
#######################################################
#' Get indices of all matches to a list
#'
#' Just a handy shortcut function
#' 
#' @param what The item to find
#' @param inlist The list to search in 
#' @return \code{matching_indices} List of the matching indices
#' @export
#' @seealso \code{\link{get_daughters}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
findall <- function(what, inlist)
	{
	TFmatches = inlist == what
	indices = 1:length(inlist)
	matching_indices = indices[TFmatches]
	return(matching_indices)
	}




#######################################################
# prflag
#######################################################
#' Utility function to conditionally print intermediate results
#'
#' Just a handy shortcut function
#' 
#' @param x What to print.
#' @param printflag If TRUE, do the printing
#' @return nothing
#' @export
#' @seealso \code{\link{get_daughters}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
prflag <- function(x, printflag=TRUE)
	{
	# A standard function to print (or not) certain variables,
	#   based on a master printflag
	# This avoids having to comment in/out various code chunks
	#   while debugging.
	if (printflag == TRUE)
		{
		# CAT instead of PRINT if it's a string or numeric
		if (is.character(x))
			{
			cat(x, "\n", sep="")
			}
		if (is.numeric(x))
			{
			cat(x, "\n", sep="")
			} else {
			print(x)
			}
		}
	else
		{
		pass="BLAH"
		}
	}




#######################################################
# get_parent
#######################################################
#' Get the direct parent node of a node
#' 
#' @param nodenum The node number to get the parent of
#' @param t An ape phylo object
#' @return \code{parent_nodenum}The parent node number
#' @export
#' @seealso \code{\link{findall}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_parent <- function(nodenum, t)
	{
	matching_edges = findall(nodenum, t$edge[,2])
	parent_nodenum = t$edge[,1][matching_edges][1]
	return(parent_nodenum)
	}


#######################################################
# get_level
#######################################################
#' Get a node's level in the tree
#'
#' Finds how many nodes deep a node is.
#' 
#' @param nodenum The node number to get the parent of
#' @param t An ape phylo object
#' @param tmplevel A starting level (the function is recursive)
#' @return \code{tmplevel} The level of the node.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_level <- function(nodenum, t, tmplevel=0)
	{
	parent_nodenum = get_parent(nodenum, t)
	if (is.na(parent_nodenum))
		{
		#tmplevel = 0
		return(tmplevel)
		}
	else
		{
		#print(paste("parent_nodenum: ", parent_nodenum, " level: ", tmplevel, sep=""))
		tmplevel = tmplevel + 1
		tmplevel = get_level(parent_nodenum, t, tmplevel)
		return(tmplevel)
		}
	# If an error occurs
	return(NA)
	}


#######################################################
# get_TF_tips
#######################################################
#'  Get TRUE/FALSE for nodes being tips
#'
#' A utility function
#' 
#' @param obj An ape phylo object
#' @return \code{TF_tips} The TRUE/FALSE list for each tip.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_TF_tips <- function(obj)
	{
	# Get TF for nodes being tips
	
	# BIG CHANGE?
	#TF_tips = match_list1_in_list2(1:length(dists_from_root), obj$tip.label)
	TF_tips = match_list1_in_list2(1:length(obj$edge), 1:length(obj$tip.label))
	#TF_tips = obj$tip.label[TF_tips_indices]
	return(TF_tips)
	}


#######################################################
# get_node_ages_of_tips
#######################################################
#' Get the ages of each tip above the root
#'
#' A utility function.
#' 
#' @param obj An ape phylo object
#' @return \code{TF_tips} The age (from the root) of each tip.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_node_ages_of_tips <- function(obj)
	{
	TF_tips = get_TF_tips(obj)
	root_node_num = get_nodenum_structural_root(obj)
	dists_from_root = dist.nodes(obj)[root_node_num, ]
	node_ages_of_tips = dists_from_root[TF_tips]
	return(node_ages_of_tips)
	}


#######################################################
# get_all_node_ages
#######################################################
#' Get the ages of all the nodes in the tree (above the root)
#'
#' A utility function. Use of \code{\link[ape]{dist.nodes}} may be slow.
#' 
#' @param obj An ape phylo object
#' @return \code{TF_tips} The age (from the root) of each node.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_all_node_ages <- function(obj)
	{
	node_ages = dist.nodes(obj)[get_nodenum_structural_root(obj), ]
	return(node_ages)
	}


#######################################################
# get_max_height_tree
#######################################################
#' Get the maximum age of all the nodes (above the root)
#'
#' I.e., the distance of the highest node above the root.  A utility function. 
#' Use of \code{\link[ape]{dist.nodes}} may be slow.
#' 
#' @param obj An ape phylo object
#' @return \code{max_height} The age (from the root) of the highest node.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_max_height_tree <- function(obj)
	{
	max_height = max(get_node_ages_of_tips(obj))
	return(max_height)
	}



#######################################################
# get_edge_times_before_present
#######################################################
#' Get the times of the top and bottom of each edge
#'
#' A utility function. 
#' 
#' @param t An ape phylo object
#' @return \code{edge_times_bp} A 2-column matrix with the age (from the present) of the top and bottom of each edge.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_edge_times_before_present <- function(t)
	{
	#height above root
	hts_at_end_of_branches_aka_at_nodes = t$edge.length
	hts_at_end_of_branches_aka_at_nodes = get_all_node_ages(t)
	h = hts_at_end_of_branches_aka_at_nodes

	# times before present, below (ultrametric!) tips
	# numbers are positive, i.e. in millions of years before present
	#                       i.e. mybp, Ma
	times_before_present = get_max_height_tree(t) - h

	
	# fill in the ages of each node for the edges
	edge_ages = t$edge
	edge_ages[,1] = h[t$edge[,1]]	# bottom of branch
	edge_ages[,2] = h[t$edge[,2]]	# top of branch

	# fill in the times before present of each node for the edges
	edge_times_bp = t$edge
	edge_times_bp[,1] = times_before_present[t$edge[,1]]	# bottom of branch
	edge_times_bp[,2] = times_before_present[t$edge[,2]]	# top of branch
	
	return(edge_times_bp)
	}


#######################################################
# is.not.na
#######################################################
#' Check for not NA
#'
#' A utility function. 
#' 
#' @param x Thing to check for NA
#' @return \code{TRUE} or \code{FALSE}
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
is.not.na <- function(x)
	{
	return(is.na(x) == FALSE)
	}
