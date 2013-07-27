# source("/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_simulate_v1.R")

require("ape")
require("rexpokit")
require("cladoRcpp")


# Qmat contains the d & e probs
# COO_probs_columnar contains the simulation probs


#######################################################
# simulate_biogeog_history
#######################################################
#' Simulate a biogeographical history, given a transition matrix and cladogenesis model
#' 
#' This function simulates a biogeographical history, given a Q transition matrix, a cladogenesis
#' model giving the relative probability of different range inheritance scenarios, a phylogeny, 
#' and a 0-based index value deciding the starting state (which could be randomly generated
#' according to a prior distribution of states).
#' 
#' @param  phy An R \code{phylo} object.
#' @param Qmat A (square, dense) Q transition matrix.  Using a sparse matrix would require writing another function.
#' @param COO_probs_columnar A speciation/cladogenesis transition matrix, in COO-like form, as produced by \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}}.
#' @param index_Qmat_0based_of_starting_state An integer index value, between 0 and \code{(numstates-1)}, which specifies what state will be the starting point for the simulation.
#' @return \code{simulated_states_by_node} A numeric matrix, giving the 0-based index of the state at each node and tip in the simulated history.  Getting a more detailed 
#' history would require a version of stochastic mapping (\cite{Huelsenbeck_etal_2003_stochastic_mapping}, \cite{Bollback_2005}, \cite{Bollback_2006_SIMMAP}), but
#' customized for the nonreversible and cladogenic aspects of biogeographical range evolution models.
#' @export
#' @seealso \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'   @cite Huelsenbeck_etal_2003_stochastic_mapping
#'   @cite Bollback_2005
#'   @cite Bollback_2006_SIMMAP
#' @examples
#' testval=1
#' 
simulate_biogeog_history <- function(phy, Qmat, COO_probs_columnar, index_Qmat_0based_of_starting_state)
	{
	defaults='
	index_Qmat_0based_of_starting_state = 1
	'

	example='
	# Simulate a history
	simulated_states_by_node = simulate_biogeog_history(phy, Qmat, COO_probs_columnar, index_Qmat_0based_of_starting_state=15)
	
	# Save to a PHYLIP-formatted LAGRANGE-type file
	out_geogfn = np(simulated_indexes_to_tipranges_file(simulated_states_by_node, areas_list, states_list, trfn, out_geogfn="lagrange_area_data_file.data"))
	moref(out_geogfn)
	'
	
	# Reorder the phylogeny to pruning-wise order...
	phy2 <- reorder(phy, "pruningwise")
	
	numedges = nrow(phy2$edge)
	ntips = length(phy2$tip.label)
	num_internal_nodes = phy2$Nnode
	numnodes = ntips + num_internal_nodes
	nodenums = c(1:numnodes)
	
	# Number of states in the full model
	numstates = nrow(Qmat)
	
	# Get the ancestral node
	# (The last 2 node in the ancestor column in a pruningwise edge matrix
	#  are the ancestor)
	anc_nodenum = phy2$edge[numedges, 1]

	# Set up the list of simulated states
	simulated_states_by_node = rep(NA, numnodes)
	simulated_states_by_node[anc_nodenum] = index_Qmat_0based_of_starting_state
	
	# simulated_states_by_node_LR_after_speciation
	simulated_states_by_node_LR_after_speciation = matrix(data=NA, nrow=numnodes, ncol=2)
	
	# Label the edges in reverse pruningwise order
	edges_to_visit_j = seq(from=numedges, by=-2, length.out=num_internal_nodes)
	edges_to_visit_i = edges_to_visit_j - 1
	
	# 
	for (i in 1:length(edges_to_visit_i))
		{
		
		# Do left descendant
		left_edge_to_visit = edges_to_visit_i[i]
		right_edge_to_visit = edges_to_visit_j[i]
		starting_nodenum = phy2$edge[left_edge_to_visit,1]
		
		# left and right descendants
		left_desc_nodenum = phy2$edge[left_edge_to_visit,2]
		right_desc_nodenum = phy2$edge[right_edge_to_visit,2]
		
		# Get the starting state for the branches
		starting_state_Qmat_0index = simulated_states_by_node[starting_nodenum]
		
		# Get the descendent states just after speciation
		simulated_states_by_node_LR_after_speciation[starting_nodenum, ] = given_a_starting_state_simulate_split(index_Qmat_0based_of_starting_state=starting_state_Qmat_0index, COO_probs_columnar=COO_probs_columnar, numstates=numstates)
		
		# Get the descendant state and the ends of the two branches
		
		# Simulate end of left branch
		simulated_states_by_node[left_desc_nodenum] = given_a_starting_state_simulate_branch_end(index_Qmat_0based_of_starting_state=simulated_states_by_node_LR_after_speciation[starting_nodenum,1], Qmat=Qmat, branchlength=phy2$edge.length[left_edge_to_visit], all_tips_living=TRUE)
		
		# Simulate end of right branch
		simulated_states_by_node[right_desc_nodenum] = given_a_starting_state_simulate_branch_end(index_Qmat_0based_of_starting_state=simulated_states_by_node_LR_after_speciation[starting_nodenum,2], Qmat=Qmat, branchlength=phy2$edge.length[right_edge_to_visit], all_tips_living=TRUE)		
		}
	
	return(simulated_states_by_node)
	}




#######################################################
# given_a_starting_state_simulate_branch_end
#######################################################
#' Given the state at the start of a branch, simulate the state at the end of the branch
#' 
#' This function simulates a biogeographical history, given a Q transition matrix, a starting state, 
#' and a branch length.  All this involves is exponentiating the Q transition matrix, producing a 
#' P transition probability matrix, and then producing a random draw from this P matrix, conditional on the ancestor.
#'
#' This could be sped up in various ways, if needed.
#' 
#' @param index_Qmat_0based_of_starting_state An integer index value, between 0 and \code{(numstates-1)}, which specifies what state is the
#' starting point for the branch.
#' @param Qmat A (square, dense) Q transition matrix.  Using a sparse matrix would require writing another function.
#' @param branchlength The length of the branch, or branch segment if you are dealing with a stratified phylogeny.
#' @param all_tips_living Currently this is the only assumption.  If, hypothetically, you had a phylogeny with extinct tips (representing
#' the ends of the ranges of fossil taxa), you might want to treat them differently, IF you think that the time-invariant geographic range
#' addition/subtraction process is the same one that made lineages go extinct (it could be something else, e.g. mass extinction).  False attribution
#' of extinctions to the range loss process will dramatically elevate the rate of range loss, and also range expansion to compensate, and the
#' resulting high rates can substantially degrade inference (\cite{Matzke_Maguire_2011_SVP}).
#' @return \code{state_desc} 0-based index of the descendant state (just before cladogenesis, if below a node).
#' @export
#' @seealso \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'   @cite Matzke_Maguire_2011_SVP
#' @examples
#' testval=1
#' 
given_a_starting_state_simulate_branch_end <- function(index_Qmat_0based_of_starting_state=1, Qmat, branchlength=1, all_tips_living=TRUE)
	{
	defaults='
	# This is the 0-based index for the full Qmat; however, once we reduce the 16x16 Qmat
	# to a 15x15 Pmat (because the null range can be neither a starting, nor ending, state here)
	# then index_Qmat_0based_of_starting_state can be used without modification
	index_Qmat_0based_of_starting_state = 1
	all_tips_living = TRUE
	'
	# Exponentiate the rate matrix Q
	Pmat = expokit_dgpadm_Qmat2(Qmat=Qmat, times=branchlength, transpose_needed=TRUE)

	# If we assume that all the tips are living, we have to do a correction for impossible
	# results
	if (all_tips_living == TRUE)
		{
		# Remove row 1 / col 1
		tmpPmat = Pmat[-1,-1]
		
		# Divide all the cells in a row by the sum of the row
		#rowsums_Pmat = rowSums(tmpPmat)
		# Pmat2 = apply(X=tmpPmat, MARGIN=2, FUN="/", rowsums_Pmat)
		
		# Only need to do the sum operation on the 1 row descendending
		# from the ancestor
		descendent_probs_list = tmpPmat[index_Qmat_0based_of_starting_state,] / sum(tmpPmat[index_Qmat_0based_of_starting_state,])
		
		prob_starts = rep(0, length(descendent_probs_list))
		prob_ends = rep(0, length(descendent_probs_list))
		
		current_startprob = 0
		current_endprob = descendent_probs_list[1]
		for (i in 1:length(prob_starts))
			{
			prob_starts[i] = current_startprob
			current_endprob = current_startprob + descendent_probs_list[i]
			prob_ends[i] = current_endprob
			current_startprob = current_endprob
			}
			
		prob_bins = adf(cbind(prob_starts, prob_ends))
		#print(prob_bins)
	
	
	
		#######################################################
		# Simulate the descendent split
		#######################################################
		random_draw = runif(n=1, min=0, max=1)
		
		random_val_above_bin_mins_TF = random_draw > prob_bins$prob_starts
		random_val_below_bin_maxes_TF = random_draw <= prob_bins$prob_ends
		bin_match_TF = (random_val_above_bin_mins_TF + random_val_below_bin_maxes_TF) == 2
		
		# Descendent state, in Qmat 0-based index form
		possible_states_Qmat_0based_indexes = 1:length(descendent_probs_list)
		state_desc = possible_states_Qmat_0based_indexes[bin_match_TF]
		
		return(state_desc)
		
		} else {
		stop("ERROR: Need to implement simulation with extinct tips")
		}
	
	return(NA)
	}



#######################################################
# given_a_starting_state_simulate_split
#######################################################
#' Given the state just below a node, simulate the states after speciation
#' 
#' This function simulates a biogeographical history during a speciation/cladogenesis range inheritance event,
#' given a cladogenesis probability transition matrix and a starting state.
#' 
#' @param index_Qmat_0based_of_starting_state An integer index value, between 0 and \code{(numstates-1)}, which specifies what state is the
#' starting point for the branch.
#' @param COO_probs_columnar A speciation/cladogenesis transition matrix, in COO-like form, as produced by \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}}.
#' @param numstates The number of states/geographic ranges.
#' @return \code{split_desc} 0-based indices of the descendant states in the two daughters.
#' @export
#' @seealso \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}}, \code{\link{rcpp_calc_rowsums_for_COOweights_columnar}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'   @cite Matzke_Maguire_2011_SVP
#' @examples
#' testval=1
#' 
given_a_starting_state_simulate_split <- function(index_Qmat_0based_of_starting_state=1, COO_probs_columnar, numstates)
	{
	defaults='
	# Note that the speciation matrix is always missing the original state 0 (null range);
	# Thus the 0th state is actually original state 1 in the Qmat (starting with 0)
	index_Qmat_0based_of_starting_state = 1
	
	for (i in 1:1000)
		{
		given_a_starting_state_simulate_split(index_Qmat_0based_of_starting_state=1, COO_probs_columnar, numstates=16)
		}
	'
	
	# If there are 16 states, there are 15 non-null rowSums
	Rsp_rowsums = rcpp_calc_rowsums_for_COOweights_columnar(COO_probs_columnar, numstates=numstates)

	# Load up the speciation matrix
	# These indexes are 0-14, i.e. 1-15 but 0-based
	RCOO_probs_columnar_anc_i_list = COO_probs_columnar[[1]]
	RCOO_left_i_list = COO_probs_columnar[[2]]
	RCOO_right_j_list = COO_probs_columnar[[3]]
	RCOO_probs_list = COO_probs_columnar[[4]]

	# Number of nonzero cells (i.e. in the speciation matrix
	num_nonzero_cells_in_sp_matrix = length(RCOO_probs_columnar_anc_i_list)
	
	# Get the ancestors that match the input ancestor
	# COO_probs_columnar have 0-based indexes
	anc_match_TF = RCOO_probs_columnar_anc_i_list == (index_Qmat_0based_of_starting_state-1)
	num_nonzero_descendent_splits = sum(anc_match_TF)
	
	#print(num_nonzero_descendent_splits)
	
	# To refer to the correct index in Rsp_rowsums, add 1 to the 0-based index
	# of 15 non-null possible ancestral states
	descendent_probs_list = RCOO_probs_list[anc_match_TF] / Rsp_rowsums[index_Qmat_0based_of_starting_state-1+1]
	
	

	#######################################################
	# Recursively num the probability bin starts and ends
	#######################################################
	prob_starts = rep(0, num_nonzero_descendent_splits)
	prob_ends = rep(0, num_nonzero_descendent_splits)
	
	current_startprob = 0
	current_endprob = descendent_probs_list[1]
	for (i in 1:length(prob_starts))
		{
		prob_starts[i] = current_startprob
		current_endprob = current_startprob + descendent_probs_list[i]
		prob_ends[i] = current_endprob
		current_startprob = current_endprob
		}
		
	prob_bins = adf(cbind(prob_starts, prob_ends))
	#print(prob_bins)



	#######################################################
	# Simulate the descendent split
	#######################################################
	random_draw = runif(n=1, min=0, max=1)
	
	random_val_above_bin_mins_TF = random_draw > prob_bins$prob_starts
	random_val_below_bin_maxes_TF = random_draw <= prob_bins$prob_ends
	bin_match_TF = (random_val_above_bin_mins_TF + random_val_below_bin_maxes_TF) == 2
	
	
	#######################################################
	# Find the descendents of the split
	#######################################################
	left_desc_index_sp_0based = RCOO_left_i_list[anc_match_TF][bin_match_TF]
	right_desc_index_sp_0based = RCOO_right_j_list[anc_match_TF][bin_match_TF]
	
	# Convert spmat 0-based indexing to Qmat 0-based indexing
	left_desc_index_0based = left_desc_index_sp_0based + 1 
	right_desc_index_0based = right_desc_index_sp_0based + 1 
	
	split_desc = c(left_desc_index_0based, right_desc_index_0based)
	#cat(split_desc)
	#cat("\n")
	
	return(split_desc)
	
	}




# Convert simulated Qmat 0-based indexes to a tipranges object
#######################################################
# simulated_indexes_to_tipranges_object
#######################################################
#' Convert simulated Qmat 0-based indexes to a tipranges object
#' 
#' This function takes simulated state indices (ranging from 0 to numstates-1, i.e. number of 
#' possible geographic ranges-1) and converts them to a tipranges object.  This can then be
#' converted into a C++-\code{LAGRANGE}-style PHYLIP geographic ranges file.
#'
#' @param simulated_states_by_node The simulated states/geographic ranges, in 0-based index form, ordered as the tips & nodes are ordered
#' in a \code{pruningwise}-ordered \code{phylo} object in \code{APE}.
#' @param areas_list A list of the desired area names/abbreviations/letters.
#' @param states_list A list of the possible states/geographic ranges, in 0-based index form.
#' @param trfn The filename of the source Newick tree.
#' @return \code{tipranges_object} An object of class \code{tipranges}.
#' @export
#' @seealso \code{\link{define_tipranges_object}}, \code{\link{getareas_from_tipranges_object}}, \code{\link{simulated_indexes_to_tipranges_file}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#'	 @cite SmithRee2010_CPPversion
#' @examples
#' testval=1
#' 
simulated_indexes_to_tipranges_object <- function(simulated_states_by_node, areas_list, states_list, trfn)
	{
	phy = read.tree(trfn)
	phy2 <- reorder(phy, "pruningwise")
	
	tiplabels = phy2$tip.label
	ntips = length(tiplabels)
	
	tipstates = simulated_states_by_node[1:ntips]
	
	tipranges_table = matrix(data=0, nrow=ntips, ncol=length(areas_list))
	
	for (i in 1:nrow(tipranges_table))
		{
		state_Qmat_0index = tipstates[i]
		areas_occupied_1index = 1 + states_list[[state_Qmat_0index+1]]
		
		tipranges_table[i, areas_occupied_1index] = 1
		}
	
	tipranges_table = as.data.frame(tipranges_table)
	
	rownames(tipranges_table) = phy2$tip.label
	colnames(tipranges_table) = areas_list
	
	tipranges_object = define_tipranges_object()
	tipranges_object@df = tipranges_table
	
	return(tipranges_object)
	}



# Convert simulated Qmat 0-based indexes to a tipranges file
#######################################################
# simulated_indexes_to_tipranges_file
#######################################################
#' Convert simulated Qmat 0-based indexes to a tipranges file
#' 
#' This function takes simulated state indices (ranging from 0 to numstates-1, i.e. number of 
#' possible geographic ranges-1) and converts them to a C++-\code{LAGRANGE}-style PHYLIP geographic ranges file.
#'
#' @param simulated_states_by_node The simulated states/geographic ranges, in 0-based index form, ordered as the tips & nodes are ordered
#' in a \code{pruningwise}-ordered \code{phylo} object in \code{APE}.
#' @param areas_list A list of the desired area names/abbreviations/letters.
#' @param states_list A list of the possible states/geographic ranges, in 0-based index form.
#' @param trfn The filename of the source Newick tree.
#' @param out_geogfn The output filename.
#' @return \code{out_geogfn} The output filename.
#' @export
#' @seealso \code{\link{define_tipranges_object}}, \code{\link{getareas_from_tipranges_object}}, \code{\link{simulated_indexes_to_tipranges_object}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#'	 @cite SmithRee2010_CPPversion
#' @examples
#' testval=1
#' 
simulated_indexes_to_tipranges_file <- function(simulated_states_by_node, areas_list, states_list, trfn, out_geogfn="lagrange_area_data_file.data")
	{
	
	tipranges_object = simulated_indexes_to_tipranges_object(simulated_states_by_node, areas_list, states_list, trfn)
	save_tipranges_to_LagrangePHYLIP(tipranges_object, lgdata_fn=out_geogfn)
	
	return(out_geogfn)
	}










#######################################################
# Processing simulations
#######################################################

# Get ML states
#######################################################
# get_ML_states
#######################################################
#' Get ML states from a BioGeoBEARS model results list
#' 
#' This function extracts the ML states from the results list produced by \code{\link{bears_2param_standard_fast}}
#' or a similar ML search function.
#'
#' Currently, the scaled conditional probabilities are used to determine the optimum states.  However, this is 
#' not strictly correct, as these use only tips-down information (\cite{Felsenstein2004}; see also this post by Revell: \url{http://blog.phytools.org/2013/03/marginal-ancestral-state-reconstruction.html}).  This is what \code{LAGRANGE} seems to do
#' when reporting ancestral states, also (personal observation, perhaps imperfect, especially if the scaled conditional 
#' likelihoods and the marginal ancestral state probabilities turn out to be
#' very close). What is desired is the marginal ancestral state
#' reconstructions.  Most authors discuss ML ancestral state reconstruction as being a matter of re-rooting the tree at each node, 
#' yielding the marginal estimate for that node, conditional on the rest of the tree.  However, this procedure assumes a 
#' time-reversible model on both branches and cladogenesis events, and we have neither in biogeography.  Probably, the solution is just
#' an up-pass from the root, calculating the probabilities on the forward model and multiplying by likelihoods from the downpass.  
#' However, this has not yet been implemented.
#'
#' @param relprobs_matrix A relative probabilities matrix returned by \code{\link{bears_2param_standard_fast}} or a similar function. 
#' The user should specify WHICH matrix in the results_object -- i.e., scaled conditional likelihoods on downpass or uppass, or 
#' actual marginal probabilities of ancestral states.  (The latter is the main thing of interest.)  This specification 
#' is done via e.g. \code{relprobs_matrix = results_object$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS}.
#' @param unlist_TF Unlist the output? Default TRUE.
#' @return \code{inf_statesvec} The inferred vector of states.
#' @export
#' @seealso \code{\link{get_ML_probs}}, \code{\link{bears_2param_standard_fast}}, \code{\link{get_ML_state_indices}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://blog.phytools.org/2013/03/marginal-ancestral-state-reconstruction.html} 
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' 
get_ML_states <- function(relprobs_matrix, unlist_TF=TRUE)
	{
	inf_statesvec = as.list(rep(-1, times=nrow(relprobs_matrix)))
	
	state_indexes_0based = seq(0, ncol(relprobs_matrix)-1, 1)
	
	#state_indexes_0based_matrix = matrix(data=state_indexes_0based, nrow=nrow(relprobs_matrix), ncol=ncol(relprobs_matrix), byrow=TRUE)
	
	
	maxprobs = apply(X=relprobs_matrix, MARGIN=1, FUN=max)
	maxprobs
	
	# Which match max
	for (rownum in 1:nrow(relprobs_matrix))
		{
		match_max_TF = relprobs_matrix[rownum,] == maxprobs[rownum]
		
		# Fix NA
		#if (is.na(match_max_TF))
		if (any(sapply(X=match_max_TF, FUN=is.na)))
			{
			inf_statesvec[[rownum]] = NA
			next()
			}
		
		# Number of matches
		nummatches = sum(match_max_TF)
		nummatches
		#print(nummatches)
		
		if (nummatches == 1)
			{
			#ML_probs_vec[[rownum]] = c(relprobs_matrix[rownum,][match_max_TF])
			inf_statesvec[[rownum]] = c(state_indexes_0based[match_max_TF])
			} 
		
		if (nummatches > 1)
			{
			cat("\nNote: multiple states tied\n")
			# unlist_TF
			if (unlist_TF == TRUE)
				{
				cat("\nNote: picking the first state in the tie; use unlist_TF=FALSE to see all states.\n")
				inf_statesvec[[rownum]] = c(state_indexes_0based[match_max_TF][1])
				} else {
				inf_statesvec[[rownum]] = c(state_indexes_0based[match_max_TF])
				}
			} 
		
		}
	
	if (unlist_TF == TRUE)
		{
		inf_statesvec = unlist(inf_statesvec)
		}
	
	
	# Return the list of states
	# (as 0-based indices)
	return(inf_statesvec)
	}








# Get ML probs of ML states
#######################################################
# get_ML_probs
#######################################################
#' Get the probability of the ML state for each node, from a BioGeoBEARS model results list
#' 
#' This function extracts the probability of the ML states from the results list produced by \code{\link{bears_2param_standard_fast}}
#' or a similar ML search function.
#'
#' This is useful for displaying e.g. pie charts of the probability of the ML ancestral state at each node.
#'
#' Note, though, that it is somewhat peculiar and arbitrary to focus on the ancestral states just at nodes, particularly in the context of
#' fossils with time ranges and geographic ranges.
#'
#' @param relprobs_matrix A relative probabilities matrix returned by \code{\link{bears_2param_standard_fast}} or a similar function. 
#' The user should specify WHICH matrix in the results_object -- i.e., scaled conditional likelihoods on downpass or uppass, or 
#' actual marginal probabilities of ancestral states.  (The latter is the main thing of interest.)  This specification 
#' is done via e.g. \code{relprobs_matrix = results_object$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS}..
#' @param unlist_TF Unlist the output? Default TRUE.
#' @return \code{inf_probsvec} The inferred vector of probabilities of ML states.
#' @export
#' @seealso \code{\link{get_ML_probs}}, \code{\link{bears_2param_standard_fast}}, \code{\link{get_ML_state_indices}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://blog.phytools.org/2013/03/marginal-ancestral-state-reconstruction.html} 
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' 
get_ML_probs <- function(relprobs_matrix, unlist_TF=TRUE)
	{
	inf_probsvec = as.list(rep(-1, times=nrow(relprobs_matrix)))
	
	#state_indexes_0based = seq(0, ncol(relprobs_matrix)-1, 1)
	
	#state_indexes_0based_matrix = matrix(data=state_indexes_0based, nrow=nrow(relprobs_matrix), ncol=ncol(relprobs_matrix), byrow=TRUE)
	
	
	maxprobs = apply(X=relprobs_matrix, MARGIN=1, FUN=max)
	maxprobs
	
	# Which match max
	for (rownum in 1:nrow(relprobs_matrix))
		{
		match_max_TF = relprobs_matrix[rownum,] == maxprobs[rownum]

		# Fix NA
		#if (is.na(match_max_TF))
		if (any(sapply(X=match_max_TF, FUN=is.na)))
			{
			inf_probsvec[[rownum]] = NA
			next()
			}
	
		# Number of matches
		nummatches = sum(match_max_TF)
		nummatches

		if (nummatches == 1)
			{
			#ML_probs_vec[[rownum]] = c(relprobs_matrix[rownum,][match_max_TF])
			inf_probsvec[[rownum]] = c(relprobs_matrix[rownum,][match_max_TF])
			} 
		
		if (nummatches > 1)
			{
			cat("\nNOTE: multiple states tied\n")
			# unlist_TF
			if (unlist_TF == TRUE)
				{
				cat("\nNote: in get_ML_probs(), picking the first state in the tie; use unlist_TF=FALSE to see all states.\n")
				inf_probsvec[[rownum]] = c(relprobs_matrix[rownum,][match_max_TF][1])
				} else {
				inf_probsvec[[rownum]] = c(relprobs_matrix[rownum,][match_max_TF])
				}
			} 

		}
	
	if (unlist_TF == TRUE)
		{
		inf_probsvec = unlist(inf_probsvec)
		}
	
	# Return the list of states
	return(inf_probsvec)
	}

# Load the simulation information
# These simstates are 0-based
#######################################################
# get_simstates
#######################################################
#' Load the simulation information from an underscore delimited text string.
#' 
#' If the simulated states are stored in a big text file, it can be useful to store
#' them as a single string in a single cell per row, so that the number of columns
#' doesn't have to change with each different-sized tree. This function extracts the 
#' simulated states from this format.
#' 
#' @param simhist_row A row from a table, which must have a column named \code{simulated_states_by_node_txt}.
#' @return \code{simulated_states_by_node} A numeric vector of 0-based state indices.
#' @export
#' @seealso \code{\link[utils]{read.table}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' 
get_simstates <- function(simhist_row)
	{
	simulated_states_by_node_txt = simhist_row$simulated_states_by_node_txt
	simulated_states_by_node = strsplit(simulated_states_by_node_txt, split="_")[[1]]
	simulated_states_by_node = as.integer(simulated_states_by_node)
	return(simulated_states_by_node)	
	}

#######################################################
# infprobs_to_probs_of_each_area
#######################################################
#' Convert probabilities of each state, to the probabilities of presence in each area
#' 
#' Biogeographic inference in LAGRANGE and DIVA has focused heavily on inference of the 
#' exact ancestral state/geographic range.  However, when the state space is large, there
#' is often considerable uncertainty in the exact ancestral range.  Even the ancestral state
#' that confers the maximum likelihood on the data, and thus is the most probable ancestor, 
#' may have less than 50% probability, or even less (25%, 5%...), depending on the size of 
#' the state space.  This function converts the probability of specific states/geographic 
#' ranges into the probability of presence/absence in each area.  This can typically be 
#' inferred with much higher confidence.
#' 
#' @param relprobs_matrix A relative probabilities matrix returned by \code{\link{bears_2param_standard_fast}} or a similar function. 
#' The user should specify WHICH matrix in the results_object -- i.e., scaled conditional likelihoods on downpass or uppass, or 
#' actual marginal probabilities of ancestral states.  (The latter is the main thing of interest.)  This specification 
#' is done via e.g. \code{relprobs_matrix = results_object$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS}..
#' @param states_list A list of the possible states/geographic ranges, in 0-based index form.
#' @return \code{area_probs} The probability of presence in each area.
#' @export
#' @seealso \code{\link{bears_2param_standard_fast}}, \code{\link{get_ML_states}}, \code{\link{get_ML_probs}}, \code{\link{infprobs_to_probs_of_each_area_from_relprobs}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' 
infprobs_to_probs_of_each_area <- function(relprobs_matrix, states_list)
	{
	# Get the areas from the states list
	areas = unique(unlist(states_list))
	# Remove null range
	areas = areas[!is.na(areas)]
	
	# Area probabilities table for each node
	area_probs = matrix(0, nrow=nrow(relprobs_matrix), ncol=length(areas))

	# Go through the states
	for (i in 1:length(states_list))
		{
		if (is.na(states_list[[i]] && (length(states_list[[i]]==1)) ))
			{
			next()
			} else {
			# Convert 0-based states to 1-based states
			areas_in_this_state = states_list[[i]] + 1

			# Go through the rows (the ancestral nodes)
			for (rownum in 1:nrow(relprobs_matrix))
				{
				# Prob of a particular state
				tmpprob = relprobs_matrix[rownum,i]
				
				# Every area in this state gets this probability
				area_probs[rownum,areas_in_this_state] = area_probs[rownum,areas_in_this_state] + tmpprob
				}
			}
		}
	
	return(area_probs)
	}



#######################################################
# infprobs_to_probs_of_each_area_from_relprobs
#######################################################
#' Convert relative probabilities matrix to the probabilities of presence in each area
#' 
#' Biogeographic inference in LAGRANGE and DIVA has focused heavily on inference of the 
#' exact ancestral state/geographic range.  However, when the state space is large, there
#' is often considerable uncertainty in the exact ancestral range.  Even the ancestral state
#' that confers the maximum likelihood on the data, and thus is the most probable ancestor, 
#' may have less than 50% probability, or even less (25%, 5%...), depending on the size of 
#' the state space.  This function converts the probability of specific states/geographic 
#' ranges into the probability of presence/absence in each area.  This can typically be 
#' inferred with much higher confidence.
#' 
#' @param relprobs_matrix A matrix with nrows for nodes and columns for states, with each
#' cell holding the relative probability of that state at that node.
#' @param states_list A list of the possible states/geographic ranges, in 0-based index form.
#' @return \code{area_probs} The probability of presence in each area.
#' @export
#' @seealso \code{\link{bears_2param_standard_fast}}, \code{\link{get_ML_states}}, \code{\link{get_ML_probs}}, \code{\link{infprobs_to_probs_of_each_area}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' 
infprobs_to_probs_of_each_area_from_relprobs <- function(relprobs_matrix, states_list)
	{
	# Get the areas from the states list
	areas = unique(unlist(states_list))
	# Remove null range
	areas = areas[!is.na(areas)]
	
	# Area probabilities table for each node
	area_probs = matrix(0, nrow=nrow(relprobs_matrix), ncol=length(areas))

	# Go through the states
	for (i in 1:length(states_list))
		{
		if (is.na(states_list[[i]] && (length(states_list[[i]]==1)) ))
			{
			next()
			} else {
			# Convert 0-based states to 1-based states
			areas_in_this_state = states_list[[i]] + 1

			# Go through the rows (the ancestral nodes)
			for (rownum in 1:nrow(relprobs_matrix))
				{
				# Prob of a particular state
				tmpprob = relprobs_matrix[rownum,i]
				
				# Every area in this state gets this probability
				area_probs[rownum,areas_in_this_state] = area_probs[rownum,areas_in_this_state] + tmpprob
				}
			}
		}
	
	return(area_probs)
	}



#######################################################
# simstates_to_probs_of_each_area
#######################################################
#' Convert simulated states to probabilities of each area
#' 
#' Basically this function assigns probability 1 to occupied areas according to the simulated state for a node, and 
#' probability 0 for the other areas.  These data -- the simulated truth -- can then be compared to the inferred 
#' probabilities of presence in each area, from \code{\link{infprobs_to_probs_of_each_area}}.
#' 
#' @param simulated_states_by_node The simulated states by node (0-based indices).
#' @param states_list A list of the possible states/geographic ranges, in 0-based index form.
#' @param relprobs_matrix A relative probabilities matrix returned by \code{\link{bears_2param_standard_fast}} or a similar function. 
#' The user should specify WHICH matrix in the results_object -- i.e., scaled conditional likelihoods on downpass or uppass, or 
#' actual marginal probabilities of ancestral states.  (The latter is the main thing of interest.)  This specification 
#' is done via e.g. \code{relprobs_matrix = results_object$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS}..
#' @return \code{area_probs} The probability of presence in each area.
#' @export
#' @seealso \code{\link{simulate_biogeog_history}}, \code{\link{infprobs_to_probs_of_each_area}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#'
simstates_to_probs_of_each_area <- function(simulated_states_by_node, states_list, relprobs_matrix)
	{
	# Convert simulated states from 0-based to 1-based
	simulated_states_by_node_1based = simulated_states_by_node + 1

	# Get the areas from the states list
	areas = unique(unlist(states_list))
	# Remove null range
	areas = areas[!is.na(areas)]
	
	# Area probabilities table for each node
	area_probs = matrix(0, nrow=nrow(relprobs_matrix), ncol=length(areas))
	
	# Go through the nodes (rows)
	for (rownum in 1:length(simulated_states_by_node_1based))
		{
		# Convert 0-based states to 1-based states
		tmpstate = simulated_states_by_node_1based[rownum]
		tmp_states_list = states_list[[tmpstate]] + 1
		
		if (is.na(tmp_states_list) && (length(tmp_states_list)==1) )
			{
			next()
			} else {
			areas_in_this_state = tmp_states_list
			# tmpprob
			tmpprob = 1
			
			# Every area in this state gets this probability
			area_probs[rownum,areas_in_this_state] = area_probs[rownum,areas_in_this_state] + tmpprob
			}
		}
	
	return(area_probs)
	}



# Get the probabilities of the true (simulated) states
#######################################################
# get_infprobs_of_simstates
#######################################################
#' Get the probabilities of the true (simulated) states
#' 
#' Basically this function assigns probability 1 to the simulated state/geographic range, and 
#' probability 0 for the other states/geographic ranges.  These data -- the simulated truth -- can then be compared to the inferred 
#' probabilities for the states, from e.g. \code{\link{get_ML_probs}}.
#' 
#' @param relprobs_matrix A relative probabilities matrix returned by \code{\link{bears_2param_standard_fast}} or a similar function. 
#' The user should specify WHICH matrix in the results_object -- i.e., scaled conditional likelihoods on downpass or uppass, or 
#' actual marginal probabilities of ancestral states.  (The latter is the main thing of interest.)  This specification 
#' is done via e.g. \code{relprobs_matrix = results_object$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS}..
#' @param simhist_row A row from a table, which must have a column named \code{simulated_states_by_node_txt}.
#' @return \code{infprobs_of_simstates} The probability of each state at each node (all 1s and 0s).
#' @export
#' @seealso \code{\link{simulate_biogeog_history}}, \code{\link{infprobs_to_probs_of_each_area}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#'
get_infprobs_of_simstates <- function(relprobs_matrix, simhist_row)
	{
	infprobs_of_simstates = as.list(rep(-1, times=nrow(relprobs_matrix)))
	
	# Get true simulated states
	simulated_states_by_node = get_simstates(simhist_row)
	
	# Get the probs of these states
	# Which match these states
	for (rownum in 1:nrow(relprobs_matrix))
		{
		infprobs_of_simstates[[rownum]] = relprobs_matrix[rownum,(1+simulated_states_by_node[rownum])]
		}
	
	return(infprobs_of_simstates)
	}


# Get the inferred parameters
#######################################################
# get_infparams_optimx
#######################################################
#' Get the inferred parameters from an ML optimization
#' 
#' This function extracts the ML parameter values, and associated statistics and codes, from the 
#' \code{relprobs_matrix} returned by \code{\link{bears_2param_standard_fast}} and similar functions.
#'
#' The function has subroutines for recognizing a variety of currently-implemented models, assuming they
#' used \code{\link[optimx]{optimx}} internally to do the ML search.  New models 
#' would require addition of new subroutines.
#'
#' \code{\link{get_infparams_optimx}} and \code{\link{get_infparams_optimx_nosim}} differ only in the format of the filenames.
#' 
#' @param results_object The results returned by \code{\link{bears_2param_standard_fast}} or a similar function.
#' @param inffn The filename holding the results_object, which specifies which model was run.
#' @return \code{infparams} The vector of inferred parameters.
#' @export
#' @seealso \code{\link{get_infparams_optimx_nosim}}, \code{\link{bears_2param_standard_fast}}, \code{\link{get_inf_LgL_etc_optimx}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#'
get_infparams_optimx <- function(results_object, inffn)
	{
	parvec = results_object$optim_result$par$par
	# function_name = simhist_row$run
	
	numpars = length(parvec)
	
	infparams = list()
	infparams$d = parvec[1]
	infparams$e = parvec[2]

	if (grepl(pattern="2param_standard_fast_symOnly_sim", x=inffn))
		{
		infparams$j = 0
		
		# Calculate the actual weights
		ysv = 1-infparams$j
		v = ysv * 0.0
		ys = ysv - v

		infparams$v = v
		infparams$ys = ys
		infparams$maxent01 = 1
		infparams$maxent01v = 1
		}

	if (grepl(pattern="2param_standard_fast_sim", x=inffn))
		{
		infparams$j = 0
		
		# Calculate the actual weights
		ysv = 1-infparams$j
		v = ysv * 0.5
		ys = ysv - v

		infparams$v = v
		infparams$ys = ys
		infparams$maxent01 = 0.0001
		infparams$maxent01v = 0.0001
		}
	
	if (grepl(pattern="3param_standard_fast_sim", x=inffn))
		{
		infparams$j = parvec[3]
		
		# Calculate the actual weights
		ysv = 1-infparams$j
		v = ysv * 0.5
		ys = ysv - v

		infparams$v = v
		infparams$ys = ys
		infparams$maxent01 = 0.0001
		infparams$maxent01v = 0.0001
		}

	if (grepl(pattern="4param_standard_fast_sim", x=inffn))
		{
		infparams$j = parvec[3]
		
		# Calculate the actual weights
		ysv = 1-infparams$j
		v = ysv * parvec[4]
		ys = ysv - v

		infparams$v = v
		infparams$ys = ys
		infparams$maxent01 = 0.0001
		infparams$maxent01v = 0.0001
		}

	# unique(simhist_records2$run)
	if (grepl(pattern="5param_standard_fast_v_sim", x=inffn))
		{
		infparams$j = parvec[3]
		
		# Calculate the actual weights
		ysv = 1-infparams$j
		v = ysv * parvec[4]
		ys = ysv - v

		infparams$v = v
		infparams$ys = ys
		infparams$maxent01 = 0.0001
		infparams$maxent01v = parvec[5]
		}

	if (grepl(pattern="5param_standard_fast_sim", x=inffn))
		{
		infparams$j = parvec[3]
		
		# Calculate the actual weights
		ysv = 1-infparams$j
		v = ysv * parvec[4]
		ys = ysv - v

		infparams$v = v
		infparams$ys = ys
		infparams$maxent01 = parvec[5]
		infparams$maxent01v = 0.0001
		}

	if (grepl(pattern="6param_standard_fast_sim", x=inffn))
		{
		infparams$j = parvec[3]
		
		# Calculate the actual weights
		ysv = 1-infparams$j
		v = ysv * parvec[4]
		ys = ysv - v

		infparams$v = v
		infparams$ys = ys
		infparams$maxent01 = parvec[5]
		infparams$maxent01v = parvec[6]
		}
	
	return(infparams)
	}



# Get the inferred parameters
#######################################################
# get_inf_LgL_etc_optimx
#######################################################
#' Get the inferred parameters from a results object (utility function)
#' 
#' This function extracts the ML parameter values from the 
#' \code{results_object} returned by \code{\link{bears_2param_standard_fast}} and similar functions.
#'
#' This is primarily a utility function for \code{\link{get_infparams_optimx}}.
#'
#' @param results_object The results returned by \code{\link{bears_2param_standard_fast}} or a similar function.
#' @return \code{infparams} The vector of inferred parameters.
#' @export
#' @seealso \code{\link{bears_2param_standard_fast}}, \code{\link{get_infparams_optimx}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#'
get_inf_LgL_etc_optimx <- function(results_object)
	{
	inference_stats = results_object$optim_result[2:length(results_object$optim_result)]
	
	# Rename "fvalues" to "LnL"
	names(inference_stats)[names(inference_stats)=="fvalues"] = "LnL"
	
	return(inference_stats)
	}



# Get the inferred parameters
#######################################################
# get_infparams_optimx_nosim
#######################################################
#' Get the inferred parameters from an ML optimization (different filenames)
#' 
#' Like \code{\link{get_infparams_optimx}}, this function extracts the ML parameter values, and associated statistics and codes, from the 
#' \code{results_object} returned by \code{\link{bears_2param_standard_fast}} and similar functions.
#'
#' The function has subroutines for recognizing a variety of currently-implemented models, assuming they
#' used \code{\link[optimx]{optimx}} internally to do the ML search.  New models 
#' would require addition of new subroutines.
#'
#' \code{\link{get_infparams_optimx}} and \code{\link{get_infparams_optimx_nosim}} differ only in the format of the filenames.
#'
#' @param results_object The results returned by \code{\link{bears_2param_standard_fast}} or a similar function.
#' @param inffn The filename holding the results_object, which specifies which model was run.
#' @return \code{infparams} The vector of inferred parameters.
#' @export
#' @seealso \code{\link{get_infparams_optimx}}, \code{\link{bears_2param_standard_fast}}, \code{\link{get_inf_LgL_etc_optimx}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#'
get_infparams_optimx_nosim <- function(results_object, inffn)
	{
	parvec = results_object$optim_result$par$par
	
	numpars = length(parvec)
	
	infparams = list()
	infparams$d = parvec[1]
	infparams$e = parvec[2]

	if (grepl(pattern="2param_standard_fast_symOnly", x=inffn))
		{
		infparams$j = 0
		
		# Calculate the actual weights
		ysv = 1-infparams$j
		v = ysv * 0.0
		ys = ysv - v

		infparams$v = v
		infparams$ys = ys
		infparams$maxent01 = 1
		infparams$maxent01v = NA
		return(infparams)
		}

	if (grepl(pattern="2param_standard_fast", x=inffn))
		{
		infparams$j = 0
		
		# Calculate the actual weights
		ysv = 1-infparams$j
		v = ysv * 0.5
		ys = ysv - v

		infparams$v = v
		infparams$ys = ys
		infparams$maxent01 = 0.0001
		infparams$maxent01v = 0.0001
		return(infparams)
		}

	if (grepl(pattern="2param_DIVA_fast", x=inffn))
		{
		infparams$j = 0
		
		# Calculate the actual weights
		ysv = 1-infparams$j
		v = ysv * 0.5
		ys = ysv - v

		infparams$v = v
		infparams$ys = ys
		infparams$maxent01 = 0.0001
		infparams$maxent01v = 0.5
		return(infparams)
		}

	
	if (grepl(pattern="3param_standard_fast", x=inffn))
		{
		infparams$j = parvec[3]
		
		# Calculate the actual weights
		ysv = 1-infparams$j
		v = ysv * 0.5
		ys = ysv - v

		infparams$v = v
		infparams$ys = ys
		infparams$maxent01 = 0.0001
		infparams$maxent01v = 0.0001
		return(infparams)
		}

	if (grepl(pattern="4param_standard_fast", x=inffn))
		{
		infparams$j = parvec[3]
		
		# Calculate the actual weights
		ysv = 1-infparams$j
		v = ysv * parvec[4]
		ys = ysv - v

		infparams$v = v
		infparams$ys = ys
		infparams$maxent01 = 0.0001
		infparams$maxent01v = 0.0001
		return(infparams)
		}

	# unique(simhist_records2$run)
	if (grepl(pattern="5param_standard_fast_v", x=inffn))
		{
		infparams$j = parvec[3]
		
		# Calculate the actual weights
		ysv = 1-infparams$j
		v = ysv * parvec[4]
		ys = ysv - v

		infparams$v = v
		infparams$ys = ys
		infparams$maxent01 = 0.0001
		infparams$maxent01v = parvec[5]
		return(infparams)
		}

	if (grepl(pattern="5param_standard_fast", x=inffn))
		{
		infparams$j = parvec[3]
		
		# Calculate the actual weights
		ysv = 1-infparams$j
		v = ysv * parvec[4]
		ys = ysv - v

		infparams$v = v
		infparams$ys = ys
		infparams$maxent01 = parvec[5]
		infparams$maxent01v = 0.0001
		return(infparams)
		}

	if (grepl(pattern="6param_standard_fast", x=inffn))
		{
		infparams$j = parvec[3]
		
		# Calculate the actual weights
		ysv = 1-infparams$j
		v = ysv * parvec[4]
		ys = ysv - v

		infparams$v = v
		infparams$ys = ys
		infparams$maxent01 = parvec[5]
		infparams$maxent01v = parvec[6]
		return(infparams)
		}
	
	return(infparams)
	}






# Get the simulated model parameters
#######################################################
# get_simparams
#######################################################
#' Get the simulated model parameters from the row of a table
#' 
#' Basically this function assigns probability 1 to the simulated state/geographic range, and 
#' probability 0 for the other states/geographic ranges.  These data -- the simulated truth -- can then be compared to the inferred 
#' probabilities for the states, from e.g. \code{\link{get_ML_probs}}.
#' 
#' @param simhist_row A row from a table, which must have a column named \code{simulated_states_by_node_txt}.
#' @return \code{simparams} A list of the parameter values.
#' @export
#' @seealso \code{\link{simulate_biogeog_history}}, \code{\link{infprobs_to_probs_of_each_area}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#'
get_simparams <- function(simhist_row)
	{
	simparams = list()
	simparams$d = as.numeric(simhist_row$d)
	simparams$e = as.numeric(simhist_row$e)
	simparams$j = as.numeric(simhist_row$j)
	simparams$v = as.numeric(simhist_row$v)
	simparams$ys = as.numeric(simhist_row$ys)
	simparams$maxent01 = as.numeric(simhist_row$maxent01)
	simparams$maxent01v = as.numeric(simhist_row$maxent01v)
	
	return(simparams)
	}









#######################################################
# Process the ancestral state probs
#######################################################
#######################################################
# get_ML_states_from_relprobs
#######################################################
#' Extract the ML states at each node, from a table of relative probabilities -- old version
#' 
#' Given a table with the rows representing nodes, and the columns representing the relative
#' probabilities of each state, this function finds the ML (maximum likelihood) state(s)
#' for each node.
#' 
#' If possible, the input matrix should be the actual ML estimate of the state probabilities at each node,
#' rather than just the scaled conditional likelihoods at each node. The latter reflect only the 
#' tips-down information, whereas the former (the marginal ancestral state reconstruction) uses all of the information, and the probabilities of the states
#' at the root and in the outgroup(s) can influence the estimates in the ingroups.  This would not likely be particularly
#' important in a pure continuous-time model, but in a model with cladogenesis it could matter quite a bit.
#' 
#' See \url{http://blog.phytools.org/2013/03/marginal-ancestral-state-reconstruction.html} for more discussion of marginal
#' ancestral state reconstructions, versus mere scaled conditional likelihoods.
#' 
#' Revell and other sources (\cite{Felsenstein2004}) advocate the "re-rooting" method for obtaining the marginal ancestral state
#' reconstructions; however, re-rooting requires a time-reversible model and a tree with no root.  In biogeography we have
#' a \emph{non}-reversible model, and typically a time-scaled chronogram.  However, the same result can be obtained by
#' modifying the scaled conditional likelihoods obtained from a downpass from the tips, via an doing an up-pass from the root
#' scaled conditional likelihoods, being careful to transfer probabilities via the time-forward version of the 
#' Q-matrix and cladogenesis/speciation matrix.
#' 
#' \bold{Note:} further notes as this is implemented (required!)
#'
#' @param relprobs A numeric matrix of relative probabilities
#' @param statenames The names of the states/geographic ranges (e.g., A, AB, CDE, ABD, etc...)
#' @param returnwhat If "indices", return the 0-based indices of the states. If "states", return the name of the state, based on statenames.
#' @param if_ties What to do with ties. Currently, the only option is to take the first (this will be
#' shown in e.g. a pie chart, of course).
#' @return \code{ML_states} or \code{ML_states_indices}, depending on \code{returnwhat}.
#' @export
#' @seealso \code{\link{get_ML_state_indices}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://blog.phytools.org/2013/03/marginal-ancestral-state-reconstruction.html} 
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Felsenstein2004
#' @examples
#' testval=1
#' 
get_ML_states_from_relprobs <- function(relprobs, statenames, returnwhat="states", if_ties="takefirst")
	{
	# Get the maximum probability values for each row of the relative probs
	maxprobs = apply(relprobs, 1, max)
	
	# Indices, 1 thru number of states
	nums = 1:ncol(relprobs)
	
	ML_states = as.list(rep(NA,nrow(relprobs)))
	ML_states_indices = as.list(rep(NA,nrow(relprobs)))
	
	# get index (col #) of the ML state(s)
	# 
	for (i in 1:nrow(relprobs))
		{
		relprobs_row = relprobs[i,]
		maxprob = maxprobs[i]
		ML_states_indices[[i]] = get_ML_state_indices(relprobs_row, nums, maxprob, if_ties="takefirst")
		ML_states[[i]] = get_ML_state_indices(relprobs_row, nums, maxprob, if_ties="takefirst")
		}

	# return values
	if (returnwhat == "states")
		{
		foo = function(MLstate, statenames)
			{
			statenames[MLstate]
			}
		ML_states2 = unlist(sapply(X=ML_states, FUN=foo, statenames))
		ML_states2
		
		return(ML_states2)
		}
	if (returnwhat == "indices")
		{
		return(ML_states_indices)
		}
		
	return("get_ML_states error, no identified 'returnwhat'")
	}

#######################################################
# get_ML_state_indices
#######################################################
#' Extract the indices for the ML states at each node, given a row of relative probabilities
#' 
#' Given a table with the rows representing nodes, and the columns representing the relative
#' probabilities of each state, this function finds the ML (maximum likelihood) state(s)
#' for each node; \code{\link{get_ML_state_indices}} does this for a row, \code{\link{get_ML_states}}
#' iterates over all the rows.
#'
#' @param relprobs_row A row from a \code{relprobs}, a numeric matrix of relative probabilities
#' @param nums Numbers indexing the states from 1 to numstates
#' @param maxprob The value of the maximum probability for the row.
#' @param if_ties What to do with ties. Currently, the only option is to take the first (this will be
#' shown in e.g. a pie chart, of course).
#' @return \code{index_of_ML_state_s} 
#' @export
#' @seealso \code{\link{get_ML_states}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' 
get_ML_state_indices <- function(relprobs_row, nums, maxprob, if_ties="takefirst")
	{
	# get index (col #) of the ML state(s)
	# This could be more than one state, obviously, if several states have
	# equal probability
	index_of_ML_state_s = nums[relprobs_row == maxprob]
	
	if (length(index_of_ML_state_s) > 1)
		{
		index_of_ML_state_s = index_of_ML_state_s[1]
		}
	
	return(index_of_ML_state_s)
	}




