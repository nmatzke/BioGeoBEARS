#######################################################
# The v1 try is rather bogus; instead, read Bollback 
# 2007.  Key thing: use downpass likelihoods, then
# sample the root, then calc probs at the above nodes 
# conditional on this sample, multiple times likelihoods
# and sample, etc. This gives a VALID JOINT SAMPLE
#######################################################

# source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_extract_Qmat_COOmat_v1.R')
# source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_stochastic_mapping_v2.R')


#######################################################
# Stochastic mapping in a non-stratified analysis
#######################################################

# This function assumes that you start with the results of 
# an optimization run.  Alternatively, a user could 
# input their own desired parameter values (or values 
# sampled during a step of an MCMC search) and run
# calc_loglike_sp().
#
# We will sample "joint" histories, i.e. the state at each 
# node will be sampled, *CONDITIONAL* on the states already
# sampled at previously-sampled nodes.  These histories 
# will be valid draws from the distribution of histories
# specified by the model parameters and known tip states 
# (and any other constraints).
# 
# The well-known stochastic mapping papers by Nielsen (2002) 
# and Huelsenbeck et al. (2003) seem to use a different method:
# 
# 1. Calculate downpass conditional likelihoods at each node
#    (note that these are NOT full ancestral state probabilities,
#     which require an uppass and then multiplication; in 
#     stochastic mapping, the uppass calculation is replaced by
#     an uppass simulation).
#
# 2. Sample a state for the root node, from the downpass
#    likelihoods (which equal the ancestral state probabilities
#    at the root)
# 
# 3. Calculate the probability of a node descending from the 
#    root, now assuming the ancestral state is known.
#
# 4. Sample from this
#
# 5. Repeat for all nodes up the tree
#
# (e.g., Nielsen 2002, p. 731)
# 
# # Start by:
# # Assuming you've run this 
#   res = bears_optim_run(BioGeoBEARS_run_object)
#   res 
# 
# 
# stochastic_map <- function(res, rootedge=TRUE, statenum_bottom_root_branch_1based=NULL, printlevel=1, stratified=FALSE)
# 	{
# 
# 	# Cluster, if desired
# 	cluster_already_open = res$inputs$cluster_already_open
# 
# 	# Get the number of states
# 	numstates = ncol(res$ML_marginal_prob_each_state_at_branch_top_AT_node)
# 
# 	# Get the tree
# 	tr = read.tree(res$inputs$trfn)
# 	phy2 = reorder(tr, "pruningwise") # Do this, 
# 
# 	# Basic tree info
# 	ntips = length(phy2$tip.label)
# 	num_internal_nodes = phy2$Nnode
# 	tipnums = 1:ntips
# 	root_nodenum = ntips+1
# 	nodenums = root_nodenum:(ntips+num_internal_nodes)
# 	
# 	# Make a table holding the states etc.
# 	trtable = prt(phy2)
# 	trtable
# 	
# 	# Add a column for the sampled node states
# 	sampled_states_AT_nodes = rep(NA, nrow(trtable))
# 	sampled_states_AT_brbots = rep(NA, nrow(trtable))
# 	trtable = cbind(trtable, sampled_states_AT_nodes, sampled_states_AT_brbots)
# 	trtable
# 
# 	# Add the right and left descendant node numbers
# 	leftright_nodes_matrix = get_leftright_nodes_matrix_from_results(phy2)
# 	
# 	left_desc_nodes = rep(NA, nrow(trtable))
# 	right_desc_nodes = rep(NA, nrow(trtable))
# 	
# 	# dcorner = descendant corner (i.e. right after speciation)
# 	samp_LEFT_dcorner = rep(NA, nrow(trtable))
# 	samp_RIGHT_dcorner = rep(NA, nrow(trtable))
# 	
# 	trtable = cbind(trtable, left_desc_nodes, right_desc_nodes, samp_LEFT_dcorner, samp_RIGHT_dcorner)
# 	trtable$left_desc_nodes[nodenums] = leftright_nodes_matrix$left
# 	trtable$right_desc_nodes[nodenums] = leftright_nodes_matrix$right
# 	trtable[nodenums,]
# 
# 
# 		
# 	#returned_mats1 = get_Qmat_COOmat_from_BioGeoBEARS_run_object(default_BioGeoBEARS_run_object)
# 	#returned_mats1
# 	
# 	returned_mats2 = get_Qmat_COOmat_from_res(res)
# 	returned_mats2
# 	
# 	# Extract output
# 	Qmat = returned_mats2$Qmat
# 	COO_weights_columnar = returned_mats2$COO_weights_columnar
# 	Rsp_rowsums = returned_mats2$Rsp_rowsums
# 	
# 	# Calculate the likelihood P((left_state,right_state)|anc_state)
# 	# for each scenario (unconstrained)
# 	# Note:
# 	# COO_weights_columnar indices are 0-based, with no null_range
# 	# So, add 2 to get it so that e.g. state 0 = state 2 = Kauai
# 	#
# 	# Or, add 1 to get the 1based state indexes INSIDE COO_weights_columnar
# 	# 
# 	# COO_weights_columnar =
# 	# ancestral index, left index, right index, conditional
# 	# probability given ancestral states. (assuming likelihood
# 	# of descendants is 1)
# 	# Probabilities of each range-inheritance scenario, conditional
# 	# on ancestral state (without constraints on Left Branch state)
# 	like_LeftRight_given_AncState = COO_weights_columnar[[4]] / (Rsp_rowsums[1+COO_weights_columnar[[1]]])
# 	like_LeftRight_given_AncState
# 
# 	
# 	# Calculate the total number of range-inheritance scenarios
# 	# under the model
# 	# (this is the number of scenarios with weight > 0)
# 	# (weight per event/ sum(weights) = prob. per event)
# 	num_scenarios = length(COO_weights_columnar[[1]])
# 
# 
# 	##########################################################
# 	# 1. Sample a state at the root
# 	# 2. Given the root state, calculate uppass probabilities to the corners
# 	# 3. Multiply the corner uppass probs by the downpass probs
# 	# 
# 	##########################################################
# 	
# 	# Sample a state at the root
# 	root_stateprobs = res$ML_marginal_prob_each_state_at_branch_top_AT_node[root_nodenum,]
# 	statenums = 1:numstates
# 	statenum_1based = sample(x=statenums, size=1, replace=TRUE, prob=root_stateprobs)
# 	statenum_1based
# 	
# 	# Calculate the probability of each range inheritance scenario, 
# 	# given the chosen root state
# 	index_Qmat_0based_of_starting_state = statenum_1based - 1
# 	
# 	RCOO_weights_list_given_ancestor = given_a_starting_state_get_prob_of_each_split_scenario(index_Qmat_0based_of_starting_state, COO_weights_columnar, numstates=1+max(sapply(X=COO_weights_columnar, FUN=max)[1:3]), include_null_range=TRUE)
# 	
# 	uppass_probs_of_scenarios_given_root_state = RCOO_weights_list_given_ancestor
# 	uppass_probs_of_scenarios_given_root_state
# 
# 
# 	} # end stochastic_mapping()



given_a_starting_state_get_prob_of_each_split_scenario <- function(index_Qmat_0based_of_starting_state=1, COO_weights_columnar, numstates=numstates, include_null_range=TRUE)
	{
	defaults='
	# Note that the speciation matrix is always missing the original state 0 (null range);
	# Thus the 0th state is actually original state 1 in the Qmat (starting with 0)
	index_Qmat_0based_of_starting_state = 1
	
	for (i in 1:1000)
		{
		given_a_starting_state_simulate_split(index_Qmat_0based_of_starting_state=1, COO_probs_columnar, numstates=16)
		}
	' # END defaults
	
	# Necessary setting to avoid getting numbers etc. in the stochastic mapping output
	options(stringsAsFactors=FALSE)

	if (include_null_range == TRUE)
		{
		numstates_during_cladogenesis = numstates - 1
		} else {
		numstates_during_cladogenesis = numstates - 0
		}

	
	# Error check; include_null_range=TRUE is default, even for old models
	if ( (length(include_null_range) == 0) || (is.na(include_null_range)) || (is.null(include_null_range))  ) 
		{
		include_null_range = TRUE
		}

	
	if (include_null_range == TRUE)
		{
		index_shift_from_Qmat_to_COOmat = -1
		} else {
		index_shift_from_Qmat_to_COOmat = 0		
		}
	
	
	# If there are 16 states, there are 15 non-null rowSums
# 	COO_weights_columnar_1based = COO_weights_columnar
# 	COO_weights_columnar_1based[[1]] = 1 + COO_weights_columnar_1based[[1]]
# 	COO_weights_columnar_1based[[2]] = 1 + COO_weights_columnar_1based[[2]]
# 	COO_weights_columnar_1based[[3]] = 1 + COO_weights_columnar_1based[[3]]
# 	Rsp_rowsums = rcpp_calc_rowsums_for_COOweights_columnar(COO_weights_columnar_1based, numstates=numstates-1)
# 	Rsp_rowsums
# 	# Get the conditional probs_of_each_scenario
# 	condprobs_each_split_scenario = COO_weights_columnar_1based[[4]]
# 	for (ii in 1:length(COO_weights_columnar_1based[[1]]))
# 		{
# 		condprobs_each_split_scenario[ii] = condprobs_each_split_scenario[ii] / 
# 	
		
	# Original
	
	Rsp_rowsums = rcpp_calc_rowsums_for_COOweights_columnar(COO_weights_columnar, numstates=numstates_during_cladogenesis)
	Rsp_rowsums

	# Load up the speciation matrix
	# These indexes are 0-14, i.e. 1-15 but 0-based
	RCOO_weights_columnar_anc_i_list = COO_weights_columnar[[1]]
	RCOO_left_i_list = COO_weights_columnar[[2]]
	RCOO_right_j_list = COO_weights_columnar[[3]]
	RCOO_weights_list = COO_weights_columnar[[4]]

	RCOO_condprobs_list = rep(0, length(RCOO_weights_list))
	for (i in 1:length(Rsp_rowsums))
		{
		# Need to do something here?
		}

	# Number of nonzero cells (i.e. in the speciation matrix
	num_nonzero_cells_in_sp_matrix = length(RCOO_weights_columnar_anc_i_list)
	
	# Get the ancestors that match the input ancestor
	# COO_probs_columnar have 0-based indexes, 
	# and NEVER have null range
	anc_match_TF = RCOO_weights_columnar_anc_i_list == (index_Qmat_0based_of_starting_state+index_shift_from_Qmat_to_COOmat)
	num_nonzero_descendent_splits = sum(anc_match_TF)
	num_nonzero_descendent_splits
	
	# Get a list of the relative probabilities of each inheritance scenario,
	# by setting the non-found ancestor probabilities to 0
	RCOO_relProbs_list_given_ancestor = RCOO_weights_list
	RCOO_relProbs_list_given_ancestor[anc_match_TF == FALSE] = 0
	
	# Since we only have 1 ancestral state, we can get the
	# actual conditional probabilities of each scenario by
	# dividing by the rowsums; but the rowsums for each 
	# ancestor state have already been calculated for us
	
	# To refer to the correct index in Rsp_rowsums, add 1 to the 0-based index
	# of 15 non-null possible ancestral states
	#RCOO_relProbs_list_given_ancestor[anc_match_TF] = RCOO_relProbs_list_given_ancestor[anc_match_TF] / Rsp_rowsums[index_Qmat_0based_of_starting_state+index_shift_from_Qmat_to_COOmat+1]
	
	#RCOO_weights_list_given_ancestor = RCOO_relProbs_list_given_ancestor
	RCOO_weights_list_given_ancestor = RCOO_relProbs_list_given_ancestor / sum(RCOO_relProbs_list_given_ancestor)
	
	# This is the probabilities for all scenarios (will include many 0s)
	return(RCOO_weights_list_given_ancestor)
	}





# Go through each inheritance scenario, multiply the uppass conditional probs
# (conditional on the sampled node state) by the downpass probs of the 
# assumed descendants
# Note: just a product function, mapply-ed

sample_split_scenario2 <- function(COO_weights_columnar, probs_ancstate, left_branch_downpass_likes, right_branch_downpass_likes, sample_which="both", return_prob_each_split_scenario=FALSE, include_null_range=TRUE, Rsp_rowsums=NULL, numstates_wo_null=NULL, printflag=FALSE)
	{
	defaults='
	probs_ancstate = c(0.000, 0.907, 0.000, 0.093)
	return_prob_each_split_scenario=TRUE
	include_null_range=TRUE
	'
	
	# Get the rowsums of the cladogenesis matrix, for calculating
	# conditional probabilities of each scenario
	if (is.null(Rsp_rowsums))
		{
		if (is.null(numstates_wo_null))
			{
			numstates_during_cladogenesis = 1 + max(sapply(X = COO_weights_columnar, 
    FUN = max)[1:3])
			} # END if (is.null(numstates_during_cladogenesis))
		# Calculate the Rsp_rowsums
		Rsp_rowsums = rcpp_calc_rowsums_for_COOweights_columnar(COO_weights_columnar=COO_weights_columnar, numstates=numstates_during_cladogenesis)
		} # END if (is.null(Rsp_rowsums))

	
	# Error check; include_null_range=TRUE is default, even for old models
	if ( (length(include_null_range) == 0) || (is.na(include_null_range)) || (is.null(include_null_range))  ) 
		{
		include_null_range = TRUE
		}
	
	# If include_null_range==TRUE, to translate from R_1based states
	# 1-16 in the Qmat & downpass likes requires +2 in the 
	# 0-based COOmat
	
	if (include_null_range == TRUE)
		{
		COOmat_0based_to_Qmat_1based = 2
		} else {
		COOmat_0based_to_Qmat_1based = 1
		}
	

	# Calculate the probability of each scenario, given 
	# ancprobs, 
	# left downpass likes
	# right downpass likes
	# the COO_weights_columnar cladogenesis model
	
	# Make the columns first before cbinding
	ancprobs = probs_ancstate[ COO_weights_columnar[[1]] + COOmat_0based_to_Qmat_1based]
	Lprobs = left_branch_downpass_likes[ COO_weights_columnar[[2]] + COOmat_0based_to_Qmat_1based]
	Rprobs = right_branch_downpass_likes[ COO_weights_columnar[[3]] + COOmat_0based_to_Qmat_1based]
	# Weights divided by the sum of the weights for that row
	scenario_condprob = COO_weights_columnar[[4]] / Rsp_rowsums[ COO_weights_columnar[[1]] + 1 ]
	
	
	# 
	input_probs_each_split_scenario = cbind(ancprobs, Lprobs, Rprobs, scenario_condprob)
	input_probs_each_split_scenario

	relprob_each_split_scenario = apply(X=input_probs_each_split_scenario, MARGIN=1, FUN=prod)
	relprob_each_split_scenario
	
	prob_each_split_scenario = relprob_each_split_scenario / sum(relprob_each_split_scenario)
	
	if (printflag)
		{
		print(round(prob_each_split_scenario, 3))
		}
	
	# R_combine_uppass_splitprobs_w_downpass_condlikes
	# Make a table
	# Everything, including anc, left, right state indices
	#cbind(COO_weights_columnar[[1]], COO_weights_columnar[[2]], COO_weights_columnar[[3]], uppass_probs_of_scenarios_given_root_state, left_branch_downpass_likes[1+COO_weights_columnar[[2]]], right_branch_downpass_likes[1+COO_weights_columnar[[3]]])
	# Just the probs
	
	# ***** double-check 2 and 3 here...NOPE, THEY ARE CORRECT
	# ***** 2014-12-29_NJM
	
	# Number the scenarios
	split_scenario_nums = 1:length(COO_weights_columnar[[1]])
	# Sample one of them
	split_scenario_num = sample(x=split_scenario_nums, size=1, replace=TRUE, prob=prob_each_split_scenario)
	
	# Get the left and right descendant states
	# State 0 in COO_weights_columnar[[2]] is state 2 (K) in normal...
	left_decstate_1based = COOmat_0based_to_Qmat_1based + COO_weights_columnar[[2]][split_scenario_num]
	right_decstate_1based = COOmat_0based_to_Qmat_1based + COO_weights_columnar[[3]][split_scenario_num]
		
	sampled_split_descendants = NULL
	sampled_split_descendants$left_decstate_1based = left_decstate_1based
	sampled_split_descendants$right_decstate_1based = right_decstate_1based
	
	if (return_prob_each_split_scenario == TRUE)
		{
		sampled_split_descendants$prob_each_split_scenario = prob_each_split_scenario
		}
	
	return(sampled_split_descendants)
	}



sample_split_scenario <- function(COO_weights_columnar, uppass_probs_of_scenarios_given_root_state, left_branch_downpass_likes, right_branch_downpass_likes, sample_which="both", return_prob_each_split_scenario=FALSE, include_null_range=TRUE)
	{
	defaults='
	return_prob_each_split_scenario=FALSE
	include_null_range=TRUE
	'
	
	# Error check; include_null_range=TRUE is default, even for old models
	if ( (length(include_null_range) == 0) || (is.na(include_null_range)) || (is.null(include_null_range))  ) 
		{
		include_null_range = TRUE
		}
	
	# If include_null_range==TRUE, to translate from R_1based states
	# 1-16 in the Qmat & downpass likes requires +2 in the 
	# 0-based COOmat
	
	if (include_null_range == TRUE)
		{
		COOmat_0based_to_Qmat_1based = 2
		} else {
		COOmat_0based_to_Qmat_1based = 1
		}
	
	
	# R_combine_uppass_splitprobs_w_downpass_condlikes
	# Make a table
	# Everything, including anc, left, right state indices
	#cbind(COO_weights_columnar[[1]], COO_weights_columnar[[2]], COO_weights_columnar[[3]], uppass_probs_of_scenarios_given_root_state, left_branch_downpass_likes[1+COO_weights_columnar[[2]]], right_branch_downpass_likes[1+COO_weights_columnar[[3]]])
	# Just the probs
	
	# ***** double-check 2 and 3 here...NOPE, THEY ARE CORRECT
	# ***** 2014-12-29_NJM
	
	# Check joint versus individual sampling
	if (sample_which == "both")
		{
		# COO_weights_columnar range from 0-14
		# ...which corresponds to anagenetic 1-based states 2-16
		input_probs_each_split_scenario = cbind(uppass_probs_of_scenarios_given_root_state, left_branch_downpass_likes[ COOmat_0based_to_Qmat_1based + COO_weights_columnar[[2]] ], right_branch_downpass_likes[ COOmat_0based_to_Qmat_1based + COO_weights_columnar[[3]] ])
		input_probs_each_split_scenario
		relprob_each_split_scenario = apply(X=input_probs_each_split_scenario, MARGIN=1, FUN=prod)
		relprob_each_split_scenario
		prob_each_split_scenario = relprob_each_split_scenario / sum(relprob_each_split_scenario)
		prob_each_split_scenario
	
		# Number the scenarios
		split_scenario_nums = 1:length(COO_weights_columnar[[1]])
		# Sample one of them
		split_scenario_num = sample(x=split_scenario_nums, size=1, replace=TRUE, prob=prob_each_split_scenario)
	
		# Get the left and right descendant states
		# State 0 in COO_weights_columnar[[2]] is state 2 (K) in normal...
		left_decstate_1based = COOmat_0based_to_Qmat_1based + COO_weights_columnar[[2]][split_scenario_num]
		right_decstate_1based = COOmat_0based_to_Qmat_1based + COO_weights_columnar[[3]][split_scenario_num]
		} # END if (sample_which == "both")

	if (sample_which == "left")
		{
		# COO_weights_columnar range from 0-14
		# ...which corresponds to anagenetic 1-based states 2-16
		input_probs_each_split_scenario = cbind(uppass_probs_of_scenarios_given_root_state, left_branch_downpass_likes[ COOmat_0based_to_Qmat_1based + COO_weights_columnar[[2]] ], right_branch_downpass_likes[ COOmat_0based_to_Qmat_1based + COO_weights_columnar[[3]] ])
		input_probs_each_split_scenario
		relprob_each_split_scenario = apply(X=input_probs_each_split_scenario, MARGIN=1, FUN=prod)
		relprob_each_split_scenario
		prob_each_split_scenario = relprob_each_split_scenario / sum(relprob_each_split_scenario)
		prob_each_split_scenario
	
		# Number the scenarios
		split_scenario_nums = 1:length(COO_weights_columnar[[1]])
		# Sample one of them
		split_scenario_num = sample(x=split_scenario_nums, size=1, replace=TRUE, prob=prob_each_split_scenario)
	
		# Get the left and right descendant states
		# State 0 in COO_weights_columnar[[2]] is state 2 (K) in normal...
		left_decstate_1based = COOmat_0based_to_Qmat_1based + COO_weights_columnar[[2]][split_scenario_num]
		right_decstate_1based = COOmat_0based_to_Qmat_1based + COO_weights_columnar[[3]][split_scenario_num]
		} # END if (sample_which == "both")

		
	sampled_split_descendants = NULL
	sampled_split_descendants$left_decstate_1based = left_decstate_1based
	sampled_split_descendants$right_decstate_1based = right_decstate_1based
	
	if (return_prob_each_split_scenario == TRUE)
		{
		sampled_split_descendants$prob_each_split_scenario = prob_each_split_scenario
		}
	
	return(sampled_split_descendants)
	} # END sample_split_scenario <- function(COO_weights_columnar, uppass_probs_of_scenarios_given_root_state, left_branch_downpass_likes, right_branch_downpass_likes, sample_which="both", return_prob_each_split_scenario=FALSE, include_null_range=TRUE)





stochastic_map_branch <- function(nodenum_at_top_of_branch, trtable, Qmat, state_indices_0based, ranges_list, areas, single_branch=FALSE, stratified=FALSE, maxtries=40000, manual_history_for_difficult_branches=TRUE)
	{
	#2016-05-07
	bugfix_on_Indicatoraceae='
	seedval = 27579383
	timeperiod_i = 1
	piecenum = 10
	print(seedval)
	print(timeperiod_i)
	print(piecenum)

	nodenum_at_top_of_branch=subtable_rownum; 
	trtable=master_table_timeperiod_i; 
	Qmat; 
	state_indices_0based; 
	ranges_list; 
	areas; 
	single_branch=single_branch; 
	stratified=stratified; 
	maxtries=maxtries
	manual_history_for_difficult_branches=TRUE
	'


	#######################################################
	# Simulate events on branches
	#######################################################
	defaults='
	nodenum_at_top_of_branch = 31
	'
	
	states_list = state_indices_0based
	numstates = length(state_indices_0based)
	
	# Setup
	statenum_1based_at_branch_bottom = trtable$sampled_states_AT_brbots[nodenum_at_top_of_branch]
	statenum_1based_at_branch_top = trtable$sampled_states_AT_nodes[nodenum_at_top_of_branch]
	
	#print("Trying to get from:")
	#print(statenum_1based_at_branch_bottom)
	#print(statenum_1based_at_branch_top)
	
	names_in_trtable = names(trtable)
	
	# Stratified analyses have SUBedge.length columns, however these may not have
	# been updated accurately during tree chainsawing. The below independently
	# calculates the branch length, which should not be larger than the 
	# time bin width, which is $reltimept
	if ( ("SUBedge.length" %in% names_in_trtable) == FALSE)
		{
		brlen = trtable$edge.length[nodenum_at_top_of_branch]
		if (brlen <= 0)
			{
			stop("STOP ERROR00 in stochastic_map_branch(): brlen <= 0")
			} # END if (brlen_in_section <= 0)

		} else {
		
		# 2014-05-25_NJM: check if reltimept smaller?
		if (is.na(trtable$SUBedge.length[nodenum_at_top_of_branch]))
			{
			# If it's a sub-branch, set edges of length NA to 0
			trtable$SUBedge.length[nodenum_at_top_of_branch] = 0
			}
		
		# Check if sub-edge length and edge length are the same:
		subedge_length_equals_edge_length_WORRY = FALSE
		if (trtable$SUBedge.length[nodenum_at_top_of_branch] > trtable$reltimept[nodenum_at_top_of_branch])
			{
			subedge_length_equals_edge_length_WORRY = TRUE

			# We can fix sub-branches easily:
			if (trtable$piececlass[nodenum_at_top_of_branch] == "subbranch")
				{
				# Fix to reltimept IF it's *NOT* a fossil:
				if (  (is.na(trtable$fossils[nodenum_at_top_of_branch]) == TRUE) || (trtable$fossils[nodenum_at_top_of_branch] == FALSE) )
					{
					brlen_in_section = trtable$reltimept[nodenum_at_top_of_branch]
					subedge_length_equals_edge_length_WORRY = FALSE
					if (brlen_in_section <= 0)
						{
						stop("STOP ERROR0a in stochastic_map_branch(): brlen_in_section <= 0")
						} # END if (brlen_in_section <= 0)
					} else {
					# It *IS* a fossil, it's brlen is based on time_bp
					brlen_in_section = trtable$time_bot[nodenum_at_top_of_branch] - trtable$time_bp[nodenum_at_top_of_branch]
					subedge_length_equals_edge_length_WORRY = FALSE
					if (brlen_in_section <= 0)
						{
						print(trtable[nodenum_at_top_of_branch,])
						print(brlen_in_section)
						stop("STOP ERROR0b in stochastic_map_branch(): brlen_in_section <= 0")
						} # END if (brlen_in_section <= 0)
					}
				} # END check for fossils
			} else {
			brlen_in_section = trtable$SUBedge.length[nodenum_at_top_of_branch]
			if (brlen_in_section <= 0)
				{
				stop("STOP ERROR0c in stochastic_map_branch(): brlen_in_section <= 0")
				} # END if (brlen_in_section <= 0)
			}
		if (brlen_in_section <= 0)
			{
			stop("STOP ERROR1 in stochastic_map_branch(): brlen_in_section <= 0")
			} # END if (brlen_in_section <= 0)
		
		
		# Alternatively, if it's a root branch of a subtree, the branch length is just
		# the node time_bp minus the timeperiod bottom
		if ((trtable$SUBnode.type[nodenum_at_top_of_branch] == "root") && (trtable$piececlass[nodenum_at_top_of_branch] == "subtree"))
			{
			brlen_in_section = trtable$time_bot[nodenum_at_top_of_branch] - trtable$time_bp[nodenum_at_top_of_branch]
			subedge_length_equals_edge_length_WORRY = FALSE
			if (brlen_in_section <= 0)
				{
				stop("STOP ERROR2 in stochastic_map_branch(): brlen_in_section <= 0")
				} # END if (brlen_in_section <= 0)
			}
		
		
		if (subedge_length_equals_edge_length_WORRY == TRUE)
			{
			errortxt = paste("\n\nError in stochastic_map_branch(): your master_tree table, at row 'nodenum_at_top_of_branch'=", nodenum_at_top_of_branch, "\nhas an SUBedge.length > reltimept and was not corrected, as it's not a subbranch.\n\n", sep="")
			cat(errortxt)

			cat("\n\n")
			print("nodenum_at_top_of_branch:")
			cat("\n\n")
			print(nodenum_at_top_of_branch)
			cat("\n\n")
			print("trtable[nodenum_at_top_of_branch,]:")
			cat("\n\n")
			print(trtable[nodenum_at_top_of_branch,])
			cat("\n\n")
			}
	
		# 2014-05-25_NJM new:
		brlen = brlen_in_section
		if (brlen_in_section <= 0)
			{
			stop("STOP ERROR3 in stochastic_map_branch(): brlen_in_section <= 0")
			} # END if (brlen_in_section <= 0)
		
		# Original:
		#brlen = trtable$SUBedge.length[nodenum_at_top_of_branch]
		}
	#edgenum = trtable$parent_br[nodenum_at_top_of_branch]
	
	# Get the rates from the Qmat
	# Source: https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/R/simulation.R?root=proteinevoutk20
	# Rate of an event, conditional on each possible state
	rates_de = -diag(Qmat) #exponential rate of waiting, diagonal of Qmat



	# Simulate up branch
	# Number of tries
	trynum = 0
	#maxtries = 40000


	# Run this while loop until you hit a successful simulation
	while (1)
		{
		# current time starts at 0 (bottom of branch)
		# Actually, start it at the time of the branch bottom
		# And make it time before present (time_bp)
		if (stratified == TRUE)
			{
			# Stratified analysis:
			# Take the time at the top of the time bin, add the time to the top of the node, then add the brlen before
			abs_time_at_branch_bottom = trtable$time_top[nodenum_at_top_of_branch] + trtable$SUBtime_bp[nodenum_at_top_of_branch] + brlen
			} else {
			# Non-stratified analysis: 
			# Take the time at the top of the node, then add the brlen before
			abs_time_at_branch_bottom = trtable$time_bp[nodenum_at_top_of_branch] + brlen
			}
		curr_time = 0
		time_stop = brlen	# (subtracting since we are going up towards 0)
		current_rangenum_1based = statenum_1based_at_branch_bottom

		trynum = trynum + 1
		# Keep track of events
		event_time = NULL
		event_type = NULL
		event_txt = NULL
		dispersal_to = NULL
		extirpation_from = NULL
		events_table_for_branch = NULL
	
		# Repeat until one of the break statements is reached
		repeatnum = 1
		maxrepeats = 10000
		repeat
			{
		
			# Exponentially-distributed waiting time
			# based on current rate
			rate = rates_de[current_rangenum_1based]
			
			# If you hit the null state, you'll get a rate of 0,
			# scotch that simulation
			if (rate == 0)
				{
				break
				}
			
			waiting_time = rexp(n=1, rate=rate)
			curr_time = curr_time + waiting_time

			# Check if you've passed the stop of the branch
			time_stop_hit = FALSE

			#print(dt)
			#print(edge.zone[alive_TF]==1)
			#print(edge.zone[alive_TF]==2)
			# Counting time down to 0
			#cat("\n", curr_time, " <= ", time_stop, "	", (curr_time <= time_stop), sep="")
			#cat("\n", waiting_time, " < ", brlen, "	", (waiting_time < brlen), sep="")
			
			
			#if (curr_time > time_stop)
			#if (waiting_time  brlen)
			#	{
			#	stop(cat("\n\n\nREPEATNUM: ", repeatnum, "\n\n\n"))
			#	}
			#print(trynum)
			#print(brlen)
			if (curr_time > brlen)
				{
				curr_time = brlen
				time_stop_hit = TRUE
				break
				}
			#print("B")
			# Don't allow more than 100000 events on a branch (!)
			if (repeatnum >= maxrepeats)
				{
				break
				} else {
				repeatnum = repeatnum + 1
				}
			#print("C")
		
			# Anagenetic range expansion/contraction event
			# (can result in extinction if ranges of size 1 drop to 0)
			# (such simulated branches would fail, however)
			probs_of_new_ranges = Qmat[current_rangenum_1based, ]
			# Zero out the diagonal, which is negative (and you 
			# know something changed, since you drew this)
			probs_of_new_ranges[probs_of_new_ranges < 0] = 0
			probs_of_new_ranges = probs_of_new_ranges / sum(probs_of_new_ranges)

			#cat("\ncurrent_rangenum_1based:", sep="")
			#cat("\n", current_rangenum_1based, sep="")

			#cat("\nprobs_of_new_ranges:", sep="")
			#cat("\n", probs_of_new_ranges, sep="")

			new_rangenum_1based = sample(x=1:numstates, size=1, replace=FALSE, prob=probs_of_new_ranges)
			#print(new_rangenum_1based)
			
			if (isblank_TF(new_rangenum_1based) == TRUE)
				{
				txt = paste0("STOP ERROR_line741 in stochastic_map_branch(): new_rangenum_1based is BLANK. new_rangenum_1based='", new_rangenum_1based, "'")
				cat("\n\n")
				cat(txt)
				cat("\n\n")
				
				print("probs_of_new_ranges:")
				print(probs_of_new_ranges)
				
				stop(txt)
				}
			
			
			
			#cat("\nnew_rangenum_1based:", sep="")
			#cat("\n", new_rangenum_1based, sep="")
			
			# Save data on event
			# The absolute event time is the age of the bottom of 
			# the branch, minus the current time
			abs_event_time = abs_time_at_branch_bottom - curr_time
			event_time = curr_time
		
			##################
			# Get event type
			##################
			# "dispersal" (range expansion)
			if ( length(states_list[[new_rangenum_1based]]) > length(states_list[[current_rangenum_1based]]) )
				{
				event_type = "d"
			
				# Identify the new area
				TF = states_list[[new_rangenum_1based]] %in% states_list[[current_rangenum_1based]]
				new_area_TF = TF == FALSE
				new_area_num_0based = states_list[[new_rangenum_1based]][new_area_TF]
				new_area_num_1based = new_area_num_0based + 1
				lost_area_num_1based = "-"
				dispersal_to = areas[new_area_num_1based]
				extirpation_from = "-"
				}

			# "extinction" (range contraction / local extirpation)
			if ( length(states_list[[new_rangenum_1based]]) < length(states_list[[current_rangenum_1based]]) )
				{
				event_type = "e"

				# Identify the area lost
				TF = states_list[[current_rangenum_1based]] %in% states_list[[new_rangenum_1based]]
				lost_area_TF = TF == FALSE
				lost_area_num_0based = states_list[[current_rangenum_1based]][lost_area_TF]
				lost_area_num_1based = lost_area_num_0based + 1
				new_area_num_1based = "-"
				dispersal_to = "-"
				extirpation_from = areas[lost_area_num_1based]
				}
		
			# If the new and old states are the same size, check further
			if ( length(states_list[[new_rangenum_1based]]) == length(states_list[[current_rangenum_1based]]) )
				{
				# Contraction will result in the same length, check for this
				if (ranges_list[[new_rangenum_1based]] == "_")
					{
					event_type = "e"

					# Identify the area lost
					TF = states_list[[current_rangenum_1based]] %in% states_list[[new_rangenum_1based]]
					lost_area_TF = TF == FALSE
					lost_area_num_0based = states_list[[current_rangenum_1based]][lost_area_TF]
					lost_area_num_1based = lost_area_num_0based + 1
					new_area_num_1based = "-"
					dispersal_to = "-"
					extirpation_from = areas[lost_area_num_1based]
					} else {
					event_type = "a" # anagenetic range-switching
				
					# In a sense, we have joint dispersal and extinction
					#print(ranges_list)
					#print(states_list[[new_rangenum_1based]])
					
					new_area_num_1based = 1+states_list[[new_rangenum_1based]]
					dispersal_to = areas[new_area_num_1based]

					lost_area_num_1based = 1+states_list[[current_rangenum_1based]]
					extirpation_from = areas[lost_area_num_1based]
					}
				} # END if ( length(states_list[[new_rangenum_1based]]) == 
				  #          length(states_list[[current_rangenum_1based]]) )
			
			##################
			# ENDING Get event type
			##################
		
			# event_txt
			event_txt = paste(ranges_list[[current_rangenum_1based]], "->", ranges_list[[new_rangenum_1based]], sep="")
		
			# Make a row for a data.table
			current_rangetxt = ranges_list[[current_rangenum_1based]]
			new_rangetxt = ranges_list[[new_rangenum_1based]]
			tmprow = c(nodenum_at_top_of_branch, trynum, brlen, current_rangenum_1based, new_rangenum_1based, current_rangetxt, new_rangetxt, abs_event_time, event_time, event_type, event_txt, new_area_num_1based, lost_area_num_1based, dispersal_to, extirpation_from)
		
			events_table_for_branch = rbind(events_table_for_branch, tmprow)
		
		
			# Update the current branch statenum
			current_rangenum_1based = new_rangenum_1based
			} # END repeat
	
		# Look for successful hit; if not, repeat
		# If, at the end of the simulation, current_rangenum_1based equals
		# the known endpoint, then break out of the loop
		#print(c(current_rangenum_1based, statenum_1based_at_branch_top, current_rangenum_1based == statenum_1based_at_branch_top))
		# Check if you are in the right state, at the end of the branch simulation
		if (current_rangenum_1based == statenum_1based_at_branch_top)
			{
			# Successful simulation of branch!
			#print(paste0("branch below node ", nodenum_at_top_of_branch, " successful!"))
			error_check_Psychotria_all_tips_size1 = FALSE
			if (error_check_Psychotria_all_tips_size1)
				{
				if (trtable$time_top[nodenum_at_top_of_branch] == 0)
					{
					cat("\n\nSuccessful simulation of branch...but check if the simulated state matches the sampled state! \n\n")
			
					print(statenum_1based_at_branch_top)
			
					print(current_rangenum_1based)
			
					print(nodenum_at_top_of_branch)
			
					print(events_table_for_branch)
					} # END if (trtable$time_top[nodenum_at_top_of_branch] == 0)
				} # END if (error_check_Psychotria_all_tips_size1)
			
			# Break out of loop, simulation was successful			
			break()
			} # END repeat
		
		
		
		# print(c(trynum, maxtries))
		if (trynum > maxtries)
			{
			node_in_master_tree_w_error = trtable$node[nodenum_at_top_of_branch]
			errortxt = paste("\n\nError_in_stochastic_simulation; no success on branch below node ", node_in_master_tree_w_error, " after ", maxtries, " tries. Printing starting/ending states, and events table for last attempt.", sep="")
			cat(errortxt)
			
			cat("\n\nStarting state (1-based): ", current_rangenum_1based, sep="")
			cat("\n\nEnding state (1-based): ", new_rangenum_1based, sep="")
			cat("\nGoal state (1-based): ", statenum_1based_at_branch_top, sep="")
			
			cat("\n\n")
			print("Last attempt at events_table_for_branch:")
			cat("\n\n")
			print(events_table_for_branch)
			cat("\n\n")
			
			#events_table_for_branch = rep(
			#stop("Stopping...")
			
			if (manual_history_for_difficult_branches == FALSE)
				{
				cat("\nReturning object of class 'try-error'.\n")
			
				error_msg_for_stop = paste("Error in stochastic_map_branch(): Stochastic mapping failed after maxtries=", maxtries, " tries on branch below node=", nodenum_at_top_of_branch, sep="")
			
				attr(error_msg_for_stop, which="class") = "try-error"
				attr(error_msg_for_stop, which="condition") = error_msg_for_stop
				stop(error_msg_for_stop)
				} # END if (manual_history_for_difficult_branches == FALSE)

			if (manual_history_for_difficult_branches == TRUE)
				{
				error_msg_for_stop = paste("Error in stochastic_map_branch(): Stochastic mapping failed after maxtries=", maxtries, " tries on branch below node=", nodenum_at_top_of_branch, sep="")
				cat("\n")
				cat(error_msg_for_stop)
				
				cat("\nAs manual_history_for_difficult_branches==TRUE, we are manually devising a history to force-fit a difficult branch (e.g. AB -> ACEDFGH)... \n")

				cat("\nNOTE: This fix will run, but will NOT work correctly for pure '+a' models (as of 2016-05-07) where d and/or e equals 0, since the manual fix forces a random series of d and e events, on the theory that the events must be super-rare since maxtries failed, so we might as well order them randomly (subject to the constraint that the range is never NULL). For +a models, try increasing maxtries instead. \n")

				# Manual histories are denoted by maxtries+2
				trynum = maxtries + 2
				
				# Look at the starting and ending ranges:

				range_at_branch_bottom = ranges_list[[statenum_1based_at_branch_bottom]]
				range_at_branch_top = ranges_list[[statenum_1based_at_branch_top]]
				
				cat("Starting range: ", range_at_branch_bottom, "\n")
				cat("Desired ending range: ", range_at_branch_top, "\n")
				
				# Force a series of range/loss events at uniform random intervals
				starting_ranges = strsplit(range_at_branch_bottom, split="")[[1]]
				ending_ranges = strsplit(range_at_branch_top, split="")[[1]]
				
				area_indices_0based_branch_bottom = state_indices_0based[[statenum_1based_at_branch_bottom]]
				area_indices_0based_branch_top = state_indices_0based[[statenum_1based_at_branch_top]]

				
				
				# Ranges to add:
				addTF = (area_indices_0based_branch_top %in% area_indices_0based_branch_bottom) == FALSE
				area_indices_0based_to_add = area_indices_0based_branch_top[addTF]
				# Ranges to subtract:
				subTF = (area_indices_0based_branch_bottom %in% area_indices_0based_branch_top) == FALSE
				area_indices_0based_to_subtract = area_indices_0based_branch_bottom[subTF]
				
				# Make a table of add/loss events
				add_cols = cbind(area_indices_0based_to_add, rep("add", times=length(area_indices_0based_to_add)))
				sub_cols = cbind(area_indices_0based_to_subtract, rep("sub", times=length(area_indices_0based_to_subtract)))
				manual_table = rbind(add_cols, sub_cols)
				manual_table = as.data.frame(manual_table, stringsAsFactors=FALSE)
				names(manual_table) = c("area_0index", "addsub")
				manual_table$area_0index = as.numeric(manual_table$area_0index)
				
				
				# Randomly sort table, make sure you never pass through a null range
				# Repeat until one of the break statements is reached
				while (1)
					{
					rownums = 1:nrow(manual_table)
					rownums = sample(x=rownums, size=length(rownums), replace=FALSE)
					sorted_manual_table = manual_table[rownums,]
					
					tmp_area_indices = area_indices_0based_branch_bottom
					
					break_condition = TRUE
					for (mm in 1:(nrow(sorted_manual_table)) )
						{
						if (sorted_manual_table[mm,2] == "add")
							{
							tmp_area_indices = sort(c(tmp_area_indices, sorted_manual_table[mm,1]))
							} else {
							tmp_area_indices = tmp_area_indices[tmp_area_indices != sorted_manual_table[mm,1]]
							}
						#print(tmp_area_indices)	
							
						if (length(tmp_area_indices) == 0)
							{
							break_condition = FALSE
							# Double-check (null range OK at the end of the events, if end state is null)
							if ( is.na(area_indices_0based_branch_top) && (mm == nrow(sorted_manual_table)) )
								{
								break_condition = TRUE
								}
							}
						
						} # END for (mm in 1:(nrow(sorted_manual_table)-1) )
					
					if (break_condition == TRUE)
						{
						break()
						}
					} # END while (1)
				
				# Now check if the events go through a null state (except the end state, which 
				# could be null)
				curr_time = 0
				time_stop = brlen	# (subtracting since we are going up towards 0)

				# Event times: uniformly distributed under this manual scenario
				# 2016-05-31_bug: sometimes this didn't produce enough 
				#                 event times, producing NAs
				# (I don't know why I put length() in the first place, nrow() is
				#  what makes sense and produces the correct number of events)
				# (Probably this means stochastic maps that made frequent use of
				#  the manual-histories backups should be re-run, because if the 
				#  list of times was too long, only the first chunk would be used, 
				#  so they might not follow the actual uniform distribution on the
				#  manual-history branch.)
				#event_times = sort(runif(n=length(sorted_manual_table), min=curr_time, max=time_stop))
				# 2016-05-31_fix:
				event_times = sort(runif(n=nrow(sorted_manual_table), min=curr_time, max=time_stop))
				
				current_rangetxt = range_at_branch_bottom
				current_area_indices_0based = area_indices_0based_branch_bottom
				events_table_for_branch = NULL
				
				# 2016-05-07: taking the -1 out of nrow(sorted_manual_table)-1
				for (mm in 1:(nrow(sorted_manual_table)) )
					{
					curr_time = event_times[mm]
					
					# Bug check
					if (is.na(curr_time) == TRUE)
						{
						print("NA error in stochastic_map_branch()!")
						print("The program is in the manually-sorted events section, reserved for rare cases where a successful history could not be simulated.")
						print("sorted_manual_table:")
						print(sorted_manual_table)
						print("dim(sorted_manual_table)")
						print(dim(sorted_manual_table))
						print("mm")
						print(mm)
						print("event_times:")
						print(event_times)
						print("curr_time:")
						print(curr_time)
						
						print("current_rangetxt:")
						print(current_rangetxt)
						print("current_area_indices_0based:")
						print(current_area_indices_0based)
						stop("Stopped here 12345")
						}
					
					# Keep track of events
					event_time = NULL
					event_type = NULL
					event_txt = NULL
					dispersal_to = NULL
					extirpation_from = NULL
				
					abs_event_time = abs_time_at_branch_bottom - curr_time
					event_time = curr_time
					
					# Range gain/loss events
					if (sorted_manual_table[mm,2] == "add")
						{
						event_type = "d" # range expansion dispersal
						
						new_area_num_1based = 1+sorted_manual_table[mm,1]
						dispersal_to = areas[new_area_num_1based]
						lost_area_num_1based = "-"
						extirpation_from = "-"
						
						current_area_indices_0based = sort(c(current_area_indices_0based, sorted_manual_table[mm,1]))
						# 2016-05-07_bugfix: add 1 to 0-based areas!!
						new_rangetxt = paste0(areas[current_area_indices_0based+1], collapse="")
						
						# Get the range numbers
						new_range_TF = (new_rangetxt == ranges_list)
						new_rangenum_1based = (1:length(ranges_list))[new_range_TF]
						current_range_TF = (current_rangetxt == ranges_list)
						current_rangenum_1based = (1:length(ranges_list))[current_range_TF]
						}

					if (sorted_manual_table[mm,2] == "sub")
						{
						event_type = "e" # range contraction dispersal
						
						new_area_num_1based = "-"
						dispersal_to = "-"
						lost_area_num_1based = 1+sorted_manual_table[mm,1]
						extirpation_from = areas[lost_area_num_1based]
						
						current_area_indices_0based = current_area_indices_0based[current_area_indices_0based != sorted_manual_table[mm,1]]
						# 2016-05-07_bugfix: add 1 to 0-based areas!!
						new_rangetxt = paste0(areas[current_area_indices_0based+1], collapse="")
						
						# Get the range numbers
						current_range_TF = (current_rangetxt == ranges_list)
						current_rangenum_1based = (1:length(ranges_list))[current_range_TF]
						new_range_TF = (new_rangetxt == ranges_list)
						new_rangenum_1based = (1:length(ranges_list))[new_range_TF]
						}

					# Make a row for a data.table
				
					# event_txt
					event_txt = paste(ranges_list[[current_rangenum_1based]], "->", ranges_list[[new_rangenum_1based]], sep="")
		
					tmprow = c(nodenum_at_top_of_branch, trynum, brlen, current_rangenum_1based, new_rangenum_1based, current_rangetxt, new_rangetxt, abs_event_time, event_time, event_type, event_txt, new_area_num_1based, lost_area_num_1based, dispersal_to, extirpation_from)

					events_table_for_branch = rbind(events_table_for_branch, tmprow)
		
		
					# Update the current branch statenum
					current_rangenum_1based = new_rangenum_1based
					current_rangetxt = new_rangetxt
					} # END for (mm in 1:(nrow(sorted_manual_table)-1) )

				cat("\n\n")
				print("stochastic_map_branch() found this 'manual' (required range add/subtract events, uniform random distribution on branch) events_table_for_branch:")
				cat("\n\n")
				print(events_table_for_branch)
				cat("\n\n")

				# Break out of loop, simulation was successful			
				break()
				} # END if (manual_history_for_difficult_branches == FALSE)
				
			} # END if (trynum > maxtries)
	
		# Otherwise, try again
		#txt = paste("\nSimulation #", trynum, "/", maxtries, " produced ending state ", ranges_list[[current_rangenum_1based]], " (goal: ", ranges_list[[statenum_1based_at_branch_top]], "); trying again.", sep="")
		#cat(txt)
	
		} # END while

	events_table_for_branch = as.data.frame(events_table_for_branch, stringsAsFactors=FALSE)
	
	
	if (length(events_table_for_branch) > 0)
		{
		# Format the events table
		# 2014-05-27_NJM: This had an error as abs_event_time was calculated as NULL in non-stratified analyses
		# due to default NULL state for abs. time branch bottom; fixeds
		#print(events_table_for_branch)
		
		names(events_table_for_branch) = c("nodenum_at_top_of_branch", "trynum", "brlen", "current_rangenum_1based", "new_rangenum_1based", "current_rangetxt", "new_rangetxt", "abs_event_time", "event_time", "event_type", "event_txt", "new_area_num_1based", "lost_area_num_1based", "dispersal_to", "extirpation_from")
		row.names(events_table_for_branch) = NULL
		events_table_for_branch
		
		# ERROR CHECK
		# This is an error check to run on e.g. Psychotria data
		# which all have size 1
		# so if you get a tip with a larger size
		# throw an error.
		error_check_Psychotria_all_tips_size1 = FALSE
		if (error_check_Psychotria_all_tips_size1)
			{
			if (stratified == TRUE)
				{
				# The time top should be zero AND the treepiece should be an orig_tip
				TF1 = trtable$time_top[nodenum_at_top_of_branch] == 0
				TF2 = trtable$piececlass[nodenum_at_top_of_branch] == "orig_tip"
				TF = (TF1 + TF2) == 2
				if (TF == TRUE)
					{
					if (events_table_for_branch[4] > 5)
						{
						print(trtable[nodenum_at_top_of_branch,])
						print(events_table_for_branch)
						errortxt = paste("\n\nERROR!  You have a stochastic mapping tip of state 6 or above, but the Psychotria dataset tips are all state 1-4 (size 1 area).\n\n", sep="")
						cat(errortxt)
				
						stop("Stopping on error.")
						} # END if (events_table_for_branch[4] > 5)
					} # END if (events_table_for_branch[4] > 5)
				} else {
				# The time top should be zero AND the treepiece should be an orig_tip
				TF = trtable$time_bp[nodenum_at_top_of_branch] == 0
				if (TF == TRUE)
					{
					if (events_table_for_branch[4] > 5)
						{
						print(trtable[nodenum_at_top_of_branch,])
						print(events_table_for_branch)
						errortxt = paste("\n\nERROR!  You have a stochastic mapping tip of state 6 or above, but the Psychotria dataset tips are all state 1-4 (size 1 area).\n\n", sep="")
						cat(errortxt)
				
						stop("Stopping on error.")
						} # END if (events_table_for_branch[4] > 5)
					} # END if (events_table_for_branch[4] > 5)				
				}
			} # END if (error_check_Psychotria_all_tips_size1)

		} else {
		# If there were no branch events, put 0
		events_table_for_branch = NA
		} # END if (length(events_table_for_branch) > 0)
	

	return(events_table_for_branch)
	}

# Convert a row of events_table_for_branch into a 
# text string, for later storage in trtable
events_table_row_into_txt <- function(tmprow)
	{
	event_desc_txt = paste("nodenum_at_top_of_branch:", tmprow[1], ",trynum:", tmprow[2], ",brlen:", tmprow[3], ",current_rangenum_1based:", tmprow[4], ",new_rangenum_1based:", tmprow[5], ",current_rangetxt:", tmprow[6], ",new_rangetxt:", tmprow[7], ",abs_event_time:", tmprow[8], ",event_time:", tmprow[9], ",event_type:", tmprow[10], ",event_txt:", tmprow[11], ",new_area_num_1based:", tmprow[12], ",lost_area_num_1based:", tmprow[13], ",dispersal_to:", tmprow[14], ",extirpation_from:", tmprow[15], sep="")
	
	return(event_desc_txt)
	}

events_table_into_txt <- function(events_table_for_branch)
	{
	if (is.null(nrow(events_table_for_branch)))
		{
		branch_events_txt = "none"
		} else {
		# Collapse into a list of lists, ;-delimited
		words = apply(X=events_table_for_branch, MARGIN=1, FUN=events_table_row_into_txt)
		branch_events_txt = paste(words, sep="", collapse=";")
		}
	return(branch_events_txt)
	}

events_txt_list_into_events_table <- function(events_txt_list, trtable=NULL, recalc_abs_ages=TRUE)
	{
	defaults='
	events_txt_list=master_table_w_stochastic_maps$anagenetic_events_txt_below_node
	trtable=master_table_w_stochastic_maps
	recalc_abs_ages=TRUE
	'
	
	if ((length(trtable) > 0) && (class(trtable) == "list"))
		{
		txt = "STOP ERROR in events_txt_list_into_events_table(). Input 'trtable' was a list, but should be a data.frame table, or empty. Printing 'trtable':"
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		print(trtable)
		
		cat("\n\n")
		stop(txt)		
		}
	
	# 2014-05-27_NJM: Note that if you have NO events in the WHOLE TREE, you will
	# have a "NULL" in events_txt_list.
	# 
	# Solution: test for NULL. Really, you should pre-allocate this so that 
	# it's never NULL.
	#
	# Actually, this wasn't the problem since I already did pre-allocation. The problem
	# was looking for BioGeoBEARS_run_object$states_list when it should have been 
	# looking for res$inputs$states_list
	if (is.null(events_txt_list))
		{
		errortxt = paste("\nWARNING in events_txt_list_into_events_table(): your events_txt_list has NO events!\n\nThis means your tree has NO d/e/a events across the whole tree.\nThis is *expected* e.g. if you inferred d=e=0 under DEC+J. Input a list of '' or NA to avoid this error.\n\n", sep="")
		cat(errortxt)
		errortxt2 = paste("events_txt_list_into_events_table() is returning NULL which will might cause issues later.\n\n", sep="")
		cat(errortxt2)
		return(NULL)
		}
	
	# Convert NAs to "none"
	events_txt_list[is.na(events_txt_list)] = "none"
	
	# Remove lines with no events or NA:
	noneTF = events_txt_list == "none"
	keepTF = (noneTF == FALSE)
	events_txt_list = events_txt_list[keepTF]
	
	#print(events_txt_list)
	
	# If no anagenetic events, return NULL
	if (length(events_txt_list) == 0)
		{
		events_table = NULL
		return(events_table)
		}

	#print("here3")	

	# Include the trtable, if that is input
	if (length(trtable) > 0)
		{
		trtable_subset = NULL
		}

	#print("here4")	

	
	# Convert the events text back into a table:
	tmptable = NULL
	for (i in 1:length(events_txt_list))
		{
		#print(events_txt_list)
		tmptable_rows = events_txt_into_events_table(events_txt_list[i])
		# NJM 2015-06-08
		# NJM 2016-05-05 bug fix: add "as.numeric"
		rownums_in_trtable = as.numeric(tmptable_rows$nodenum_at_top_of_branch)
		#print(tmptable_rows)
		num_newrows = nrow(tmptable_rows)
		tmptable = rbind(tmptable, tmptable_rows)

		# Include the trtable, if that is input
		if (length(trtable) > 0)
		#if ( (is.null(trtable) == FALSE) && (trtable != list()) )
			{
			for (nnr in 1:num_newrows)
				{
				#print(trtable)
				#print(keepTF)
				#print(trtable)
				# NJM 2015-04-05
				#trtable_subset = rbind(trtable_subset, trtable[keepTF,][nnr,])
				# NJM 2015-06-08
				trtable_subset = rbind(trtable_subset, trtable[rownums_in_trtable[nnr],])
				} # END for (nnr in 1:num_newrows)
			} # END if (length(trtable) > 0)
		} # END for (i in 1:length(events_txt_list))
	events_table = dfnums_to_numeric(adf2(tmptable))
	names(events_table) = c("nodenum_at_top_of_branch", "trynum", "brlen", "current_rangenum_1based", "new_rangenum_1based", "current_rangetxt", "new_rangetxt", "abs_event_time", "event_time", "event_type", "event_txt", "new_area_num_1based", "lost_area_num_1based", "dispersal_to", "extirpation_from")

	# Include the trtable, if that is input
	if (length(trtable) > 0)
		{
		# Remove the event txt column
		trtable_subset_col_keepTF = names(trtable_subset) != "anagenetic_events_txt_below_node"
		trtable_subset = trtable_subset[,trtable_subset_col_keepTF]
		events_table = cbind(trtable_subset, events_table)
		
		# Recalculate the absolute ages of anagenetic events
		# (necessary for stratified analysis master_table
		# (And, you can only do when SUBedge.length exists as a 
		#  table column)
		if ((recalc_abs_ages == TRUE) && ( ("SUBedge.length" %in% names(events_table))==TRUE ) )
			{
			
			# 2014-05-25_NJM: checking for SUBedge.length too long on 
			# tip branches
			# Is reltimept smaller?
			#naTF = is.na(events_table$SUBedge.length)
			#events_table$SUBedge.length[naTF] = 0

			#SUBedge_too_big_TF = events_table$reltimept < events_table$SUBedge.length
			#brlen_in_section = events_table$SUBedge.length
			
			#brlen_in_section[SUBedge_too_big_TF] = events_table$reltimept[SUBedge_too_big_TF]
			
			# Absolute event age = node age at top of branch, PLUS
			# the length of the SUBedge.length below, MINUS the relative
			# event time above the bottom
			#old_abs_event_time = events_table$abs_event_time
			#events_table$abs_event_time = events_table$time_bp + events_table$SUBedge.length - events_table$event_time
			#events_table$abs_event_time = events_table$time_top + events_table$SUBtime_bp + brlen - events_table$event_time
			
			#old_abs_event_time
			events_table$abs_event_time = events_table$abs_event_time 
			} # END if doing recalculation
		} # END if (is.null(trtable) == FALSE)
	
	
	return(events_table)
	}



events_txt_into_events_table <- function(branch_events_txt)
	{
	# Split on semi-colon
	branch_events_txt = as.character(branch_events_txt)
	words = strsplit(branch_events_txt, split=";")[[1]]
	
	events_table_for_branch = t(sapply(X=words, FUN=event_txt_into_events_row))
	row.names(events_table_for_branch) = NULL
	events_table_for_branch
	
	events_table_for_branch = adf2(events_table_for_branch)
	events_table_for_branch
	names(events_table_for_branch) = c("nodenum_at_top_of_branch", "trynum", "brlen", "current_rangenum_1based", "new_rangenum_1based", "current_rangetxt", "new_rangetxt", "abs_event_time", "event_time", "event_type", "event_txt", "new_area_num_1based", "lost_area_num_1based", "dispersal_to", "extirpation_from")
	
	return(events_table_for_branch)
	}

event_txt_into_events_row <- function(word)
	{
	words2 = strsplit(word, split=",")[[1]]
	output = sapply(X=words2, FUN=split_key_item)
	tmprow = matrix(data=output[2,], nrow=1)
	return(tmprow)
	}


split_key_item <- function(word2)
	{
	output_pair = c("", "")
	words3 = strsplit(word2, split=":")[[1]]
	numwords = length(words3)
	
	output_pair[1:numwords] = words3[1:numwords]
	
	return(output_pair)
	}


add_cladogenetic_events_to_trtable <- function(trtable, BioGeoBEARS_run_object, tipranges, stratified=FALSE, piecenum=NULL)
	{
	#######################################################
	# Label the cladogenetic events
	#######################################################
	
	if (stratified == FALSE)
		{
		ntips = nrow(tipranges@df)
		nodenums = (ntips+1):(ntips+(ntips-1))
		
		# Need this to generalize for both stratified/non-stratified
		just_this_subtree_table_TF = rep(TRUE, (ntips+length(nodenums)))
		
		} else {
		# Count the number of subtree tip nodes
		# (counting tips doesn't work, since some of those are single branches)
		just_this_subtree_table_TF = trtable$piecenum == piecenum
		just_this_subtree_table = trtable[just_this_subtree_table_TF,]
		tip_nodes_TF = just_this_subtree_table$SUBnode.type == "tip"
		ntips = sum(tip_nodes_TF)
		nodenums = (ntips+1):(ntips+(ntips-1))
		}
	
	# Get the list of geographic areas
	#tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=BioGeoBEARS_run_object$geogfn)
	areas = getareas_from_tipranges_object(tipranges)
	areas_list = seq(0, length(areas)-1, 1)		# 0-base indexes
	areanames = areas
	areas
	areas_list

	# Initialize key information about the state space
	include_null_range = BioGeoBEARS_run_object$include_null_range
	maxareas = BioGeoBEARS_run_object$max_range_size


	# Get the 0-based states list, if needed	
	if (is.null(BioGeoBEARS_run_object$states_list))
		{
		state_indices_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=BioGeoBEARS_run_object$max_range_size, include_null_range=include_null_range)
		state_indices_0based
		states_list = state_indices_0based

		# Get the ranges
		ranges_list = areas_list_to_states_list_new(areas=areas, maxareas=maxareas,
		include_null_range=include_null_range, split_ABC=FALSE)
		ranges_list
		ranges = unlist(ranges_list)
		ranges
		} else {
		
		state_indices_0based = BioGeoBEARS_run_object$states_list
		states_list = state_indices_0based
		
		# Get the list of text ranges
		ranges_list = list()
		for (i in 1:length(states_list))
			{
			if ( is.na(states_list[[i]]) )
				{
				ranges_list[[i]] = "_"
				} else {
				ranges_list[[i]] = paste(areas[ 1+states_list[[i]] ], sep="", collapse="")
				}
			}
		ranges_list
		ranges = unlist(ranges_list)
		ranges		
		}


	# Build the output table
	cols_to_translate_cladogenesis_events = c("sampled_states_AT_nodes", "samp_LEFT_dcorner", "samp_RIGHT_dcorner")
	table_of_cladogenetic_events_translated = trtable[just_this_subtree_table_TF,][nodenums, cols_to_translate_cladogenesis_events]
	names(table_of_cladogenetic_events_translated) = c("parent_range", "Left_state", "Right_state")
	table_of_cladogenetic_events_translated

	# Label/classify each cladogenetic event
	cladogenetic_event_labels = label_table_of_cladogenetic_events(table_of_cladogenetic_events_translated, states_list, ranges_list, track_dispersal_dest=TRUE, areanames=areanames)

	blanks = rep("", times=ntips)
	trtable$clado_event_type[just_this_subtree_table_TF] = c(blanks, cladogenetic_event_labels$event_type)
	trtable$clado_event_txt[just_this_subtree_table_TF] = c(blanks, cladogenetic_event_labels$event_txt)
	trtable$clado_dispersal_to[just_this_subtree_table_TF] = c(blanks, cladogenetic_event_labels$dispersal_to)
	trtable[just_this_subtree_table_TF,][nodenums,]

	return(trtable)
	} # END add_cladogenetic_events_to_trtable


# Get inputs for non-stratified, and stratified, Biogeographical Stochastic Mapping analysis
get_inputs_for_stochastic_mapping <- function(res, cluster_already_open=FALSE, rootedge=FALSE, statenum_bottom_root_branch_1based=NULL, printlevel=1, min_branchlength=0.000001)
	{
	defaults='
	res=res
	cluster_already_open=NULL
	rootedge=FALSE
	statenum_bottom_root_branch_1based=NULL
	printlevel=1
	min_branchlength=0.000001
	' # END defaults
	
	# Necessary setting to avoid getting numbers etc. in the stochastic mapping output
	options(stringsAsFactors=FALSE)



	# High performance computing
	# HPC using parallel package in R 2.15 or higher, which allows
	# mcmapply (multicore apply)
	# Don't use multicore if using R.app ('AQUA')
	num_cores_to_use = res$inputs$num_cores_to_use
	cluster_already_open = cluster_already_open

	cluster_was_open = FALSE
	if (.Platform$GUI != "AQUA" && ((is.na(num_cores_to_use) == TRUE) || ( (is.na(num_cores_to_use)==FALSE) && (num_cores_to_use > 1))) )
		{
		# We are doing manual, optional processing on several cores;
		# this seems to have less overhead/hassle/incompatibility issues
		# than using mcmapply, mclapply, etc...
		#require("parallel") #<- do this higher up

		num_cores_computer_has = detectCores()
		txt = paste0("Your computer has ", num_cores_computer_has, " cores.")
		cat("\n")
		cat(txt)
		cat("\n")		

		if (num_cores_to_use > num_cores_computer_has)
			{
			txt = paste0("WARNING from bears_optim_run(): You specified num_cores_to_use=", num_cores_to_use, " cores, but your computer only has ", num_cores_computer_has, ". Resetting to ", num_cores_computer_has, ".")
			cat("\n")
			cat(txt)
			cat("\n")
			warning(txt)
			num_cores_to_use = num_cores_computer_has
			}
		
		
		if (is.null(num_cores_to_use))
			{
			num_cores_to_use = num_cores_computer_has
			}

		# Don't do this, if the cluster is already open
		cat("\nYour computer has ", num_cores_computer_has, " cores. You have chosen to use:\nnum_cores_to_use = ", num_cores_to_use, " cores for the matrix exponentiations in the likelihood calculations.\n", sep="")
		
		if ( is.logical(cluster_already_open) == TRUE )
			{
			if (cluster_already_open == FALSE)
				{
				cluster_already_open = makeCluster(rep("localhost",num_cores_to_use), type = "SOCK")
				cat("Started cluster with ", num_cores_to_use, " cores.\n\n", sep="")

				# Flag so that you remember to close cluster at the end
				cluster_open=TRUE
				cluster_was_open = FALSE
				}
			} else {
			cluster_was_open = TRUE
			cat("Cluster with ", num_cores_to_use, " cores already open.\n\n", sep="")
			}
		} else {
		# You are using R.app and clusters don't work...

		num_cores_computer_has = detectCores()
		txt = paste0("Your computer has ", num_cores_computer_has, " cores.")
		cat("\n")
		cat(txt)
		cat("\n")	
		
		if (num_cores_to_use > 1)
			{
			txt = paste0("WARNING from bears_optim_run(): You specified num_cores_to_use=", num_cores_to_use, " cores, but in R.app, multicore functionality doesn't work. Resetting num_cores_to_use=1.")
			cat("\n")
			cat(txt)
			cat("\n")
			warning(txt)
			num_cores_to_use = num_cores_computer_has
			}
		
		
		cluster_already_open = NULL
		cluster_was_open = FALSE
		}

	cat("\n\nCurrently, cluster_already_open=")
	print(cluster_already_open)




	
	# Stratified or non-stratified -- is there a times filename? (timesfn)
	if (is.na(res$inputs$timesfn) == TRUE)
		{
		# Non-stratified inputs
		stratified = FALSE

		if (printlevel >= 1)
			{
			cat("\nDetected NON-stratified input. Running get_inputs_for_stochastic_mapping_from_results_object().\n")
			}

		stochastic_mapping_inputs = get_inputs_for_stochastic_mapping_from_results_object(res=res, cluster_already_open=cluster_already_open, rootedge=rootedge, statenum_bottom_root_branch_1based=statenum_bottom_root_branch_1based, printlevel=printlevel, stratified=FALSE, min_branchlength=min_branchlength, timeperiod_i=1)		
		} else {
		# Stratified inputs
		stratified = TRUE
		
		if (printlevel >= 1)
			{
			cat("\nDetected STRATIFIED input. Running get_inputs_for_stochastic_mapping_stratified().\n")
			}
		
		stochastic_mapping_inputs = get_inputs_for_stochastic_mapping_stratified(res=res, cluster_already_open=cluster_already_open, rootedge=rootedge, statenum_bottom_root_branch_1based=statenum_bottom_root_branch_1based, printlevel=printlevel, min_branchlength=min_branchlength)	
		} # END if (is.na(res$inputs$timesfn) == TRUE)
	
	return(stochastic_mapping_inputs)
	}


# Get inputs for non-stratified Biogeographical Stochastic Mapping analysis
get_inputs_for_stochastic_mapping_from_results_object <- function(res, cluster_already_open=FALSE, rootedge=FALSE, statenum_bottom_root_branch_1based=NULL, printlevel=1, stratified=FALSE, min_branchlength=0.000001, timeperiod_i=1)
	{
	defaults='
	cluster_already_open=NULL
	rootedge=FALSE
	statenum_bottom_root_branch_1based=NULL
	printlevel=1
	stratified=FALSE
	min_branchlength=0.000001
	timeperiod_i=1
	'
	# Necessary setting to avoid getting numbers etc. in the stochastic mapping output
	options(stringsAsFactors=FALSE)
		
	# Get the input run object
	#BioGeoBEARS_run_object = res$inputs
	
	# Error check
	if ((stratified == FALSE) && (timeperiod_i != 1))
		{
		errortxt = paste("\n\nERROR in get_inputs_for_stochastic_mapping_from_results_object(): stratified=FALSE, but timeperiod_i is: ", timeperiod_i, ";\ntimeperiod_i=1 is required if stratified=FALSE (since you only have 1 stratum!\n\n", sep="")
		cat(errortxt)
		
		stop("Stopping on error.")
		}
	
	
	stochastic_mapping_inputs = list()
	
	# Load defaults
	stochastic_mapping_inputs$rootedge = rootedge
	stochastic_mapping_inputs$statenum_bottom_root_branch_1based = statenum_bottom_root_branch_1based
	stochastic_mapping_inputs$printlevel = printlevel
	stochastic_mapping_inputs$stratified = stratified
	stochastic_mapping_inputs$min_branchlength = min_branchlength
	# Get likelihoods on cluster, if desired
	stochastic_mapping_inputs$cluster_already_open = cluster_already_open


	# Get the number of states
	numstates = ncol(res$ML_marginal_prob_each_state_at_branch_top_AT_node)
	
	

	

	# Get the tree
	if (stratified == FALSE)
		{
		#tr = read.tree(res$inputs$trfn)
		tr = check_trfn(trfn=res$inputs$trfn)
		phy2 = reorder(tr, "pruningwise") # Do this
		} else {
		# Get the list of tree pieces for timeperiod_i
		tree_sections_list = res$inputs$tree_sections_list[[timeperiod_i]]
		}

	# Get geographic ranges at tips
	if (res$inputs$use_detection_model == FALSE)
		{
		tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=res$inputs$geogfn)
		}
	if (res$inputs$use_detection_model == TRUE)
		{
		tipranges = tipranges_from_detects_fn(detects_fn=res$inputs$detects_fn)
		} # END if (res$inputs$use_detection_model == TRUE)



	# Get the areas and states list
	areas = getareas_from_tipranges_object(tipranges)
	areas_list = seq(0, length(areas)-1, 1)		# 0-base indexes
	areanames = areas
	areas
	areas_list

	# Get the 0-based states list, if needed
	include_null_range = res$inputs$include_null_range
	if (is.null(res$inputs$states_list))
		{
		state_indices_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=res$inputs$max_range_size, include_null_range=include_null_range)
		state_indices_0based
		states_list = state_indices_0based

		# Get the ranges
		ranges_list = areas_list_to_states_list_new(areas=areas, maxareas=res$inputs$max_range_size,
		include_null_range=include_null_range, split_ABC=FALSE)
		ranges_list
		ranges = unlist(ranges_list)
		ranges
		} else {
		state_indices_0based = res$inputs$states_list
		states_list = state_indices_0based
		
		# Get the list of text ranges
		ranges_list = list()
		for (i in 1:length(states_list))
			{
			if ( is.na(states_list[[i]]) )
				{
				ranges_list[[i]] = "_"
				} else {
				ranges_list[[i]] = paste(areas[ 1+states_list[[i]] ], sep="", collapse="")
				}
			}
		ranges_list
		ranges = unlist(ranges_list)
		ranges	
		}



	# Make the tree_table
	if (stratified == FALSE)
		{
		# Make a table holding the states etc.
		trtable = prt(phy2, printflag=FALSE)
		trtable
		} else {
		# Get the table corresponding to the tree pieces for timeperiod_i
		TF = res$inputs$master_table$stratum == timeperiod_i
		master_table_timeperiod_i = res$inputs$master_table[TF,]
		}


	
	#######################################################
	# Get the Qmat, etc., and calculate the independent 
	# likelihoods (ONCE!)
	#######################################################
	
	# timeperiod_i will be 1 for non-stratified analysis
	returned_mats2 = get_Qmat_COOmat_from_res(res, timeperiod_i=timeperiod_i, include_null_range=include_null_range)
	returned_mats2

	# Extract output
	Qmat = returned_mats2$Qmat
	COO_weights_columnar = returned_mats2$COO_weights_columnar
	Rsp_rowsums = returned_mats2$Rsp_rowsums

	######################################
	# Pre-calculate the independent likelihoods on each branch
	######################################
	
	cat("\n\nPre-calculating the independent likelihoods on each branch...\n")
	cat("Currently, cluster_already_open=")
	print(cluster_already_open)
	cat("\n\n")
	
	if (stratified == FALSE)
		{
		crash_error_check = FALSE
		if (crash_error_check == TRUE)
			{
			print("Checking error")
			print(phy2)
			print(Qmat)
			print(cluster_already_open)
			} # END if (crash_error_check = TRUE)
		
		independent_likelihoods_on_each_branch = calc_independent_likelihoods_on_each_branch(phy2, Qmat, cluster_already_open=cluster_already_open, Qmat_is_sparse=FALSE)
		} else {
		# This gets interesting with a stratified analysis!!!!
		tree_pieces = tree_sections_list$return_pieces_list
		# Go through the tree pieces
		
		independent_likelihoods_by_tree_piece = list()
		
		for (tp in 1:length(tree_pieces))
			{
			# What sort of tree piece is it?
			if (is.numeric(tree_pieces[[tp]]))
				{
				# If it's a single branch section, it's numeric (just a number)
				branch_length = tree_pieces[[tp]]
				independent_likelihoods_on_single_branch = expokit_dgpadm_Qmat2(times=branch_length, Qmat=Qmat, transpose_needed=TRUE)

				independent_likelihoods_by_tree_piece[[tp]] = independent_likelihoods_on_single_branch
				} else {
				# If it's a subtree
				phy = tree_pieces[[tp]]
				phy2 = reorder(phy, "pruningwise")
				independent_likelihoods_on_each_branch = calc_independent_likelihoods_on_each_branch(phy2, Qmat, cluster_already_open=cluster_already_open, Qmat_is_sparse=FALSE)

				independent_likelihoods_by_tree_piece[[tp]] = independent_likelihoods_on_each_branch
				} # END if (is.numeric(tree_pieces[[tp]]))
			} # END for (tp in 1:length(tree_pieces))
		} # END if (stratified == FALSE)
	
	# Store the results and model
	stochastic_mapping_inputs$res = res
	stochastic_mapping_inputs$Qmat = Qmat
	stochastic_mapping_inputs$COO_weights_columnar = COO_weights_columnar
	stochastic_mapping_inputs$Rsp_rowsums = Rsp_rowsums

	# Store the tree info
	stochastic_mapping_inputs$tipranges = tipranges
	stochastic_mapping_inputs$areas = areas
	stochastic_mapping_inputs$state_indices_0based = state_indices_0based
	stochastic_mapping_inputs$ranges_list = ranges_list
	stochastic_mapping_inputs$numstates = numstates
	
	if (stratified == FALSE)
		{
		stochastic_mapping_inputs$independent_likelihoods_on_each_branch = independent_likelihoods_on_each_branch
		stochastic_mapping_inputs$trtable = trtable
		stochastic_mapping_inputs$tr = tr
		stochastic_mapping_inputs$phy2 = phy2
		} else {
		# The independent likelihoods for the branches
		# of JUST THIS STRATUM
		stochastic_mapping_inputs$independent_likelihoods_by_tree_piece_for_timeperiod_i = independent_likelihoods_by_tree_piece
		
		# The tree section for this stratum (one or more pieces, 3 items per piece)
		stochastic_mapping_inputs$tree_sections_list = tree_sections_list
		

		# Add the necessary columns to the subtable
		# States at nodes and branch bottoms below nodes
		sampled_states_AT_nodes = rep(NA, nrow(master_table_timeperiod_i))
		sampled_states_AT_brbots = rep(NA, nrow(master_table_timeperiod_i))
	
		# Add the right and left descendant node numbers
		left_desc_nodes = rep(NA, nrow(master_table_timeperiod_i))
		right_desc_nodes = rep(NA, nrow(master_table_timeperiod_i))

		# dcorner = descendant corner (i.e. right after speciation)
		samp_LEFT_dcorner = rep(NA, nrow(master_table_timeperiod_i))
		samp_RIGHT_dcorner = rep(NA, nrow(master_table_timeperiod_i))
		
		# Cladogenesis events
		clado_event_type = rep(NA, nrow(master_table_timeperiod_i))
		clado_event_txt = rep(NA, nrow(master_table_timeperiod_i))
		clado_dispersal_to = rep(NA, nrow(master_table_timeperiod_i))

		# Events on branches
		anagenetic_events_txt_below_node = rep(NA, nrow(master_table_timeperiod_i))

		master_table_timeperiod_i = cbind(master_table_timeperiod_i, sampled_states_AT_nodes, sampled_states_AT_brbots, left_desc_nodes, right_desc_nodes, samp_LEFT_dcorner, samp_RIGHT_dcorner, clado_event_type, clado_event_txt, clado_dispersal_to, anagenetic_events_txt_below_node)

		# Get the left/right node nums for each subtree in the tree pieces
		tree_pieces = stochastic_mapping_inputs$tree_sections_list$return_pieces_list
		num_pieces = length(tree_pieces)
		for (piecenum in 1:num_pieces)
			{
			tree_piece = tree_pieces[[piecenum]]
			tree_piece_rows_TF = master_table_timeperiod_i$piecenum == piecenum
			
			# If it's just a branch section, put NA for the descendant nodes
			if (is.numeric(tree_piece))
				{
				# It's a single branch
				master_table_timeperiod_i$left_desc_nodes[tree_piece_rows_TF] = NA
				master_table_timeperiod_i$right_desc_nodes[tree_piece_rows_TF] = NA
				} else {
				# It's a phylo object 
				leftright_nodes_matrix = get_leftright_nodes_matrix_from_results(tree_piece)
				
				subtable_SUBnodenums_to_put_in_subtable = row.names(leftright_nodes_matrix)
				
				row_indices_to_use = match(x=subtable_SUBnodenums_to_put_in_subtable, table=master_table_timeperiod_i$SUBnode)
				
				master_table_timeperiod_i$left_desc_nodes[row_indices_to_use] = leftright_nodes_matrix$right
				master_table_timeperiod_i$right_desc_nodes[row_indices_to_use] = leftright_nodes_matrix$left
				master_table_timeperiod_i
				} # END if (is.numeric(tree_piece))
			} # END for (piecenum in 1:num_pieces)
		
		# Save the subtable
		stochastic_mapping_inputs$master_table_timeperiod_i = master_table_timeperiod_i

		} # END if (stratified == FALSE)

	return(stochastic_mapping_inputs)
	} # END get_inputs_for_stochastic_mapping_from_results_object <- function(res, cluster_already_open=FALSE, rootedge=FALSE, statenum_bottom_root_branch_1based=NULL, printlevel=1, stratified=FALSE, min_branchlength=0.000001, timeperiod_i=1)



#######################################################
# Stochastic mapping on a non-stratified analysis
#######################################################
# 
# This function does stochastic mapping on non-stratified analysis, or 
# calls stochastic_mapping_on_stratified() if the input is stratified.
#  
# stochastic_mapping_inputs = get_inputs_for_stochastic_mapping_from_results_object(res=resDEC)
# stochastic_mapping_results = stochastic_map_given_inputs(stochastic_mapping_inputs)
# 
# For stochastic mapping on stratified analysis:
# 
# stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping_stratified(res=resDEC)
# stochastic_mapping_results = stochastic_mapping_on_stratified(res=resDEC, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list)
# 


#######################################################
# get_inputs_for_stochastic_mapping_stratified()
# This is a separate function, because it only has
# to be run once.
#######################################################
get_inputs_for_stochastic_mapping_stratified <- function(res, cluster_already_open=FALSE, rootedge=FALSE, statenum_bottom_root_branch_1based=NULL, printlevel=1, min_branchlength=0.000001)
	{
	defaults='
	res=res
	cluster_already_open=NULL
	rootedge=FALSE
	statenum_bottom_root_branch_1based=NULL
	printlevel=1
	min_branchlength=0.000001
	
	timeperiod_i=1
	stratified=TRUE 
	
	' # END defaults
	
	# Necessary setting to avoid getting numbers etc. in the stochastic mapping output
	options(stringsAsFactors=FALSE)
	
	# Get the input run object
	#BioGeoBEARS_run_object = res$inputs
	
	# The conditional likelihoods downpass table for all nodes
	# and stratum breakpoints
	dim(res$condlikes_table)

	# Stratified analysis
	if ((is.numeric(res$inputs$timeperiods))) #&& (length(inputs$timeperiods) > 1))
		{
		stratified = TRUE
		} else {
		errortxt = paste("\n\nERROR in get_inputs_for_stochastic_mapping_stratified(): Your input results_object does not appear to be from a time-stratified analysis.\n\n", sep="")
		stop(errortxt)
		}

	# Get stochastic mapping inputs for stratified analysis
	stochastic_mapping_inputs_list = list()
	
	#number_of_timepoints_at_which_to_calculate_likelihoods = length(unique(res$inputs$stratum))
	# length(res$inputs$tree_sections_list)
	
	for (timeperiod_i in 1:length(res$inputs$tree_sections_list))
		{
		stochastic_mapping_inputs = get_inputs_for_stochastic_mapping_from_results_object(res=res, cluster_already_open=cluster_already_open, rootedge=rootedge, statenum_bottom_root_branch_1based=statenum_bottom_root_branch_1based, printlevel=printlevel, stratified=stratified, min_branchlength=min_branchlength, timeperiod_i=timeperiod_i)
	
		stochastic_mapping_inputs_list[[timeperiod_i]] = stochastic_mapping_inputs
		} # END for (timeperiod_i in 1:length(res$inputs$tree_sections_list))
	
	# Specify that this is for a stratified analysis
	#stochastic_mapping_inputs_list$stratified = TRUE
	
	return(stochastic_mapping_inputs_list)
	} # END get_inputs_for_stochastic_mapping_stratified()







# Not used in default stochastic mapping (2016-05-07)
make_BSMs <- function(res1, model_name="DEC", BSM_tables_dir="BSM_tables_M3_v1", BSM_fn_base="BSM_M3_", num_maps=5, maxtries=1000000, groupname="", seedval=as.numeric(Sys.time()))
	{
	#######################################################
	#######################################################
	# Make stochastic maps under ML model
	#######################################################
	#######################################################

	defaults='
	# Directory to hold individual stochastic map objects
	# ("master_table_cladogenetic_events")
	res1 = resDEC
	model_name = "DEC"
	
	BSM_tables_dir = "BSM_tables_M3_v1"
	BSM_fn_base = "BSM_M3_"
	num_maps = 110
	maxtries = 1000000
	groupname = ""
	'
	# Specify cluster_already_open
	cluster_already_open = res1$cluster_already_open
	
	# Error check; include_null_range=TRUE is default, even for old models
	include_null_range = res1$inputs$include_null_range
	if ( (length(include_null_range) == 0) || (is.na(include_null_range)) || (is.null(include_null_range))  ) 
		{
		include_null_range = TRUE
		res1$inputs$include_null_range = include_null_range
		}

	# Error check; include_null_range=TRUE is default, even for old models
	use_detection_model = res1$inputs$use_detection_model
	if ( (length(use_detection_model) == 0) || (is.na(use_detection_model)) || (is.null(use_detection_model))  ) 
		{
		use_detection_model = TRUE
		res1$inputs$use_detection_model = use_detection_model
		}
	
	
	
	# Is res stratified?
	if (is.na(res1$inputs$timesfn) == TRUE)
		{
		stratified = FALSE
		res1$inputs$stratified = FALSE
		stochastic_mapping_inputs = get_inputs_for_stochastic_mapping_from_results_object(res=res1, cluster_already_open=cluster_already_open)
		} else {
		stratified = TRUE
		res1$inputs$stratified = TRUE
		stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping_stratified(res=res1, cluster_already_open=cluster_already_open)

		# Get timeperiods
		#timeperiods = read_times_fn(inputs=res1$inputs)
		#timeperiods
		}
	
	
	# Make stochastic maps under ML DEC model
	
	# Get tree
	#trfn = res1$inputs$trfn
	#tr = read.tree(trfn)
	
	#getwd()
	
	# Store the filenames
	BSM_fns = NULL
	
	#seed = 12345
	#for (mapnum in 1:num_maps)
	
	# Number of maps kept
	mapnum = 1
	numBSM_tries = 0
	
	# If you have a catastrophic failure rate, fail the function
	numBSM_tries_max = num_maps * 10 
	
	while( mapnum <= num_maps )
		{
		#######################################################
		# Stochastic mapping under DEC
		#######################################################

		analysis_titletxt = paste("Stochastic map #", mapnum, "/", num_maps, " under ", model_name, " model...", sep="")
	
		cat("\n")
		cat(analysis_titletxt)
	
	
		#######################################################
		# Do stochastic mapping conditional on the ML model parameters
		#######################################################
		tempseed = seedval + numBSM_tries + mapnum
		#print(paste0("make_BSMs tempseed=", tempseed, sep=""))
		if (stratified == TRUE)
			{
			# Stratified stochastic mapping
			cmdstr = paste("stochastic_mapping_results = stochastic_mapping_on_stratified(res=res1, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxtries=", maxtries, ", seedval=", tempseed, ")", sep="")
			# (add include_null_range??? 2016-05-05 )
			} else {
			# Non-stratified stochastic mapping
			cmdstr = paste("stochastic_mapping_results = stochastic_map_given_inputs(stochastic_mapping_inputs=stochastic_mapping_inputs, maxtries=", maxtries, ", seedval=", tempseed, ", include_null_range=include_null_range)", sep="")
			}
		
		# Increment the seed and try a stochastic map
		#set.seed(seedval + numBSM_tries + mapnum)
		try_SM = try(expr=eval(parse(text=cmdstr)))

		if (class(try_SM) == "try-error")
			{
			cat("failed, not saving.")
			# Increment either way (especially for seed!)
			numBSM_tries = numBSM_tries + 1
			next()
			}
			
		stochastic_mapping_results = try_SM
		master_table_cladogenetic_events = stochastic_mapping_results$master_table_cladogenetic_events
		BSM_results = stochastic_mapping_results
	
		# Save the output
		tmpi = sprintf("%05.0f", mapnum)
		BSMfn = slashslash(paste(addslash(BSM_tables_dir), BSM_fn_base, model_name, "_", tmpi, ".Rdata", sep=""))
		cat("\nBSMfn: '", BSMfn, "'...", sep="")

		# Loads to: BSM_results
		save(BSM_results, file=BSMfn)
		cat("saved.")
		# Store the list of filenames
		BSM_fns = c(BSM_fns, BSMfn)
		mapnum = mapnum + 1
		
		# Break if it's not working at all
		# (90%+ failure rate)
		if (numBSM_tries > numBSM_tries_max)
			{
			errortxt = paste("\n\nERROR: stopping BSM on on this dataset since numBSM_tries > numBSM_tries_max\n(", numBSM_tries, " > ", numBSM_tries_max, "\n\n", sep="")
			cat(errortxt)
			break;
			}
		} # END for (mapnum in 1:num_maps)
	
	return(BSM_fns)
	} # END make_BSMs_stratified




plot_BSM_fns <- function(res1, res2=NULL, BSM_fns_res1, BSM_fns_res2=NULL, res1name="", res2name="", pdffn="default.pdf", pdf_mfrow=NULL, pdfwidth=NULL, pdfheight=NULL, groupname="", include_null_range=TRUE)
	{
	defaults = '
	include_null_range=TRUE
	areanames = c("K","O","M","H")
	max_range_size = 4
	groupname = ""
	pdffn = "stochastic_maps_DEC_vs_DECj_M0_v1.pdf"

	'
	
	# Setup
	# PDF width
	if (is.null(pdfwidth) == TRUE)
		{
		if (is.null(res2) == FALSE)
			{
			pdfwidth = 12
			} else {
			pdfwidth = 6
			}
		}	

	# PDF height
	if (is.null(pdfheight) == TRUE)
		{
		pdfheight = 6
		}	
	
	
	# Start the pdf
	pdf(file=pdffn, width=pdfwidth, height=pdfheight)
	
	# If 2 models, do 2 subplots side-by-side (default)
	if (is.null(pdf_mfrow) == TRUE)
		{
		if (is.null(res2) == FALSE)
			{
			pdf_mfrow = c(1,2)
			par(mfrow=pdf_mfrow)
			}
		} # END if (is.null(pdf_mfrow) == TRUE)
	
	# Get the tree
	trfn = res1$inputs$trfn
	#tr = read.tree(trfn)
	tr = check_trfn(trfn=trfn)
	
	# Get the tipranges, area names, and max_range_size
	# And states_list
	geogfn = res1$inputs$geogfn
	
	# Get geographic ranges at tips
	if (res1$inputs$use_detection_model == FALSE)
		{
		tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=res1$inputs$geogfn)
		}
	if (res1$inputs$use_detection_model == TRUE)
		{
		tipranges = tipranges_from_detects_fn(detects_fn=res1$inputs$detects_fn)
		} # END if (res1$inputs$use_detection_model == TRUE)

	
	
	areas = getareas_from_tipranges_object(tipranges)
	areanames = areas
	max_range_size = res1$inputs$max_range_size
	states_list_0based = res1$inputs$states_list_0based
	
	if (is.null(states_list_0based) == TRUE)
		{
		states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)
		} # END if (is.null(states_list_0based) == TRUE)
	
	
	num_maps = length(BSM_fns_res1)
	for (mapnum in 1:num_maps)
		{
		#######################################################
		# Stochastic mapping under e.g. DEC
		#######################################################
		if (groupname == "")
			{
			groupname_txt = ""
			} else {
			groupname_txt = paste(" on ", groupname, sep="")
			}
		
		analysis_titletxt = paste("Extracting stochastic map #", mapnum, "/", num_maps, " under model ", res1name, groupname_txt, sep="")
	
		cat("\n")
		cat(analysis_titletxt)

		
		#######################################################
		# Load stochastic mapping result for res1
		#######################################################
		BSM_fn = BSM_fns_res1[mapnum]
		cat("...loading '", BSM_fn, "'...", sep="")
		
		# Loads to BSM_results
		load(BSM_fn)
		
		# Is it an error?
		if (class(BSM_results) == "try-error")
			{
			cat("\n\nBSM_fn=", BSM_fn, " is of class 'try-error', skipping plot.\n")
			next()
			}
		
		# Is it stratified?
		if (class(BSM_results) == "data.frame")
			{
			# Non-stratified. Extract cladogenetic events table, 
			# calculate anagenetic events table
			BSM_strat_TF = FALSE
			stochastic_mapping_results = BSM_results
			master_table_cladogenetic_events = stochastic_mapping_results

			# Calculate anagenetic events table
			events_table = events_txt_list_into_events_table(events_txt_list=master_table_cladogenetic_events$anagenetic_events_txt_below_node, trtable=master_table_cladogenetic_events, recalc_abs_ages=TRUE)
			events_table
			table_w_anagenetic_events = events_table
		
			} else {
			# Stratified. Extract cladogenetic and anagenetic events tables
			BSM_strat_TF = TRUE
			stochastic_mapping_results = BSM_results
			master_table_cladogenetic_events = stochastic_mapping_results$master_table_cladogenetic_events
			table_w_anagenetic_events = stochastic_mapping_results$table_w_anagenetic_events
			} # END if (class(BSM_results) == "data.frame")
	

	
	
		#######################################################
		# Do stochastic mapping conditional on the ML model parameters
		#######################################################
		cat("...plotting...")
		
		#######################################################
		# Plot the states and paint the branches
		#######################################################
		# Get colors_list_for_states
		# Setup 
		#include_null_range = TRUE
		#areanames = c("K", "O", "M", "H")
		#areas = areanames
		#max_range_size = 4
		#states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)

		# Get colors
		colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size)

		# Plot the tree and states at nodes/corners
		resmod = stochastic_map_states_into_res(res=res1, master_table_cladogenetic_events, stratified=BSM_strat_TF)

		scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
		plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt=analysis_titletxt, addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=FALSE, titlecex=0.6, tipcex=0.6, tr=tr, tipranges=tipranges)

		# Paint on the branch states
		paint_stochastic_map_branches(res=resmod, master_table_cladogenetic_events, colors_list_for_states=colors_list_for_states, lwd=5, lty=par("lty"), root.edge=TRUE, stratified=BSM_strat_TF)

		plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="", addl_params=list("j"), plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=TRUE, tr=tr, tipranges=tipranges)





		#######################################################
		# Stochastic mapping under e.g. DEC+J
		#######################################################
		
		if (is.null(BSM_fns_res2) == FALSE)
			{
			analysis_titletxt = paste("Extracting stochastic map #", mapnum, "/", num_maps, " under model ", res2name, groupname_txt, sep="")
	
			cat("\n")
			cat(analysis_titletxt)

			#######################################################
			# Load stochastic mapping result for res2
			#######################################################
			BSM_fn = BSM_fns_res2[mapnum]
			cat("...loading '", BSM_fn, "'...", sep="")
		
			# Loads to BSM_results
			load(BSM_fn)
		
			# Is it an error?
			if (class(BSM_results) == "try-error")
				{
				cat("\n\nBSM_fn=", BSM_fn, " is of class 'try-error', skipping plot.\n")
				plot(1,1,pch=".")
				title("'try-error' on this stochastic map...delete")
				}
		
			# Is it stratified?
			if (class(BSM_results) == "data.frame")
				{
				# Non-stratified. Extract cladogenetic events table, 
				# calculate anagenetic events table
				BSM_strat_TF = FALSE
				stochastic_mapping_results = BSM_results
				master_table_cladogenetic_events = stochastic_mapping_results

				# ***********
				# Calculate anagenetic events table
				# (why is this 0 under DEC?)
				events_table = events_txt_list_into_events_table(events_txt_list=master_table_cladogenetic_events$anagenetic_events_txt_below_node, trtable=master_table_cladogenetic_events, recalc_abs_ages=TRUE)
				events_table
				table_w_anagenetic_events = events_table
		
				} else {
				# Stratified. Extract cladogenetic and anagenetic events tables
				BSM_strat_TF = TRUE
				stochastic_mapping_results = BSM_results
				master_table_cladogenetic_events = stochastic_mapping_results$master_table_cladogenetic_events
				table_w_anagenetic_events = stochastic_mapping_results$table_w_anagenetic_events
				} # END if (class(BSM_results) == "data.frame")
	

	
	
			#######################################################
			# Do stochastic mapping conditional on the ML model parameters
			#######################################################
			cat("...plotting...")

			#######################################################
			# Plot the states and paint the branches
			#######################################################
			# Get colors_list_for_states
			# Setup 
			#include_null_range = TRUE
			#areanames = c("K", "O", "M", "H")
			#areas = areanames
			#max_range_size = 4
			#states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)

			# Get colors
			colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size)

			# Plot the tree and states at nodes/corners
			resmod = stochastic_map_states_into_res(res=res2, master_table_cladogenetic_events, stratified=BSM_strat_TF)

			scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
			plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt=analysis_titletxt, addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=FALSE, titlecex=0.6, tipcex=0.6, tr=tr, tipranges=tipranges)

			# Paint on the branch states
			paint_stochastic_map_branches(res=resmod, master_table_cladogenetic_events, colors_list_for_states=colors_list_for_states, lwd=5, lty=par("lty"), root.edge=TRUE, stratified=BSM_strat_TF)

			plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="", addl_params=list("j"), plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=TRUE, tr=tr, tipranges=tipranges)
			} # END if (is.null(BSM_fns_res2) == FALSE)

		} # END for (mapnum in 1:num_maps)



	dev.off()
	cmdstr = paste("open ", pdffn, sep="")
	system(cmdstr)

	return(pdffn)
	} # END plot_BSM_fns
