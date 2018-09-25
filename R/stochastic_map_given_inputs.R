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


stochastic_map_given_inputs <- function(stochastic_mapping_inputs, piecenum=NULL, maxtries=40000, seedval=as.numeric(Sys.time()), include_null_range, master_nodenum_toPrint=0)
	{
	#######################################################
	# running_stochastic_mapping, given a results object
	#######################################################
	defaults='
	res
	rootedge=FALSE
	statenum_bottom_root_branch_1based=NULL
	printlevel=1
	stratified=FALSE
	min_branchlength=0.000001
	
	stochastic_mapping_inputs=stochastic_mapping_inputs_list
	piecenum=NULL
	maxtries=40000
	seedval=12346	
	include_null_range=TRUE
	'
	
	# Set the seed
	#print(paste0("stochastic_map_given_inputs seedval=", seedval, sep=""))
	if (seedval > 2147483647)
		{
		seedval = seedval %% 2147483647
		}
	set.seed(seedval)
	
	# Necessary setting to avoid getting numbers etc. in the stochastic mapping output
	options(stringsAsFactors=FALSE)
		
	# Load inputs 
	rootedge = stochastic_mapping_inputs$rootedge
	statenum_bottom_root_branch_1based = stochastic_mapping_inputs$statenum_bottom_root_branch_1based
	printlevel = stochastic_mapping_inputs$printlevel
	stratified = stochastic_mapping_inputs$stratified
	
	if (is.null(stratified) == TRUE)
		{
		stop("ERROR in stochastic_map_given_inputs(); stochastic_mapping_inputs$stratified is NULL")
		}
	
	
	min_branchlength = stochastic_mapping_inputs$min_branchlength
	cluster_already_open = stochastic_mapping_inputs$cluster_already_open
	res = stochastic_mapping_inputs$res
	Qmat = stochastic_mapping_inputs$Qmat
	COO_weights_columnar = stochastic_mapping_inputs$COO_weights_columnar
	Rsp_rowsums = stochastic_mapping_inputs$Rsp_rowsums

	tipranges = stochastic_mapping_inputs$tipranges
	areas = stochastic_mapping_inputs$areas
	state_indices_0based = stochastic_mapping_inputs$state_indices_0based
	states_list = state_indices_0based
	ranges_list = stochastic_mapping_inputs$ranges_list
	numstates = stochastic_mapping_inputs$numstates

	if (include_null_range == TRUE)
		{
		numstates_during_cladogenesis = numstates - 1
		} else {
		numstates_during_cladogenesis = numstates - 0
		}
	
	
	
	if (stratified == FALSE)
		{
		trtable = stochastic_mapping_inputs$trtable
		tr = stochastic_mapping_inputs$tr
		phy2 = stochastic_mapping_inputs$phy2
		
		independent_likelihoods_on_each_branch = stochastic_mapping_inputs$independent_likelihoods_on_each_branch

		
		} else {
		subtable_rowsTF = stochastic_mapping_inputs$master_table_timeperiod_i$piecenum == piecenum
		subtable_rownums = (1:nrow(stochastic_mapping_inputs$master_table_timeperiod_i))[subtable_rowsTF]
		
 		# Convert master_table_timeperiod_i to trtable
		trtable = stochastic_mapping_inputs$master_table_timeperiod_i[subtable_rownums,]
		# trtable parts to use/change are ONLY the parts corresponding to this particular treepiece
		#trtable_rows_correct_pieces_TF = trtable_all_pieces$piecenum == piecenum
		#trtable = trtable_all_pieces[trtable_rows_correct_pieces_TF,]
		
		independent_likelihoods_on_each_branch = stochastic_mapping_inputs$independent_likelihoods_by_tree_piece_for_timeperiod_i[[piecenum]]
		
		phy = stochastic_mapping_inputs$tree_sections_list$return_pieces_list[[piecenum]]
		phy2 = reorder(phy, "pruningwise") # Do this, 
		} # END if (stratified == FALSE)


	# Basic tree info
	ntips = length(phy2$tip.label)
	num_internal_nodes = phy2$Nnode
	tipnums = 1:ntips
	root_nodenum = ntips+1
	nodenums = root_nodenum:(ntips+num_internal_nodes)
	
	
	


	if (stratified == FALSE)
		{
		# Add a column for the sampled node states
		sampled_states_AT_nodes = rep(NA, nrow(trtable))
		sampled_states_AT_brbots = rep(NA, nrow(trtable))

		# Add the right and left descendant node numbers
		leftright_nodes_matrix = get_leftright_nodes_matrix_from_results(phy2)

		left_desc_nodes = rep(NA, nrow(trtable))
		right_desc_nodes = rep(NA, nrow(trtable))

		# dcorner = descendant corner (i.e. right after speciation)
		samp_LEFT_dcorner = rep(NA, nrow(trtable))
		samp_RIGHT_dcorner = rep(NA, nrow(trtable))

		# Events on branches
		anagenetic_events_txt_below_node = rep(NA, nrow(trtable))

		trtable = cbind(trtable, sampled_states_AT_nodes, sampled_states_AT_brbots, left_desc_nodes, right_desc_nodes, samp_LEFT_dcorner, samp_RIGHT_dcorner, anagenetic_events_txt_below_node)
		
		# This works, the reverse does not -- 2015-04-06
		# (LR for tree plots is opposite for the internal structure)
		trtable$left_desc_nodes[nodenums] = leftright_nodes_matrix$right
		trtable$right_desc_nodes[nodenums] = leftright_nodes_matrix$left
		trtable[nodenums,]
		} else {
		# Stratified analysis: 
		# You've already got all these columns
		junk = 1
		} # END if (stratified == FALSE)


	#returned_mats1 = get_Qmat_COOmat_from_BioGeoBEARS_run_object(default_BioGeoBEARS_run_object)
	#returned_mats1

# 	returned_mats2 = get_Qmat_COOmat_from_res(res)
# 	returned_mats2
# 
# 	# Extract output
# 	Qmat = returned_mats2$Qmat
# 	COO_weights_columnar = returned_mats2$COO_weights_columnar
# 	Rsp_rowsums = returned_mats2$Rsp_rowsums

	# Calculate the likelihood P((left_state,right_state)|anc_state)
	# for each scenario (unconstrained)
	# Note:
	# COO_weights_columnar indices are 0-based, with no null_range
	# So, add 2 to get it so that e.g. state 0 = state 2 = Kauai
	#
	# Or, add 1 to get the 1based state indexes INSIDE COO_weights_columnar
	# 
	# COO_weights_columnar =
	# ancestral index, left index, right index, conditional
	# probability given ancestral states. (assuming likelihood
	# of descendants is 1)
	# Probabilities of each range-inheritance scenario, conditional
	# on ancestral state (without constraints on Left Branch state)
	#like_LeftRight_given_AncState = COO_weights_columnar[[4]] / (Rsp_rowsums[1+COO_weights_columnar[[1]]])
	#like_LeftRight_given_AncState


	# Calculate the total number of range-inheritance scenarios
	# under the model
	# (this is the number of scenarios with weight > 0)
	# (weight per event/ sum(weights) = prob. per event)
	num_scenarios = length(COO_weights_columnar[[1]])



	#######################################################
	#######################################################
	# THIS IS AN UPPASS FROM THE TIPS TO THE ROOT
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

	##########################################################
	# 1. Sample a state at the root of the tree or subtree (if it doesn't exist)
	# 2. Given the root state, calculate uppass probabilities to the corners
	# 3. Multiply the corner uppass probs by the downpass probs
	# 
	##########################################################
	
	
	# This is for sampling just the root node state
	# If there is no rootedge, and if the starting state
	# at the root node is not determined
	
	# If stratified analysis, no root edge, no pre-determined root state
	TF1 = ((stratified == TRUE) && (rootedge == FALSE) && (is.na(trtable$sampled_states_AT_nodes[root_nodenum])) && (trtable$node.type[root_nodenum] == "root"))

	# If non-stratified analysis, whether or not there is a root edge, no pre-determined root state
	#TF2 = ((stratified == FALSE) && (is.na(trtable$sampled_states_AT_nodes[root_nodenum])) && (trtable$node.type[root_nodenum] == "root"))
	TF2 = (stratified == FALSE)
	#print(TF2)
	if ( TF1 || TF2 )
		{
		# Sample a state at the root
		
		# The global root nodenum will be different than the subtree root nodenum
		if (stratified == TRUE)
			{
			global_rootTF = trtable$node.type == "root"
			# Error check
			if (sum(global_rootTF) != 1)
				{
				stop("\n\nStop ERROR: no global root node in subtree.\n\n")
				}
			global_root_nodenum = trtable$node[global_rootTF]
			
			} else {
			global_root_nodenum = root_nodenum
			} # END if (stratified == TRUE)
		
		# Get the stateprobs at the global root, and sample from them
		probs_branch_top = res$ML_marginal_prob_each_state_at_branch_top_AT_node[global_root_nodenum,]
		statenums = 1:numstates
		statenum_1based = sample(x=statenums, size=1, replace=TRUE, prob=probs_branch_top)
		statenum_1based
		#print(global_root_nodenum)
		#print(round(node_stateprobs,3))
		#print(statenum_1based)
		
		# Store the sampled state at the root
		trtable$sampled_states_AT_nodes[root_nodenum] = statenum_1based
		}

	
	# If stratified analysis, with root edge, no pre-determined root state
	# Simulate up the root branch, if that exists
	# (and if its a stratified analysis)
	#print(stratified)
	if ( (stratified == TRUE) && (rootedge == TRUE))
		{
		# Find the downpass conditional likelihoods (normalized) that have
		# been pre-calculated
		#res$condlikes
		#res$inputs$master_table
		subtable_rootTF = trtable$SUBnode.type == "root"
		subtable_rownum = (1:nrow(trtable))[subtable_rootTF]
		
		TF1 = res$inputs$master_table$node == trtable$node[subtable_rownum]
		TF2 = res$inputs$master_table$stratum == trtable$stratum[subtable_rownum]
		TF3 = res$inputs$master_table$piecenum == trtable$piecenum[subtable_rownum]
		TF4 = res$inputs$master_table$piececlass == "subtree"
		TF5 = res$inputs$master_table$SUBnode.type == "root"
		TF = (TF1 + TF2 + TF3 + TF4 + TF5) == 5
		rownums = 1:nrow(res$inputs$master_table)
		rownum = rownums[TF]
		rownum
		downpass_condlikes_at_branch_top = res$condlikes[rownum,]
		downpass_relprobs_at_branch_top = downpass_condlikes_at_branch_top / sum(downpass_condlikes_at_branch_top)
		
		# Now you just need to exponentiate up, given the previous-done 
		# independent likelihoods
		starting_state_1based = trtable$sampled_states_AT_brbots[root_nodenum]
		condprobs_branch_top = rep(0, times=numstates)
		condprobs_branch_bot = rep(0, times=numstates)
		condprobs_branch_bot[starting_state_1based] = 1

		# 2017-04-06_error check
		if (sum(condprobs_branch_bot) == 0)
			{
			txt = "STOP ERROR BB2: sum(condprobs_branch_bot) == 0.  Printing starting_state_1based of condprobs_branch_bot[starting_state_1based] = 1:"
			cat("\n\n")
			print(txt)
			print("print(starting_state_1based):")
			print(starting_state_1based)
			cat("\n\n")
			stop(txt)
			}

		
		# Exponentiate up (well, sorta, exponential pre-calculated)
		#condprobs_branch_top = condprobs_branch_bot %*% independent_likelihoods_by_tree_piece_for_timeperiod_i[[piecenum]]
		branch_length = trtable$SUBedge.length[root_nodenum]
		independent_likelihoods_on_root_branch_of_subtree = expokit_dgpadm_Qmat2(times=branch_length, Qmat=Qmat, transpose_needed=TRUE)
		
		condprobs_branch_top = condprobs_branch_bot %*% independent_likelihoods_on_root_branch_of_subtree
		
		if (include_null_range == TRUE)
			{
			condprobs_branch_top[1] = 0	# zero out the NULL range, since it is impossible in a survivor
			}
		
		# State probabilities at the top of the branch
		probs_branch_top = condprobs_branch_top * downpass_relprobs_at_branch_top
		probs_branch_top = probs_branch_top / sum(probs_branch_top)
		
		

		master_nodenum = res$inputs$master_table$node[rownum]
		if (master_nodenum == master_nodenum_toPrint)
			{
			print("stochastic_map_given_inputs():")
			print("rownum:")
			print(rownum)
			print("master_nodenum_toPrint:")
			print(master_nodenum_toPrint)
			print("condprobs_branch_bot:")
			print(round(condprobs_branch_bot, 3))
			print("downpass_relprobs_at_branch_top:")
			print(round(downpass_relprobs_at_branch_top, 3))
			print("condprobs_branch_top:")
			print(round(condprobs_branch_top, 3))
			print("probs_branch_top:")
			print(round(probs_branch_top, 3))
			}
		
		
		
		
		# Sample the state
		sampled_state_branch_top_1based = sample(x=1:numstates, size=1, replace=TRUE, prob=probs_branch_top)
		sampled_state_branch_top_1based

		if (is.na(sampled_state_branch_top_1based) == TRUE)
			{
			txt = paste0("STOP ERROR_line347 in stochastic_map_given_inputs(): sampled_state_branch_top_1based is NA.")
			cat("\n\n")
			cat(txt)
			cat("\n\n")

			print("probs_branch_top:")
			print(probs_branch_top)
			stop(txt)
			}


		
		# Store the state
		trtable$sampled_states_AT_nodes[subtable_rownum] = sampled_state_branch_top_1based


		#######################################################
		# Stochastic mapping, once the states at branch bottoms
		# and branch tops have been sampled
		# Specifically for root branches of sub-trees!!
		# (left this out the first time -- 2014-05-28_NJM)
		#######################################################
		
		# Stochastic mapping of events on the subtree root branch
		events_table_for_branch_below_subtree_root_node = stochastic_map_branch(nodenum_at_top_of_branch=subtable_rownum, trtable=trtable, Qmat=Qmat, state_indices_0based=state_indices_0based, ranges_list=ranges_list, areas=areas, stratified=stratified, maxtries=maxtries)

		# Store the text representation
		# (extract to table with events_txt_into_events_table() )
		subtree_root_branch_events_txt = events_table_into_txt(events_table_for_branch_below_subtree_root_node)
	
		trtable$anagenetic_events_txt_below_node[subtable_rownum] = subtree_root_branch_events_txt
		# End stochastic mapping on the branches below subtree root
		} # END if ( (stratified == TRUE) && (rootedge == TRUE))
	

	
	#######################################################
	# UPPASS THROUGH THE NODE FROM THE ROOT NODE
	# START FROM A NODE STATE, THEN SIMULATE THE TWO NODE STATES ABOVE,
	# THEN SIMULATE THE BRANCH EVENTS
	#######################################################
	
	# Visit edges in reverse order from the downpass
	edges_to_visit_uppass = seq(from=(num_internal_nodes*2), by=-2, length.out=num_internal_nodes)
	# Since we are going backwards
	#print(edges_to_visit_uppass)
	#print(i)
	#print(j)
	#cat("\n")

	#for (i in edges_to_visit_uppass)
	#j=edges_to_visit_uppass[1]
	cat("\nBeginning stochastic mapping UPPASS (i:Leftnode,j:Rightnode,anc:Ancnode;):\n", sep="")
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
		anc <- phy2$edge[i, 1]
		# Store the node number (starting with the root)
		nodenum = anc

		cat(i, ":", left_desc_nodenum, ",", j, ":", right_desc_nodenum, ",anc:", anc, "; ", sep="")

	
		# get the correct edges
		left_edge_TF = phy2$edge[,2] == left_desc_nodenum
		right_edge_TF = phy2$edge[,2] == right_desc_nodenum
	
		# Check the branchlength of each edge
		# It's a hook if either branch is super-short
		is_leftbranch_hook_TF = phy2$edge.length[left_edge_TF] < min_branchlength
		is_rightbranch_hook_TF = phy2$edge.length[right_edge_TF] < min_branchlength
		hooknode_TF = (is_leftbranch_hook_TF + is_rightbranch_hook_TF) > 0



		# Get the state at the current node (anc)
		statenum_1based = trtable$sampled_states_AT_nodes[nodenum] 
		
		if (is.na(statenum_1based) == TRUE)
			{
			txt = paste0("STOP ERROR_line426 in stochastic_map_given_inputs(): statenum_1based is NA.")
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			print("i:")
			print(i)
			print("j:")
			print(j)
			print("hooknode_TF:")
			print(hooknode_TF)
			print("nodenum:")
			print(nodenum)
			stop(txt)
			}
		
		
		#######################################################
		# STOCHASTIC MAPPING OF CLADOGENETIC PROCESS
		#######################################################

		# FIRST, get the node nums
		# 2017-04-07_bug fix: hooknodes weren't getting these
		# The downpass probs of each state at each branch
		if (stratified == FALSE)
			{
			# stratified == FALSE
			# *************** perhaps use 
			# left_desc_nodenum and right_desc_nodenum
			left_branch_decnode = trtable$left_desc_nodes[nodenum]
			right_branch_decnode = trtable$right_desc_nodes[nodenum]

			# *************** perhaps use 
			# left_desc_nodenum and right_desc_nodenum
			# NOPE, THIS IS FINE -- 2014-12-29_NJM
			if ( (left_desc_nodenum != left_branch_decnode) | (right_desc_nodenum != right_branch_decnode) )
				{
				print("left_desc_nodenum, right_desc_nodenum")
				cat(left_desc_nodenum, right_desc_nodenum)
				
				print("left_branch_decnode, right_branch_decnode")
				cat(left_branch_decnode, right_branch_decnode)
				
				stop()
				}


			left_branch_downpass_likes = res$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[left_branch_decnode, ]
			right_branch_downpass_likes = res$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[right_branch_decnode, ]
			} else {
			# stratified==TRUE
			daughter_nodenums_global = trtable$daughter_nds[nodenum][[1]]
			# names(leftright_nodes_matrix) = c("right", "left")
			left_branch_decnode_global = daughter_nodenums_global[1]
			right_branch_decnode_global = daughter_nodenums_global[2]
			left_branch_downpass_likes = res$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[left_branch_decnode_global, ]
			right_branch_downpass_likes = res$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[right_branch_decnode_global, ]
			
			# Get the LOCAL left/right nodes, for the 
			# sub-table
			left_branch_decnode_TF = trtable$node == left_branch_decnode_global
			left_branch_decnode = (1:nrow(trtable))[left_branch_decnode_TF]

			right_branch_decnode_TF = trtable$node == right_branch_decnode_global
			right_branch_decnode = (1:nrow(trtable))[right_branch_decnode_TF]
			
			left_branch_decnode
			right_branch_decnode
			} # END if (stratified == FALSE)




		# If it's a hooknode, then copy the node state up
		# (no sampling needed)
		if (hooknode_TF == TRUE)
			{
			# Just copy the node state to the corners
			#sampled_split_descendants = list()
			#sampled_split_descendants$left_decstate_1based = statenum_1based
			#sampled_split_descendants$right_decstate_1based = statenum_1based
			
			# 2017-04-07_bug_fix
			sample_uppass_res = list()
			sample_uppass_res$left_decstate_1based = statenum_1based
			sample_uppass_res$right_decstate_1based = statenum_1based
			
			if (is.na(statenum_1based) == TRUE)
				{
				print("print(sample_uppass_res)")
				print(sample_uppass_res)
				stop("STOP ERROR_line477")
				}
			#print(statenum_1based)
			#stop("STOP ERROR_line477a")
			
			# These will get copied to the table outside of this loop
			} else {
			# If NOT a hooknode (typical), sample a pair of descendant states	
			# Calculate the probability of each range inheritance scenario, 
			# given the chosen root state
			
			# -1 regardless of whether there is a null range (1 to 0 conversion)
			index_Qmat_0based_of_starting_state = statenum_1based - 1
	
			#RCOO_probs_list_given_ancestor = given_a_starting_state_get_prob_of_each_split_scenario(index_Qmat_0based_of_starting_state, COO_weights_columnar, numstates=numstates_during_cladogenesis, include_null_range=TRUE)
	
			#uppass_probs_of_scenarios_given_root_state = RCOO_probs_list_given_ancestor
			#uppass_probs_of_scenarios_given_root_state

			#cbind(COO_weights_columnar[[1]], COO_weights_columnar[[2]], COO_weights_columnar[[3]], uppass_probs_of_scenarios_given_root_state)

				
			#sampled_split_descendants = sample_split_scenario(COO_weights_columnar, uppass_probs_of_scenarios_given_root_state, left_branch_downpass_likes, right_branch_downpass_likes)
			
			# Input the ancestral state into the state probabilities at the current node
			statenum_1based
			probs_ancstate = rep(0, length(left_branch_downpass_likes))
			probs_ancstate[statenum_1based] = 1
			probs_ancstate
			
			# NJM -- for bug checking...
			if (nodenum == 21)
				{
				printflag = TRUE
				} else {
				printflag = FALSE
				}
			
			# OLD (2014-05-ish)
			#sampled_split_descendants = sample_split_scenario2(COO_weights_columnar, probs_ancstate, left_branch_downpass_likes, right_branch_downpass_likes, sample_which="both", return_prob_each_split_scenario=FALSE, include_null_range=TRUE, Rsp_rowsums=NULL, numstates_wo_null=NULL, printflag=printflag)
			
			# NEW (2015-01-10)
			sample_uppass_res = sample_uppass_split_scenario_given_probs_ancstate(probs_ancstate=probs_ancstate, COO_weights_columnar=COO_weights_columnar, numstates=numstates, include_null_range=include_null_range, left_branch_downpass_likes=left_branch_downpass_likes, right_branch_downpass_likes=right_branch_downpass_likes, Rsp_rowsums=NULL)

			} # END if (hooknode_TF == TRUE)
		
		

		
		# OLD 
		#left_decstate_1based = sampled_split_descendants$left_decstate_1based
		#right_decstate_1based = sampled_split_descendants$right_decstate_1based
		# NEW 
		left_decstate_1based = as.numeric(sample_uppass_res$left_decstate_1based)
		right_decstate_1based = as.numeric(sample_uppass_res$right_decstate_1based)
		
		#print(c("oldLeft", "oldRight", "newLeft", "newRight"))
		#print(c(sampled_split_descendants$left_decstate_1based, sampled_split_descendants$right_decstate_1based, sample_uppass_res$left_decstate_1based, sample_uppass_res$right_decstate_1based))
		
		# Put these into the trtable
		trtable$samp_LEFT_dcorner[nodenum] = left_decstate_1based
		trtable$samp_RIGHT_dcorner[nodenum] = right_decstate_1based
	
		# And put them in as the sampled states at branch bottoms for the appropriate
		# descendant nods
		# left_branch_decnode and right_branch_decnode are LOCAL
		# for the SUBTREE
		trtable$sampled_states_AT_brbots[left_branch_decnode] = left_decstate_1based
		trtable$sampled_states_AT_brbots[right_branch_decnode] = right_decstate_1based

		if (isblank_TF(left_decstate_1based) == TRUE)
			{
			txt = paste0("STOP ERROR_line576 in stochastic_map_given_inputs(): left_decstate_1based is BLANK: left_decstate_1based='", left_decstate_1based, "'.")
			stop(txt)		
			}


		if (isblank_TF(right_decstate_1based) == TRUE)
			{
			txt = paste0("STOP ERROR_line576 in stochastic_map_given_inputs(): left_decstate_1based is BLANK: right_decstate_1based='", right_decstate_1based, "'.")
			stop(txt)		
			}

	
	
		#######################################################
		# STOCHASTIC MAPPING OF ANAGENETIC PROCESS
		#######################################################
		# Now, evolution ALONG branches
		#independent_likelihoods_on_each_branch = calc_independent_likelihoods_on_each_branch(phy2, Qmat, cluster_already_open=NULL, Qmat_is_sparse=FALSE)
		# Steps:
		# a. Given a state at a corner, calculate the conditional probabilities
		#    of states at the branch top.
		# b. Multiply these by the saved downpass probabilities
		# c. Sample from this distribution, & store at the nodes at the top
	
		# Initialize the starting probabilities at branch bottoms
		# (setting the P(known sampled state) to equal 1!!)
		condprobs_Left_branch_top = rep(0, times=numstates)
		condprobs_Right_branch_top = rep(0, times=numstates)

		condprobs_Left_branch_bot = rep(0, times=numstates)
		condprobs_Right_branch_bot = rep(0, times=numstates)
		condprobs_Left_branch_bot[left_decstate_1based] = 1
		condprobs_Right_branch_bot[right_decstate_1based] = 1
	
		# Dense matrix exponentiation, which has been done already!
		TF2 = ( (length(cluster_already_open)==1) && (cluster_already_open==FALSE) )
		if (is.null(cluster_already_open) || (TF2))
			{
			# Relative probabilities of states at the top of left branch
			condprobs_Left_branch_top = try(condprobs_Left_branch_bot %*% independent_likelihoods_on_each_branch[,,i])
			
			if (class(condprobs_Left_branch_top) == "try-error")
				{
				print(i)
				print(length(condprobs_Left_branch_bot))
				print(dim(independent_likelihoods_on_each_branch[,,i]))
				#print(condprobs_Left_branch_top)

				save(condprobs_Left_branch_bot, file="condprobs_Left_branch_bot.Rdata")
				save(independent_likelihoods_on_each_branch, file="independent_likelihoods_on_each_branch.Rdata")
				save(independent_likelihoods_on_each_branch[,,i], file="independent_likelihoods_on_each_branch_i.Rdata")

				print("STOPPING on error in 'condprobs_Left_branch_bot %*% independent_likelihoods_on_each_branch[,,i]'")
				stop("STOPPING on error in 'condprobs_Left_branch_bot %*% independent_likelihoods_on_each_branch[,,i]'")
				}
					
			# Relative probabilities of states at the top of right branch
			condprobs_Right_branch_top = try(condprobs_Right_branch_bot %*% independent_likelihoods_on_each_branch[,,j])

			if (class(condprobs_Right_branch_top) == "try-error")
				{
				print(j)
				print(length(condprobs_Right_branch_bot))
				print(dim(independent_likelihoods_on_each_branch[,,j]))
				#print(condprobs_Right_branch_top)
				
				save(condprobs_Right_branch_bot, file="condprobs_Right_branch_bot.Rdata")
				save(independent_likelihoods_on_each_branch, file="independent_likelihoods_on_each_branch.Rdata")
				save(independent_likelihoods_on_each_branch[,,j], file="independent_likelihoods_on_each_branch_i.Rdata")
				
				print("STOPPING on error in 'condprobs_Right_branch_bot %*% independent_likelihoods_on_each_branch[,,j]'")
				stop("STOPPING on error in 'condprobs_Right_branch_bot %*% independent_likelihoods_on_each_branch[,,j]'")
				}
			} else {
		
			# Here, the independent_likelihoods_on_each_branch are stored in a list of matrices
			# Relative probabilities of states at the top of left branch
			condprobs_Left_branch_top = condprobs_Left_branch_bot %*% independent_likelihoods_on_each_branch[[i]]
					
			# Relative probabilities of states at the top of right branch
			condprobs_Right_branch_top = condprobs_Right_branch_bot %*% independent_likelihoods_on_each_branch[[j]]
			} # END if (is.null(cluster_already_open))

		# zero out the NULL range, since it is impossible in a survivor
		if (include_null_range == TRUE)
			{
			condprobs_Left_branch_top[1] = 0
			condprobs_Right_branch_top[1] = 0
			} # END if (include_null_range == TRUE)


		# Get the probabilities at the branch tops for the two branches above the node
		# under consideration
		# In non-stratified -- these are just the node numbers
		# In stratified analysis:
		# Here, left_desc_nodenum & right_desc_nodenum are for the SUBTREE
		if (stratified == FALSE)
			{
			# OK, now multiply the UPPASS and DOWNPASS probabilities
			probs_Left_branch_top = condprobs_Left_branch_top * res$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,]
			probs_Right_branch_top = condprobs_Right_branch_top * res$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum,]
			
			# In case users want to trace what's going on
			left_desc_nodenum_global = left_desc_nodenum
			left_desc_nodenum_global = right_desc_nodenum
			} else {
			# When stratified == TRUE, we have to dig up the corresponding rows of res
			# to get the correct downpass condlikes
			
			# Left node
			TF1 = res$inputs$master_table$node == trtable$node[left_desc_nodenum]
			TF2 = res$inputs$master_table$stratum == trtable$stratum[left_desc_nodenum]
			TF3 = res$inputs$master_table$piecenum == trtable$piecenum[left_desc_nodenum]
			TF4 = res$inputs$master_table$piececlass == "subtree"
			left_desc_nodenum_global_TF = (TF1 + TF2 + TF3 + TF4) == 4
			left_desc_nodenum_global = (1:nrow(res$condlikes))[left_desc_nodenum_global_TF]
			
			# Right node
			TF1 = res$inputs$master_table$node == trtable$node[right_desc_nodenum]
			TF2 = res$inputs$master_table$stratum == trtable$stratum[right_desc_nodenum]
			TF3 = res$inputs$master_table$piecenum == trtable$piecenum[right_desc_nodenum]
			TF4 = res$inputs$master_table$piececlass == "subtree"
			right_desc_nodenum_global_TF = (TF1 + TF2 + TF3 + TF4) == 4
			right_desc_nodenum_global = (1:nrow(res$condlikes))[right_desc_nodenum_global_TF]

			# OK, now multiply the UPPASS and DOWNPASS probabilities
			probs_Left_branch_top = condprobs_Left_branch_top * res$condlikes[left_desc_nodenum_global,]
			probs_Right_branch_top = condprobs_Right_branch_top * res$condlikes[right_desc_nodenum_global,]
			} # END if (stratified == FALSE)
	
		# Normalize by sum so they add to 1
		probs_Left_branch_top = probs_Left_branch_top / sum(probs_Left_branch_top)
		probs_Right_branch_top = probs_Right_branch_top / sum(probs_Right_branch_top)

# 		print("left_desc_nodenum_global:")
# 		print(left_desc_nodenum_global)
# 		print("right_desc_nodenum_global:")
# 		print(right_desc_nodenum_global)
		
		
		#print("Checking:")
		#print(probs_Left_branch_top)
		#print(probs_Right_branch_top)

		#for (zzz in 1:20)
		#	{
		# Sample states at the top of the two descendant branches
		# (these are conditional on the uppass probs, conditional on the 
		#  present node, AND on the saved downpass probs
		sampled_state_Left_branch_top_1based = sample(x=1:numstates, size=1, replace=TRUE, prob=probs_Left_branch_top)

		if (isblank_TF(sampled_state_Left_branch_top_1based) == TRUE)
			{
			txt = paste0("STOP ERROR_line715 in stochastic_map_given_inputs(): sampled_state_Left_branch_top_1based is NA.")
			cat("\n\n")
			cat(txt)
			cat("\n\n")

			print("probs_Left_branch_top:")
			print(probs_Left_branch_top)
			stop(txt)
			}

		sampled_state_Right_branch_top_1based = sample(x=1:numstates, size=1, replace=TRUE, prob=probs_Right_branch_top)
		if (isblank_TF(sampled_state_Right_branch_top_1based) == TRUE)
			{
			txt = paste0("STOP ERROR_line728 in stochastic_map_given_inputs(): sampled_state_Right_branch_top_1based is NA.")
			cat("\n\n")
			cat(txt)
			cat("\n\n")

			print("probs_Right_branch_top:")
			print(probs_Right_branch_top)
			stop(txt)
			}


	
		# Store these states
		trtable$sampled_states_AT_nodes[left_branch_decnode] = sampled_state_Left_branch_top_1based

		trtable$sampled_states_AT_nodes[right_branch_decnode] = sampled_state_Right_branch_top_1based



		if (stratified == FALSE)
			{
			master_nodenum = left_desc_nodenum
			} else {
			master_nodenum = res$inputs$master_table$node[left_desc_nodenum_global]
			} # END if (stratified == FALSE)
		if (master_nodenum == master_nodenum_toPrint)
			{
			cat("\n")
			print("left_branch_decnode:")
			print(left_branch_decnode)
			print("condprobs_Left_branch_bot:")
			print(round(condprobs_Left_branch_bot,3))
			print("condprobs_Left_branch_top:")
			print(round(condprobs_Left_branch_top,3))
			print("res$condlikes[left_desc_nodenum_global,]:")
			print(round(res$condlikes[left_desc_nodenum_global,],3))
			print("probs_Left_branch_top:")
			print(round(probs_Left_branch_top,3))
			print("sampled_state_Left_branch_top_1based:")
			print(round(sampled_state_Left_branch_top_1based,3))
			sampled_stateprobs = rep(0, numstates)
			sampled_stateprobs[sampled_state_Left_branch_top_1based] = 1
			print("sampled_stateprobs:")
			print(sampled_stateprobs)
			cat("\n")
			}


		if (stratified == FALSE)
			{
			master_nodenum = right_desc_nodenum
			} else {
			master_nodenum = res$inputs$master_table$node[right_desc_nodenum_global]
			} # END if (stratified == FALSE)
		if (master_nodenum == master_nodenum_toPrint)
			{
			cat("\n")
			print("right_branch_decnode:")
			print(right_branch_decnode)
			print("condprobs_Right_branch_bot:")
			print(round(condprobs_Right_branch_bot,3))
			print("condprobs_Right_branch_top:")
			print(round(condprobs_Right_branch_top,3))
			print("res$condlikes[right_desc_nodenum_global,]:")
			print(round(res$condlikes[right_desc_nodenum_global,],3))
			print("probs_Right_branch_top:")
			print(round(probs_Right_branch_top,3))
			print("sampled_state_Right_branch_top_1based:")
			print(round(sampled_state_Right_branch_top_1based,3))
			sampled_stateprobs = rep(0, numstates)
			sampled_stateprobs[sampled_state_Right_branch_top_1based] = 1
			print("sampled_stateprobs:")
			print(sampled_stateprobs)
			cat("\n")
			}
		




		
		#txt = paste0("zzz: ", zzz, ", L: ", left_decstate_1based, "->", sampled_state_Left_branch_top_1based, "; R: ", right_decstate_1based, "->", sampled_state_Right_branch_top_1based)
		#cat("\n")
		#cat(txt)
		#} # END zzz
		#cat("\n")
	
		# CHECK LEFT-RIGHT NODE NUMBERING
		check_left_vs_right_numbering = FALSE	
		# NOTE: FOR TREE ITERATIONS, YOU HAVE TO SWITCH
		# LEFT AND RIGHT, WHICH IS DIFFERENT THAN
		# FOR GRAPHING:

		# # Add the right and left descendant node numbers
		# leftright_nodes_matrix = get_leftright_nodes_matrix_from_results(phy2)
		# 
		# left_desc_nodes = rep(NA, nrow(trtable))
		# right_desc_nodes = rep(NA, nrow(trtable))
		# 
		# # dcorner = descendant corner (i.e. right after speciation)
		# samp_LEFT_dcorner = rep(NA, nrow(trtable))
		# samp_RIGHT_dcorner = rep(NA, nrow(trtable))
		# 
		# trtable = cbind(trtable, left_desc_nodes, right_desc_nodes, samp_LEFT_dcorner, samp_RIGHT_dcorner)
		# trtable$left_desc_nodes[nodenums] = leftright_nodes_matrix$right
		# trtable$right_desc_nodes[nodenums] = leftright_nodes_matrix$left
		# trtable[nodenums,]

		if (check_left_vs_right_numbering == TRUE)
			{
			if (left_branch_decnode <= ntips)
				{
				print("Left tip:")
				print(left_branch_decnode)
				print(left_desc_nodenum)
				print(sampled_state_Left_branch_top_1based)
				print(probs_Left_branch_top)
				} # END if (left_branch_decnode <= ntips)

			if (right_branch_decnode <= ntips)
				{
				print("Right tip:")
				print(right_branch_decnode)
				print(right_desc_nodenum)
				print(sampled_state_Right_branch_top_1based)
				print(probs_Right_branch_top)
				} # END if (right_branch_decnode <= ntips)
			} # END if (check_left_vs_right_numbering == TRUE)
	
	
		#######################################################
		# Stochastic mapping, once the states at branch bottoms
		# and branch tops have been sampled
		#######################################################
		
		# Stochastic mapping of events on the left branch
		events_table_for_branch_below_Left_node = stochastic_map_branch(nodenum_at_top_of_branch=left_branch_decnode, trtable=trtable, Qmat=Qmat, state_indices_0based=state_indices_0based, ranges_list=ranges_list, areas=areas, stratified=stratified, maxtries=maxtries)

		# Stochastic mapping of events on the right branch
		events_table_for_branch_below_Right_node = stochastic_map_branch(nodenum_at_top_of_branch=right_branch_decnode, trtable=trtable, Qmat=Qmat, state_indices_0based=state_indices_0based, ranges_list=ranges_list, areas=areas, stratified=stratified, maxtries=maxtries)
	
		# Store the text representation
		# (extract to table with events_txt_into_events_table() )
		left_branch_events_txt = events_table_into_txt(events_table_for_branch_below_Left_node)
	
		right_branch_events_txt = events_table_into_txt(events_table_for_branch_below_Right_node)
	
		trtable$anagenetic_events_txt_below_node[left_branch_decnode] = left_branch_events_txt
		trtable$anagenetic_events_txt_below_node[right_branch_decnode] = right_branch_events_txt
		
# 		print(left_branch_decnode)
# 		print(left_branch_events_txt)
# 		print(right_branch_decnode)
# 		print(right_branch_events_txt)
	
		} # END for (j in edges_to_visit_uppass)
		  # (ENDING uppass loop)
	
	# Some of these are getting skipped??
	#stochastic_mapping_inputs_list[[11]]$master_table_timeperiod_i$sampled_states_AT_nodes
	#master_table_timeperiod_i$sampled_states_AT_nodes
	
	if (any(is.na(trtable$sampled_states_AT_nodes)))
		{
		txt = paste0("STOP ERROR_line923 in stochastic_map_given_inputs(): These trtable$sampled_states_AT_nodes are still NA!")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		print("Rows that are still NA:")
		rows_still_NA = (1:length(trtable$sampled_states_AT_nodes))[is.na(trtable$sampled_states_AT_nodes)]
		print( rows_still_NA )
		
		print("print(trtable[is.na(trtable$sampled_states_AT_nodes),])")
		print(trtable[is.na(trtable$sampled_states_AT_nodes),])
		
		print(trtable[trtable$SUBedge.length<min_branchlength,])
		
		stop(txt)
		}
	
	
	
	cat("\n...finished stochastic mapping UPPASS.\n", sep="")
	
#	trtable$anagenetic_events_txt_below_node
#	trtable = add_cladogenetic_events_to_trtable(trtable=trtable_good, BioGeoBEARS_run_object=res$inputs, tipranges, stratified=stratified, piecenum=NULL)
#	trtable$anagenetic_events_txt_below_node
	#print(trtable$anagenetic_events_txt_below_node)
	
	cat("Adding cladogenetic events and re-arranging columns...\n", sep="")
	# Add cladogenetic events and re-arrange columns
	if (stratified == FALSE)
		{
		trtable = add_cladogenetic_events_to_trtable(trtable, BioGeoBEARS_run_object=res$inputs, tipranges, stratified=stratified, piecenum=NULL)
		} else {
		trtable = add_cladogenetic_events_to_trtable(trtable, BioGeoBEARS_run_object=res$inputs, tipranges, stratified=stratified, piecenum=piecenum)
		}
	cat("DONE adding cladogenetic events and re-arranging columns.\n", sep="")

	#print(trtable$anagenetic_events_txt_below_node)

	cat("FINISHING stochastic_map_given_inputs().\n", sep="")
		
	# If stratified, don't rearrange
	# If not stratified, *do* rearrange
	if (stratified == FALSE)
		{
		first_colnums = 1:(ncol(trtable)-4)
		last4_colnums = (ncol(trtable)-3):(ncol(trtable))
		new_colnums = c(first_colnums, last4_colnums[c(2,3,4,1)])
		trtable = trtable[,new_colnums]
		}

	# Convert master_table_timeperiod_i to trtable
	if (stratified == TRUE)
		{
		#cat("\n\nprint(dim(trtable)):\n\n")
		#print(dim(trtable))

		#cat("\n\nprint(master_table_timeperiod_i[subtable_rownums,]):\n\n")
		#print(master_table_timeperiod_i[subtable_rownums,])

		
		# We may be modifying just PART of the subtable
		stochastic_mapping_inputs$master_table_timeperiod_i[subtable_rownums,] = trtable
		# Look at results
		stochastic_mapping_inputs$master_table_timeperiod_i
		
		return(stochastic_mapping_inputs$master_table_timeperiod_i)
		} else {
		# Look at results
		trtable[nodenums,]
		return(trtable)
		}
	
	return(stop("ERROR: you should not reach this."))
	} # END stochastic_map_given_inputs <- function(stochastic_mapping_inputs, piecenum=NULL, maxtries=40000, seedval=12345, include_null_range)















# 
# For stochastic mapping on non-stratified analysis:
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
# stochastic_mapping_on_stratified()
#######################################################
# 
# Do stochastic mapping on the results of a stratified analysis, given
# (1) a results object from bears_optim_run
# (2) a list of stochastic mapping inputs generated
#     as follows:
# # Get stochastic mapping inputs for stratified analysis
# stochastic_mapping_inputs_list = list()
# for (timeperiod_i in 1:length(res$inputs$tree_sections_list))
# 	{
# 	stochastic_mapping_inputs = get_inputs_for_stochastic_mapping_from_results_object(res=res, stratified=stratified, timeperiod_i=timeperiod_i)
# 	
# 	stochastic_mapping_inputs_list[[timeperiod_i]] = stochastic_mapping_inputs
# 	} # END for (timeperiod_i in 1:length(res$inputs$tree_sections_list))
# 
#######################################################
stochastic_mapping_on_stratified <- function(res, stochastic_mapping_inputs_list, maxtries=40000, seedval=as.numeric(Sys.time()), master_nodenum_toPrint=0)
	{
	# Necessary setting to avoid getting numbers etc. in the stochastic mapping output
	options(stringsAsFactors=FALSE)

	# Set the seed
	#print(paste0("stochastic_mapping_on_stratified seedval=", seedval, sep=""))
	if (seedval > 2147483647)
		{
		seedval = seedval %% 2147483647
		}
	set.seed(seedval)
	
	# Seed increment counter for tree pieces
	treepiece_seed_counter_increment = 90127 # (large prime number)
	treepiece_seed_counter = 1
	
	
	# Set stratified=TRUE, obviously -- put in error catch for non-stratified analysis.
	if (is.na(res$inputs$timesfn) == TRUE)
		{
		stratified=FALSE
		errortxt = "STOP ERROR in stochastic_mapping_on_stratified().\n\nThe 'res' object that you input to this function has no input times filename.\n\nSpecifically, res$inputs$timesfn equals NA.  This suggests that the 'res' input\n\nwas generated by a non-stratified analysis.  In that case, you should use\n\nstochastic_map_given_inputs() instead of stochastic_mapping_on_stratified() for \n\nBiogeographical Stochastic Mapping (BSM) in BioGeoBEARS\n\n."
		cat(errortxt)
		stop(errortxt)		
		} else {
		stratified=TRUE
		}
	
	# Get results from the likelihood analysis
	# (for input to:
	#  get_tip_likelihoods_of_subbranch_from_resCondlikes_given_master_table )
	master_table = res$inputs$master_table
	condlikes = res$condlikes_table
	relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE = res$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE
	
	# Error check
	if (is.null(relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE))
		{
		errortxt = "STOP ERROR in stochastic_mapping_on_stratified():\n\nres$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE is NULL.\n\nThis will occur on results_objects (res) that were made by bears_optim_run()\n before about May 2014. Please re-run the ML inference with updated package/sourcefiles before\nrunning Biogeographic Stochastic Mapping.\n\nOr, perhaps you are accidentally using stochastic_mapping_on_stratified() when your input ML analysis was non-stratified. In that case, you should use stochastic_map_given_inputs() instead.\n\nHave a nice day.\n\n"
		cat(errortxt)
		stop(errortxt)
		} # END if (is.null(relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE))
	
	
	# Get the root state
	num_timeperiods = length(res$inputs$timeperiods)
	stochastic_mapping_inputs = stochastic_mapping_inputs_list[[num_timeperiods]]
	subtable = stochastic_mapping_inputs$master_table_timeperiod_i
	
	
	
	rootTF = subtable$node.type == "root"
	if (sum(rootTF) != 1)
		{
		errortxt = paste("\n\nERROR: Your bottom tree section, timeperiod_i=", num_timeperiods, ", does not contain a node of node.type 'root'\nin $master_table_timeperiod_i:\n\n", sep="")
		cat(errortxt)
		} else {
		rootnode = subtable$node[rootTF]
		}
	rootnode

	numstates = ncol(res$ML_marginal_prob_each_state_at_branch_top_AT_node)
	node_stateprobs = res$ML_marginal_prob_each_state_at_branch_top_AT_node[rootnode,]
	statenums = 1:numstates
	statenum_1based = sample(x=statenums, size=1, replace=TRUE, prob=node_stateprobs)
	statenum_1based

	if (is.na(statenum_1based) == TRUE)
		{
		txt = paste0("STOP ERROR_line1013 in stochastic_map_given_inputs(): statenum_1based is NA.")
		cat("\n\n")
		cat(txt)
		cat("\n\n")

		print("node_stateprobs:")
		print(node_stateprobs)
		stop(txt)
		}
	


	# Store the root state in the master table
	stochastic_mapping_inputs_list[[num_timeperiods]]$master_table_timeperiod_i$sampled_states_AT_nodes[rootTF] = statenum_1based
	# The treepiece in the BOTTOM stratum should ALWAYS be a subtree and 
	# NEVER a root branch

	# Now, you have to walk up the time pieces (including root branches)
	for (timeperiod_i in num_timeperiods:1)
		{
		stratum = timeperiod_i
		timeperiod_i_up = timeperiod_i - 1
		
		stochastic_mapping_inputs = stochastic_mapping_inputs_list[[timeperiod_i]]
		master_table_timeperiod_i = stochastic_mapping_inputs$master_table_timeperiod_i


# 		if (any(is.na(master_table_timeperiod_i$sampled_states_AT_nodes)))
# 			{
# 			txt = paste0("STOP ERROR_line1085 in (): Some of master_table_timeperiod_i$sampled_states_AT_nodes are NA. Specifically, these:")
# 			cat("\n\n")
# 			cat(txt)
# 			print("print( (1:length(master_table_timeperiod_i$sampled_states_AT_nodes))[is.na(master_table_timeperiod_i$sampled_states_AT_nodes)] )")
# 		
# 			print( (1:length(master_table_timeperiod_i$sampled_states_AT_nodes))[is.na(master_table_timeperiod_i$sampled_states_AT_nodes)] )
# 			print("timeperiod_i")
# 			print(timeperiod_i)
# 
# 			print("piecenum")
# 			print(piecenum)
# 		
# 			cat("\n\n")
# 			stop(txt)
# 			}



		tree_sections = stochastic_mapping_inputs$tree_sections_list
	
		tree_pieces = tree_sections$return_pieces_list
		num_pieces = length(tree_pieces)
	
		# Get the independent likelihoods
		independent_likelihoods_by_tree_piece_for_timeperiod_i = stochastic_mapping_inputs$independent_likelihoods_by_tree_piece_for_timeperiod_i
	
		# Get the transition rate parameters
		Qmat = stochastic_mapping_inputs$Qmat
		COO_weights_columnar = stochastic_mapping_inputs$COO_weights_columnar
		Rsp_rowsums = stochastic_mapping_inputs$Rsp_rowsums
		state_indices_0based = stochastic_mapping_inputs$state_indices_0based
		ranges_list = stochastic_mapping_inputs$ranges_list
		areas = stochastic_mapping_inputs$areas
	
	  

		for (piecenum in 1:num_pieces)
			{
			# 2016-05-07_bugfix
# 			seedval = 24785446
# 			timeperiod_i =1
# 			piecenum = 10
			#print("Apiecenum:")
			
			tree_piece = tree_pieces[[piecenum]]

			error_check_Psychotria_all_tips_size1 = FALSE
			if (error_check_Psychotria_all_tips_size1)
				{
				cat("\n\n\n\n")
	
	
				print("timeperiod_i")
				print(timeperiod_i)
	
				print("stratum")
				print(stratum)

				print("piecenum")
				print(piecenum)
				} # END if (error_check_Psychotria_all_tips_size1)
		
			if (is.numeric(tree_piece))
				{
				#print("Apiecenum:")
				#print(tree_piece)
				
				# Piece is a single branch; stochastically map on just that branch
				single_branch = TRUE
			
				#######################################################
				# NOTE: For SINGLE BRANCHES, let's use:
				# sampled_states_AT_brbots = state at the bottom of the branch
				# 							 (previously sampled)
				# sampled_states_AT_nodes  = state at the top of the branch
				# 							 (sampled at the end of this step)
				#######################################################
			
				# Simulate up from the (pre-determined) bottom of the branch
				#master_table_timeperiod_i$sampled_states_AT_brbots[subtable_rownum] = 2
			
				# Get the pre-determined starting state
				rowTF = master_table_timeperiod_i$piecenum == piecenum
				subtable_rownum = (1:nrow(master_table_timeperiod_i))[rowTF]
				starting_state_1based = master_table_timeperiod_i$sampled_states_AT_brbots[subtable_rownum]
				starting_state_1based
				
				#print("print(timeperiod_i):")
				#print(timeperiod_i)

				#print("print(subtable_rownum):")
				#print(subtable_rownum)
				
				#print("print(master_table_timeperiod_i[subtable_rownum, ]):")
				#print(master_table_timeperiod_i[subtable_rownum, ])
				
				if (isblank_TF(starting_state_1based) == TRUE)
					{
					stop("STOP_line_1086")
					}
				
				# Find the downpass conditional likelihoods (normalized) that have
				# been pre-calculated
				# 
				# Inputs:
				# res$condlikes -- contains downpass branch-top likelihoods for
				# nodes at the tops of branches in subtrees, and for the 
				# bottoms of subbranches 

				# res$condlikes -- contains downpass branch-top likelihoods for
				# nodes at the tops of branches in subtrees, and for the 
				# bottoms of subbranches 
				
				# Contains all of the tree pieces and their references to the master tree
				# res$inputs$master_table
				
	
	
				# START CHECK IF THERE IS ONLY ONE ROW (redundant)
				TF1 = res$inputs$master_table$node == master_table_timeperiod_i$node[subtable_rownum]
				TF2 = res$inputs$master_table$stratum == master_table_timeperiod_i$stratum[subtable_rownum]
				TF3 = res$inputs$master_table$piecenum == master_table_timeperiod_i$piecenum[subtable_rownum]
				TF4 = res$inputs$master_table$piececlass != "subtree"
				TF = (TF1 + TF2 + TF3 + TF4) == 4
				rownums = 1:nrow(res$inputs$master_table)
				rownum_master_table = rownums[TF]
				rownum = rownum_master_table
				
				
				#master_nodenum = res$inputs$master_table$node[rownum_master_table]
				#print(master_nodenum)
				
				if (length(rownum) != 1)
					{
					errortxt = paste("\n\nERROR in identifying corresponding subbranch in res$inputs$master_table.\nlength(rownum) should be 1, but is actually: ", length(rownum), "\n\n", sep="")
					cat(errortxt)
					
					cat("\n\n")
					cat("timeperiod_i=", timeperiod_i, sep="")
					cat("\n")
					
					cat("piecenum=", piecenum, sep="")
					cat("\n\n")
					
					
					cat("master_table_timeperiod_i:")
					cat("\n\n")
					printall(master_table_timeperiod_i)
					
					stop("\n\nStopping on error.")
					} # END check if there is only 1 row
	
				# Inputs
				condlikes = res$condlikes
				relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE = res$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE
				
				
				# Get the tip likelihoods, including subbranches in the
				# top stratum (#1) corresponding to master_tree tips, and
				# fossil tips which will also be in the master_tree_tips
				#print("print(AC_rowsums_condlikes):")
				#print(sum(rowSums(condlikes) == 0))

				#print("print(AC_relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE):")
				#print(sum(rowSums(relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE) == 0))
				
				downpass_condlikes_at_branch_top = get_tip_likelihoods_of_subbranch_from_resCondlikes_given_master_table(stratum=stratum, piecenum=piecenum, master_table=master_table, condlikes=condlikes, relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE=relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE)
				# really these are downpass relprobs, but whatever
				
				# Normalize just to make sure
				#print("print(AC_downpass_condlikes_at_branch_top):")
				#print(downpass_condlikes_at_branch_top)

				downpass_relprobs_at_branch_top = downpass_condlikes_at_branch_top / sum(downpass_condlikes_at_branch_top)

# 				if (master_nodenum == 36)
# 					{
# 					print(downpass_relprobs_at_branch_top)
# 					}
			
				# Now you just need to exponentiate up, given the previous-done 
				# independent likelihoods
				condprobs_branch_top = rep(0, times=numstates)
				condprobs_branch_bot = rep(0, times=numstates)
				condprobs_branch_bot[starting_state_1based] = 1
				
				
				# 2017-04-06_error check
				if (sum(condprobs_branch_bot) == 0)
					{
					txt = "STOP ERROR BB1: sum(condprobs_branch_bot) == 0.  Printing starting_state_1based of condprobs_branch_bot[starting_state_1based] = 1:"
					cat("\n\n")
					print(txt)
					print("print(starting_state_1based):")
					print(starting_state_1based)
					cat("\n\n")
					stop(txt)
					}
				
				# Exponentiate up (well, sorta, exponential pre-calculated)
				condprobs_branch_top = condprobs_branch_bot %*% independent_likelihoods_by_tree_piece_for_timeperiod_i[[piecenum]]
				
				if (res$inputs$include_null_range == TRUE)
					{
					condprobs_branch_top[1] = 0	# zero out the NULL range, since it is impossible in a survivor
					} # END if (res$inputs$include_null_range == TRUE)
			
				# State probabilities at the top of the branch
				#print("print(AD_condprobs_branch_top):")
				#print(condprobs_branch_top)
				#print("print(AD_downpass_relprobs_at_branch_top):")
				#print(downpass_relprobs_at_branch_top)

				probs_branch_top = condprobs_branch_top * downpass_relprobs_at_branch_top
				probs_branch_top = probs_branch_top / sum(probs_branch_top)
			
				# Sample the state
				treepiece_seed_counter = treepiece_seed_counter + 1
				seedval = seedval + (treepiece_seed_counter * treepiece_seed_counter_increment)
				if (seedval > 2147483647)
					{
					seedval = seedval %% 2147483647
					}				
				set.seed(seed=seedval)
				
				#print("print(AD_probs_branch_top):")
				#print(probs_branch_top)
				
				sampled_state_branch_top_1based = sample(x=1:numstates, size=1, replace=TRUE, prob=probs_branch_top)
				sampled_state_branch_top_1based
				
				# ERROR CHECK
				print_big_error_check = FALSE
				if (isblank_TF(sampled_state_branch_top_1based) == TRUE)
					{
					print_big_error_check = TRUE
					}
				
				
				if (print_big_error_check == TRUE)
					{
					txt = paste0("\n\nSTOP ERROR123 AFTER 'sampled_state_branch_top_1based = sample(x=1:numstates, size=1, replace=TRUE, prob=probs_branch_top)'\n\n")
					cat(txt)
				
					print("timeperiod_i")
					print(timeperiod_i)
				
					print("stratum")
					print(stratum)

					print("piecenum")
					print(piecenum)

					print("rownum")
					print(rownum)
				
				
					print("master_table_timeperiod_i")
					printall(master_table_timeperiod_i)


					print("master_table_timeperiod_i[subtable_rownum,]")
					print(master_table_timeperiod_i[subtable_rownum,])


					print("subtable_rownum")
					print(subtable_rownum)

					print("dim(master_table_timeperiod_i)")
					print(dim(master_table_timeperiod_i))


				
					print("rownum")
					print(rownum)

					print("condprobs_branch_top")
					print(condprobs_branch_top)


					print("downpass_relprobs_at_branch_top")
					print(downpass_relprobs_at_branch_top)

					print("probs_branch_top")
					print(probs_branch_top)

				
					print("sampled_state_branch_top_1based")
					print(sampled_state_branch_top_1based)
					stop(txt)
					} # END if (print_big_error_check) == TRUE
				
				
				# Store the state
				master_table_timeperiod_i$sampled_states_AT_nodes[subtable_rownum] = sampled_state_branch_top_1based
			
				# Also, put the state at the branch bottom of the next stratum up
				# (as long as you're not in the top stratum)
				if (timeperiod_i != 1)
					{
					# 1 stratum up:
					timeperiod_i_up = timeperiod_i - 1
					stochastic_mapping_inputs_up = stochastic_mapping_inputs_list[[timeperiod_i_up]]
					master_table_timeperiod_i_up = stochastic_mapping_inputs_up$master_table_timeperiod_i
					tree_sections_up = stochastic_mapping_inputs_up$tree_sections_list
				
					# Current stratum
					node_at_top_of_branch = master_table_timeperiod_i$node[subtable_rownum]
				
					# Find the node at the branch top (which could be in the stratum above, 
					# or even higher, no worries, they all have the same node at the top)
					node_above_TF = master_table_timeperiod_i_up$node == node_at_top_of_branch
					rownum_above = (1:nrow(master_table_timeperiod_i_up))[node_above_TF]
					rownum_above
				
					# This is the row for the node at the top of the branch
					# (or equivalent branch section in the stratum above)
					# master_table_timeperiod_i_up[rownum_above,]
				
					# Store the (now known) ancestral state at the bottom of the branch
					master_table_timeperiod_i_up$sampled_states_AT_brbots[rownum_above] = sampled_state_branch_top_1based
				
					# Put these back into the master tables
				
					# Store the node in the stratum above
					stochastic_mapping_inputs_list[[timeperiod_i_up]]$master_table_timeperiod_i = master_table_timeperiod_i_up

					} # END if (timeperiod_i != 1)
			
			
				# Simulate a history along the branch, given the starting and
				# ending states
				# Stochastic map on a sub-branch
				treepiece_seed_counter = treepiece_seed_counter + 1
				seedval = seedval + (treepiece_seed_counter * treepiece_seed_counter_increment)
				# Correction for seeds over the integer max
				if (seedval > 2147483647)
					{
					seedval = seedval %% 2147483647
					}
				set.seed(seed=seedval)
				
				
				#######################################################
				# 2016-05-07_bugfix
				#######################################################
				#print(seedval)
				#print(timeperiod_i)
				#print(piecenum)
# 				seedval = 27579383
# 				timeperiod_i = 1
# 				piecenum = 10
# 				print(seedval)
# 				print(timeperiod_i)
# 				print(piecenum)
# 
# 				nodenum_at_top_of_branch=subtable_rownum; 
# 				trtable=master_table_timeperiod_i; 
# 				Qmat; 
# 				state_indices_0based; 
# 				ranges_list; 
# 				areas; 
# 				single_branch=single_branch; 
# 				stratified=stratified; 
# 				maxtries=maxtries
# 				manual_history_for_difficult_branches=TRUE
				
				# Check for branches where the tips are BELOW the time-bin!
				branch_top_full_tree = master_table_timeperiod_i$time_bp[subtable_rownum]
				branch_bot_full_tree = master_table_timeperiod_i$time_bp[subtable_rownum] - master_table_timeperiod_i$edge.length[subtable_rownum]
				branch_length_full_tree = master_table_timeperiod_i$edge.length[subtable_rownum]
				# The SUBedge.length produced by tree sectioning; MAY INCLUDE FAKE BRANCH LENGTH!!
				branch_length_subsection = master_table_timeperiod_i$SUBedge.length[subtable_rownum]

				#cat("\nABC1_stochastic_map_given_inputs():\n")
				
				if (branch_length_subsection > branch_length_full_tree)
					{
					events_table_for_branch = NULL
					} else {
					events_table_for_branch = stochastic_map_branch(nodenum_at_top_of_branch=subtable_rownum, trtable=master_table_timeperiod_i, Qmat, state_indices_0based, ranges_list, areas, single_branch=single_branch, stratified=stratified, maxtries=maxtries, manual_history_for_difficult_branches=TRUE)
					} # END if (branch_length_subsection > branch_length_full_tree)
				#cat("\nABC2_stochastic_map_given_inputs():\n")
			
				# Convert the events to text
				branch_events_txt = events_table_into_txt(events_table_for_branch)
				branch_events_txt
				#cat("\nABC3_stochastic_map_given_inputs():\n")

				# Store the node history in the current stratum
				master_table_timeperiod_i$anagenetic_events_txt_below_node[subtable_rownum] = branch_events_txt
				# Store the updated subtable
				stochastic_mapping_inputs$master_table_timeperiod_i = master_table_timeperiod_i
				# Make sure it's in the global data structure
				stochastic_mapping_inputs_list[[timeperiod_i]]$master_table_timeperiod_i = master_table_timeperiod_i
				#cat("\nABC4_stochastic_map_given_inputs():\n")
				#print("print(piecenum):")
				#print(piecenum)
				#print("print(timeperiod_i):")
				#print(timeperiod_i)
				
				# End stochastic mapping on a branch
				} else {
				#######################################################
				#######################################################
				# Stochastically map on sub-tree (not a single branch)
				#######################################################
				#######################################################
			
				# Piece is a subtree; stochastically map on that subtree, 
				# STARTING FROM THE ROOT BRANCH
				# timeperiod_i = 5
				# piecenum = 1
				# stochastic_mapping_inputs = stochastic_mapping_inputs_list[[timeperiod_i]]
				# master_table_timeperiod_i = stochastic_mapping_inputs$master_table_timeperiod_i
				# tree_sections_list = stochastic_mapping_inputs$tree_sections_list
				# tree_pieces = tree_sections_list$return_pieces_list
				# num_pieces = length(tree_pieces)
				# tree_piece = tree_pieces[[piecenum]]
				txt = paste("timeperiod_i= ", timeperiod_i, ", piecenum=", piecenum, sep="")
				#("\n")
				#cat(txt)
			
				# Find the root node of the subtree
				contains_rootTF = master_table_timeperiod_i$SUBnode.type == "root"
				sub_piecenumTF = master_table_timeperiod_i$piecenum == piecenum
				sub_rootTF = (contains_rootTF + sub_piecenumTF) == 2
			
				# Error check
				if (sum(sub_rootTF) != 1)
					{
					errortxt = paste("\n\nERROR: subtree table must have a single node of SUBnode.type=='root'.\nPrinting subtree table:\n\n", sep="")
					cat(errortxt)
					print(master_table_timeperiod_i)
					stop("\nStopping on error.\n")
					}
			
				# Get the row number in the subtable
				subtable_rownum = (1:nrow(master_table_timeperiod_i))[sub_rootTF]
			
				# Check if there is a branch below the root node
				# (unless it's the global root, in which case, ignore)
				rootedge = FALSE
				if (master_table_timeperiod_i$node.type[subtable_rownum] == "root")
					{
					# It's the root of the full tree, so *NO* root edge below
					rootedge = FALSE
					} else {
					if ( (!is.na(master_table_timeperiod_i$edge.length[subtable_rownum])) && (!is.na(master_table_timeperiod_i$sampled_states_AT_brbots[subtable_rownum])) )
						{
						rootedge = TRUE
						}
					}
				
				
		
				
				
				
				
				#print(rootedge)
				stochastic_mapping_inputs$rootedge = rootedge
			
				# Stochastically map 
				stochastic_mapping_inputs$master_table_timeperiod_i = master_table_timeperiod_i
			
				#cat("\n\nTrying stochastic mapping on subtree:\n\n")
				treepiece_seed_counter = treepiece_seed_counter + 1
				seedval = seedval + (treepiece_seed_counter * treepiece_seed_counter_increment)

				if (seedval > 2147483647)
					{
					seedval = seedval %% 2147483647
					}
				set.seed(seed=seedval)
				
				# Here, we can just use res$condlikes, since for subtrees, the node downpass likelihoods are
				# for branch tops

# 				if (any(is.na(master_table_timeperiod_i$sampled_states_AT_nodes)))
# 					{
# 					txt = paste0("STOP ERROR_line1469 in (): Some of master_table_timeperiod_i$sampled_states_AT_nodes are NA. Specifically, these:")
# 					cat("\n\n")
# 					cat(txt)
# 					print("print( (1:length(master_table_timeperiod_i$sampled_states_AT_nodes))[is.na(master_table_timeperiod_i$sampled_states_AT_nodes)] )")
# 					
# 					print( (1:length(master_table_timeperiod_i$sampled_states_AT_nodes))[is.na(master_table_timeperiod_i$sampled_states_AT_nodes)] )
# 					print("timeperiod_i")
# 					print(timeperiod_i)
# 
# 					print("piecenum")
# 					print(piecenum)
# 					
# 					cat("\n\n")
# 					stop(txt)
# 					}

				#cat("\nA_stochastic_map_given_inputs():\n")
				master_table_timeperiod_i = stochastic_map_given_inputs(stochastic_mapping_inputs, piecenum=piecenum, maxtries=maxtries, seedval=seedval, include_null_range=res$inputs$include_null_range, master_nodenum_toPrint=master_nodenum_toPrint)
				#cat("\nB_stochastic_map_given_inputs():\n")
				
# 				if (any(is.na(master_table_timeperiod_i$sampled_states_AT_nodes)))
# 					{
# 					txt = paste0("STOP ERROR_line1490 in (): Some of master_table_timeperiod_i$sampled_states_AT_nodes are NA. Specifically, these:")
# 					cat("\n\n")
# 					cat(txt)
# 					print("print( (1:length(master_table_timeperiod_i$sampled_states_AT_nodes))[is.na(master_table_timeperiod_i$sampled_states_AT_nodes)] )")
# 					
# 					print( (1:length(master_table_timeperiod_i$sampled_states_AT_nodes))[is.na(master_table_timeperiod_i$sampled_states_AT_nodes)] )
# 					print("timeperiod_i")
# 					print(timeperiod_i)
# 
# 					print("piecenum")
# 					print(piecenum)
# 					
# 					cat("\n\n")
# 					stop(txt)
# 					}
				
				
				#cat("\n\nEnding stochastic mapping on subtree:\n\n")

				# Store the subtable back in inputs
				stochastic_mapping_inputs$master_table_timeperiod_i = master_table_timeperiod_i
			
				# Save the result in the list of inputs
				stochastic_mapping_inputs_list[[timeperiod_i]] = stochastic_mapping_inputs
			
			
				#######################################################
				# AND, we HAVE to copy up the simulated states to the timeperiod_i above
				#######################################################
				if (timeperiod_i != 1)
					{
					# 1 stratum up:
					timeperiod_i_up = timeperiod_i - 1
					stochastic_mapping_inputs_up = stochastic_mapping_inputs_list[[timeperiod_i_up]]
					master_table_timeperiod_i_up = stochastic_mapping_inputs_up$master_table_timeperiod_i
			
					# Find the nodes at the branch tops of the left and right branches
					# (which could be in the stratum above, 
					# or even higher, no worries, they all have the same node at the top)
				
					rows_w_correct_piecenum_TF = master_table_timeperiod_i$piecenum == piecenum
					rows_w_tips_TF = master_table_timeperiod_i$SUBnode.type == "tip"
					rows_w_subtree_tips_TF = (rows_w_correct_piecenum_TF + rows_w_tips_TF) == 2
					rownums_subtree_tips = (1:nrow(master_table_timeperiod_i))[rows_w_subtree_tips_TF]
				
					for (subtree_tipnum in rownums_subtree_tips)
						{
						#print(subtree_tipnum)
						# Get the global node number of this subtree tip
						desc_node_node_master_tree = master_table_timeperiod_i$node[subtree_tipnum]
					
						# Get the rownum of this global node number in the 
						# stratum up
						row_timeperiod_i_up_TF = master_table_timeperiod_i_up$node == desc_node_node_master_tree
						rownum_timeperiod_i_up = (1:nrow(master_table_timeperiod_i_up))[row_timeperiod_i_up_TF]
					
						# Store the known descendant state at the branch bottom 
						# in the next stratum up
						sampled_state_AT_brbot = master_table_timeperiod_i$sampled_states_AT_nodes[subtree_tipnum]
						if (isblank_TF(sampled_state_AT_brbot))
							{
							txt = paste0("STOP ERROR_line1514 in (): sampled_state_AT_brbot is blank: '", sampled_state_AT_brbot,  "'\n\nrownum_timeperiod_i_up=", rownum_timeperiod_i_up, ".\n\nsubtree_tipnum=", subtree_tipnum, ".")
							stop(txt)
							}
						
						master_table_timeperiod_i_up$sampled_states_AT_brbots[rownum_timeperiod_i_up] = sampled_state_AT_brbot
						
						if (isblank_TF(sampled_state_AT_brbot))
							{
							txt = paste0("STOP ERROR_line1520 in (): master_table_timeperiod_i_up$sampled_states_AT_brbots[rownum_timeperiod_i_up] is blank: '", master_table_timeperiod_i_up$sampled_states_AT_brbots[rownum_timeperiod_i_up], "'\n\nrownum_timeperiod_i_up=", rownum_timeperiod_i_up, ".\n\nsubtree_tipnum=", subtree_tipnum, ".")
							stop(txt)
							}
						
						} # END for (subtree_tipnum in rownums_subtree_tips)
				
					# This is the row for the node at the top of the branch
					# (or equivalent branch section in the stratum above)
					# master_table_timeperiod_i_up[rownum_above,]
				
					# Store the (now known) ancestral state at the bottom of the branches
					# in the stratum above
				
					# Store the node in the stratum above
					stochastic_mapping_inputs_list[[timeperiod_i_up]]$master_table_timeperiod_i = master_table_timeperiod_i_up 

					} # END if (timeperiod_i != 1)	
				} # END if (is.numeric(tree_piece))
			} # END for (piecenum in 1:num_pieces)
		} # END for (timeperiod_i in num_timeperiods:1)
	
	cat("\nFINISHED_stochastic_mapping_on_stratified()\n")

	# Display the stochastic mapping results
	master_table_w_stochastic_maps = NULL
	for (i in 1:length(stochastic_mapping_inputs_list))
		{
		master_table_w_stochastic_maps = rbind(master_table_w_stochastic_maps, stochastic_mapping_inputs_list[[i]]$master_table_timeperiod_i)
		}
	#printall(master_table_w_stochastic_maps[,-ncol(master_table_w_stochastic_maps)])
	dim(master_table_w_stochastic_maps)


	# Get the anagenetic events in a nice text form
	#print("print(master_table_w_stochastic_maps$anagenetic_events_txt_below_node):")
	#print(master_table_w_stochastic_maps$anagenetic_events_txt_below_node)
	events_table = events_txt_list_into_events_table(events_txt_list=master_table_w_stochastic_maps$anagenetic_events_txt_below_node, trtable=master_table_w_stochastic_maps)
	#print(events_table)

	stochastic_mapping_results = NULL
	stochastic_mapping_results$master_table_cladogenetic_events = master_table_w_stochastic_maps

#print(master_table_w_stochastic_maps)

	stochastic_mapping_results$table_w_anagenetic_events = events_table

	return(stochastic_mapping_results)
	} # END stochastic_mapping_on_stratified()









