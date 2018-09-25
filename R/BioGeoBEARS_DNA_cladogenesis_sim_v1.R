


# Set up the tip conditional likelihoods for DNA
tip_simstate_nums_to_tip_condlikes_of_data_on_each_state <- function(simulated_states_by_node, numtips, numstates=NULL)
	{
	if (is.null(numstates))
		{
		numstates = length(unique(sort(simulated_states_by_node)))
		}
	
	tip_condlikes_of_data_on_each_state = matrix(0, nrow=numtips, ncol=numstates)
	
	for (i in 1:numtips)
		{
		statenum_1based = as.numeric(simulated_states_by_node[i])
		tip_condlikes_of_data_on_each_state[i, statenum_1based] = 1
		}
	return(tip_condlikes_of_data_on_each_state)
	}

calc_loglike_for_optim_DNA_4param <- function(params, phy, Qmat_txt, cladogenesis_txt, tip_condlikes_of_data_on_each_state)
	{
	defaults='
	phy = yuletree
	params = c(a1, a2, c1, c2)
	params_lower = c(0, 0, 0, 0)
	params_upper = c(2, 2, 2, 2)
	'
	
	# Extract the parameters
	a1 = params[1]
	a2 = params[2]
	c1 = params[3]
	c2 = params[4]
	
	# Calculate the Q matrix
	Qmat = DNA_anagenesis_to_Qmat(a1, a2, Qmat_txt)
	Qmat
	
	# Calculate the cladogenesis matrix
	cladogenesis_probs = DNA_cladogenesis_to_inheritance_condprobs(c1, c2, cladogenesis_txt)
	cladogenesis_probs

	total_loglikelihood = calc_DNA_loglike(tip_condlikes_of_data_on_each_state, phy, Qmat, cladogenesis_probs, returnwhat="loglike")
	total_loglikelihood
	
	tmprow = matrix(data=c(a1, a2, c1, c2, total_loglikelihood), nrow=1)
	result = adf2(tmprow)
	names(result) = c("a1", "a2", "c1", "c2", "LnL")
	
	print (result)
	
	return(total_loglikelihood)
	}




calc_loglike_for_optim_DNA_2ana <- function(params, phy, Qmat_txt, cladogenesis_txt, tip_condlikes_of_data_on_each_state)
	{
	defaults='
	phy = yuletree
	params = c(a1, a2)
	params_lower = c(0, 0)
	params_upper = c(2, 2)
	'
	
	# Extract the parameters
	a1 = params[1]
	a2 = params[2]
	c1 = 0
	c2 = 0
	
	# Calculate the Q matrix
	Qmat = DNA_anagenesis_to_Qmat(a1, a2, Qmat_txt)
	Qmat
	
	# Calculate the cladogenesis matrix
	cladogenesis_probs = DNA_cladogenesis_to_inheritance_condprobs(c1, c2, cladogenesis_txt)
	cladogenesis_probs

	total_loglikelihood = calc_DNA_loglike(tip_condlikes_of_data_on_each_state, phy, Qmat, cladogenesis_probs, returnwhat="loglike")
	total_loglikelihood
	
	tmprow = matrix(data=c(a1, a2, c1, c2, total_loglikelihood), nrow=1)
	result = adf2(tmprow)
	names(result) = c("a1", "a2", "c1", "c2", "LnL")
	
	print (result)
	
	return(total_loglikelihood)
	}



calc_loglike_for_optim_DNA_2clada <- function(params, phy, Qmat_txt, cladogenesis_txt, tip_condlikes_of_data_on_each_state)
	{
	defaults='
	phy = yuletree
	params = c(c1, c2)
	params_lower = c(0, 0)
	params_upper = c(2, 2)
	'
	
	# Extract the parameters
	a1 = 0
	a2 = 0
	c1 = params[1]
	c2 = params[2]
	
	# Calculate the Q matrix
	Qmat = DNA_anagenesis_to_Qmat(a1, a2, Qmat_txt)
	Qmat
	
	# Calculate the cladogenesis matrix
	cladogenesis_probs = DNA_cladogenesis_to_inheritance_condprobs(c1, c2, cladogenesis_txt)
	cladogenesis_probs

	total_loglikelihood = calc_DNA_loglike(tip_condlikes_of_data_on_each_state, phy, Qmat, cladogenesis_probs, returnwhat="loglike")
	total_loglikelihood
	
	tmprow = matrix(data=c(a1, a2, c1, c2, total_loglikelihood), nrow=1)
	result = adf2(tmprow)
	names(result) = c("a1", "a2", "c1", "c2", "LnL")
	
	print (result)
	
	return(total_loglikelihood)
	}




# Calculate likelihood given anagenetic and cladogenetic model parameters
calc_DNA_loglike <- function(tip_condlikes_of_data_on_each_state, phy, Qmat, cladogenesis_probs, returnwhat="loglike", numstates=4)
	{
	defaults='
	phy=yuletree
	returnwhat="loglike"
	numstates = 4
	'
	

	numstates = ncol(tip_condlikes_of_data_on_each_state)
	edgelengths = phy$edge.length
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
	computed_likelihoods_at_each_node[tipnums] = rowSums(tip_condlikes_of_data_on_each_state)

	condlikes_of_each_state <- matrix(data=0, nrow=numnodes, ncol=numstates)
	relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS <- matrix(data=0, nrow=numnodes, ncol=numstates)


	relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[tipnums, ] = tip_condlikes_of_data_on_each_state #/ rowSums(tip_condlikes_of_data_on_each_state)
	
	# BUT, DO NOT DIVIDE THIS BY rowSums(tip_condlikes_of_data_on_each_state), it forces normalization which prevents e.g.
	# passing down likelihoods during stratification
	condlikes_of_each_state[tipnums, ] = tip_condlikes_of_data_on_each_state

	relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS <- matrix(data=NA, nrow=numnodes, ncol=numstates)
	ML_marginal_prob_each_state_at_branch_bottom_below_node <- matrix(data=NA, nrow=numnodes, ncol=numstates)
	ML_marginal_prob_each_state_at_branch_top_AT_node <- matrix(data=NA, nrow=numnodes, ncol=numstates)
	ML_marginal_prob_each_split_at_branch_top_AT_node = list()
	
	
	
	
	
	independent_likelihoods_on_each_branch = vector("list", length(phy2$edge.length))
	tmpmatrix = matrix(data=0, nrow=nrow(Qmat), ncol=ncol(Qmat))
	for (m in 1:length(phy2$edge.length))
		{
		independent_likelihoods_on_each_branch[[m]] = tmpmatrix
		}
		
	independent_likelihoods_on_each_branch = mapply_likelihoods(Qmat, phy2, transpose_needed=TRUE)
	
	i = 1
	edges_to_visit = seq(from=1, by=2, length.out=num_internal_nodes)
	
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
		
		# And for the ancestor edge (i or j shouldn't matter, should produce the same result!!!)
		anc <- phy2$edge[i, 1]
		
		txt = paste("anc:", anc, " left:", left_desc_nodenum, " right:", right_desc_nodenum, sep="")
		#print(txt)
		
		# Conditional likelihoods of states at the bottom of left branch
		condlikes_Left = independent_likelihoods_on_each_branch[,,i] %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,]
					
		# Conditional likelihoods of states at the bottom of right branch
		condlikes_Right = independent_likelihoods_on_each_branch[,,j] %*% relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum,]

		# Every node (except maybe the root) has a branch below it, and there is also a 
		# relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS at the bottom of this branch
		relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[left_desc_nodenum,] = condlikes_Left / sum(condlikes_Left)
		relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[right_desc_nodenum,] = condlikes_Right / sum(condlikes_Right)
		
		# get the correct edge
		left_edge_TF = phy2$edge[,2] == left_desc_nodenum
		right_edge_TF = phy2$edge[,2] == right_desc_nodenum
		
		# for each ancestral state, get prob of branch pairs 
		outmat = matrix(0, nrow=nrow(cladogenesis_probs), ncol=numstates)
		
		anc_condLikes_for_this_node = rep(0, numstates)
		
		# cn = cladogenesis matrix rownumber
		for (cn in 1:nrow(cladogenesis_probs))
			{
			startprob = anc_condLikes_for_this_node[cladogenesis_probs$anc_ind[cn]]
			downpass_prob = condlikes_Left[cladogenesis_probs$L_ind[cn]] * condlikes_Right[cladogenesis_probs$R_ind[cn]] * cladogenesis_probs$probs[cn]
			
			anc_condLikes_for_this_node[cladogenesis_probs$anc_ind[cn]] = startprob + downpass_prob
			}
		
		# Store results
		node_likelihood_with_speciation = anc_condLikes_for_this_node
		node_likelihood = node_likelihood_with_speciation
		total_likelihood_for_node = sum(node_likelihood)
		
		# Total likelihoods
		computed_likelihoods_at_each_node[anc] = total_likelihood_for_node
		relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[anc, ] = node_likelihood / total_likelihood_for_node
		condlikes_of_each_state[anc, ] = node_likelihood
		} # end downpass
	#######################################################
	# END DOWNPASS FROM THE TIPS TO THE ROOT
	#######################################################
	
	
	relative_probs_of_each_state_at_bottom_of_root_branch = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[anc, ]

	total_loglikelihood = sum(log(computed_likelihoods_at_each_node))



	calc_loglike_sp_results = list()
	calc_loglike_sp_results$computed_likelihoods_at_each_node = computed_likelihoods_at_each_node
	calc_loglike_sp_results$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS
	calc_loglike_sp_results$condlikes_of_each_state = condlikes_of_each_state
	calc_loglike_sp_results$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS
	
	calc_loglike_sp_results$relative_probs_of_each_state_at_bottom_of_root_branch = relative_probs_of_each_state_at_bottom_of_root_branch
	calc_loglike_sp_results$total_loglikelihood = total_loglikelihood
	class(calc_loglike_sp_results) = "calc_loglike_sp_results"
	
	if (returnwhat == "all")
		{
		return(calc_loglike_sp_results)
		}

	if (returnwhat == "loglike")
		{
		return(total_loglikelihood)
		}
	}



# Simulate DNA on the tree, given Qmat and cladogenesis_probs
sim_DNA_clado <- function(phy=yuletree, Qmat, cladogenesis_probs, starting_state=1)
	{
	defaults='
	phy=yuletree
	Qmat
	cladogenesis_probs
	starting_state=1	# start in A
	'
	
	# Reorder the phylogeny to pruning-wise order...
	phy2 <- reorder(phy, "pruningwise")
	
	numedges = nrow(phy2$edge)
	ntips = length(phy2$tip.label)
	num_internal_nodes = phy2$Nnode
	numnodes = ntips + num_internal_nodes
	nodenums = c(1:numnodes)
	internal_nodenums = (ntips+1):numnodes
	
	# Number of states in the full model
	numstates = nrow(Qmat)
	
	# Get the ancestral node
	# (The last 2 node in the ancestor column in a pruningwise edge matrix
	#  are the ancestor)
	anc_nodenum = phy2$edge[numedges, 1]

	# Set up the list of simulated states
	simulated_states_by_node = rep(NA, numnodes)
	simulated_states_by_node[anc_nodenum] = starting_state
	
	# simulated_states_by_node_LR_after_speciation
	simulated_states_by_node_LR_after_speciation = matrix(data=NA, nrow=numnodes, ncol=2)
	
	# Label the edges in reverse pruningwise order
	edges_to_visit_j = seq(from=numedges, by=-2, length.out=num_internal_nodes)
	edges_to_visit_i = edges_to_visit_j - 1
	
	# Travel up the tree and simulate cladogenesis and anagenesis events
	for (i in 1:length(edges_to_visit_i))
		{
		#######################################################
		# cladogenetic range-inheritance
		#######################################################
		# Do left descendant
		left_edge_to_visit = edges_to_visit_i[i]
		right_edge_to_visit = edges_to_visit_j[i]
		starting_nodenum = phy2$edge[left_edge_to_visit,1]
		
		# left and right descendants
		left_desc_nodenum = phy2$edge[left_edge_to_visit,2]
		right_desc_nodenum = phy2$edge[right_edge_to_visit,2]
		
		# Get the starting state for the branches
		starting_state_Qmat_1index = simulated_states_by_node[starting_nodenum]
		
		# Get the descendent states just after speciation
		TF = cladogenesis_probs$anc_ind == starting_state_Qmat_1index
		transition_matrix_conditional_on_ancstate = cladogenesis_probs[TF,]
		scenario_nums = 1:nrow(transition_matrix_conditional_on_ancstate)
		scenario_chosen = sample(x=scenario_nums, size=1, replace=FALSE, prob=transition_matrix_conditional_on_ancstate$probs)
		
		simulated_states_by_node_LR_after_speciation[starting_nodenum, ] = unlist(transition_matrix_conditional_on_ancstate[scenario_chosen, c("L_ind","R_ind")])
		
		simulated_states_by_node_LR_after_speciation
		
		#######################################################
		# anagenetic range-inheritance
		#######################################################
		# Get the descendant state and the ends of the two branches
		
		# Simulate end of LEFT branch
		branchlength = phy2$edge.length[left_edge_to_visit]
		branch_bot_state_1based = simulated_states_by_node_LR_after_speciation[starting_nodenum, 1] - 0
		branch_bot_state_1based
		
		# Exponentiate forward
		# (True/false doesn't matter for time-reversible model)
		Pmat = expokit_dgpadm_Qmat2(Qmat=Qmat, times = branchlength, 
        transpose_needed = TRUE)

		desc_state_probs = Pmat[branch_bot_state_1based,]
		desc_state_1based = sample(x=1:4, size=1, replace=FALSE, prob=desc_state_probs)
		simulated_states_by_node[left_desc_nodenum] = desc_state_1based
		
		
		# Simulate end of RIGHT branch
		branchlength = phy2$edge.length[right_edge_to_visit]
		branch_bot_state_1based = simulated_states_by_node_LR_after_speciation[starting_nodenum, 2] - 0
		branch_bot_state_1based
		
		# Exponentiate forward
		# (True/false doesn't matter for time-reversible model)
		Pmat = expokit_dgpadm_Qmat2(Qmat=Qmat, times = branchlength, 
        transpose_needed = TRUE)

		desc_state_probs = Pmat[branch_bot_state_1based,]
		desc_state_1based = sample(x=1:4, size=1, replace=FALSE, prob=desc_state_probs)
		simulated_states_by_node[right_desc_nodenum] = desc_state_1based
		} # end loop through branches
	# End simulation
	return(simulated_states_by_node)
	}



# Get the Qmat for anagenesis, given param values
DNA_anagenesis_to_Qmat <- function(a1, a2, Qmat_txt)
	{
	defaults='
	a1=0.25
	a2=0.25
	'
	
	Qmat = matrix(data=NA, nrow=nrow(Qmat_txt), ncol=ncol(Qmat_txt))
	TF = Qmat_txt == "a1"
	Qmat[TF] = a1
	TF = Qmat_txt == "a2"
	Qmat[TF] = a2
	
	diag(Qmat) = -rowSums(Qmat, na.rm=TRUE)
	
	return(Qmat)
	}



# Get the probabilities of each cladogenetic transition, 
# given param values
DNA_cladogenesis_to_inheritance_condprobs <- function(c1, c2, cladogenesis_txt)
	{
	defaults='
	# JC69 (Jukes-Cantor)
	c1 = 1
	c2 = 1
	
	# K2p (Kimura 1980 2-parameter)
	c1 = 2
	c2 = 0.5
	
	# Standard DNA model (no cladogenesis process)
	c1 = 0
	c2 = 0
	
	' # end defaults
	
	cladogenesis_probs = matrix(data=NA, ncol=6, nrow=nrow(cladogenesis_txt))
	
	# Code A as index 1, C as index 2, etc.
	TF = cladogenesis_txt[,1:3] == "A"
	cladogenesis_probs[,1:3][TF] = 1
	TF = cladogenesis_txt[,1:3] == "C"
	cladogenesis_probs[,1:3][TF] = 2
	TF = cladogenesis_txt[,1:3] == "G"
	cladogenesis_probs[,1:3][TF] = 3
	TF = cladogenesis_txt[,1:3] == "T"
	cladogenesis_probs[,1:3][TF] = 4
	cladogenesis_probs

	# Assign "-" a weight of 1
	TF = cladogenesis_txt == "-"
	cladogenesis_probs[TF] = 1

	# Assign parameter values
	TF = cladogenesis_txt == "c1"
	cladogenesis_probs[TF] = c1
	TF = cladogenesis_txt == "c2"
	cladogenesis_probs[TF] = c2

	# Calculate sum of weights for each of the 4 ancestral states
	for (ind in 1:4)
		{
		# Ancestral node index
		TF = cladogenesis_probs[,1] == ind
		sum_of_weights = sum(cladogenesis_probs[TF,4])
		cladogenesis_probs[TF, 5] = sum_of_weights
		}
	
	# Calculate the conditional probabilities of each scenario
	cladogenesis_probs[,6] = cladogenesis_probs[,4] / cladogenesis_probs[,5]
	
	cladogenesis_probs = adf2(cladogenesis_probs)
	names(cladogenesis_probs) = c("anc_ind", "L_ind", "R_ind", "weights", "rowsums", "probs")
	cladogenesis_probs
	
	# The sum of the conditional probabilities for 4 different ancestors should be 4
	sum(cladogenesis_probs$probs)
	
	return(cladogenesis_probs)
	}











