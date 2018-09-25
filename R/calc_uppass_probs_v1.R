

sample_uppass_split_scenario_given_probs_ancstate <- function(probs_ancstate, COO_weights_columnar, numstates, include_null_range, left_branch_downpass_likes, right_branch_downpass_likes, Rsp_rowsums=NULL)
	{
	options(stringsAsFactors=FALSE)
	
	# Converting relative scenario probabilities to absolute scenario probabilities

	# Setup
	if (is.null(Rsp_rowsums))
		{
		Rsp_rowsums = rcpp_calc_rowsums_for_COOweights_columnar(COO_weights_columnar=COO_weights_columnar, numstates=numstates, printmat=FALSE)
		} # END if (is.null(Rsp_rowsums))
	if (include_null_range == TRUE)
		{
		COOmat_0based_to_Qmat_1based = 2
		} else {
		COOmat_0based_to_Qmat_1based = 1
		} # END if (include_null_range == TRUE)

	# Converting relative scenario probabilities to absolute scenario probabilities
	# Uppass scenario probabilities INCLUDING input downpass likelihoods
	condprob_each_split_scenario_df2 = calc_uppass_scenario_probs_new2(probs_ancstate, COO_weights_columnar, numstates, include_null_range=include_null_range, left_branch_downpass_likes=left_branch_downpass_likes, right_branch_downpass_likes=right_branch_downpass_likes, Rsp_rowsums=Rsp_rowsums)
	condprob_each_split_scenario_df2 = matrix(data=c(condprob_each_split_scenario_df2), nrow=length(COO_weights_columnar[[1]]), byrow=FALSE)
	condprob_each_split_scenario_df2
	#round(apply(X=condprob_each_split_scenario_df2, MARGIN=1, mean),3)

	# Converting relative scenario probabilities to absolute scenario probabilities
	#relprobs_each_split_scenario = apply(X=condprob_each_split_scenario_df2, MARGIN=1, sum)
	relprobs_each_split_scenario = condprob_each_split_scenario_df2
	probs_each_split_scenario = relprobs_each_split_scenario / sum(relprobs_each_split_scenario)
	round(probs_each_split_scenario, 3)
	
	# Sample a split scenario
	split_scenario_nums = 1:length(COO_weights_columnar[[1]])
	split_scenario_num = sample(x=split_scenario_nums, size=1, replace=TRUE, prob=probs_each_split_scenario)
	
	# Get the ancestor state, Left corner, Right corner
	ancstate_1based = COOmat_0based_to_Qmat_1based + COO_weights_columnar[[1]][split_scenario_num]
	left_decstate_1based = COOmat_0based_to_Qmat_1based + COO_weights_columnar[[2]][split_scenario_num]
	right_decstate_1based = COOmat_0based_to_Qmat_1based + COO_weights_columnar[[3]][split_scenario_num]


	
	# Outputs
	sample_uppass_res = list()
	sample_uppass_res$probs_ancstate = probs_ancstate
	sample_uppass_res$condprob_each_split_scenario_df2 = condprob_each_split_scenario_df2
	sample_uppass_res$probs_each_split_scenario = probs_each_split_scenario
	sample_uppass_res$split_scenario_num = split_scenario_num
	sample_uppass_res$ancstate_1based = ancstate_1based
	sample_uppass_res$left_decstate_1based = left_decstate_1based
	sample_uppass_res$right_decstate_1based = right_decstate_1based
	
	extract='
	probs_ancstate = sample_uppass_res$probs_ancstate
	condprob_each_split_scenario_df2 = sample_uppass_res$condprob_each_split_scenario_df2
	probs_each_split_scenario = sample_uppass_res$probs_each_split_scenario
	split_scenario_num = sample_uppass_res$split_scenario_num
	ancstate_1based = sample_uppass_res$ancstate_1based
	left_decstate_1based = sample_uppass_res$left_decstate_1based
	right_decstate_1based = sample_uppass_res$right_decstate_1based	
	' # END extraact
	
	return(sample_uppass_res)
	} # END sample_uppass_split_scenario_given_probs_ancstate()




calc_uppass_scenario_probs_new2 <- function(probs_ancstate, COO_weights_columnar, numstates, include_null_range=include_null_range, left_branch_downpass_likes=NULL, right_branch_downpass_likes=NULL, Rsp_rowsums=NULL)
	{
	defaults='
	probs_ancstate = relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc,]
	' # END defaults
	
	# Setup
	
	# If no left_branch_downpass_likes, right_branch_downpass_likes,
	# then set them to 0
	if (is.null(left_branch_downpass_likes))
		{
		left_branch_downpass_likes = rep(1, numstates)
		} # END if (is.null(left_branch_downpass_likes))
	if (is.null(right_branch_downpass_likes))
		{
		right_branch_downpass_likes = rep(1, numstates)
		} # END if (is.null(right_branch_downpass_likes))
	if (is.null(Rsp_rowsums))
		{
		Rsp_rowsums = rcpp_calc_rowsums_for_COOweights_columnar(COO_weights_columnar=COO_weights_columnar, numstates=numstates, printmat=FALSE)
		} # END if (is.null(Rsp_rowsums))
	if (include_null_range == TRUE)
		{
		COOmat_0based_to_Qmat_1based = 2
		} else {
		COOmat_0based_to_Qmat_1based = 1
		} # END if (include_null_range == TRUE)
	
	# Uppass probabilities
	condprob_each_split_scenario_df = NULL

	# Calculate probabilities, assuming each possible 
	# ancestral state at the node (we will multiply
	# by the node probabilities later)
	#for (s in 1:numstates)
	#	{
	#sampled_state = s
	#probs_ancstate_tmp = rep(0, length(probs_ancstate))
	#probs_ancstate_tmp[sampled_state] = 1
	
	# For sampling
	ancprobs = probs_ancstate[ COO_weights_columnar[[1]] + COOmat_0based_to_Qmat_1based]
	# For average probs
	#ancprobs = probs_ancstate[ COO_weights_columnar[[1]] + COOmat_0based_to_Qmat_1based]

	Lprobs = left_branch_downpass_likes[ COO_weights_columnar[[2]] + COOmat_0based_to_Qmat_1based]
	Rprobs = right_branch_downpass_likes[ COO_weights_columnar[[3]] + COOmat_0based_to_Qmat_1based]
	# Weights divided by the sum of the weights for that row
	# COO_weights_columnar[[1]] is the ancestral node, 0-based, cladogenesis-only
	# So, you will always add 1 to the index to get the rowsum conditional
	# on a particular ancestor (3 non-null states means 3 Rsp_rowsums)
	scenario_condprob = COO_weights_columnar[[4]] / Rsp_rowsums[ COO_weights_columnar[[1]] + 1 ]

	# Bind the probabilities of the ancestor states, and 
	# downpass likelihoods from each descendant, and the scenario
	# probabilities
	input_probs_each_split_scenario = cbind(ancprobs, Lprobs, Rprobs, scenario_condprob)
	round(input_probs_each_split_scenario, 3)
	
	
	
	# Multiply the probabilities through
	relprob_each_split_scenario = apply(X=input_probs_each_split_scenario, MARGIN=1, FUN=prod)
	round(relprob_each_split_scenario,3)
	
	# Error check
	sum_relprob_each_split_scenario = sum(relprob_each_split_scenario)

	#cat("Printing sum(relprob_each_split_scenario):")
	#print(sum(relprob_each_split_scenario))
	
	if (any(is.na(sum_relprob_each_split_scenario)))
		{
		txt = paste0("STOP ERROR in calc_uppass_scenario_probs_new2(): NAs found in sum_relprob_each_split_scenario")
		cat("\n\n")
		cat(txt)

		cat("\n\n")
		cat("Printing relprob_each_split_scenario:")
		print(relprob_each_split_scenario)


		cat("\n\n")
		cat("Printing sum(relprob_each_split_scenario):")
		print(sum(relprob_each_split_scenario))

		# Print for this bug check:
		# Uppass starting for marginal ancestral states estimation!
		#  Error in if (sum(relprob_each_split_scenario) > 0) { : 
		#   missing value where TRUE/FALSE needed
		cat("\n\n")
		cat("Printing COO_weights_columnar:")
		print(COO_weights_columnar)
		
		cat("\n\n")
		cat("Printing probs_ancstate:")
		print(probs_ancstate)
		
		cat("\n\n")
		cat("Printing left_branch_downpass_likes:")
		print(left_branch_downpass_likes)
		
		cat("\n\n")
		cat("Printing right_branch_downpass_likes:")
		print(right_branch_downpass_likes)
		
		cat("\n\n")
		cat("Printing Lprobs:")
		print(Lprobs)

		cat("\n\n")
		cat("Printing Rprobs:")
		print(Rprobs)

		cat("\n\n")
		cat("Printing input_probs_each_split_scenario:")
		print(input_probs_each_split_scenario)
		
		
		stop(txt)
		}
	
	if (sum(relprob_each_split_scenario) > 0)
		{
		# Old 2015
		#prob_each_split_scenario = relprob_each_split_scenario / sum(relprob_each_split_scenario)
		# New 2015
		prob_each_split_scenario = relprob_each_split_scenario #* probs_ancstate[s]
		round(prob_each_split_scenario, 3)

		#split_scenario_nums = 1:length(COO_weights_columnar[[1]])
		# Sample one of them
		#split_scenario_num = sample(x=split_scenario_nums, size=1, replace=TRUE, prob=prob_each_split_scenario)
		} else {
		prob_each_split_scenario = rep(0, length(ancprobs))
		}
	
	# condprob_split_starting_from_1based_state1 = prob_each_split_scenario
	#condprobs_name = paste0("condprob_split_starting_from_1based_state", s)
	#cmdstr = paste0(condprobs_name, " = prob_each_split_scenario")
	#eval(parse(text=cmdstr))

	#condprob_each_split_scenario_df = cbind(condprob_each_split_scenario_df, prob_each_split_scenario)
	#cmdstr = paste0("condprob_each_split_scenario_df = cbind(condprob_each_split_scenario_df, ", condprobs_name, ")")
	#eval(parse(text=cmdstr))
		
	#	} # END for (s in 1:numstates)
	#condprob_each_split_scenario_df

	#head(condprob_each_split_scenario_df)
	#round(apply(X=condprob_each_split_scenario_df, MARGIN=2, mean),3)

	# Copy the split scenario probabilities, conditional on each 
	# possible ancestral state
	#condprob_each_split_scenario_df2 = condprob_each_split_scenario_df

	# Sum, weighting by the probability of the ancestral state
	#for (r in 1:length(probs_ancstate))
	#	{
		#condprob_each_split_scenario_df2[,r] = condprob_each_split_scenario_df[,r] * probs_ancstate[r]
	#	}
	
	#condprob_each_split_scenario_df2
	#round(rowSums(condprob_each_split_scenario_df2),3)
	#sum(rowSums(condprob_each_split_scenario_df2))
	
	# Convert to relative and absolute scenario probabilities:
	conversion='
	relprobs_each_split_scenario = apply(X=condprob_each_split_scenario_df2, MARGIN=1, sum)
	probs_each_split_scenario = relprobs_each_split_scenario / sum(relprobs_each_split_scenario)
	round(probs_each_split_scenario, 3)
	' # END conversion
	
	prob_each_split_scenario = prob_each_split_scenario / sum(prob_each_split_scenario)
	
	return(prob_each_split_scenario)
	} # END calc_uppass_scenario_probs_new2()


prob_each_split_scenario_to_uppass_probs <- function(prob_each_split_scenario, COO_weights_columnar, numstates, include_null_range=TRUE)
	{
	ancprobs = rep(0, numstates)
	Lprobs = rep(0, numstates)
	Rprobs = rep(0, numstates)

	if (include_null_range == TRUE)
		{
		COOmat_0based_to_Qmat_1based = 2
		} else {
		COOmat_0based_to_Qmat_1based = 1
		} # END if (include_null_range == TRUE)

	# Add up the probabilities over all scenarios
	for (nn in 1:length(prob_each_split_scenario))
		{
		ancstate_tmp = (COO_weights_columnar[[1]][nn] + COOmat_0based_to_Qmat_1based)
		ancprobs[ancstate_tmp] = ancprobs[ancstate_tmp] + prob_each_split_scenario[nn]
		
		Lstate_tmp = (COO_weights_columnar[[2]][nn] + COOmat_0based_to_Qmat_1based)
		Lprobs[Lstate_tmp] = Lprobs[Lstate_tmp] + prob_each_split_scenario[nn]

		Rstate_tmp = (COO_weights_columnar[[3]][nn] + COOmat_0based_to_Qmat_1based)
		Rprobs[Rstate_tmp] = Rprobs[Rstate_tmp] + prob_each_split_scenario[nn]
		}
	
	clado_uppass_probs = NULL
	clado_uppass_probs$ancprobs = ancprobs
	clado_uppass_probs$Lprobs = Lprobs
	clado_uppass_probs$Rprobs = Rprobs
	
	return(clado_uppass_probs)
	} # END prob_each_split_scenario_to_uppass_probs <- function(prob_each_split_scenario, COO_weights_columnar, numstates, include_null_range=TRUE)



calc_uppass_probs_new2 <- function(probs_ancstate, COO_weights_columnar, numstates, include_null_range=include_null_range, left_branch_downpass_likes=NULL, right_branch_downpass_likes=NULL, Rsp_rowsums=NULL)
	{
# 	print("COO_weights_columnar")
# 	print(COO_weights_columnar)
# 	print("probs_ancstate")
# 	print(probs_ancstate)
# 	
# 	print("left_branch_downpass_likes")
# 	print(left_branch_downpass_likes)
# 	print("right_branch_downpass_likes")
# 	print(right_branch_downpass_likes)
# 	
# 	print("numstates")
# 	print(numstates)

	# Error check
	if (FALSE)
		{
		txt = paste0("ERROR CHECK in calc_uppass_probs_new2()")
		cat("\n\n")
		cat(txt)
		cat("\n\n")

		# Print for this bug check:
		# Uppass starting for marginal ancestral states estimation!
		#  Error in if (sum(relprob_each_split_scenario) > 0) { : 
		#   missing value where TRUE/FALSE needed
		cat("\n\n")
		cat("Printing COO_weights_columnar:")
		print(COO_weights_columnar)
		
		cat("\n\n")
		cat("Printing probs_ancstate:")
		print(probs_ancstate)
		
		cat("\n\n")
		cat("Printing left_branch_downpass_likes:")
		print(left_branch_downpass_likes)
		
		cat("\n\n")
		cat("Printing right_branch_downpass_likes:")
		print(right_branch_downpass_likes)
		
		#stop(txt)
		}


	
	# Uppass scenario probabilities
	condprob_each_split_scenario_df2 = calc_uppass_scenario_probs_new2(probs_ancstate, COO_weights_columnar, numstates, include_null_range, left_branch_downpass_likes, right_branch_downpass_likes, Rsp_rowsums)
	condprob_each_split_scenario_df2 = matrix(data=c(condprob_each_split_scenario_df2), nrow=length(COO_weights_columnar[[1]]), byrow=FALSE)
	
	# Converting relative scenario probabilities to absolute scenario probabilities
	#relprobs_each_split_scenario = apply(X=condprob_each_split_scenario_df2, MARGIN=1, sum)
	probs_each_split_scenario = condprob_each_split_scenario_df2 / sum(condprob_each_split_scenario_df2)
	round(probs_each_split_scenario, 3)

	# Setup
	Lprobs = rep(0, numstates)
	Lstates_1based = 1:numstates
	Rprobs = rep(0, numstates)
	Rstates_1based = 1:numstates
	uppass_Lprobs = rep(0, numstates)
	uppass_Rprobs = rep(0, numstates)

	if (include_null_range == TRUE)
		{
		COOmat_0based_to_Qmat_1based = 2
		} else {
		COOmat_0based_to_Qmat_1based = 1
		} # END if (include_null_range == TRUE)

	# Add up the probabilites of the left states, and right states
	for (n in 1:numstates)
		{
		scenarios_that_match_TF = Lstates_1based[n] == (COO_weights_columnar[[2]] + COOmat_0based_to_Qmat_1based)
	
		tmp_Lprobs = probs_each_split_scenario[scenarios_that_match_TF]
		uppass_Lprobs[n] = sum(tmp_Lprobs)

		scenarios_that_match_TF = Rstates_1based[n] == (COO_weights_columnar[[3]] + COOmat_0based_to_Qmat_1based)
	
		tmp_Rprobs = probs_each_split_scenario[scenarios_that_match_TF]
		uppass_Rprobs[n] = sum(tmp_Rprobs)
		}
	round(uppass_Lprobs, 3)
	round(uppass_Rprobs, 3)

	relprobs_just_after_speciation = list()
	relprobs_just_after_speciation$condprob_each_split_scenario_df2 = condprob_each_split_scenario_df2
	relprobs_just_after_speciation$relprobs_just_after_speciation_UPPASS_Left = uppass_Lprobs
	relprobs_just_after_speciation$relprobs_just_after_speciation_UPPASS_Right = uppass_Rprobs
	
	extract='
	relprobs_just_after_speciation = calc_uppass_probs_orig(anc, relative_probs_of_each_state_at_branch_top_AT_node_UPPASS, COO_weights_columnar, numstates, include_null_range)
	condprob_each_split_scenario_df2 = relprobs_just_after_speciation$condprob_each_split_scenario_df2
	relprobs_just_after_speciation_UPPASS_Left = relprobs_just_after_speciation$relprobs_just_after_speciation_UPPASS_Left
	relprobs_just_after_speciation_UPPASS_Right = relprobs_just_after_speciation$relprobs_just_after_speciation_UPPASS_Right
	' # END extract
	
	return(relprobs_just_after_speciation)
	} # END calc_uppass_probs_new





