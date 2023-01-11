

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
	
	# Error check
	TF1 = is.na(split_scenario_nums)
	TF2 = is.nan(split_scenario_nums)
	TF = (TF1 + TF2) > 0
	if (sum(TF) > 0)
		{
		txt = paste0("STOP ERROR in: sample_uppass_split_scenario_given_probs_ancstate() -- NA or NaN in split_scenario_nums. Printing 'COO_weights_columnar' and 'split_scenario_nums':")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		print("COO_weights_columnar")
		print(COO_weights_columnar)
		cat("\n\n")
		print("split_scenario_nums")
		print(split_scenario_nums)
		stop(txt)
		}
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










#######################################################
# EXAMPLE TO CHECK calc_uppass_probs_new2
#######################################################
calc_uppass_probs_new2_example <- function()
	{
# dput(resDEC)

resDEC = structure(list(computed_likelihoods_at_each_node = c(1, 1, 1, 
1, 0.36502788215878, 0.227783575654047, 0.136167332530107), relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS = structure(c(0, 
0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 6.42693382782259e-14, 7.46190547654875e-14, 
3.13642570191625e-13, 1, 0, 1, 1, 0.757828232249181, 0.344749222130095, 
3.13642570191625e-13, 0, 0, 0, 0, 0.242171767750754, 0.65525077786983, 
0.999999999999373), dim = c(7L, 4L)), condlikes_of_each_state = structure(c(0, 
0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2.3460100439447e-14, 1.69969951064079e-14, 
4.27078721508805e-14, 1, 0, 1, 1, 0.276628434658051, 0.0785282105207443, 
4.27078721508805e-14, 0, 0, 0, 0, 0.0883994475007057, 0.149255365133286, 
0.136167332530022), dim = c(7L, 4L)), relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS = structure(c(0, 
0, 0, 0, NA, 0, 0, 5.22737626427221e-14, 0.999999999998895, 2.16444459101218e-13, 
5.04411011774348e-13, NA, 0.0576312779227834, 0.0806194459381856, 
0.999999999998895, 5.22737626427221e-14, 0.999999999997567, 0.999999999995991, 
NA, 0.342775485350427, 0.0806194459381856, 1.05227375864476e-12, 
1.05227375864476e-12, 2.21644445110031e-12, 3.5044109997644e-12, 
NA, 0.599593236726789, 0.838761108123629), dim = c(7L, 4L)), 
    relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS = structure(c(0, 
    0, 0, 0, NA, 0, 0, 0.477617887308721, 0.0266167777965505, 
    0.308323346796742, 0.292600516988572, NA, 0.125000000000644, 
    0.046615618818038, 0.261191056346013, 0.946766444406928, 
    0.649573122214116, 0.616448394398207, NA, 0.749999999999466, 
    0.90676876236405, 0.261191056345266, 0.0266167777965214, 
    0.0421035309891421, 0.090951088613221, NA, 0.12499999999989, 
    0.0466156188179121), dim = c(7L, 4L)), relative_probs_of_each_state_at_branch_top_AT_node_UPPASS = structure(c(7.02704868858637e-13, 
    9.25815984473407e-13, 1.73460095194895e-12, 2.35254814927231e-12, 
    0.25, 8.32240548158005e-13, 9.06794445714384e-13, 0.431710549603533, 
    0.0240584452060588, 0.251901392427722, 0.216078384422621, 
    0.25, 0.112985338561828, 0.0421350517952136, 0.236086079443338, 
    0.855765818076806, 0.530703807120855, 0.455232186573681, 
    0.25, 0.67791203136619, 0.819612604898638, 0.332203370953129, 
    0.120175736717135, 0.217394800451423, 0.328689429003699, 
    0.25, 0.209102630071982, 0.138252343306148), dim = c(7L, 
    4L)), ML_marginal_prob_each_state_at_branch_bottom_below_node = structure(c(0, 
    0, 0, 0, NA, 0, 0, 9.55885872371469e-14, 0.999999999997088, 
    1.02736516865682e-13, 2.3942137600093e-13, NA, 0.0212357703981079, 
    0.0324086154040224, 0.999999999998852, 1.85939277741373e-12, 
    0.999999999999754, 0.999999999999244, NA, 0.757828224601766, 
    0.630413600097194, 1.05227375864171e-12, 1.05227375864171e-12, 
    1.43663791560109e-13, 5.17042461743951e-13, NA, 0.220936005000126, 
    0.337177784498784), dim = c(7L, 4L)), ML_marginal_prob_each_state_at_branch_top_AT_node = structure(c(0, 
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 6.42693382782259e-14, 2.27415872607374e-14, 
    9.55885855108166e-14, 1, 0, 1, 1, 0.757828232249181, 0.630413602214152, 
    1.85939274383416e-12, 0, 0, 0, 0, 0.242171767750754, 0.369586397785825, 
    0.999999999998045), dim = c(7L, 4L)), relative_probs_of_each_state_at_bottom_of_root_branch = c(0, 
    6.42693382782259e-14, 0.757828232249181, 0.242171767750754
    ), total_loglikelihood = -4.48101163256014, inputs = list(
        geogfn = "geog.data", trfn = "tree.newick", abbr = "default", 
        description = "defaults", BioGeoBEARS_model_object = new("BioGeoBEARS_model", 
            params_table = structure(list(type = c("free", "free", 
            "fixed", "fixed", "fixed", "fixed", "fixed", "fixed", 
            "fixed", "3-j", "ysv*2/3", "ysv*1/3", "ysv*1/3", 
            "ysv*1/3", "fixed", "mx01", "mx01", "mx01", "mx01", 
            "fixed", "fixed", "fixed", "fixed"), init = c(0.101055674351941, 
            1e-12, 0, 1, 0, 0, 1, 0, 0, 2.99999, 1.99999, 1, 
            1, 1, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 0.5, 0.1, 
            1, 0), min = c(1e-12, 1e-12, 1e-12, 1e-12, -2.5, 
            -10, -10, -10, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 
            1e-05, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 
            0.005, 0.005, 0.005), max = c(4.999999999999, 4.999999999999, 
            4.999999999999, 0.999999999999, 2.5, 10, 10, 10, 
            2.99999, 3, 2, 1, 1, 1, 0.9999, 0.9999, 0.9999, 0.9999, 
            0.9999, 0.9999, 0.995, 0.995, 0.995), est = c(0.101055674351941, 
            1e-12, 0, 1, 0, 0, 1, 0, 0, 3, 2, 1, 1, 1, 1e-04, 
            1e-04, 1e-04, 1e-04, 1e-04, 0.5, 0.1, 1, 0), note = c("works", 
            "works", "works", "non-stratified only", "works", 
            "works", "works", "works", "works", "works", "works", 
            "works", "works", "works", "works", "works", "works", 
            "works", "works", "no", "yes", "yes", "yes"), desc = c("anagenesis: rate of 'dispersal' (range expansion)", 
            "anagenesis: rate of 'extinction' (range contraction)", 
            "anagenesis: rate of range-switching (i.e. for a standard char.)", 
            "anagenesis: exponent on branch lengths", "exponent on distance (modifies d, j, a)", 
            "exponent on environmental distance (modifies d, j, a)", 
            "exponent on manual dispersal multipliers (modifies d, j, a)", 
            "anagenesis: exponent on extinction risk with area (modifies e)", 
            "cladogenesis: relative per-event weight of jump dispersal", 
            "cladogenesis: y+s+v", "cladogenesis: y+s", "cladogenesis: relative per-event weight of sympatry (range-copying)", 
            "cladogenesis: relative per-event weight of subset speciation", 
            "cladogenesis: relative per-event weight of vicariant speciation", 
            "cladogenesis: controls range size of smaller daughter", 
            "cladogenesis: controls range size of smaller daughter", 
            "cladogenesis: controls range size of smaller daughter", 
            "cladogenesis: controls range size of smaller daughter", 
            "cladogenesis: controls range size of smaller daughter", 
            "root: controls range size probabilities of root", 
            "mean frequency of truly sampling OTU of interest", 
            "detection probability per true sample of OTU of interest", 
            "false detection of OTU probability per true taphonomic control sample"
            )), row.names = c("d", "e", "a", "b", "x", "n", "w", 
            "u", "j", "ysv", "ys", "y", "s", "v", "mx01", "mx01j", 
            "mx01y", "mx01s", "mx01v", "mx01r", "mf", "dp", "fdp"
            ), class = "data.frame")), timesfn = NA, distsfn = NA, 
        dispersal_multipliers_fn = NA, area_of_areas_fn = NA, 
        areas_allowed_fn = NA, areas_adjacency_fn = NA, detects_fn = NA, 
        controls_fn = NA, max_range_size = 2, force_sparse = FALSE, 
        use_detection_model = FALSE, print_optim = TRUE, printlevel = 0, 
        on_NaN_error = -1e+50, wd = "/GitHub/PhyBEARS.jl/test/apes_SSE", 
        num_cores_to_use = 1, cluster_already_open = FALSE, use_optimx = TRUE, 
        rescale_params = FALSE, return_condlikes_table = TRUE, 
        calc_TTL_loglike_from_condlikes_table = TRUE, calc_ancprobs = TRUE, 
        speedup = TRUE, include_null_range = TRUE, useAmbiguities = FALSE, 
        min_branchlength = 1e-06, allow_null_tips = FALSE, all_geog_states_list_usually_inferred_from_areas_maxareas = list(
            NA, 0L, 1L, 0:1)), outputs = new("BioGeoBEARS_model", 
        params_table = structure(list(type = c("free", "free", 
        "fixed", "fixed", "fixed", "fixed", "fixed", "fixed", 
        "fixed", "3-j", "ysv*2/3", "ysv*1/3", "ysv*1/3", "ysv*1/3", 
        "fixed", "mx01", "mx01", "mx01", "mx01", "fixed", "fixed", 
        "fixed", "fixed"), init = c(0.101055674351941, 1e-12, 
        0, 1, 0, 0, 1, 0, 0, 2.99999, 1.99999, 1, 1, 1, 1e-04, 
        1e-04, 1e-04, 1e-04, 1e-04, 0.5, 0.1, 1, 0), min = c(1e-12, 
        1e-12, 1e-12, 1e-12, -2.5, -10, -10, -10, 1e-05, 1e-05, 
        1e-05, 1e-05, 1e-05, 1e-05, 1e-04, 1e-04, 1e-04, 1e-04, 
        1e-04, 1e-04, 0.005, 0.005, 0.005), max = c(4.999999999999, 
        4.999999999999, 4.999999999999, 0.999999999999, 2.5, 
        10, 10, 10, 2.99999, 3, 2, 1, 1, 1, 0.9999, 0.9999, 0.9999, 
        0.9999, 0.9999, 0.9999, 0.995, 0.995, 0.995), est = c(0.101055674351941, 
        1e-12, 0, 1, 0, 0, 1, 0, 0, 3, 2, 1, 1, 1, 1e-04, 1e-04, 
        1e-04, 1e-04, 1e-04, 0.5, 0.1, 1, 0), note = c("works", 
        "works", "works", "non-stratified only", "works", "works", 
        "works", "works", "works", "works", "works", "works", 
        "works", "works", "works", "works", "works", "works", 
        "works", "no", "yes", "yes", "yes"), desc = c("anagenesis: rate of 'dispersal' (range expansion)", 
        "anagenesis: rate of 'extinction' (range contraction)", 
        "anagenesis: rate of range-switching (i.e. for a standard char.)", 
        "anagenesis: exponent on branch lengths", "exponent on distance (modifies d, j, a)", 
        "exponent on environmental distance (modifies d, j, a)", 
        "exponent on manual dispersal multipliers (modifies d, j, a)", 
        "anagenesis: exponent on extinction risk with area (modifies e)", 
        "cladogenesis: relative per-event weight of jump dispersal", 
        "cladogenesis: y+s+v", "cladogenesis: y+s", "cladogenesis: relative per-event weight of sympatry (range-copying)", 
        "cladogenesis: relative per-event weight of subset speciation", 
        "cladogenesis: relative per-event weight of vicariant speciation", 
        "cladogenesis: controls range size of smaller daughter", 
        "cladogenesis: controls range size of smaller daughter", 
        "cladogenesis: controls range size of smaller daughter", 
        "cladogenesis: controls range size of smaller daughter", 
        "cladogenesis: controls range size of smaller daughter", 
        "root: controls range size probabilities of root", "mean frequency of truly sampling OTU of interest", 
        "detection probability per true sample of OTU of interest", 
        "false detection of OTU probability per true taphonomic control sample"
        )), row.names = c("d", "e", "a", "b", "x", "n", "w", 
        "u", "j", "ysv", "ys", "y", "s", "v", "mx01", "mx01j", 
        "mx01y", "mx01s", "mx01v", "mx01r", "mf", "dp", "fdp"
        ), class = "data.frame")), optim_result = structure(list(
        p1 = 0.101055674351941, p2 = 1e-12, value = -4.48101163256014, 
        fevals = 39, gevals = NA_real_, niter = NA_real_, convcode = 0, 
        kkt1 = FALSE, kkt2 = TRUE, xtime = 0.354), row.names = "bobyqa", class = c("optimx", 
    "data.frame"), details = structure(list("bobyqa", c(0.0313402484276759, 
    0), structure(c(-65.7496196774338, 0, 0, 0), dim = c(2L, 
    2L)), c(-65.7496196774338, 0), "none"), dim = c(1L, 5L), dimnames = list(
        "bobyqa", c("method", "ngatend", "nhatend", "hev", "message"
        ))), maximize = TRUE, npar = 2L, follow.on = FALSE)), class = "calc_loglike_sp_results")


numstates = 4
tmpres = get_Qmat_COOmat_from_res(resDEC, numstates=numstates, include_null_range=TRUE, max_range_size=NULL, timeperiod_i=1)

probs_ancstate = rep(0.25, numstates)
COO_weights_columnar = tmpres$COO_weights_columnar
include_null_range = TRUE
left_branch_downpass_likes = rep(1, numstates)
right_branch_downpass_likes = rep(1, numstates)
Rsp_rowsums = tmpres$Rsp_rowsums


clado_table = cbind(tmpres$COO_weights_columnar[[1]]+1, tmpres$COO_weights_columnar[[2]]+1, tmpres$COO_weights_columnar[[3]]+1, tmpres$COO_weights_columnar[[4]])

probs = rep(0, nrow(clado_table))
for (i in 1:length(probs))
	probs[i] = clado_table[i,4] / Rsp_rowsums[clado_table[i,1]]
end
clado_table = cbind(clado_table, probs)
clado_table_df = adf2(clado_table)
names(clado_table_df) = c("i", "j", "k", "wt", "prob")
clado_table_df

calc_uppass_probs_new2(probs_ancstate, COO_weights_columnar, numstates, include_null_range = include_null_range, 
    left_branch_downpass_likes = NULL, right_branch_downpass_likes = NULL, 
    Rsp_rowsums = NULL)

	} # END calc_uppass_probs_new2_example <- function()
	
	
	
	



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





