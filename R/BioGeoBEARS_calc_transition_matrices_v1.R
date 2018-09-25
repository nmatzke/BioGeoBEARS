spPmat_inputs_from_BioGeoBEARS_model_object <- function(BioGeoBEARS_run_object, states_list, dispersal_multipliers_matrix=NULL)
	{
	BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
	
	# Get the dispersal_multipliers_matrix, if needed
	if (is.null(dispersal_multipliers_matrix))
		{
		dispersal_multipliers_matrix = dispersal_multipliers_matrix_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object)
		}
	
	
	#######################################################
	# Cladogenic model
	#######################################################
	j = BioGeoBEARS_model_object@params_table["j","est"]
	ysv = BioGeoBEARS_model_object@params_table["ysv","est"]
	ys = BioGeoBEARS_model_object@params_table["ys","est"]
	v = BioGeoBEARS_model_object@params_table["v","est"]
	y = BioGeoBEARS_model_object@params_table["y","est"]
	s = BioGeoBEARS_model_object@params_table["s","est"]
	sum_SPweights = y + s + j + v


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

	# Note that this gets the dispersal multipliers matrix, which is applied to 
	# e.g. the j events, NOT the dmat_times_d above which is d*dispersal_multipliers_matrix
	dmat = dispersal_multipliers_matrix
	spPmat_inputs$dmat = dmat

	states_indices = states_list
	
	if (BioGeoBEARS_run_object$include_null_range == TRUE)
		{
		states_indices[1] = NULL	# shorten the states_indices by 1 (cutting the null range state from the speciation matrix)
		}
	spPmat_inputs$l = states_indices
	spPmat_inputs$s = s
	spPmat_inputs$v = v
	spPmat_inputs$j = j
	spPmat_inputs$y = y
	spPmat_inputs$maxent01s_param = maxent01s_param
	spPmat_inputs$maxent01v_param = maxent01v_param
	spPmat_inputs$maxent01j_param = maxent01j_param
	spPmat_inputs$maxent01y_param = maxent01y_param
	
	return(spPmat_inputs)
	}


dispersal_multipliers_matrix_from_BioGeoBEARS_run_object <- function(BioGeoBEARS_run_object)
	{
	BioGeoBEARS_model_object = BioGeoBEARS_run_object
	
	# Get geographic ranges at tips
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=np(BioGeoBEARS_run_object$geogfn))
	# Get the list of geographic areas
	areas = getareas_from_tipranges_object(tipranges)
	areas_list = seq(0, length(areas)-1, 1)		# 0-base indexes

	# Put the parameters into the BioGeoBEARS_model_object, so that they can be universally read out
	# into any function
	BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
	
	# Update linked parameters
	BioGeoBEARS_model_object = calc_linked_params_BioGeoBEARS_model_object(BioGeoBEARS_model_object)
	
	# Set the dispersal and extinction rate
	d = BioGeoBEARS_model_object@params_table["d","est"]
	e = BioGeoBEARS_model_object@params_table["e","est"]
	a = BioGeoBEARS_model_object@params_table["a","est"]


	#######################################################
	#######################################################
	# Do distance-dependence and dispersal multipliers matrix
	#######################################################
	#######################################################
	# Equal dispersal in all directions (unconstrained)
	# Equal extinction probability for all areas
	areas = areas_list
	
	# If there is a distance matrix, use the first one 
	# (non-stratified analysis, here)
	if ( (is.null(BioGeoBEARS_run_object$list_of_distances_mats) == FALSE))
		{
		distances_mat = BioGeoBEARS_run_object$list_of_distances_mats[[1]]
		} else {
		# Default is all areas effectively equidistant
		distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))
		}

	# Get the exponent on distance, apply to distances matrix
	x = BioGeoBEARS_model_object@params_table["x","est"]
	dispersal_multipliers_matrix = distances_mat ^ x

	
	# Environmental distances
	if ( (is.null(BioGeoBEARS_run_object$list_of_envdistances_mats) == FALSE))
		{
		envdistances_mat = BioGeoBEARS_run_object$list_of_envdistances_mats[[1]]
		} else {
		# Default is all areas effectively equidistant
		envdistances_mat = matrix(1, nrow=length(areas), ncol=length(areas))
		}

	# Get the exponent on environmental distance, apply to distances matrix
	n = BioGeoBEARS_model_object@params_table["n","est"]
	dispersal_multipliers_matrix = dispersal_multipliers_matrix * (envdistances_mat ^ n)


	# Apply manual dispersal multipliers, if any
	# If there is a manual dispersal multipliers matrix, use the first one 
	# (non-stratified analysis, here)
	if ( (is.null(BioGeoBEARS_run_object$list_of_dispersal_multipliers_mats) == FALSE))
		{
		manual_dispersal_multipliers_matrix = BioGeoBEARS_run_object$list_of_dispersal_multipliers_mats[[1]]
		} else {
		# Default is all areas effectively equidistant
		manual_dispersal_multipliers_matrix = matrix(1, nrow=length(areas), ncol=length(areas))
		}

	# Get the exponent on manual dispersal multipliers
	w = BioGeoBEARS_model_object@params_table["w","est"]

	# Apply element-wise
	dispersal_multipliers_matrix = dispersal_multipliers_matrix * manual_dispersal_multipliers_matrix ^ w

	return(dispersal_multipliers_matrix)
	}


spPmat_inputs_to_COO_weights_columnar <- function(spPmat_inputs, cppSpMethod=3, numstates_in_cladogenesis_matrix, printmat=FALSE, printlevel=0, m=NULL, include_null_range=TRUE, jts_matrix=NULL)
	{
	defaults='
	cppSpMethod=3
	printmat=FALSE
	printlevel=0
	'
	
	
	# Fix "l" (which is states_indices, i.e. a list of lists of state_indices)
	# so that there is no "null geographic range", i.e. no "_", no "-", no "NA"
	#if ( is.null(spPmat_inputs)==FALSE )
	#	{
	spPmat_inputs$l[spPmat_inputs$l == c("_")] = NULL
	spPmat_inputs$l[spPmat_inputs$l == c("-")] = NULL
	spPmat_inputs$l[spPmat_inputs$l == c("-1")] = NULL
	#spPmat_inputs$l[spPmat_inputs$l == c(-1)] = NULL
	#	}

	
	#if ( is.null(spPmat_inputs)==FALSE )
	#	{
	# Calculate the rowsums (for input into rcpp_calc_anclikes_sp()
	# Actually, just do this ONCE

	# (above) Fix "l" (which is states_indices, i.e. a list of lists of state_indices)
	# so that there is no "null geographic range", i.e. no "_", no "-", no "NA"
	l = spPmat_inputs$l		# states_indices
	
	states_indices = spPmat_inputs$l
	
	s = spPmat_inputs$s
	v = spPmat_inputs$v
	j = spPmat_inputs$j
	y = spPmat_inputs$y

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


	tmpca_1 = rep(1, numstates_in_cladogenesis_matrix)
	tmpcb_1 = rep(1, numstates_in_cladogenesis_matrix)

	# Print the matrix to screen from C++
	printmat = FALSE

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
			txt = "STOP ERROR in spPmat_inputs_to_COO_weights_columnar(): if m is not NULL, only option cppSpMethod==3 can be used. You have cppSpMethod==2."
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			stop(txt)
			} # END if (is.null(m) == FALSE)		
				
		COO_probs_list = rcpp_calc_anclikes_sp_COOprobs(Rcpp_leftprobs=tmpca_1, Rcpp_rightprobs=tmpcb_1, l=l, s=s, v=v, j=j, y=y, dmat=dmat, maxent01s=maxent01s, maxent01v=maxent01v, maxent01j=maxent01j, maxent01y=maxent01y, max_minsize_as_function_of_ancsize=max_minsize_as_function_of_ancsize, printmat=printmat)

		# Sum the probabilities (list [[3]]) in each row ([[3]][[1]] list of probs
		# through [[3]][[15]] list of probs)
		Rsp_rowsums = sapply(X=COO_probs_list[[3]], FUN=sum)

		COOmat_Rsp_rowsums = NULL
		COOmat_Rsp_rowsums$COO_probs_list = COO_probs_list
		COOmat_Rsp_rowsums$Rsp_rowsums = Rsp_rowsums
		} # END if (cppSpMethod == 2)

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
		
		# NJM MOD HERE: , m=m
		# mod done: 2017-08-15

		# m_null_range Is the null range included in the state space in the general 
		#' analysis? (The function needs to know this, when there are traits, to index
		#' the state space correctly.)
		
		Rcpp_leftprobs=tmpca_1
		Rcpp_rightprobs=tmpcb_1
		m_null_range=include_null_range
		
		
		COO_weights_columnar = rcpp_calc_anclikes_sp_COOweights_faster(Rcpp_leftprobs=tmpca_1, Rcpp_rightprobs=tmpcb_1, l=l, s=s, v=v, j=j, y=y, dmat=dmat, maxent01s=maxent01s, maxent01v=maxent01v, maxent01j=maxent01j, maxent01y=maxent01y, max_minsize_as_function_of_ancsize=max_minsize_as_function_of_ancsize, printmat=printmat, m=m, m_null_range=include_null_range, jts_matrix=jts_matrix)
	
		# combine with C++ function
		# This causes an error with spPmat=NULL; spPmat_inputs=NULL; use_cpp=TRUE; sparse=FALSE
		# i.e. gives 16 states with a 0 on the end, rather than 15 states
		#Rsp_rowsums = rcpp_calc_rowsums_for_COOweights_columnar(COO_weights_columnar=COO_weights_columnar, numstates=numstates)
	
		# This gives 15 states
		Rsp_rowsums = rcpp_calc_rowsums_for_COOweights_columnar(COO_weights_columnar=COO_weights_columnar, numstates=numstates_in_cladogenesis_matrix)

		COOmat_Rsp_rowsums = NULL
		COOmat_Rsp_rowsums$COO_weights_columnar = COO_weights_columnar
		COOmat_Rsp_rowsums$Rsp_rowsums = Rsp_rowsums
		} # END if (cppSpMethod == 3)
	#	}
	return(COOmat_Rsp_rowsums)
	}