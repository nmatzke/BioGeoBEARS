
#######################################################
# Create a little fake function
#######################################################
little_Qmats_to_big_Qmat <- function(params=NULL, BioGeoBEARS_run_object)
	{
	
	# Input the current parameters
	if (is.null(params))
		{
		# do nothing; use current BioGeoBEARS_run_object
		} else {
		BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
		BioGeoBEARS_model_object = params_into_BioGeoBEARS_model_object(BioGeoBEARS_model_object=BioGeoBEARS_model_object, params=params)
	
		# Update linked parameters
		BioGeoBEARS_model_object = calc_linked_params_BioGeoBEARS_model_object(BioGeoBEARS_model_object)
		
		# Store back in the run object
		BioGeoBEARS_run_object$BioGeoBEARS_model_object = BioGeoBEARS_model_object
		} # END if (is.null(params))
	
	
	if (is.null(BioGeoBEARS_run_object$numstates))
		{
		txt = "STOP ERROR in the function 'Mk_Qmat', typically specified as a custom Qmat-generating function in BioGeoBEARS_run_object$custom_Qmat_fn_text. This function requires that the user specify BioGeoBEARS_run_object$numstates ahead of time."
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		} else {
		numstates = BioGeoBEARS_run_object$numstates
		}# END default_Qmat <- function(BioGeoBEARS_run_object)
	
	#states_list = rcpp_areas_list_to_states_list(areas=areanames, maxareas=max_range_size, include_null_range=BioGeoBEARS_run_object$include_null_range)
	states_list = as.list(1:numstates)
	states_list
	
	
	# Just put in "a"
	a = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","est"]
	
	Qmat = matrix(data=a, nrow=numstates, ncol=numstates)
	diag(Qmat) = 0
	diag(Qmat) = -rowSums(Qmat)
	
	return(Qmat)
	}




#######################################################
# Create a little fake function
#######################################################
Mk_Qmat <- function(params=NULL, BioGeoBEARS_run_object)
	{
	
	# Input the current parameters
	if (is.null(params))
		{
		# do nothing; use current BioGeoBEARS_run_object
		} else {
		BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
		BioGeoBEARS_model_object = params_into_BioGeoBEARS_model_object(BioGeoBEARS_model_object=BioGeoBEARS_model_object, params=params)
	
		# Update linked parameters
		BioGeoBEARS_model_object = calc_linked_params_BioGeoBEARS_model_object(BioGeoBEARS_model_object)
		
		# Store back in the run object
		BioGeoBEARS_run_object$BioGeoBEARS_model_object = BioGeoBEARS_model_object
		} # END if (is.null(params))
	
	
	if (is.null(BioGeoBEARS_run_object$numstates))
		{
		txt = "STOP ERROR in the function 'Mk_Qmat', typically specified as a custom Qmat-generating function in BioGeoBEARS_run_object$custom_Qmat_fn_text. This function requires that the user specify BioGeoBEARS_run_object$numstates ahead of time."
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		} else {
		numstates = BioGeoBEARS_run_object$numstates
		}# END default_Qmat <- function(BioGeoBEARS_run_object)
	
	#states_list = rcpp_areas_list_to_states_list(areas=areanames, maxareas=max_range_size, include_null_range=BioGeoBEARS_run_object$include_null_range)
	states_list = as.list(1:numstates)
	states_list
	
	
	# Just put in "a"
	a = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","est"]
	
	Qmat = matrix(data=a, nrow=numstates, ncol=numstates)
	diag(Qmat) = 0
	diag(Qmat) = -rowSums(Qmat)
	
	return(Qmat)
	}


# Make the chromosome add/loss Q matrix
Qmat_chromo_add_subtract <- function(BioGeoBEARS_model_object=define_ChromoBEARS_model_object(), maxnum_chromosomes=12, minnum_chromosomes=2, chrcount_Qmat=NULL)
	{
	defaults='
	BioGeoBEARS_model_object=define_ChromoBEARS_model_object()
	maxnum_chromosomes = 12
	minnum_chromosomes = 2
	'
	states_list_chrcount = seq(minnum_chromosomes, maxnum_chromosomes, 1)
	numstates_chrcount = length(states_list_chrcount)

	if (is.null(chrcount_Qmat))
		{
		chrcount_Qmat = matrix(0, nrow=numstates_chrcount, ncol=numstates_chrcount)
		colnames(chrcount_Qmat) = states_list_chrcount
		rownames(chrcount_Qmat) = states_list_chrcount
		} # END if (is.null(chrcount_Qmat))

	# lamda: rate of adding 1 chromosome
	anc_statenums = 1:numstates_chrcount
	dec_statenums = 1 + anc_statenums
	keepTFanc = ((anc_statenums > 0) + (anc_statenums <= numstates_chrcount)) == 2
	keepTFdec = ((dec_statenums > 0) + (dec_statenums <= numstates_chrcount)) == 2
	keepTF = (keepTFanc + keepTFdec) == 2
	anc_statenums1 = anc_statenums[keepTF]
	dec_statenums1 = dec_statenums[keepTF]
	anc_statenums1
	dec_statenums1

	# mu: rate of subtracting 1 chromosome
	anc_statenums2 = dec_statenums[keepTF]
	dec_statenums2 = anc_statenums[keepTF]

	
	BioGeoBEARS_model_object
	lambda = BioGeoBEARS_model_object@params_table["lambda","est"]
	mu = BioGeoBEARS_model_object@params_table["mu","est"]

	for (i in 1:length(anc_statenums1))
		{
		chrcount_Qmat[anc_statenums1[i], dec_statenums1[i]] = lambda + chrcount_Qmat[anc_statenums1[i], dec_statenums1[i]]
		chrcount_Qmat[anc_statenums2[i], dec_statenums2[i]] = mu + chrcount_Qmat[anc_statenums2[i], dec_statenums2[i]]
		} # END for (i in 1:length(anc_statenums))

	diag(chrcount_Qmat) = 0

	return(chrcount_Qmat)
	} # END Qmat_chromo_add_subtract <- function(maxnum_chromosomes=12, minnum_chromosomes=2)



# Make the odd bridge increase, +1x to ploidy
Qmat_chromo_alpha <- function(BioGeoBEARS_model_object=define_ChromoBEARS_model_object(), ploidy=2, maxnum_chromosomes=12, minnum_chromosomes=2, chrcount_Qmat_alpha=NULL)
	{
	defaults='
	BioGeoBEARS_model_object=define_ChromoBEARS_model_object()
	ploidy=2
	maxnum_chromosomes = 12
	minnum_chromosomes = 2
	'
	states_list_chrcount = seq(minnum_chromosomes, maxnum_chromosomes, 1)
	numstates_chrcount = length(states_list_chrcount)

	if (is.null(chrcount_Qmat_alpha))
		{
		chrcount_Qmat_alpha = matrix(0, nrow=numstates_chrcount, ncol=numstates_chrcount)
		colnames(chrcount_Qmat_alpha) = states_list_chrcount
		rownames(chrcount_Qmat_alpha) = states_list_chrcount
		} # END if (is.null(chrcount_Qmat_alpha))

	# Figure out the number of chromosomes in the 1x ploidy state
	# Assumed ploidy: ploidy
	# If that is the ploidy, then 1x is:
	numchrom_in_1x_ploidy = floor( (1/ploidy) * states_list_chrcount )
	numchrom_in_1x_ploidy

	dec_states_list_chrcount_add1x = states_list_chrcount + numchrom_in_1x_ploidy
	
	##########################################
	# Adding 1x (odd bridge increase)
	##########################################
	# What states do these counts correspond to?
	anc_statenums = 1:numstates_chrcount
	dec_statenums = match(x=dec_states_list_chrcount_add1x, table=states_list_chrcount)
	
	# Filter
	keepTF = is.na(dec_statenums) == FALSE
	anc_statenums_a1x = anc_statenums[keepTF]
	dec_statenums_a1x = dec_statenums[keepTF]


	# alpha: rate of odd bridge increase (+1x)
	alpha = BioGeoBEARS_model_object@params_table["alpha","est"]

	for (i in 1:length(anc_statenums_a1x))
		{
		chrcount_Qmat_alpha[anc_statenums_a1x[i], dec_statenums_a1x[i]] = alpha + chrcount_Qmat_alpha[anc_statenums_a1x[i], dec_statenums_a1x[i]]
		} # END for (i in 1:length(anc_statenums))

	diag(chrcount_Qmat_alpha) = 0

	return(chrcount_Qmat_alpha)
	} # END Qmat_chromo_alpha <- function(BioGeoBEARS_model_object=define_ChromoBEARS_model_object(), ploidy=2, maxnum_chromosomes=12, minnum_chromosomes=2, chrcount_Qmat_alpha=NULL)









# Make the odd bridge decrease, -1x to ploidy
Qmat_chromo_beta <- function(BioGeoBEARS_model_object=define_ChromoBEARS_model_object(), ploidy=2, maxnum_chromosomes=12, minnum_chromosomes=2, chrcount_Qmat_beta=NULL)
	{
	defaults='
	BioGeoBEARS_model_object=define_ChromoBEARS_model_object()
	ploidy=2
	maxnum_chromosomes = 12
	minnum_chromosomes = 2
	'
	states_list_chrcount = seq(minnum_chromosomes, maxnum_chromosomes, 1)
	numstates_chrcount = length(states_list_chrcount)

	if (is.null(chrcount_Qmat_beta))
		{
		chrcount_Qmat_beta = matrix(0, nrow=numstates_chrcount, ncol=numstates_chrcount)
		colnames(chrcount_Qmat_beta) = states_list_chrcount
		rownames(chrcount_Qmat_beta) = states_list_chrcount
		} # END if (is.null(chrcount_Qmat_beta))

	# Figure out the number of chromosomes in the 1x ploidy state
	# Assumed ploidy: ploidy
	# If that is the ploidy, then 1x is:
	numchrom_in_1x_ploidy = floor( (1/ploidy) * states_list_chrcount )
	numchrom_in_1x_ploidy

	dec_states_list_chrcount_sub1x = states_list_chrcount - numchrom_in_1x_ploidy
	
	##########################################
	# Subtracting 1x (odd bridge decrease)
	##########################################
	# What states do these counts correspond to?
	anc_statenums = 1:numstates_chrcount
	dec_statenums = match(x=dec_states_list_chrcount_sub1x, table=states_list_chrcount)
	
	# Filter
	keepTF = is.na(dec_statenums) == FALSE
	anc_statenums_s1x = anc_statenums[keepTF]
	dec_statenums_s1x = dec_statenums[keepTF]


	# beta: rate of odd bridge decrease (+1x)
	beta = BioGeoBEARS_model_object@params_table["beta","est"]

	for (i in 1:length(anc_statenums_s1x))
		{
		chrcount_Qmat_beta[anc_statenums_s1x[i], dec_statenums_s1x[i]] = beta + chrcount_Qmat_beta[anc_statenums_s1x[i], dec_statenums_s1x[i]]
		} # END for (i in 1:length(anc_statenums))

	diag(chrcount_Qmat_beta) = 0

	return(chrcount_Qmat_beta)
	} # END Qmat_chromo_beta <- function(BioGeoBEARS_model_object=define_ChromoBEARS_model_object(), ploidy=2, maxnum_chromosomes=12, minnum_chromosomes=2, chrcount_Qmat_beta=NULL)




# Make the even bridge increase, +2x to ploidy
Qmat_chromo_ep <- function(BioGeoBEARS_model_object=define_ChromoBEARS_model_object(), ploidy=2, maxnum_chromosomes=12, minnum_chromosomes=2, chrcount_Qmat_ep=NULL)
	{
	defaults='
	BioGeoBEARS_model_object=define_ChromoBEARS_model_object()
	ploidy=2
	maxnum_chromosomes = 12
	minnum_chromosomes = 2
	'
	states_list_chrcount = seq(minnum_chromosomes, maxnum_chromosomes, 1)
	numstates_chrcount = length(states_list_chrcount)

	if (is.null(chrcount_Qmat_ep))
		{
		chrcount_Qmat_ep = matrix(0, nrow=numstates_chrcount, ncol=numstates_chrcount)
		colnames(chrcount_Qmat_ep) = states_list_chrcount
		rownames(chrcount_Qmat_ep) = states_list_chrcount
		} # END if (is.null(chrcount_Qmat_ep))

	# Figure out the number of chromosomes in the 2x ploidy state
	# Assumed ploidy: ploidy
	# If that is the ploidy, then 2x is:
	numchrom_in_2x_ploidy = floor( (2/ploidy) * states_list_chrcount )
	numchrom_in_2x_ploidy

	dec_states_list_chrcount_add2x = states_list_chrcount + numchrom_in_2x_ploidy
	
	##########################################
	# Adding 2x (even bridge increase)
	##########################################
	# What states do these counts correspond to?
	anc_statenums = 1:numstates_chrcount
	dec_statenums = match(x=dec_states_list_chrcount_add2x, table=states_list_chrcount)
	
	# Filter
	keepTF = is.na(dec_statenums) == FALSE
	anc_statenums_a2x = anc_statenums[keepTF]
	dec_statenums_a2x = dec_statenums[keepTF]


	# epsilon: rate of even bridge increase (+2x)
	epsilon = BioGeoBEARS_model_object@params_table["epsilon","est"]

	for (i in 1:length(anc_statenums_a2x))
		{
		chrcount_Qmat_ep[anc_statenums_a2x[i], dec_statenums_a2x[i]] = epsilon + chrcount_Qmat_ep[anc_statenums_a2x[i], dec_statenums_a2x[i]]
		} # END for (i in 1:length(anc_statenums))

	diag(chrcount_Qmat_ep) = 0

	return(chrcount_Qmat_ep)
	} # END Qmat_chromo_ep <- function(BioGeoBEARS_model_object=define_ChromoBEARS_model_object(), ploidy=2, maxnum_chromosomes=12, minnum_chromosomes=2, chrcount_Qmat_ep=NULL)



# Make the rho ("pure duplication") matrix
Qmat_chromo_rho <- function(BioGeoBEARS_model_object=define_ChromoBEARS_model_object(), maxnum_chromosomes=12, minnum_chromosomes=2, chrcount_Qmat_rho=NULL)
	{
	defaults='
	BioGeoBEARS_model_object=define_ChromoBEARS_model_object()
	maxnum_chromosomes = 12
	minnum_chromosomes = 2
	'
	states_list_chrcount = seq(minnum_chromosomes, maxnum_chromosomes, 1)
	numstates_chrcount = length(states_list_chrcount)

	if (is.null(chrcount_Qmat_rho))
		{
		chrcount_Qmat_rho = matrix(0, nrow=numstates_chrcount, ncol=numstates_chrcount)
		colnames(chrcount_Qmat_rho) = states_list_chrcount
		rownames(chrcount_Qmat_rho) = states_list_chrcount
		} # END if (is.null(chrcount_Qmat_rho))

	# rho: rate of pure duplication
	# states resulting from duplication
	dec_states_list_chrcount_dup = states_list_chrcount * 2
	
	# corresponding statenums
	dec_statenums = match(x=dec_states_list_chrcount_dup, table=states_list_chrcount)
	keepTFanc = ((anc_statenums > 0) + (anc_statenums <= numstates_chrcount)) == 2
	keepTFdec = ((dec_statenums > 0) + (dec_statenums <= numstates_chrcount)) == 2
	keepTF = (keepTFanc + keepTFdec) == 2
	anc_statenums1 = anc_statenums[keepTF]
	dec_statenums1 = dec_statenums[keepTF]

	rho = BioGeoBEARS_model_object@params_table["rho","est"]

	for (i in 1:length(anc_statenums1))
		{
		chrcount_Qmat_rho[anc_statenums1[i], dec_statenums1[i]] = rho + chrcount_Qmat_rho[anc_statenums1[i], dec_statenums1[i]]
		} # END for (i in 1:length(anc_statenums1))
	
	diag(chrcount_Qmat_rho) = 0
	return(chrcount_Qmat_rho)
	} # END chrcount_Qmat_rho <- function(BioGeoBEARS_model_object=define_ChromoBEARS_model_object(), maxnum_chromosomes=12, minnum_chromosomes=2, chrcount_Qmat_rho=NULL)







# Make the delta matrix (diploidization, return to 2x, no change in chrcount)
Qmat_chromo_delta <- function(BioGeoBEARS_model_object=define_ChromoBEARS_model_object(), maxnum_chromosomes=12, minnum_chromosomes=2, chrcount_Qmat_delta=NULL)
	{
	defaults='
	BioGeoBEARS_model_object=define_ChromoBEARS_model_object()
	maxnum_chromosomes = 12
	minnum_chromosomes = 2
	'
	states_list_chrcount = seq(minnum_chromosomes, maxnum_chromosomes, 1)
	numstates_chrcount = length(states_list_chrcount)

	if (is.null(chrcount_Qmat_delta))
		{
		chrcount_Qmat_delta = matrix(0, nrow=numstates_chrcount, ncol=numstates_chrcount)
		colnames(chrcount_Qmat_delta) = states_list_chrcount
		rownames(chrcount_Qmat_delta) = states_list_chrcount
		} # END if (is.null(chrcount_Qmat_delta))

	# The chromosome counts stay the same, so just set the diagonal to 
	# delta
	delta = BioGeoBEARS_model_object@params_table["delta","est"]
	
	diag(chrcount_Qmat_delta) = delta
	
	return(chrcount_Qmat_delta)
	} # END Qmat_chromo_delta <- function(BioGeoBEARS_model_object=define_ChromoBEARS_model_object(), maxnum_chromosomes=12, minnum_chromosomes=2, chrcount_Qmat_delta=NULL)









Qmat_ploidy <- function(BioGeoBEARS_model_object=define_ChromoBEARS_model_object(), maxnum_ploidy=4, minnum_ploidy=2, ploidy_Qmat=NULL)
	{
	states_list_ploidy = seq(minnum_ploidy, maxnum_ploidy, 1)
	numstates_ploidy = length(states_list_ploidy)

	if (is.null(ploidy_Qmat))
		{
		ploidy_Qmat = matrix(0, nrow=numstates_ploidy, ncol=numstates_ploidy)
		colnames(ploidy_Qmat) = states_list_ploidy
		rownames(ploidy_Qmat) = states_list_ploidy
		} # END if (is.null(ploidy_Qmat))

	# Load the current parameter values
	alpha = BioGeoBEARS_model_object@params_table["alpha","est"]
	beta = BioGeoBEARS_model_object@params_table["beta","est"]
	epsilon = BioGeoBEARS_model_object@params_table["epsilon","est"]
	delta = BioGeoBEARS_model_object@params_table["delta","est"]
	rho = BioGeoBEARS_model_object@params_table["rho","est"]

	# alpha: +1x ploidy
	anc_statenums = 1:numstates_ploidy
	dec_statenums = 1 + anc_statenums
	keepTFanc = ((anc_statenums > 0) + (anc_statenums <= numstates_ploidy)) == 2
	keepTFdec = ((dec_statenums > 0) + (dec_statenums <= numstates_ploidy)) == 2
	keepTF = (keepTFanc + keepTFdec) == 2
	anc_statenums_alpha = anc_statenums[keepTF]
	dec_statenums_alpha = dec_statenums[keepTF]

	for (i in 1:length(anc_statenums_alpha))
		{
		ploidy_Qmat[anc_statenums_alpha[i], dec_statenums_alpha[i]] = alpha + ploidy_Qmat[anc_statenums_alpha[i], dec_statenums_alpha[i]]
		} # END for (i in 1:length(anc_statenums_alpha))
	#print(ploidy_Qmat)
	
	# beta: -1x ploidy
	anc_statenums = 1:numstates_ploidy
	dec_statenums = anc_statenums - 1
	keepTFanc = ((anc_statenums > 0) + (anc_statenums <= numstates_ploidy)) == 2
	keepTFdec = ((dec_statenums > 0) + (dec_statenums <= numstates_ploidy)) == 2
	keepTF = (keepTFanc + keepTFdec) == 2
	anc_statenums_beta = anc_statenums[keepTF]
	dec_statenums_beta = dec_statenums[keepTF]
	
	for (i in 1:length(anc_statenums_beta))
		{
		ploidy_Qmat[anc_statenums_beta[i], dec_statenums_beta[i]] = beta + ploidy_Qmat[anc_statenums_beta[i], dec_statenums_beta[i]]
		} # END for (i in 1:length(anc_statenums_beta))
	#print(ploidy_Qmat)

	# epsilon: +2x ploidy
	anc_statenums = 1:numstates_ploidy
	dec_statenums = 2 + anc_statenums
	keepTFanc = ((anc_statenums > 0) + (anc_statenums <= numstates_ploidy)) == 2
	keepTFdec = ((dec_statenums > 0) + (dec_statenums <= numstates_ploidy)) == 2
	keepTF = (keepTFanc + keepTFdec) == 2
	anc_statenums_epsilon = anc_statenums[keepTF]
	dec_statenums_epsilon = dec_statenums[keepTF]
	
	for (i in 1:length(anc_statenums_epsilon))
		{
		ploidy_Qmat[anc_statenums_epsilon[i], dec_statenums_epsilon[i]] = epsilon + ploidy_Qmat[anc_statenums_epsilon[i], dec_statenums_epsilon[i]]
		} # END for (i in 1:length(anc_statenums_epsilon))
	#print(ploidy_Qmat)

	# delta: diploidization (return to 2x)
	anc_states_list_ploidy = states_list_ploidy
	dec_states_list_ploidy = rep(2, length(anc_states_list_ploidy))
	anc_statenums = 1:numstates_ploidy
	dec_statenums = match(x=dec_states_list_ploidy, table=states_list_ploidy)
	dec_statenums = dec_statenums[is.na(dec_statenums) == FALSE]
	anc_statenums = anc_statenums[is.na(dec_statenums) == FALSE]
	keepTFanc = ((anc_statenums > 0) + (anc_statenums <= numstates_ploidy)) == 2
	keepTFdec = ((dec_statenums > 0) + (dec_statenums <= numstates_ploidy)) == 2
	keepTF = (keepTFanc + keepTFdec) == 2
	anc_statenums_delta = anc_statenums[keepTF]
	dec_statenums_delta = dec_statenums[keepTF]
	
	for (i in 1:length(anc_statenums_delta))
		{
		ploidy_Qmat[anc_statenums_delta[i], dec_statenums_delta[i]] = delta + ploidy_Qmat[anc_statenums_delta[i], dec_statenums_delta[i]]
		} # END for (i in 1:length(anc_statenums_delta))
	#print(ploidy_Qmat)

	# rho: x2 ploidy
	anc_states_list_ploidy = states_list_ploidy
	dec_states_list_ploidy = 2 * anc_states_list_ploidy
	anc_statenums = 1:numstates_ploidy
	dec_statenums = match(x=dec_states_list_ploidy, table=states_list_ploidy)
	dec_statenums = dec_statenums[is.na(dec_statenums) == FALSE]
	anc_statenums = anc_statenums[is.na(dec_statenums) == FALSE]
	keepTFanc = ((anc_statenums > 0) + (anc_statenums <= numstates_ploidy)) == 2
	keepTFdec = ((dec_statenums > 0) + (dec_statenums <= numstates_ploidy)) == 2
	keepTF = (keepTFanc + keepTFdec) == 2
	anc_statenums_rho = anc_statenums[keepTF]
	dec_statenums_rho = dec_statenums[keepTF]
	
	for (i in 1:length(anc_statenums_rho))
		{
		ploidy_Qmat[anc_statenums_rho[i], dec_statenums_rho[i]] = rho + ploidy_Qmat[anc_statenums_rho[i], dec_statenums_rho[i]]
		} # END for (i in 1:length(anc_statenums_rho))
	#print(ploidy_Qmat)
	
	diag(ploidy_Qmat) = 0
	
	return(ploidy_Qmat)
	} # END Qmat_ploidy <- function(BioGeoBEARS_model_object=define_ChromoBEARS_model_object(), maxnum_ploidy=4, minnum_ploidy=2, ploidy_Qmat=NULL)




Qmat_chromoploid <- function(params=NULL, BioGeoBEARS_model_object=define_ChromoBEARS_model_object(), maxnum_ploidy=4, minnum_ploidy=2, minnum_chromosomes=2, maxnum_chromosomes=12)
	{
	# Setup
	defaults='
	params=NULL
	BioGeoBEARS_model_object=define_ChromoBEARS_model_object()
	minnum_ploidy=2
	maxnum_ploidy=4
	minnum_chromosomes=2
	maxnum_chromosomes=12
	'

	# Input the current parameters
	if (is.null(params))
		{
		# do nothing; use current BioGeoBEARS_run_object
		} else {
		BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
		BioGeoBEARS_model_object = params_into_BioGeoBEARS_model_object(BioGeoBEARS_model_object=BioGeoBEARS_model_object, params=params)
	
		# Update linked parameters
		BioGeoBEARS_model_object = calc_linked_params_BioGeoBEARS_model_object(BioGeoBEARS_model_object)
		
		# Store back in the run object
		BioGeoBEARS_run_object$BioGeoBEARS_model_object = BioGeoBEARS_model_object
		} # END if (is.null(params))
	
	# Setup
	states_list_ploidy = seq(minnum_ploidy, maxnum_ploidy, 1)
	numstates_ploidy = length(states_list_ploidy)
	states_list_chrcount = seq(minnum_chromosomes, maxnum_chromosomes, 1)
	numstates_chrcount = length(states_list_chrcount)

	numstates = numstates_ploidy * numstates_chrcount
	Qmat = matrix(data=0, nrow=numstates, ncol=numstates)

	# Pre-define the ones that are not ancploidy-specific
	# Single-chromosome add-subtract (lambda/mu)
	chrcount_Qmat = Qmat_chromo_add_subtract(BioGeoBEARS_model_object, maxnum_chromosomes, minnum_chromosomes, chrcount_Qmat=NULL)
	# diploidization, return to 2x, no change in chrcount
	chrcount_Qmat_delta = Qmat_chromo_delta(BioGeoBEARS_model_object, maxnum_chromosomes, minnum_chromosomes, chrcount_Qmat_delta=NULL)
	# pure duplication, double ploidy
	chrcount_Qmat_rho = Qmat_chromo_rho(BioGeoBEARS_model_object, maxnum_chromosomes, minnum_chromosomes, chrcount_Qmat_rho=NULL)


	# Iterate through the smaller ploidy matrix
	# i = from
	# j = to
	# Pick a row (ancestral state), then go through the cells
	for (i in 1:numstates_ploidy)
		{
		# Go through the cells (columns of descendant states)
		for (j in 1:numstates_ploidy)
			{
			rowstart = 1 + ((i-1) * numstates_chrcount)
			rowend = 0 + ((i-0) * numstates_chrcount)
			rows = rowstart:rowend
		
			colstart = 1 + ((j-1) * numstates_chrcount)
			colend = 0 + ((j-0) * numstates_chrcount)
			cols = colstart:colend
		
			# Figure out which ploidy transition you are in
			ancploidy = states_list_ploidy[i]
			decploidy = states_list_ploidy[j]
		
			# Where there is no transition in ploidy, you
			# can still have the add one chromosome, subtract
			# one chromosome system
			if (i == j)
				{
				# Add the single-chromosome-change events to
				# the Qmat
				Qmat[rows,cols] = Qmat[rows,cols] + chrcount_Qmat
			
				# Since none of the other events apply, go to 
				# next iteration
				next()
				}
		
			# alpha: +1x ploidy
			if (decploidy == ancploidy+1)
				{
				# Calculate the transitions under this ploidy
				chrcount_Qmat_alpha = Qmat_chromo_alpha(BioGeoBEARS_model_object, ploidy=ancploidy, maxnum_chromosomes, minnum_chromosomes, chrcount_Qmat_alpha=NULL)
			
				# Add to Qmat
				Qmat[rows,cols] = Qmat[rows,cols] + chrcount_Qmat_alpha
				} # END if (decploidy == ancploidy+1)

			# beta: -1x ploidy
			if (decploidy == ancploidy-1)
				{
				# Calculate the transitions under this ploidy
				chrcount_Qmat_beta = Qmat_chromo_beta(BioGeoBEARS_model_object, ploidy=ancploidy, maxnum_chromosomes, minnum_chromosomes, chrcount_Qmat_beta=NULL)
			
				# Add to Qmat
				Qmat[rows,cols] = Qmat[rows,cols] + chrcount_Qmat_beta
				} # END if (decploidy == ancploidy-1)

			# epsilon: +2x ploidy
			if (decploidy == ancploidy+2)
				{
				# Calculate the transitions under this ploidy
				chrcount_Qmat_ep = Qmat_chromo_ep(BioGeoBEARS_model_object, ploidy=ancploidy, maxnum_chromosomes, minnum_chromosomes, chrcount_Qmat_ep=NULL)
			
				# Add to Qmat
				Qmat[rows,cols] = Qmat[rows,cols] + chrcount_Qmat_ep
				} # END if (decploidy == ancploidy+2)

			# delta: diploidization (return to 2x)
			if (decploidy == 2)
				{
				# Add to Qmat
				Qmat[rows,cols] = Qmat[rows,cols] + chrcount_Qmat_delta
				} # END if (decploidy == 2)

			# rho - pure duplication, double ploidy
			if (decploidy == ancploidy*2)
				{
				# Add to Qmat
				Qmat[rows,cols] = Qmat[rows,cols] + chrcount_Qmat_rho
				} # END if (decploidy == ancploidy*2)
		
			} # END for (i in 1:numstates_ploidy)
		} # END for (j in 1:numstates_ploidy)

	# 33x33 matrix in example
	dim(Qmat)
	# Diagonal is already 0, but enforce that
	diag(Qmat) = 0
	# Set the diagonal to -rowSums to make official Qmat
	diag(Qmat) = -rowSums(Qmat)

	# Write to text, if desired
	# write.table(x=Qmat, file="Qmat_ploidy2_4_chrcount2_12.txt", sep="\t")
	
	return(Qmat)
	} # END Qmat_chromoploid <- function(BioGeoBEARS_model_object=define_ChromoBEARS_model_object(), maxnum_ploidy=4, minnum_ploidy=2, minnum_chromosomes=2, maxnum_chromosomes=12)











define_ChromoBEARS_model_object <- function(minval_anagenesis=1e-12, minval_cladogenesis=1e-5, maxval=5)
	{
	# Define the BioGeoBEARS_model class;

	#BioGeoBEARS_model_object = new("BioGeoBEARS_model", df=tmpdf)
	ChromoBEARS_model_object = new("BioGeoBEARS_model", params_table=ChromoBEARS_model_defaults(minval_anagenesis, minval_cladogenesis, maxval))
	ChromoBEARS_model_object
	
	# you can get the dataframe with
	# BioGeoBEARS_model_object@df

	return(ChromoBEARS_model_object)
	}




ChromoBEARS_model_defaults <- function(minval_anagenesis=1e-12, minval_cladogenesis=1e-5, maxval=5)
	{
	defaults='
	minval_anagenesis = 1e-12
	minval_cladogenesis = 1e-5
	maxval = 5
	ChromoBEARS_model_defaults()
	'
	
	param_data_starter = as.data.frame(matrix(data=0, nrow=1, ncol=7), stringsAsFactors=FALSE)
	names(param_data_starter) = c("type", "init", "min", "max", "est", "note", "desc")
	param_table = NULL
	param_names = NULL
	
	# For a parameter on a branch (anagenesis)
	
	# rate of adding 1 chromosome, no change in ploidy
	param_name = "lambda"
	param_data = param_data_starter
	param_data$type = "free"
	param_data$init = 0.11
	param_data$min = 0 + minval_anagenesis
	param_data$max = maxval - minval_anagenesis
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "rate of adding 1 chromosome"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table
	
	# rate of subtracting 1 chromosome, no change in ploidy
	param_name = "mu"
	param_data = param_data_starter
	param_data$type = "free"
	param_data$init = 0.10
	param_data$min = 0 + minval_anagenesis
	param_data$max = maxval - minval_anagenesis
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "rate of subtracting 1 chromosome"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table

	# odd bridge increase, +1x ploidy
	param_name = "alpha"
	param_data = param_data_starter
	param_data$type = "fixed"
	param_data$init = 0.011
	param_data$min = 0 + minval_anagenesis
	param_data$max = maxval - minval_anagenesis
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "odd bridge increase, +1x ploidy"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table

	# odd bridge decrease, -1x ploidy
	param_name = "beta"
	param_data = param_data_starter
	param_data$type = "fixed"
	param_data$init = 0.010
	param_data$min = 0 + minval_anagenesis
	param_data$max = maxval - minval_anagenesis
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "odd bridge decrease, -1x ploidy"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table

	# even bridge increase, +2x ploidy
	param_name = "epsilon"
	param_data = param_data_starter
	param_data$type = "fixed"
	param_data$init = 0.020
	param_data$min = 0 + minval_anagenesis
	param_data$max = maxval - minval_anagenesis
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "even bridge increase, +2x ploidy"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table

	# diploidization, return to 2x
	param_name = "delta"
	param_data = param_data_starter
	param_data$type = "fixed"
	param_data$init = 0.05
	param_data$min = 0 + minval_anagenesis
	param_data$max = maxval - minval_anagenesis
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "diploidization, return to 2x"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table

	# pure duplication
	param_name = "rho"
	param_data = param_data_starter
	param_data$type = "fixed"
	param_data$init = 0.04
	param_data$min = 0 + minval_anagenesis
	param_data$max = maxval - minval_anagenesis
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "pure duplication"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table
	
	return(param_table)
	}




#######################################################
# Crucial: in chromoploid_tiplikes:
# First row: num. species, then numploidies
# Later rows:
# 1. Species name
# 2. Ploidy likelihoods (1 if possible, 0 if impossible given data)
# 3. Chromosome count (e.g. "4") or range (e.g. "4-6") if there is a polymorphism
#######################################################

chromoploid_tiplikes <- function(fn, maxnum_ploidy=4, minnum_ploidy=2, minnum_chromosomes=2, maxnum_chromosomes=12, max_numstates=1500)
	{
	# Setup
	defaults='
	fn = "/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/custom_Qmat/ploidy_chrom_data.txt"
	minnum_ploidy=2
	maxnum_ploidy=4
	minnum_chromosomes=2
	maxnum_chromosomes=12
	max_numstates=1500
	'
	
	# Setup
	states_list_ploidy = seq(minnum_ploidy, maxnum_ploidy, 1)
	numstates_ploidy = length(states_list_ploidy)
	states_list_chrcount = seq(minnum_chromosomes, maxnum_chromosomes, 1)
	numstates_chrcount = length(states_list_chrcount)
	numstates = numstates_ploidy * numstates_chrcount
	
	txt = paste0("Your inputs specify ", numstates_ploidy, " ploidy states, and ", numstates_chrcount, " chromosome count states. This means ", numstates_ploidy, "x", numstates_chrcount, "=", numstates, " total states in your state space, and thus a ", numstates, "x", numstates, " Q rate matrix.")
	cat("\n\n")
	cat(txt)
	cat("\n\n")
	
	if (numstates > max_numstates)
		{
		txt = paste0("STOP ERROR in chromoploid_tiplikes(): Your numstates=", numstates, ", which is bigger than max_numstates=", max_numstates, ", the recommended maximum. Matrix exponentiation and the likelihood calculations get exceedingly slow above 1500 states. Fast computers and parallel processing will help a little but I've never run more than ~2200 states. Raise 'max_numstates' if you want to try it anyway despite the speed issue.")
		
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		
		stop(txt)
		} # END if (numstates > max_numstates)
	
	
	# Read the 1st line, split on whitespaces
	firstline = scan(fn, what="character", nlines=1)
	
	# Parse the firstline
	ntips = as.numeric(firstline[1])
	nploidies = as.numeric(firstline[2])
	
	if (nploidies != numstates_ploidy)
		{
		txt = paste0("STOP ERROR in chromoploid_tiplikes(). Your input file '", fn, "' says that the number of ploidy states is ", nploidies, ", but your input to chromoploid_tiplikes() says minnum_ploidy=", minnum_ploidy, ", maxnum_ploidy=", maxnum_ploidy, ", meaning numstates_ploidy=", numstates_ploidy, ". This could be a failure to change the default inputs or an error in '", fn, "'.")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		
		stop(txt)
		} # END if (nploidies != numstates_ploidy)
	
	
	# Make the data.frame
	tiplikes = matrix(data=NA, nrow=ntips, ncol=numstates)
	
	# Parse the remaining rows
	tmplines = scan(fn, what="character", sep="\t", skip=1)
	tmplines = matrix(data=tmplines, ncol=3, byrow=TRUE)
	tmplines

	# Go through each tip and calculate the starting likelihood
	for (i in 1:ntips)
		{
		# Ploidy likes
		ploidy_likes = as.numeric(strsplit(tmplines[i,2], split="")[[1]])
		tmpmat = matrix(data=NA, nrow=numstates_chrcount, ncol=nploidies)

		# For each ploidy
		for (j in 1:length(ploidy_likes))
			{
			tmpmat[,j] = ploidy_likes[j]
			}
		ploidy_likes_expanded = c(tmpmat)
		
		
		# Chrcount_likes
		tmpmat = matrix(0, nrow=numstates_chrcount, ncol=nploidies)
		
		chrcount_txt = tmplines[i,3]
		simpleTF = grepl(pattern="-", x=chrcount_txt) == FALSE
		if (simpleTF == TRUE)
			{
			chrcount = as.numeric(chrcount_txt)
			chrcount_stateindex_1based = match(x=chrcount, table=states_list_chrcount)
			tmpmat[chrcount_stateindex_1based,] = 1
			} else {
			# You have a range
			chrcounts_nums = strsplit(chrcount_txt, split="-")[[1]]
			chrcount_start = as.numeric(chrcounts_nums[1])
			chrcount_end = as.numeric(chrcounts_nums[2])
			index_start = match(x=chrcount_start, table=states_list_chrcount)
			index_end = match(x=chrcount_end, table=states_list_chrcount)
			tmpmat[index_start:index_end,] = 1
			} # END if (simpleTF == TRUE)
		chrcount_likes_expanded = c(tmpmat)
		
		# Now, multiply the two
		tiplikes_row = ploidy_likes_expanded * chrcount_likes_expanded
		
		# Store 'em
		tiplikes[i,] = tiplikes_row
		} # END for (i in 1:ntips)
	
	# In the example file, we have some ambiguity/polymorphism
	if ( any(rowSums(tiplikes) == 0) )
		{
		txt = paste0("STOP ERROR in chromoploid_tiplikes(). Some of your tips have likelihood 0 for ALL states! This is a fatal error, revise your input file.")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		
		cat("rowSums(tiplikes):\n\n")
		print(rowSums(tiplikes))
		
		stop(txt)
		}
	tail(tiplikes)
	
	return(tiplikes)
	} # END chromoploid_tiplikes <- function(fn, maxnum_ploidy=4, minnum_ploidy=2, minnum_chromosomes=2, maxnum_chromosomes=12)


