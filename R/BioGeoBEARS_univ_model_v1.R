require("ape")
require("rexpokit")
require("cladoRcpp")



# Negative version of loglike, for optimization
calc_loglike_for_optim_neg <- function(params, BioGeoBEARS_run_object, phy, tip_condlikes_of_data_on_each_state, print_optim=TRUE, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, return_what="loglike", calc_ancprobs=FALSE)
	{
	logLike = calc_loglike_for_optim(params, BioGeoBEARS_run_object, phy, tip_condlikes_of_data_on_each_state, print_optim=BioGeoBEARS_run_object$print_optim, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, return_what="loglike", calc_ancprobs=FALSE)
	negLogLike = -1 * logLike
	return(negLogLike)
	}


#######################################################
# Set up the function for optimization
#######################################################	
# params are a list of the values of the FREE parameters; but everything is contained in the 
# BioGeoBEARS_model object at all times


#######################################################
# calc_loglike_for_optim
#######################################################
#' Take model parameters and the data and calculate the log-likelihood
#' 
#' This function is an input to optim or optimx, the ML estimation routines.
#' 
#' @param params A vector of parameters for optimization.
#' @param BioGeoBEARS_run_object Object containing the run parameters and the model.
#' @param phy An ape tree object
#' @param tip_condlikes_of_data_on_each_state A numeric matrix with rows representing tips, and columns representing states/geographic ranges.  The cells
#' give the likelihood of the observation data under the assumption that the tip has that state; typically this means that the known geographic range gets a 
#' '1' and all other states get a 0.
#' @param force_sparse Should sparse matrix exponentiation be used?
#' @param print_optim If TRUE (default), print the optimization steps as ML estimation progresses.
#' @param areas_list A list of the desired area names/abbreviations/letters (?).
#' @param states_list A list of the possible states/geographic ranges, in 0-based index form.
#' @param cluster_already_open If the user wants to distribute the matrix exponentiation calculations from all the branches across a number of processors/nodes on 
#' a cluster, specify the cluster here.  E.g. \code{cluster_already_open = makeCluster(rep("localhost",num_cores_to_use), type = "SOCK")}.  Note: this will work on 
#' most platforms, including Macs running R from command line, but will NOT work on Macs running the R GUI \code{R.app}, because parallel processing functions like
#' \code{MakeCluster} from e.g. \code{library(parallel)} for some reason crash R.app.  The program runs a check for R.app and will just run on 1 node if found. 
#' @param return_what What should be returned to the user? Options are "loglike" (the log-likelihood of the data under the tree, model, and model parameters), 
#' "nodelikes" (the scaled conditional likelihoods at the nodes), "rootprobs" (the relative probability of the geographic ranges/states at the root), or "all"
#' (all of the above in a list).  Typically the user will only want to return "loglike" while doing ML optimization, but then return "all" once the ML parameter
#' values have been found.
#' @param calc_ancprobs Just use this function once, return the anc probs of states.
#' @return \code{ttl_loglike} The log-likelihood of the data under the input model and parameters.
#' @export
#' @seealso \code{\link{prune_states_list}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
calc_loglike_for_optim <- function(params, BioGeoBEARS_run_object, phy, tip_condlikes_of_data_on_each_state, print_optim=TRUE, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, return_what="loglike", calc_ancprobs=FALSE)
	{
	defaults='
	print_optim=TRUE; areas_list=areas_list; states_list=states_list; force_sparse=force_sparse; cluster_already_open=cluster_already_open; return_what="loglike"; calc_ancprobs=TRUE
	'
	
	if (is.null(BioGeoBEARS_run_object$printlevel))
		{
		BioGeoBEARS_run_object$printlevel = 0
		}
	printlevel = BioGeoBEARS_run_object$printlevel
	

	# Is this a traits-based analysis?
	traitTF = is.null(BioGeoBEARS_run_object$trait) == FALSE
	
	# Initialize m
	m = NULL

	# Initialize jts_matrix, matrix of t12, t23, etc., during a j event
	jts_matrix = NULL

	
	# Put the parameters into the BioGeoBEARS_model_object, so that they can be universally read out
	# into any function
	BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
	#print(params)
	#print(BioGeoBEARS_model_object)
	BioGeoBEARS_model_object = params_into_BioGeoBEARS_model_object(BioGeoBEARS_model_object=BioGeoBEARS_model_object, params=params)

	######################################################
	# 2016-03-23_NJM: adding rescaling
	# (unscale params, if they were used before)
	######################################################
	if (BioGeoBEARS_run_object$rescale_params == TRUE)
		{
		BioGeoBEARS_model_object@params_table = unscale_BGB_params(scaled_params_table=BioGeoBEARS_model_object@params_table)
		
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table = BioGeoBEARS_model_object@params_table
		}

	
	# Update linked parameters
	BioGeoBEARS_model_object = calc_linked_params_BioGeoBEARS_model_object(BioGeoBEARS_model_object)
	
	# Update to the run object, just to be SURE
	BioGeoBEARS_run_object$BioGeoBEARS_model_object = BioGeoBEARS_model_object
	
	# Set the dispersal and extinction rate
	d = BioGeoBEARS_model_object@params_table["d","est"]
	e = BioGeoBEARS_model_object@params_table["e","est"]
	a = BioGeoBEARS_model_object@params_table["a","est"]


	#######################################################
	#######################################################
	# Do branch-length exponentiation if desired
	#######################################################
	#######################################################
	b = BioGeoBEARS_model_object@params_table["b","est"]
	# Modify the edge.lengths
	phy$edge.length = phy$edge.length ^ b
	# Make sure this doesn't duplicate a previous "^b", e.g.
	# the summarization step in bears_optim_run

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
	dispersal_multipliers_matrix = dispersal_multipliers_matrix * (manual_dispersal_multipliers_matrix ^ w)

	#######################################################
	# ALTERNATIVE DISTANCE FUNCTIONS
	#
	# The exponential function of distance is heavy-tailed.
	# Other functions might be empirically better
	# 
	# NOTE: in order to keep the base dispersal rate
	# (base dispersal rate = rate at distance=0) comparable
	# across models, the probability density function should 
	# be normalized (/divided by) the pdf evaluated at zero.
	#
	# In other words, the dispersal multiplier should always
	# come out as 1, when distance=0. I.e. pdf(x=0) = 1.
	# 
	# These will each be identified by optional parameters
	# in the model:
	# WALD
	# HNORM
	#
	#######################################################
	
	#######################################################
	#	WALD distribution
	#######################################################
	#
	# If one has a model where dispersal rate is modified as 
	# a function of the WALD distribution (inverse gaussian),
	# this will require parameters:
	# WALD_mu - mu, the mean in the WALD probability density function
	# WALD_lambda - lambda, the shape parameter
	# See: https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution
	#######################################################
	
	# Alternative distance model abbreviations
	alt_distance_models_abbr = c("WALD", "HNORM")
	tmp_param_names = row.names(BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table)

	# Calculate additional dispersal modifiers on distance
	# WALD distribution
	TF = grepl(pattern="WALD", x=tmp_param_names)
	if (sum(TF) > 0)
		{
		require(statmod)	# for dinvgauss
		tmpname = "WALD_mu"
		TF = grepl(pattern=tmpname, x=tmp_param_names)
		WALD_mu = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table[tmpname, "est"]

		tmpname = "WALD_lambda"
		TF = grepl(pattern=tmpname, x=tmp_param_names)
		WALD_lambda = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table[tmpname, "est"]
		
		# Calculate the multipliers under this model
		tmp_multipliers = dinvgauss(x=distances_mat, mean=WALD_mu, shape=WALD_lambda)
		normalizer = dinvgauss(x=0, mean=WALD_mu, shape=WALD_lambda)
		tmp_multipliers = tmp_multipliers / normalizer
		
		# Multiply element-wise
		dispersal_multipliers_matrix = dispersal_multipliers_matrix * tmp_multipliers
		} # END if (sum(TF) > 0)


	# Half-normal distribution
	# https://en.wikipedia.org/wiki/Half-normal_distribution
	# Theta is the scaled precision (inverse of the variance)
	# Theta = sqrt(pi) / (sigma * sqrt(2))
	TF = grepl(pattern="HNORM", x=tmp_param_names)
	if (sum(TF) > 0)
		{
		require(fdrtool)	# for dhalfnorm
		tmpname = "HNORM_theta"
		TF = grepl(pattern=tmpname, x=tmp_param_names)
		HNORM_theta = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table[tmpname, "est"]


		# Calculate the multipliers under this model
		tmp_multipliers = dhalfnorm(x=distances_mat, theta=HNORM_theta)
		normalizer = halfnorm(x=0, theta=HNORM_theta)
		tmp_multipliers = tmp_multipliers / normalizer
		
		# Multiply element-wise
		dispersal_multipliers_matrix = dispersal_multipliers_matrix * tmp_multipliers
		} # END if (sum(TF) > 0)



	#######################################################
	# multiply parameter d by dispersal_multipliers_matrix
	#######################################################
	dmat_times_d = dispersal_multipliers_matrix * matrix(d, nrow=length(areas), ncol=length(areas))
	amat = dispersal_multipliers_matrix * matrix(a, nrow=length(areas), ncol=length(areas))




	#######################################################
	#######################################################
	# Do area-dependence and extinction multipliers list
	#######################################################
	#######################################################
	if ( (is.null(BioGeoBEARS_run_object$list_of_area_of_areas) == FALSE))
		{
		area_of_areas = BioGeoBEARS_run_object$list_of_area_of_areas[[1]]
		} else {
		# Default is all areas effectively equidistant
		area_of_areas = rep(1, length(areas))
		}
		
	# Get the exponent on extinction, apply to extinction modifiers	
	u = BioGeoBEARS_model_object@params_table["u","est"]
	extinction_modifier_list = area_of_areas ^ (1 * u)
	
	# Apply to extinction rate
	elist = extinction_modifier_list * rep(e, length(areas))
	
	# Set up the instantaneous rate matrix (Q matrix)
	# someday we'll have to put "a" (anagenic range-switching) in here
	# (this was been done in 2014 - NJM)
	
	# Substitute here (for starters) if you want to have a custom Qmat
	if (is.null(BioGeoBEARS_run_object$custom_Qmat_fn_text) == TRUE)
		{
		# Standard analysis, no traits
		if (traitTF == FALSE)
			{
			Qmat = rcpp_states_list_to_DEmat(areas_list=areas_list, states_list=states_list, dmat=dmat_times_d, elist=elist, amat=amat, include_null_range=BioGeoBEARS_run_object$include_null_range, normalize_TF=TRUE, makeCOO_TF=force_sparse)
			}
		
		# Analysis with a trait modifying dispersal rate
		if (traitTF == TRUE)
			{
			numstates_geogtrait = ncol(tip_condlikes_of_data_on_each_state)
			
			res = modify_Qmat_with_trait(Qmat=NULL, BioGeoBEARS_run_object, numstates_geogtrait=numstates_geogtrait, areas_list=areas_list, states_list=states_list, dispersal_multipliers_matrix=dispersal_multipliers_matrix, elist=elist, force_sparse=force_sparse)
			Qmat = res$Qmat
			m = res$m

			# If the trait can change during jump events
			if (is.null(BioGeoBEARS_run_object$jts_txt_matrix) == FALSE)
				{
				jts_txt_matrix = BioGeoBEARS_run_object$jts_txt_matrix
				jts_matrix = matrix(data=0, nrow=nrow(jts_txt_matrix), ncol=ncol(jts_txt_matrix))
				TF_matrix = matrix(data=TRUE, nrow=nrow(jts_txt_matrix), ncol=ncol(jts_txt_matrix))
				diag(TF_matrix) = FALSE
				jts_txt_params = c(jts_txt_matrix[TF_matrix])
				jts_txt_params
			
				# Populate the numeric jts_matrix
				for (jts_i in 1:nrow(jts_txt_matrix))
					{
					diag_val = 1
					for (jts_j in 1:ncol(jts_txt_matrix))
						{
						if (jts_i == jts_j)
							{
							next()
							}
						jts_txt = jts_txt_matrix[jts_i,jts_j]
						newval = as.numeric(BioGeoBEARS_model_object@params_table[jts_txt, "est"])
						jts_matrix[jts_i,jts_j] = newval
						diag_val = 1-newval
						}
					# Populate the diagonal
					jts_matrix[jts_i,jts_i] = diag_val
					} # END for (jts_i in 1:nrow(jts_txt_matrix))
				} # END if (is.null(BioGeoBEARS_run_object$jts_txt_matrix) == FALSE)
			} # END if (if (traitTF == TRUE))
		
		} else {
		cat("\n\nNOTE: BioGeoBEARS is using a custom Qmat-generating function.\n\n")
		# Evaluates to "Qmat"
		eval(parse(text=BioGeoBEARS_run_object$custom_Qmat_fn_text))
		} # END if (is.null(BioGeoBEARS_run_object$custom_Qmat_fn_text) == TRUE)
	# Print the Qmat, if desired
#	print(round(Qmat, 4))
	
	
	
	#######################################################
	# Cladogenic model
	#######################################################
	spPmat_inputs = get_spPmat_inputs_from_BGB(BioGeoBEARS_run_object=BioGeoBEARS_run_object, states_list=states_list, dispersal_multipliers_matrix=dispersal_multipliers_matrix)

	#######################################################
	# Detection model
	#######################################################
	# Use detection model to generate tip likelihoods if desired; or take
	# pre-specified tip likelihoods. Otherwise, use those input from bears_optim_run
	if (is.null(BioGeoBEARS_run_object$tip_condlikes_of_data_on_each_state) == TRUE)
		{
		# Usual strategy for calculating tip-likelihoods
		# Get the detection model
		if (BioGeoBEARS_run_object$use_detection_model == TRUE)
			{
			# Calculate the initial tip likelihoods, using the detection model
			# Assumes correct order, double-check this
			numareas = length(areas)
			detects_df = BioGeoBEARS_run_object$detects_df
			controls_df = BioGeoBEARS_run_object$controls_df
			mean_frequency = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mf", "init"]
			dp = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["dp", "init"]
			fdp = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["fdp", "init"]
		
			# return_LnLs=TRUE ensures no under-flow
			tip_condlikes_of_data_on_each_state = tiplikes_wDetectionModel(states_list_0based_index=states_list, phy=phy, numareas=numareas, detects_df=detects_df, controls_df=controls_df, mean_frequency=mean_frequency, dp=dp, fdp=fdp, null_range_gets_0_like=TRUE, return_LnLs=TRUE, relative_LnLs=TRUE, exp_LnLs=TRUE, error_check=TRUE)
			}
		#print(tip_condlikes_of_data_on_each_state)

		} else {
		# Or, use pre-specified tip conditional likelihoods
		# Pre-specified (custom) tip-likelihoods
		tip_condlikes_of_data_on_each_state = BioGeoBEARS_run_object$tip_condlikes_of_data_on_each_state
		} # END if (is.null(BioGeoBEARS_run_object$tip_condlikes_of_data_on_each_state) == FALSE)



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


	
	if (calc_ancprobs == FALSE)
		{
		# E.g., during optimx(), you don't need the ancestral
		# states, nor the uppass/downpass stuff
		# NOTE: We should, though, include
		# fixlikes when calc_ancprobs = TRUE
		fixnode = BioGeoBEARS_run_object$fixnode
		fixlikes = BioGeoBEARS_run_object$fixlikes
		
		# Calculate the log-likelihood of the data, given the model parameters during this iteration
		#print(jts_matrix)
		#print("BioGeoBEARS_run_object$printlevel #1")
		#print(BioGeoBEARS_run_object$printlevel)
		ttl_loglike = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="loglike", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, cppSpMethod=3, printlevel=BioGeoBEARS_run_object$printlevel, cluster_already_open=cluster_already_open, calc_ancprobs=FALSE, fixnode=fixnode, fixlikes=fixlikes, include_null_range=BioGeoBEARS_run_object$include_null_range, m=m, jts_matrix=jts_matrix, BioGeoBEARS_model_object=BioGeoBEARS_run_object$BioGeoBEARS_model_object, on_NaN_error=BioGeoBEARS_run_object$on_NaN_error)
		ttl_loglike

		if (print_optim == TRUE)
			{
			LnL = ttl_loglike
			# If the log likelihood is successful, print it
			outvars = adf(t(c(BioGeoBEARS_model_object@params_table$est, LnL)))
			#outvars = cbind(outvars, LnL)
			
			#print("HERE #1!!!")
			names(outvars) = c(rownames(BioGeoBEARS_model_object@params_table), "LnL")
			print(round(outvars,3))
	
			#cat(ttl_loglike, "\n", sep="")
			}
	
		return(ttl_loglike)
		} else {
		# E.g., after optimx(), you *DO* usually want the ancestral
		# states, *AND* the uppass/downpass stuff
		# NOTE: We should, though, include
		# fixlikes when calc_ancprobs = TRUE
		
		# Fixing ancestral nodes
		fixnode = BioGeoBEARS_run_object$fixnode
		fixlikes = BioGeoBEARS_run_object$fixlikes
		
		# Print m
		# (m1, m2, etc.)
#		print("m:")
#		print(m)
		
		# Calculate EVERYTHING!
		#print(jts_matrix)
		model_results = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="all", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=BioGeoBEARS_run_object$printlevel, cluster_already_open=cluster_already_open, calc_ancprobs=TRUE, fixnode=fixnode, fixlikes=fixlikes, include_null_range=BioGeoBEARS_run_object$include_null_range, m=m, jts_matrix=jts_matrix, BioGeoBEARS_model_object=BioGeoBEARS_run_object$BioGeoBEARS_model_object, on_NaN_error=BioGeoBEARS_run_object$on_NaN_error)
		return(model_results)
		}
	} # END calc_loglike_for_optim









# Any parameter model, adding j, v (vicariance proportion), maxent_constraint_01 (for non-vicariant subsets), maxent_constraint_01v (weighting for size of smaller offspring)
#######################################################
# bears_optim_run
#######################################################
#' Run ML search from \code{BioGeoBEARS_run} object
#' 
#' Uses a BioGeoBEARS_run_object to simplify input.
#'
#' @param BioGeoBEARS_run_object Contains all inputs
#' @param skip_optim If TRUE, just calculate the starting 
#' likelihood, and skip the optimization (mostly for timing). 
#' @param skip_optim_option Default "return_loglike" returns the log-likelihood (default 
#' skip_optim=TRUE behavior). The other option, "return_all", will return everything
#' for the starting parameters (ancestral state probabilities, etc.)
#' Default FALSE.
#' @return \code{bears_output} A list of outputs.  bears_output$optim_result
#' @export
#' @seealso \code{\link{readfiles_BioGeoBEARS_run}}, \code{\link{bears_2param_standard_fast}}, \code{\link[cladoRcpp]{numstates_from_numareas}}, \code{\link{getranges_from_LagrangePHYLIP}}, \code{\link[ape]{read.tree}}, \code{\link{calc_loglike_sp}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' Felsenstein, Joe.  The Newick tree format.  \url{http://evolution.genetics.washington.edu/phylip/newicktree.html}
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#'	 @cite Ree2009configurator
#'	 @cite SmithRee2010_CPPversion
#'	 @cite Landis_Matzke_etal_2013_BayArea
#' @examples
#' test=1
#' 
#' # Get the example files directory
#' extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
#' # tmp hard code: 
#' # extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"
#'
#' # Set the filenames (Hawaiian Psychotria from Ree & Smith 2008)
#' trfn = np(paste(extdata_dir, "/Psychotria_5.2.newick", sep=""))
#' tr = read.tree(file=trfn)
#' 
#' geogfn = np(paste(extdata_dir, "/Psychotria_geog.data", sep=""))
#' 
#' # Look at the tree and ranges, for kicks
#' getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
#' tr
#' 
#' \dontrun{
#' # Run the ML search
#' bears_output = bears_optim_run(trfn=trfn, geogfn=geogfn)
#' bears_output
#' }
#'
bears_optim_run <- function(BioGeoBEARS_run_object = define_BioGeoBEARS_run(), skip_optim=FALSE, skip_optim_option="return_loglike")
	{
	defaults='	
	BioGeoBEARS_run_object = define_BioGeoBEARS_run()
	BioGeoBEARS_run_object
	skip_optim=FALSE
	
	skip_optim=TRUE
	skip_optim_option="return_all"
	'
	
	require(cladoRcpp)
	require(rexpokit)
	
	
	# Wipe out any old/previous warnings()
	assign("last.warning", NULL, envir = baseenv())
	
	
	if (is.null(BioGeoBEARS_run_object$include_null_range))
		{
		BioGeoBEARS_run_object$include_null_range = TRUE
		}	
	include_null_range = BioGeoBEARS_run_object$include_null_range
	
	if (is.null(BioGeoBEARS_run_object$allow_null_tips))
		{
		BioGeoBEARS_run_object$allow_null_tips = FALSE
		}
	
	# 2017-11-29
	# Error check for tr if not loaded elsewhere
	if (exists("tr") == FALSE)
		{
		tr = read.tree(BioGeoBEARS_run_object$trfn)
		}

	
	#######################################################
	# Check for traits and trait model
	#   - Need BOTH, or NEITHER
	#######################################################
	traitTF = is.null(BioGeoBEARS_run_object$trait) == FALSE

	# Initialize m, if needed
	m = NULL
	
	if (is.null(BioGeoBEARS_run_object$min_branchlength) == FALSE)
		{
		min_branchlength = BioGeoBEARS_run_object$min_branchlength
		} else {
		min_branchlength = 0.000001
		}
	
	# ERROR CHECKS FOR TRAIT MODEL
	if (traitTF)
		{
		if (is.na(BioGeoBEARS_run_object$timesfn) == FALSE)
			{
			txt = "WARNING: you have loaded a BioGeoBEARS_run_object$trait, but you have a timesfn, indicate a time-stratified analysis. Traits-based dispersal has now been implemented for time-stratified analyses, but is still experimental."
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			warning(txt)			
			}
		
		# Check for trait_Pmat_txt
		if (is.null(BioGeoBEARS_run_object$trait_Pmat_txt) == TRUE)
			{
			txt = "STOP ERROR: you have loaded a BioGeoBEARS_run_object$trait, but you are missing a BioGeoBEARS_run_object$trait_Pmat_txt"
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			stop(txt)
			}
		
		# Check for trait transition rates
		# Check for t12, t21, etc.
		trait_transition_rates_TF = grepl(pattern="trait transition rate", x=BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$desc)
		if (sum(trait_transition_rates_TF) < 1)
			{
			txt = "STOP ERROR: you have loaded a BioGeoBEARS_run_object$trait, but you need one or more 'trait transition rates' in  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table"
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			stop(txt)
			}
		
		# Check for trait-based dispersal rate multipliers		
		# Check for m1, m2, etc.
		numtraitstates = ncol(BioGeoBEARS_run_object$trait@df)
		traitbased_dispersal_Ms_TF = grepl(pattern="trait-based dispersal rate multiplier", x=BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$desc)

		if (sum(traitbased_dispersal_Ms_TF) != numtraitstates)
			{
			txt = paste0("STOP ERROR: you have loaded a BioGeoBEARS_run_object$trait, and it has ", numtraitstates, " states, so you need to have ", numtraitstates, " multipliers ('m1', 'm2', etc.) with 'desc' field 'trait-based dispersal rate multipliers...' in  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table. Instead, you have only this many: ", sum(traitbased_dispersal_Ms_TF))
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			stop(txt)
			}
		} # END if (traitTF) # ERROR CHECK

	
	# Load the trait as a (another) tipranges-class object
	if (traitTF == TRUE)
		{
		trait = BioGeoBEARS_run_object$trait
		trait_as_tip_condlikes = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges=trait, phy=tr, states_list=NULL, maxareas=1, include_null_range=FALSE, useAmbiguities=BioGeoBEARS_run_object$useAmbiguities, trait_as_tip_condlikes=NULL)
		
		# Number of traits
		ntrait_states = ncol(trait_as_tip_condlikes)
		
		# Extract these submatrices, just for dimensions, names etc.
		# ALWAYS extract parameter values from the main model_object
		# Trait modeling effect on dispersal
		# (m1, m2, etc.)
		BGB_trait_model_params_table = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table[traitbased_dispersal_Ms_TF,]
		
		# Parameters of transition matrix for the trait
		# (t12, t21, etc.)
		BGB_trait_Pmat_params_table = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table[trait_transition_rates_TF,]

		# Text transition matrix for the trait
		trait_Pmat_txt = BioGeoBEARS_run_object$trait_Pmat_txt
		} else {
		# No trait; set to NULL
		trait_as_tip_condlikes = NULL
		ntrait_states = NULL
		BGB_trait_model_params_table = NULL
		BGB_trait_Pmat_params_table = NULL
		} # END if (traitTF == TRUE)
	
	

	#######################################################
	# Load the model object
	#######################################################
	#inputs = BioGeoBEARS_run_object
	BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
	
	# Should the optim run be printed?
	print_optim = BioGeoBEARS_run_object$print_optim


	# Get geographic ranges at tips
	if (BioGeoBEARS_run_object$use_detection_model == FALSE)
		{
		tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=np(BioGeoBEARS_run_object$geogfn))
		}
	if (BioGeoBEARS_run_object$use_detection_model == TRUE)
		{
		if (BioGeoBEARS_run_object$use_detection_model == TRUE)
			{
			tipranges = tipranges_from_detects_fn(detects_fn=BioGeoBEARS_run_object$detects_fn)
			} # END if (inputs$use_detection_model == TRUE)
		} # END if (BioGeoBEARS_run_object$use_detection_model == TRUE)
	
	
	# Should we do optimx speedup?
	speedup = BioGeoBEARS_run_object$speedup
	

	# Get the list of geographic areas
#	print("print(tipranges):")
#	print(tipranges)
	areas = getareas_from_tipranges_object(tipranges)
	areas_list = seq(0, length(areas)-1, 1)		# 0-base indexes

	# Change the names to tipranges@df:
	# This converts the tipranges names to 0-based index
	# this doesn't make sense if areas_list is 0-based indexes
	# XXX - check at some point
	# REMOVED: 2017-03-14
	#names(tipranges@df) = areas_list

	#######################################################
	# Set the maximum range size (can be thought of as
	# a free parameter
	#######################################################
	if (is.na(BioGeoBEARS_run_object$max_range_size))
		{
		if (is.null(BioGeoBEARS_run_object$states_list))
			{
			# Maximum range size is all areas
			max_range_size = length(areas)
			} else {
			# If not NA
			# Get max rangesize from states list
			max_range_size = max(sapply(X=BioGeoBEARS_run_object$states_list, FUN=length), na.rm=TRUE)
			}
		} else {
		# Maximum range size hard-coded
		max_range_size = BioGeoBEARS_run_object$max_range_size
		}
	max_numareas = max_range_size
	
	#######################################################
	# Check that no tips have larger ranges than you allowed
	#######################################################
	#print("Here")
	#print(tipranges@df)
	# The dfnums_to_numeric fails, if you re-labeled the area names to 0, 1, 2, etc...
	tipranges_df_tmp = tipranges@df
	names(tipranges_df_tmp) = paste0("col", names(tipranges_df_tmp))
	tipranges_df_tmp[tipranges_df_tmp=="?"] = 0
	TF = (rowSums(dfnums_to_numeric(tipranges_df_tmp))) > max_range_size
	if (sum(TF, na.rm=TRUE) > 0)
		{
		cat("\n\nERROR: Tips with ranges too big:\n", sep="")
		print(dfnums_to_numeric(tipranges_df_tmp)[TF, ])
		cat("\n\nCheck your input geography file!\n", sep="")
		txt = paste("ERROR: Some tips (listed above) have range sizes larger than ", max_range_size, sep="")
		stop(txt)
		}

	
	
	
	
	
	# Old/slow way of getting the list of states and speciation matrix (slow)
	# old_states_list = areas_list_to_states_list(areas, include_null_range=BioGeoBEARS_run_object$include_null_range)
	# old_states_list
	# spmat = make_relprob_matrix_bi(old_states_list)
	# spmat
	
	# max_numareas = max(sapply(X=old_states_list, FUN=length), na.rm=TRUE)
	# max_numareas


	# Take the list of areas, and get list of possible states
	# (the user can manually input states if they like)
	if (is.null(BioGeoBEARS_run_object$states_list))
		{
		states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=BioGeoBEARS_run_object$include_null_range)
		states_list
		#BioGeoBEARS_run_object$states_list = states_list
		#inputs$states_list = states_list
		} else {
		states_list = BioGeoBEARS_run_object$states_list
		#BioGeoBEARS_run_object$states_list = states_list
		#inputs$states_list = states_list
		}
	
	# Used if time-changing stratified states list
	all_states_list = states_list
	
	#######################################################
	# NON-STRATIFIED: Modify the states_list if needed
	#######################################################
	if ( is.numeric(BioGeoBEARS_run_object$timeperiods) == FALSE )
		{
		#######################################################
		# If needed, modify the states_list by areas_allowed_mat
		#######################################################
		if ( (is.null(BioGeoBEARS_run_object$list_of_areas_allowed_mats) == FALSE))
			{
			# Take the first areas_allowed matrix (non-stratified)
			areas_allowed_mat = BioGeoBEARS_run_object$list_of_areas_allowed_mats[[1]]
		
			# Cut down the states accordingly (hopefully not slow!)
			original_states_list = states_list
			states_list = prune_states_list(states_list_0based_index=states_list, areas_allowed_mat=areas_allowed_mat)
			BioGeoBEARS_run_object$states_list = states_list

			print("Limiting original_states_list using an areas_allowed matrix")
			print("original_states_list")
			print(original_states_list)
			cat("\nlength(original_states_list) = ", length(original_states_list), " states/ranges.\n")
			cat("\n")

			print("states_list")
			print(states_list)
			cat("\nlength(original_states_list) = ", length(original_states_list), " states/ranges.")
			cat("\nlength(states_list) = ", length(states_list), " states/ranges.\n")


			} else {
			# Make no change
			pass = 1
			# states_list = states_list
			}

		#######################################################
		# If needed, modify the states_list by areas_adjacency_mat
		#######################################################
		if ( (is.null(BioGeoBEARS_run_object$list_of_areas_adjacency_mats) == FALSE))
			{
			# Take the first areas_adjacency matrix (non-stratified)
			areas_adjacency_mat = BioGeoBEARS_run_object$list_of_areas_adjacency_mats[[1]]
		
			# Cut down the states accordingly (hopefully not slow!)
			original_states_list = states_list
			states_list = prune_states_list_by_adjacency(states_list_0based_index=states_list, areas_adjacency_mat=areas_adjacency_mat)
			BioGeoBEARS_run_object$states_list = states_list
			
			print("Limiting original_states_list using an area adjacency matrix")
			print("original_states_list")
			print(original_states_list)
			print(length(original_states_list))
			cat("\n")

			print("states_list")
			print(states_list)
			print("length(states_list)")
			print(length(states_list))
			
			} else {
			# Make no change
			pass = 1
			# states_list = states_list
			}
		} # END if ( is.numeric(BioGeoBEARS_run_object$timeperiods) == FALSE )
	
	# Change the states_list by traits, if needed
	# (non-stratified)
	if (traitTF == TRUE)
		{
		states_list_ORIG = states_list
		#states_list_wTrait = 
		
		#trait_as_tip_condlikes
		}
	


	#######################################################
	# STRATIFIED: Modify the states_list if needed
	# (this is ONLY if the state-space is changing in the
	#  different time-slices)
	#######################################################
	# Will the state space be changing?
	TF1 = (is.null(BioGeoBEARS_run_object$list_of_areas_allowed_mats) == FALSE)
	TF2 = (is.null(BioGeoBEARS_run_object$list_of_areas_adjacency_mats) == FALSE)
	state_space_changing_TF = (TF1 + TF2) > 0
	need_to_print_list_of_states_list = TRUE
	master_states_list = states_list	# store the master list of all states;
										# check that this includes all, at some point
										# if not, warn user to change it manually
	
	if ( (is.numeric(BioGeoBEARS_run_object$timeperiods) == TRUE) && (state_space_changing_TF == TRUE) && (is.null(BioGeoBEARS_run_object$lists_of_states_lists_0based) == TRUE) )
		{
		need_to_print_list_of_states_list = FALSE
		ntimes = length(BioGeoBEARS_run_object$timeperiods)
		lists_of_states_lists_0based = list()
		
		# Go through each time bin, and make the state space different in each time bin
		for (ti in 1:ntimes)
			{
			# Initialize
			states_list_for_this_stratum = states_list
			
			
			#######################################################
			# If needed, modify the states_list by areas_allowed_mat
			#######################################################
			# Areas allowed matrix
			if (TF1 == TRUE)
				{
				# Take the first areas_allowed matrix (non-stratified)
				areas_allowed_mat = BioGeoBEARS_run_object$list_of_areas_allowed_mats[[ti]]
		
				# Cut down the states accordingly (hopefully not slow!)
				states_list_for_this_stratum = prune_states_list(states_list_0based_index=states_list_for_this_stratum, areas_allowed_mat=areas_allowed_mat)
				} else {
				# Make no change
				pass = 1
				# states_list = states_list
				}

			# Message to user
			timeslice_num = ti
			if (timeslice_num == 1)
				{
				toptime = 0
				} else {
				toptime = BioGeoBEARS_run_object$timeperiods[ti-1]
				}
			if (timeslice_num == ntimes)
				{
				bottime = BioGeoBEARS_run_object$timeperiods[ti]
				catend = "\n\n"
				} else {
				bottime = BioGeoBEARS_run_object$timeperiods[ti]
				catend = ""
				}
			txt = paste0("bears_optim_run() note: overall states_list has ", length(master_states_list), " states/ranges. In stratum #", ti, " (", toptime, "-", bottime, " mya), states_list_for_this_stratum has ", length(states_list_for_this_stratum), " states/ranges, due to areas_allowed and/or areas_adjacency matrices. See BioGeoBEARS_run_object$lists_of_states_lists_0based.")
			cat("\n")
			cat(txt)
			cat(catend)


			#######################################################
			# If needed, modify the states_list by areas_adjacency_mat
			#######################################################
			# Areas adjacency matrix
			if (TF2 == TRUE)
				{
				# Take the first areas_adjacency matrix (non-stratified)
				areas_adjacency_mat = BioGeoBEARS_run_object$list_of_areas_adjacency_mats[[ti]]
		
				# Cut down the states accordingly (hopefully not slow!)
				states_list_for_this_stratum = prune_states_list_by_adjacency(states_list_0based_index=states_list_for_this_stratum, areas_adjacency_mat=areas_adjacency_mat)
				} else {
				# Make no change
				pass = 1
				# states_list = states_list
				}

			# Store in the list of states_lists
			lists_of_states_lists_0based[[ti]] = states_list_for_this_stratum
			} # END for (ti in 1:ntimes)
		
		# Store the time-stratified list of states_lists in the BioGeoBEARS_run_object
		BioGeoBEARS_run_object$lists_of_states_lists_0based = lists_of_states_lists_0based
		
		} # END if ( (is.numeric(BioGeoBEARS_run_object$timeperiods) == TRUE) && (state_space_changing_TF == TRUE) )


	# Or, if the time-stratified stats list is pre-specified
	if (is.null(BioGeoBEARS_run_object$lists_of_states_lists_0based) == FALSE)
		{
		ntimes = length(BioGeoBEARS_run_object$timeperiods)
		
		states_allowed_TF1 = rep(TRUE, times=length(all_states_list))
		states_allowed_TF2 = rep(TRUE, times=length(all_states_list))
		states_allowed_TF3 = rep(TRUE, times=length(all_states_list))
		
		for (ntimes_i in 1:ntimes)
			{
			# Combine the 3 ways of changing states lists		
			# Areas allowed in this time bin
			if ( (is.null(BioGeoBEARS_run_object$list_of_areas_allowed_mats) == FALSE))
				{
				areas_allowed_mat = BioGeoBEARS_run_object$list_of_areas_allowed_mats[[ntimes_i]]

				states_allowed_TF1 = sapply(X=all_states_list, FUN=check_if_state_is_allowed, areas_allowed_mat)
				#states_to_use_TF = all_states_list %in% tmp_states_list
		
				if (include_null_range == TRUE)
					{
					states_allowed_TF1[1] = TRUE
					}
				# NO; use all areas for this
				# states_to_use_TF = states_allowed_TF
				} # END if ( (is.null(BioGeoBEARS_run_object$list_of_areas_allowed_mats) == FALSE))
		
			# Areas adjacency
			if ( (is.null(BioGeoBEARS_run_object$list_of_areas_adjacency_mats) == FALSE))
				{
				areas_adjacency_mat = BioGeoBEARS_run_object$list_of_areas_adjacency_mats[[ntimes_i]]

				states_allowed_TF2 = sapply(X=all_states_list, FUN=check_if_state_is_allowed_by_adjacency, areas_adjacency_mat)
				#states_to_use_TF = all_states_list %in% tmp_states_list
		
				if (include_null_range == TRUE)
					{
					states_allowed_TF2[1] = TRUE
					}
				# NO; use all areas for this
				# states_to_use_TF = states_allowed_TF
				} # END if ( (is.null(BioGeoBEARS_run_object$list_of_areas_adjacency_mats) == FALSE))

			# Manual list of allowed states
			if ( (is.null(BioGeoBEARS_run_object$lists_of_states_lists_0based) == FALSE))
				{
				states_allowed_TF3 = all_states_list %in% BioGeoBEARS_run_object$lists_of_states_lists_0based[[ntimes_i]]
			
				if (include_null_range == TRUE)
					{
					states_allowed_TF3[1] = TRUE
					}
				} # END if ( (is.null(BioGeoBEARS_run_object$lists_of_states_lists_0based) == FALSE))

			# Combine the 3 (areas_allowed, areas_adjacency, lists_of_states_lists_0based)
			states_allowed_TF = ((states_allowed_TF1 + states_allowed_TF2 + states_allowed_TF3) == 3)

			# CHANGE the inputs here, so that it can be used easily in BSM
			BioGeoBEARS_run_object$lists_of_states_lists_0based[[ntimes_i]] = all_states_list[states_allowed_TF]
			} # END for (ntimes_i in 1:ntimes)
		
		txt = paste0("bears_optim_run() note: BioGeoBEARS_run_object$lists_of_states_lists_0based has been specified. This means there is a different state space in each timebin / stratum / epoch.")
		cat("\n")
		cat(txt)
		cat("\n")
		
		# Check that number of lists of states matches the number of timebins
		number_of_lists_of_states = length(BioGeoBEARS_run_object$lists_of_states_lists_0based)
		if (ntimes == number_of_lists_of_states)
			{
			txt = paste0("bears_optim_run() note: BioGeoBEARS_run_object has ", ntimes, " timebins and ", number_of_lists_of_states, " lists of states ranges. Check passed.")
			cat("\n")
			cat(txt)
			cat("\n")
			} else {
			txt = paste0("bears_optim_run() STOP ERROR: BioGeoBEARS_run_object has ", ntimes, " timebins and ", number_of_lists_of_states, " lists of states ranges. Check FAILED.")
			cat("\n")
			cat(txt)
			cat("\n")
			stop(txt)
			} # END if (ntimes = number_of_lists_of_states)

		
		
		# Go through each time bin, and make the state 
		# space different in each time bin
		if (need_to_print_list_of_states_list == TRUE)
			{
			for (ti in 1:ntimes)
				{
				# Extract the states list in this time-stratum
				states_list_for_this_stratum = BioGeoBEARS_run_object$lists_of_states_lists_0based[[ti]]
				
				# Message to user
				timeslice_num = ti
				if (timeslice_num == 1)
					{
					toptime = 0
					} else {
					toptime = BioGeoBEARS_run_object$timeperiods[ti-1]
					}
				if (timeslice_num == ntimes)
					{
					bottime = BioGeoBEARS_run_object$timeperiods[ti]
					catend = "\n\n"
					} else {
					bottime = BioGeoBEARS_run_object$timeperiods[ti]
					catend = ""
					} # END if (timeslice_num == ntimes)				
					
					
				txt = paste0("bears_optim_run() note: overall states_list has ", length(master_states_list), " states/ranges. In stratum #", ti, " (", toptime, "-", bottime, " mya), states_list_for_this_stratum has ", length(states_list_for_this_stratum), " states/ranges, due to user-specified states_lists. See BioGeoBEARS_run_object$lists_of_states_lists_0based.")
				cat("\n")
				cat(txt)
				cat(catend)
				} # END for (ti in 1:ntimes)
			} # END if (need_to_print_list_of_states_list == TRUE)
		} # END if (is.null(BioGeoBEARS_run_object$lists_of_states_lists_0based) == TRUE)
		# END printing user-specified list of states_lists

	
	
	#######################################################
	# Sparse matrix exponentiation, if desired (dubious)
	#######################################################
	if (is.na(BioGeoBEARS_run_object$force_sparse))
		{
		if (length(states_list) > 128)
			{
			force_sparse = TRUE
			cat("\nNote: force_sparse being set to TRUE, as length(states_list) > 128\n", sep="")
			} else {
			force_sparse = FALSE
			}
		} else {
		force_sparse = BioGeoBEARS_run_object$force_sparse
		}


	
	#######################################################
	# Load the phylogenetic tree
	#######################################################
	trfn = np(BioGeoBEARS_run_object$trfn)
	#phy = read.tree(file=trfn)
	phy = check_trfn(trfn=trfn)

	# The likelihood of each state at the tips
	# Change this, if you have observations instead of presence/absence at the tips

	# Options:
	# 1. Use tipranges_to_tip_condlikes_of_data_on_each_state ()
	# 2. Use detection model to generate tip likelihoods if desired; or 
	# 3. Take pre-specified tip likelihoods. 
	if (is.null(BioGeoBEARS_run_object$tip_condlikes_of_data_on_each_state) == TRUE)
		{
		if (BioGeoBEARS_run_object$use_detection_model == FALSE)
			{
			#print("here2")
			#print(states_list)

			tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges, phy, states_list=states_list, maxareas=max_numareas, include_null_range=BioGeoBEARS_run_object$include_null_range, useAmbiguities=BioGeoBEARS_run_object$useAmbiguities, trait_as_tip_condlikes=trait_as_tip_condlikes, allow_null_tips=BioGeoBEARS_run_object$allow_null_tips)
			} else {
			# Calculate the initial tip likelihoods, using the detection model
			# Assumes correct order, double-check this
			numareas = length(areas)
			detects_df = BioGeoBEARS_run_object$detects_df
			controls_df = BioGeoBEARS_run_object$controls_df
			mean_frequency = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mf", "init"]
			dp = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["dp", "init"]
			fdp = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["fdp", "init"]
		
			# return_LnLs=TRUE ensures no under-flow
			tip_condlikes_of_data_on_each_state = tiplikes_wDetectionModel(states_list_0based_index=states_list, phy=phy, numareas=numareas, detects_df=detects_df, controls_df=controls_df, mean_frequency=mean_frequency, dp=dp, fdp=fdp, null_range_gets_0_like=TRUE, return_LnLs=TRUE, relative_LnLs=TRUE, exp_LnLs=TRUE, error_check=TRUE)
			}
		} else {
		# Or, use pre-specified tip conditional likelihoods
		# Pre-specified (custom) tip-likelihoods
		tip_condlikes_of_data_on_each_state = BioGeoBEARS_run_object$tip_condlikes_of_data_on_each_state
		} # END if (is.null(BioGeoBEARS_run_object$tip_condlikes_of_data_on_each_state) == FALSE)

	numstates = ncol(tip_condlikes_of_data_on_each_state)

	#print(tip_condlikes_of_data_on_each_state)
	
	if (is.null(BioGeoBEARS_run_object$printlevel))
		{
		BioGeoBEARS_run_object$printlevel = 0
		}
	printlevel = BioGeoBEARS_run_object$printlevel
	
	
	
	
	#######################################################
	# Read the stratification/distances input files, if any
	#######################################################
	#inputs = readfiles_BioGeoBEARS_run(inputs=BioGeoBEARS_run_object)

	#######################################################
	# Check for problems in the input files; will throw stop() if there are problems
	#######################################################
	#check_result = check_BioGeoBEARS_run(inputs=BioGeoBEARS_run_object)
	#check_result
	

	#######################################################
	# Set up the function for optimization
	#######################################################	
	# params are a list of the values of the FREE parameters; but everything is contained in the 
	# BioGeoBEARS_model object at all times
	# (moved to separate function)

	
	# defaults for optimization
	# We are using "L-BFGS-B", which is:
	#####################################################################################################
	# Method "L-BFGS-B" is that of Byrd et. al. (1995) which allows box constraints, that is each
	# variable can be given a lower and/or upper bound. The initial value must satisfy the constraints.
	# This uses a limited-memory modification of the BFGS quasi-Newton method. If non-trivial bounds
	# are supplied, this method will be selected, with a warning.
	# 
	# [...]
	#
	# Byrd, R. H., Lu, P., Nocedal, J. and Zhu, C. (1995) A limited memory algorithm for bound constrained
	# optimization. SIAM J. Scientific Computing, 16, 1190–1208.
	#####################################################################################################
	#
	# "BGFS" refers to: 4 articles, Broyden, Fletcher, Goldfarb and Shanno (1970).
	
	# Run check, before rescaling
	check_BioGeoBEARS_run(BioGeoBEARS_run_object)
	
	#######################################################
	# 2016-03-23_NJM: adding rescaling
	#######################################################
	if (BioGeoBEARS_run_object$rescale_params == TRUE)
		{
		BioGeoBEARS_model_object@params_table = scale_BGB_params(orig_params_table=BioGeoBEARS_model_object@params_table, add_smin=0, add_smax=1)
		
		BioGeoBEARS_run_object$BioGeoBEARS_model_object = BioGeoBEARS_model_object
		}
	
	
	
	params = BioGeoBEARS_model_object_to_init_params(BioGeoBEARS_model_object)
	minj = 1e-05
	# start on Lagrange results
	#params = c(3.11882,  2.51741)
	lower = BioGeoBEARS_model_object_to_params_lower(BioGeoBEARS_model_object)
	upper = BioGeoBEARS_model_object_to_params_upper(BioGeoBEARS_model_object)
	
	# High performance computing
	# HPC using parallel package in R 2.15 or higher, which allows
	# mcmapply (multicore apply)
	# Don't use multicore if using R.app ('AQUA')
	num_cores_to_use = BioGeoBEARS_run_object$num_cores_to_use
	cluster_already_open = BioGeoBEARS_run_object$cluster_already_open

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

		if (is.null(num_cores_to_use))
			{
			num_cores_to_use = num_cores_computer_has
			}

		if (num_cores_to_use > num_cores_computer_has)
			{
			txt = paste0("WARNING from bears_optim_run(): You specified num_cores_to_use=", num_cores_to_use, " cores, but your computer only has ", num_cores_computer_has, ". Resetting to ", num_cores_computer_has, ".")
			cat("\n")
			cat(txt)
			cat("\n")
			warning(txt)
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


	if (force_sparse == TRUE)
		{
		cat("\nNote: force_sparse is set to TRUE; length(states_list)=", length(states_list), "\n", sep="")

		txt = paste0("\n\nNote: sparse matrix exponentiation is being used. When on exponentiation on a branch is completed, 'L',  'R', 'S', or 'U' will print to screen  (for left and right branches, S segments in time-stratified analyses, U for uppass on a segment/branch). This will help you judge the time this analysis will take.  An ML search takes (at least) 100+ downpass calculations of the log-likelihood (lnL) of the tip data, given on left branch, given the tree, model, and parameters. Each downpass requires a matrix exponentiation on each branch of the tree. Your tree has ", length(tr$tip.label), " tips, thus ", length(tr$tip.label)+length(tr$tip.label)-1, " branches. The transition matrix has ", numstates, " states (states=possible geographic ranges), so it would be of size ", numstates, "x", numstates, " if it were dense, but you are using the sparse matrix routine to speed up calculations. Starting now...\n")
		cat(txt)
		} # END if (force_sparse == TRUE)
	
	
	
	#######################################################
	# Check if there are multiple time periods
	#######################################################
	# i.e., timeperiods must exist (not be null and be numeric) and must be of length > 1
	if ( is.numeric(BioGeoBEARS_run_object$timeperiods) ) #&& (length(BioGeoBEARS_run_object$timeperiods) > 1))
		{
		#######################################################
		#######################################################
		# STRATIFIED analysis
		#######################################################
		#######################################################
		# Run optimization on a STRATIFIED tree
		allareas = areas_list
		all_states_list = states_list
		
		use_optimx = BioGeoBEARS_run_object$use_optimx
		
		# USING OPTIM
		if ( (use_optimx == FALSE) || (use_optimx == "optim") )
			{

			
			cat("\n\nNOTE: Before running optim(), here is a test calculation of the data likelihood\nusing calc_loglike_for_optim_stratified() on initial parameter values...\nif this crashes, the error messages are more helpful\nthan those from inside optim().\n", sep="")
			
			inputs = BioGeoBEARS_run_object
			loglike = calc_loglike_for_optim_stratified(params=params, BioGeoBEARS_run_object=inputs, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open)

			cat("\ncalc_loglike_for_optim_stratified() on initial parameters loglike=", loglike, "\n\n\n\nCalculation of likelihood on initial parameters: successful.\n\nNow starting Maximum Likelihood (ML) parameter optimization with optim()...\n\n", sep="")
			
			if (skip_optim == TRUE)
				{
				# Skip optimization
				
				# Skip the optimization, just calculate the log-likelihood from the input parameters
				if (skip_optim_option == "return_loglike")
					{
					cat("Skipping ML search as skip_optim==TRUE. Returning only the log-likelihood.\n\n", sep="")
					return(loglike)
					} 

				# Skip the optimization, just calculate the log-likelihood AND EVERYTHING ELSE from the input parameters
				# (Do this, *if* the skip_optim_option is 'res' (BGB results) or another list)
				if (skip_optim_option == "return_all")
					{
					inputs = BioGeoBEARS_run_object
					cat("Skipping ML search as skip_optim==TRUE. Returning everything (ancestral probabilities etc.) conditional on 'est' parameters given in 'BioGeoBEARS_model_object' .\n\n", sep="")
					optim_result2 = put_params_into_optim_or_optimx_result(BioGeoBEARS_model_object=BioGeoBEARS_run_object$BioGeoBEARS_model_object, total_loglikelihood=loglike, use_optimx=BioGeoBEARS_run_object$use_optimx)
					}
				} else {
				inputs = BioGeoBEARS_run_object
				optim_result2 = optim(par=params, fn=calc_loglike_for_optim_stratified, BioGeoBEARS_run_object=inputs, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=print_optim, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
				} # END if (skip_optim == TRUE)
		
		#optim_result2 = nlminb(start=params, objective=calc_loglike_for_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, lower=lower, upper=upper, control=list(iter.max=50, trace=1, abs.tol=0.001))# method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
			} # END if ( (use_optimx == FALSE) || (use_optimx == "optim") )
		
		# USING OPTIMX
		if ( (use_optimx == TRUE) || (use_optimx == "optimx") )
			{
			# Compare methods with optimx
			#require(optimx)
			
			print_optim = BioGeoBEARS_run_object$print_optim
			
			# For optimx
			# Speedup if desired, using 
			# lower # of generations on optimx
			#speedup = TRUE
			if (speedup)
				{
				# use itnmax, not maxit, for optimx
				
				# IN OPTIM ONLY: default reltol: 
				# sqrt(.Machine$double.eps) = 1.490116e-08 ;
				# this should be the amount of LnL at which it stops
				
				# IN OPTIMX, L-BFGS-B method:
				# factr = controls the convergence of the 
				# "L-BFGS-B" method. Convergence occurs when the 
				# reduction in the objective is within this
				# factor of the machine tolerance. Default is
				# 1e7, that is a tolerance of about 1e-8.
				
				# IN OPTIMX, bobyqa method:
				# no control on tolerance
				
				control_list = list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE)
				# old:
				# factr=0.0001)#, reltol=0.001)#, maxit=100)
				
				# Bogus note (NJM):
				# This causes pathology: reltol=0.001
				# Actually, this was fine, it was 
				# force_sparse = TRUE that was the problem
				# (leads to different results!!  probably rounding errors)
				
				# Limit the number of iterations so it 
				# doesn't go on forever
				num_free_params = sum(BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$"type" == "free")
				num_free_params
				itnmax = 50 * num_free_params
				} else {
				control_list = list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE)
				itnmax = NULL
				} # END if (speedup)
			# For error check, on stratified analysis, just calculate the log-likelihood for hard-coded parameters
			#params = c(0.037, 0.0000000000001)
			#params = c(0.03645000, 4.49500e-08)
			
			
			cat("\n\nNOTE: Before running optimx(), here is a test calculation of the data likelihood\nusing calc_loglike_for_optim_stratified() on initial parameter values...\nif this crashes, the error messages are more helpful\nthan those from inside optimx().\n\n", sep="")
			
			inputs = BioGeoBEARS_run_object
			loglike = calc_loglike_for_optim_stratified(params=params, BioGeoBEARS_run_object=inputs, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open)

			cat("\ncalc_loglike_for_optim_stratified() on initial parameters loglike=", loglike, "\n\n\n\nCalculation of likelihood on initial parameters: successful.\n\nNow starting Maximum Likelihood (ML) parameter optimization with optimx()...\n\n", sep="")

			cat("\n\nPrinting any warnings() that occurred during calc_loglike_for_optim_stratified():\n\n")
			print(warnings())
			
			if (skip_optim == TRUE)
				{
				# Skip optimization
				
				# Skip the optimization, just calculate the log-likelihood from the input parameters
				if (skip_optim_option == "return_loglike")
					{
					cat("Skipping ML search as skip_optim==TRUE. Returning only the log-likelihood.\n\n", sep="")
					return(loglike)
					} 

				# Skip the optimization, just calculate the log-likelihood AND EVERYTHING ELSE from the input parameters
				# (Do this, *if* the skip_optim_option is 'res' (BGB results) or another list)
				if (skip_optim_option == "return_all")
					{
					inputs = BioGeoBEARS_run_object
					cat("Skipping ML search as skip_optim==TRUE. Returning everything (ancestral probabilities etc.) conditional on 'est' parameters given in 'BioGeoBEARS_model_object' .\n\n", sep="")
					optim_result2 = put_params_into_optim_or_optimx_result(BioGeoBEARS_model_object=BioGeoBEARS_run_object$BioGeoBEARS_model_object, total_loglikelihood=loglike, use_optimx=BioGeoBEARS_run_object$use_optimx)
					}
				} else {
				inputs = BioGeoBEARS_run_object
			
				# Run optimx scalecheck
				scalecheck_results = optimx:::scalecheck(par=params, lower=lower, upper=upper)
			
				cat("\n\nResults of optimx:::scalecheck() below. Note: sometimes rescaling parameters may be helpful for ML searches, when the parameters have much different absolute sizes. This can be attempted by setting BioGeoBEARS_run_object$rescale_params = TRUE.\n\n")
				print(scalecheck_results)

				minqa_TF = is.element("minqa", installed.packages()[,1])
				if (minqa_TF == FALSE)
					{
					if (packageVersion("optimx") > 2017)
						{
						txt = "Warning in bears_optim_run(): optimx version 2018.7.10 requires package 'minqa' to do optimx ML optimization with the 'bobyqa' method (optimization with mix/max limits on parameters). However, optimx 2018.7.10 doesn't load 'minqa' automatically, so you may have to do:\n\ninstall.packages('minqa')\n\n...and re-run, to get rid of this warning, and/or the error where optimx returns NA for the parameter inferences after one step, and crashes the resulting uppass calculations."
						cat("\n\n")
						cat(txt)
						cat("\n\n")
						warning(txt)
						require(minqa)
						} # END if (packageVersion("optimx") > 2017)
					} # END if (minqa_TF == FALSE)

				optim_result2 = optimx(par=params, fn=calc_loglike_for_optim_stratified, lower=lower, upper=upper, itnmax=itnmax, method=c("bobyqa"), control=control_list, BioGeoBEARS_run_object=inputs, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=print_optim, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open)# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))
				} # END if (skip_optim == TRUE)

			
			
			
			# print(condlikes_table)
			
			# Run with all methods, for testing:
			# optim_result2 = optimx(par=params, fn=calc_loglike_for_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))
	
			#######################################################
			# Compare optimization routines
			#######################################################
			
			# BEARS_results_7areas_2param
			#                        par   fvalues   method fns grs itns conv  KKT1 KKT2  xtimes
			# 6 0.010165217, 0.009422923 -57.81254   bobyqa  30  NA NULL    0 FALSE TRUE  29.761
			# 1 0.010180679, 0.009492254 -57.81255 L-BFGS-B  33  33 NULL    0 FALSE TRUE 151.137
			# 4   0.01017895, 0.00970706 -57.81263   Rcgmin  32   7 NULL    0 FALSE TRUE  42.461
			# 2 0.010242504, 0.009822486 -57.81284   nlminb  40   6    3    1 FALSE TRUE  43.399
			# 3   0.01017850, 0.01001962 -57.81293      spg  68  NA   47    0 FALSE TRUE 150.932
			# 5               0.01, 0.01 -57.81366   Rvmmin  20   1 NULL    0 FALSE TRUE  21.456
	
	
			#return (optim_result2)
			} # END if ( (use_optimx == TRUE) || (use_optimx == "optimx") )
			
		
		# USING GenSA
		if (use_optimx == "GenSA")
			{
			require(GenSA)

			cat("\n\nNOTE: You are optimizing with GenSA() ('Generalized Simulated Annealing') instead of optimx() or optim(). GenSA may be better for more complex problems (4+ parameters, wildly different scalings), but has not been extensively tested for BioGeoBEARS yet. And it may be slower.")

			
			cat("\n\nNOTE: Before running GenSA(), here is a test calculation of the data likelihood\nusing calc_loglike_for_optim_stratified() on initial parameter values...\nif this crashes, the error messages are more helpful\nthan those from inside GenSA().\n", sep="")
			
			inputs = BioGeoBEARS_run_object
			loglike = calc_loglike_for_optim_stratified(params=params, BioGeoBEARS_run_object=inputs, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=print_optim, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open)

			cat("\ncalc_loglike_for_optim_stratified() on initial parameters loglike=", loglike, "\n\n\n\nCalculation of likelihood on initial parameters: successful.\n\nNow starting Maximum Likelihood (ML) parameter optimization with GenSA()...\n\n", sep="")
			
			if (skip_optim == TRUE)
				{
				# Skip optimization
				
				# Skip the optimization, just calculate the log-likelihood from the input parameters
				if (skip_optim_option == "return_loglike")
					{
					cat("Skipping ML search as skip_optim==TRUE. Returning only the log-likelihood.\n\n", sep="")
					return(loglike)
					} 

				# Skip the optimization, just calculate the log-likelihood AND EVERYTHING ELSE from the input parameters
				# (Do this, *if* the skip_optim_option is 'res' (BGB results) or another list)
				if (skip_optim_option == "return_all")
					{
					inputs = BioGeoBEARS_run_object
					cat("Skipping ML search as skip_optim==TRUE. Returning everything (ancestral probabilities etc.) conditional on 'est' parameters given in 'BioGeoBEARS_model_object' .\n\n", sep="")
					optim_result2 = put_params_into_optim_or_optimx_result(BioGeoBEARS_model_object=BioGeoBEARS_run_object$BioGeoBEARS_model_object, total_loglikelihood=loglike, use_optimx=BioGeoBEARS_run_object$use_optimx)
					}
				} else {
				control_list = list(nb.stop.improvement=50, simple.function=TRUE, trace.mat=TRUE)			

				# NJM: I am assuming that the functions are fairly smooth in BioGeoBEARS analyses
				if (is.null(BioGeoBEARS_run_object$temperature) == FALSE)
					{
					temperature = BioGeoBEARS_run_object$temperature
					control_list = c(control_list, list(temperature=temperature))
					}
				if (is.null(BioGeoBEARS_run_object$max.call) == FALSE)
					{
					max.call = BioGeoBEARS_run_object$max.call
					control_list = c(control_list, list(max.call=max.call))
					} else {
					max.call = length(params) * 250
					control_list = c(control_list, list(max.call=max.call))
					}
			
				inputs = BioGeoBEARS_run_object
				optim_result2 = GenSA(par=params, fn=calc_loglike_for_optim_stratified_neg, BioGeoBEARS_run_object=inputs, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=print_optim, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, lower=lower, upper=upper, control=control_list)
		
			#optim_result2 = nlminb(start=params, objective=calc_loglike_for_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, lower=lower, upper=upper, control=list(iter.max=50, trace=1, abs.tol=0.001))# method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
				} # END if (skip_optim == TRUE)
			} # END if (use_optimx == "GenSA")
		

		
			
			
			
		} else {
		#######################################################
		#######################################################
		# NON-stratified analysis
		#######################################################
		#######################################################
		# Run optimization on a SINGLE tree
		use_optimx = BioGeoBEARS_run_object$use_optimx
		if ( (use_optimx == FALSE) || (use_optimx == "optim") )
			{

			# Un-comment only for error checking, then re-comment!!!!!!!!!!!!!!
			cat("\n\nNOTE: Before running optim(), here is a test calculation of the data likelihood\nusing calc_loglike_for_optim() on initial parameter values...\nif this crashes, the error messages are more helpful\nthan those from inside optim().\n\n", sep="")
			
			loglike = calc_loglike_for_optim(params, BioGeoBEARS_run_object=BioGeoBEARS_run_object, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=BioGeoBEARS_run_object$print_optim, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, return_what="loglike", calc_ancprobs=FALSE)
			
			cat("\ncalc_loglike_for_optim() on initial parameters loglike=", loglike, "\n\n\n\nCalculation of likelihood on initial parameters: successful.\n\nNow starting Maximum Likelihood (ML) parameter optimization with optim()...\n\n", sep="")			

			if (skip_optim == TRUE)
				{
				# Skip optimization
				
				# Skip the optimization, just calculate the log-likelihood from the input parameters
				if (skip_optim_option == "return_loglike")
					{
					cat("Skipping ML search as skip_optim==TRUE. Returning only the log-likelihood.\n\n", sep="")
					return(loglike)
					} 

				# Skip the optimization, just calculate the log-likelihood AND EVERYTHING ELSE from the input parameters
				# (Do this, *if* the skip_optim_option is 'res' (BGB results) or another list)
				if (skip_optim_option == "return_all")
					{
					inputs = BioGeoBEARS_run_object
					cat("Skipping ML search as skip_optim==TRUE. Returning everything (ancestral probabilities etc.) conditional on 'est' parameters given in 'BioGeoBEARS_model_object' .\n\n", sep="")
					optim_result2 = put_params_into_optim_or_optimx_result(BioGeoBEARS_model_object=BioGeoBEARS_run_object$BioGeoBEARS_model_object, total_loglikelihood=loglike, use_optimx=BioGeoBEARS_run_object$use_optimx)
					}
				} else {
				# Try parscale:
				# parscale: A vector of scaling values for the parameters. 
				# Optimization is performed on par/parscale and these should 
				# be comparable in the sense that a unit change in any element 
				# produces about a unit change in the scaled value.For optim.
				# https://www.mail-archive.com/r-help@r-project.org/msg152890.html
				# "(optimx includes parscale on all methods)."
				parscale = (upper - lower) / min(upper - lower)
				print("parscale:")
				print(parscale)

				optim_result2 = optim(par=params, fn=calc_loglike_for_optim, BioGeoBEARS_run_object=BioGeoBEARS_run_object, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=print_optim, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, return_what="loglike", calc_ancprobs=FALSE, method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, parscale=parscale))
		
			#optim_result2 = nlminb(start=params, objective=calc_loglike_for_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, lower=lower, upper=upper, control=list(iter.max=50, trace=1, abs.tol=0.001))# method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
				} # END if (skip_optim == TRUE)
			} # END if ( (use_optimx == FALSE) || (use_optimx == "optim") )
		
		if ( (use_optimx == TRUE) || (use_optimx == "optimx") )
			{
			# Compare methods with optimx
			#require(optimx)

			# For optimx
			# Speedup if desired, using 
			# lower # of generations on optimx
			#speedup = TRUE
			if (speedup)
				{
				# use itnmax, not maxit, for optimx
				
				# IN OPTIM ONLY: default reltol: 
				# sqrt(.Machine$double.eps) = 1.490116e-08 ;
				# this should be the amount of LnL at which it stops
				
				# IN OPTIMX, L-BFGS-B method:
				# factr = controls the convergence of the 
				# "L-BFGS-B" method. Convergence occurs when the 
				# reduction in the objective is within this
				# factor of the machine tolerance. Default is
				# 1e7, that is a tolerance of about 1e-8.
				
				# IN OPTIMX, bobyqa method:
				# no control on tolerance
				
				# Try parscale:
				# parscale: A vector of scaling values for the parameters. 
				# Optimization is performed on par/parscale and these should 
				# be comparable in the sense that a unit change in any element 
				# produces about a unit change in the scaled value.For optim.
				# https://www.mail-archive.com/r-help@r-project.org/msg152890.html
				# "(optimx includes parscale on all methods)."
				parscale = (upper - lower) / min(upper - lower)
				print("parscale:")
				print(parscale)
				control_list = list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE)
				# old:
				# factr=0.0001)#, reltol=0.001)#, maxit=100)
				
				# Bogus note (NJM):
				# This causes pathology: reltol=0.001
				# Actually, this was fine, it was 
				# force_sparse = TRUE that was the problem
				# (leads to different results!!  probably rounding errors)
				
				# Limit the number of iterations so it 
				# doesn't go on forever
				num_free_params = sum(BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$"type" == "free")
				num_free_params
				itnmax = 50 * num_free_params
				} else {
				
				# Try parscale:
				# parscale: A vector of scaling values for the parameters. 
				# Optimization is performed on par/parscale and these should 
				# be comparable in the sense that a unit change in any element 
				# produces about a unit change in the scaled value.For optim.
				# https://www.mail-archive.com/r-help@r-project.org/msg152890.html
				# "(optimx includes parscale on all methods)."
				parscale = (upper - lower) / min(upper - lower)
				print("parscale:")
				print(parscale)
				control_list = list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE)
				itnmax = NULL
				} # END if (speedup)

			# Un-comment only for error checking, then re-comment!!!!!!!!!!!!!!
			cat("\n\nNOTE: Before running optimx(), here is a test calculation of the data likelihood\nusing calc_loglike_for_optim() on initial parameter values...\nif this crashes, the error messages are more helpful\nthan those from inside optimx().\n\n", sep="")

			
			
			loglike = calc_loglike_for_optim(params, BioGeoBEARS_run_object=BioGeoBEARS_run_object, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=BioGeoBEARS_run_object$print_optim, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, return_what="loglike", calc_ancprobs=FALSE)
			
			cat("\ncalc_loglike_for_optim() on initial parameters loglike=", loglike, "\n\n\n\nCalculation of likelihood on initial parameters: successful.\n\nNow starting Maximum Likelihood (ML) parameter optimization with optimx()...\n\n", sep="")

			cat("\n\nPrinting any warnings() that occurred during calc_loglike_for_optim():\n\n")
			print(warnings())


			if (skip_optim == TRUE)
				{
				# Skip optimization
				
				# Skip the optimization, just calculate the log-likelihood from the input parameters
				if (skip_optim_option == "return_loglike")
					{
					cat("Skipping ML search as skip_optim==TRUE. Returning only the log-likelihood.\n\n", sep="")
					return(loglike)
					} 

				# Skip the optimization, just calculate the log-likelihood AND EVERYTHING ELSE from the input parameters
				# (Do this, *if* the skip_optim_option is 'res' (BGB results) or another list)
				if (skip_optim_option == "return_all")
					{
					inputs = BioGeoBEARS_run_object
					cat("Skipping ML search as skip_optim==TRUE. Returning everything (ancestral probabilities etc.) conditional on 'est' parameters given in 'BioGeoBEARS_model_object' .\n\n", sep="")
					optim_result2 = put_params_into_optim_or_optimx_result(BioGeoBEARS_model_object=BioGeoBEARS_run_object$BioGeoBEARS_model_object, total_loglikelihood=loglike, use_optimx=BioGeoBEARS_run_object$use_optimx)
					}
				} else {
				# optimx 2012 versus 2013
				if (packageVersion("optimx") < 2013)
					{
					# optimx 2012
					optim_result2 = optimx(par=params, fn=calc_loglike_for_optim, gr=NULL, hess=NULL, lower=lower, upper=upper, method=c("bobyqa"), itnmax=itnmax, hessian=NULL, control=control_list, BioGeoBEARS_run_object=BioGeoBEARS_run_object, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=print_optim, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, return_what="loglike", calc_ancprobs=FALSE)
					# old:
					# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))
					} else {
					# Run optimx scalecheck
					scalecheck_results = optimx:::scalecheck(par=params, lower=lower, upper=upper)
			
					cat("\n\nResults of optimx:::scalecheck() below. Note: sometimes rescaling parameters may be helpful for ML searches, when the parameters have much different absolute sizes. This can be attempted by setting BioGeoBEARS_run_object$rescale_params = TRUE.\n\n")
					print(scalecheck_results)
					
					# Check if minqa is installed, for the newest optimx (needed for optimx with 'bobyqa' optimizer)
					minqa_TF = is.element("minqa", installed.packages()[,1])
					if (minqa_TF == FALSE)
						{
						if (packageVersion("optimx") > 2017)
							{
							txt = "Warning in bears_optim_run(): optimx version 2018.7.10 requires package 'minqa' to do optimx ML optimization with the 'bobyqa' method (optimization with mix/max limits on parameters). However, optimx 2018.7.10 doesn't load 'minqa' automatically, so you may have to do:\n\ninstall.packages('minqa')\n\n...and re-run, to get rid of this warning, and/or the error where optimx returns NA for the parameter inferences after one step, and crashes the resulting uppass calculations."
							cat("\n\n")
							cat(txt)
							cat("\n\n")
							warning(txt)
							require(minqa)
							} # END if (packageVersion("optimx") > 2017)
						} # END if (minqa_TF == FALSE)

					# optimx 2013
					optim_result2 = optimx(par=params, fn=calc_loglike_for_optim, gr=NULL, hess=NULL, lower=lower, upper=upper, method=c("bobyqa"), itnmax=itnmax, hessian=FALSE, control=control_list, BioGeoBEARS_run_object=BioGeoBEARS_run_object, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=print_optim, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, return_what="loglike", calc_ancprobs=FALSE)
					# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))				
					} # end packageVersion
				} # END if (skip_optim == TRUE)
			
	

	
			# Run with all methods, for testing:
			# optim_result2 = optimx(par=params, fn=calc_loglike_for_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))
	
			#######################################################
			# Compare optimization routines
			#######################################################
			
			# BEARS_results_7areas_2param
			#                        par   fvalues   method fns grs itns conv  KKT1 KKT2  xtimes
			# 6 0.010165217, 0.009422923 -57.81254   bobyqa  30  NA NULL    0 FALSE TRUE  29.761
			# 1 0.010180679, 0.009492254 -57.81255 L-BFGS-B  33  33 NULL    0 FALSE TRUE 151.137
			# 4   0.01017895, 0.00970706 -57.81263   Rcgmin  32   7 NULL    0 FALSE TRUE  42.461
			# 2 0.010242504, 0.009822486 -57.81284   nlminb  40   6    3    1 FALSE TRUE  43.399
			# 3   0.01017850, 0.01001962 -57.81293      spg  68  NA   47    0 FALSE TRUE 150.932
			# 5               0.01, 0.01 -57.81366   Rvmmin  20   1 NULL    0 FALSE TRUE  21.456
	
	
			#return (optim_result2)
			} # END if ( (use_optimx == TRUE) || (use_optimx == "optimx") )
		
		
		# Try GenSA for more complex optimization problems (4+ parameters, or 
		# wildly different parameter scalings)
		if (use_optimx == "GenSA")
			{
			require(GenSA)
			
			cat("\n\nNOTE: You are optimizing with GenSA() ('Generalized Simulated Annealing') instead of optimx() or optim(). GenSA may be better for more complex problems (4+ parameters, wildly different scalings), but has not been extensively tested for BioGeoBEARS yet. And it may be slower.")
			
			# Un-comment only for error checking, then re-comment!!!!!!!!!!!!!!
			cat("\n\nNOTE: Before running GenSA(), here is a test calculation of the data likelihood\nusing calc_loglike_for_optim() on initial parameter values...\nif this crashes, the error messages are more helpful\nthan those from inside GenSA().\n\n", sep="")

			loglike = calc_loglike_for_optim(params, BioGeoBEARS_run_object=BioGeoBEARS_run_object, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=BioGeoBEARS_run_object$print_optim, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, return_what="loglike", calc_ancprobs=FALSE)
			
			cat("\ncalc_loglike_for_optim() on initial parameters loglike=", loglike, "\n\n\n\nCalculation of likelihood on initial parameters: successful.\n\nNow starting Maximum Likelihood (ML) parameter optimization with GenSA()...\n\n", sep="")

			cat("\n\nPrinting any warnings() that occurred during calc_loglike_for_optim():\n\n")
			print(warnings())

			if (skip_optim == TRUE)
				{
				# Skip optimization
				
				# Skip the optimization, just calculate the log-likelihood from the input parameters
				if (skip_optim_option == "return_loglike")
					{
					cat("Skipping ML search as skip_optim==TRUE. Returning only the log-likelihood.\n\n", sep="")
					return(loglike)
					} 

				# Skip the optimization, just calculate the log-likelihood AND EVERYTHING ELSE from the input parameters
				# (Do this, *if* the skip_optim_option is 'res' (BGB results) or another list)
				if (skip_optim_option == "return_all")
					{
					inputs = BioGeoBEARS_run_object
					cat("Skipping ML search as skip_optim==TRUE. Returning everything (ancestral probabilities etc.) conditional on 'est' parameters given in 'BioGeoBEARS_model_object' .\n\n", sep="")
					optim_result2 = put_params_into_optim_or_optimx_result(BioGeoBEARS_model_object=BioGeoBEARS_run_object$BioGeoBEARS_model_object, total_loglikelihood=loglike, use_optimx=BioGeoBEARS_run_object$use_optimx)
					}
				} else {
				control_list = list(nb.stop.improvement=50, simple.function=TRUE, trace.mat=TRUE)			

				# NJM: I am assuming that the functions are fairly smooth in BioGeoBEARS analyses
				if (is.null(BioGeoBEARS_run_object$temperature) == FALSE)
					{
					temperature = BioGeoBEARS_run_object$temperature
					control_list = c(control_list, list(temperature=temperature))
					}
				if (is.null(BioGeoBEARS_run_object$max.call) == FALSE)
					{
					max.call = BioGeoBEARS_run_object$max.call
					control_list = c(control_list, list(max.call=max.call))
					} else {
					max.call = length(params) * 250
					control_list = c(control_list, list(max.call=max.call))
					}
						
				optim_result2 = GenSA(par=params, fn=calc_loglike_for_optim_neg, lower=lower, upper=upper, control=control_list, BioGeoBEARS_run_object=BioGeoBEARS_run_object, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=print_optim, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, return_what="loglike", calc_ancprobs=FALSE)
				} # END if (skip_optim == TRUE)
			} # END if ( (use_optimx == TRUE) || (use_optimx == "optimx") )
		} # END if stratified
	

	
	#######################################################
	# Summarize results 
	#######################################################

	if ((skip_optim == TRUE) && (skip_optim_option == "return_loglike"))
		{
		# Skip optimization
		#cat("Just returning initial loglike as skip_optim==TRUE.\n\n", sep="")
		return(loglike)
		}
	
	# Update the parameter values in the output BioGeoBEARS_model_object using
	# the ML results
	optimx_result = optim_result2
	use_optimx = BioGeoBEARS_run_object$use_optimx
	
	if (printlevel >= 0)
		{
		cat("\n\nThis is the output from optim, optimx, or GenSA. Check the help on those functions to\ninterpret this output and check for convergence issues:\n\n")
		print(optimx_result)
		}
	
	if (printlevel >= 1)
		{
		cat("\n\nReading the optim/optimx/GenSA output into the BioGeoBEARS_model object:\n\nBioGeoBEARS_model_object =\n\n")
		}

	BioGeoBEARS_model_object = update_BioGeoBEARS_model_object_w_optimx_result(BioGeoBEARS_model_object, optimx_result, use_optimx)

	
	# ERROR CHECK
	if (any(is.na(BioGeoBEARS_model_object@params_table$est)) == TRUE)
		{
		txt = "STOP ERROR in bears_optim_run(). For some reason, your ML optimizer returned one or more NA / NaN values for the estimated parameters. Probably this is a version conflict with an update to one of the optimizer functions/packages (e.g., optim, optimx, minqa, GenSA. Printing BioGeoBEARS_model_object@params_table to screen, below.  Email the BioGeoBEARS Google Group if you cannot figure out the problem."
		
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		cat("BioGeoBEARS_model_object@params_table:\n\n")
		print(BioGeoBEARS_model_object@params_table)
		cat("\n\n")
		stop(txt)
		} # END if (any(is.na(BioGeoBEARS_model_object@params_table$est)) == TRUE)
	
	
	######################################################
	# 2016-03-23_NJM: adding rescaling
	# (unscale params, if they were used before)
	######################################################
	if (BioGeoBEARS_run_object$rescale_params == TRUE)
		{
		cat("\n(Because BioGeoBEARS_run_object$rescale_params == TRUE, using unscale_BGB_params() to return parameter estimates to the original scaling...\n")
		BioGeoBEARS_model_object@params_table = unscale_BGB_params(scaled_params_table=BioGeoBEARS_model_object@params_table)
		
		if (BioGeoBEARS_run_object$use_optimx == FALSE)
			{
			optim_result2$par = BioGeoBEARS_model_object@params_table$est[BioGeoBEARS_model_object@params_table$type=="free"]
			} # END if (BioGeoBEARS_run_object$use_optimx == FALSE)

		if ( (BioGeoBEARS_run_object$use_optimx == TRUE) || (BioGeoBEARS_run_object$use_optimx == "optimx") )
			{
			# optimx 2013+
			if (packageVersion("optimx") >= 2013)
				{
				param_names = names(optim_result2)
				param_1st_letter = substr(x=param_names, start=1, stop=1)
				param_TF = param_1st_letter == "p"
				param_names = param_names[param_TF]
			
				optim_result2[param_names] = BioGeoBEARS_model_object@params_table$est[BioGeoBEARS_model_object@params_table$type=="free"]
				}
			# optimx 2012
			if (packageVersion("optimx") < 2013)
				{
				optim_result2$par[[1]] = BioGeoBEARS_model_object@params_table$est[BioGeoBEARS_model_object@params_table$type=="free"]
				}
			} # END if (BioGeoBEARS_run_object$use_optimx == TRUE)

		#cat("...done.)\n\n")
		} # END if (BioGeoBEARS_run_object$rescale_params == TRUE)


	if (printlevel >= 1)
		{
		print(BioGeoBEARS_model_object)
		}	

	# Update the output
	BioGeoBEARS_run_object$BioGeoBEARS_model_object = BioGeoBEARS_model_object
	
	# Set the dispersal and extinction rate
# 	d = BioGeoBEARS_model_object@params_table["d","est"]
# 	e = BioGeoBEARS_model_object@params_table["e","est"]
# 	a = BioGeoBEARS_model_object@params_table["a","est"]
# 	
	# Set the branch length exponent 
	# NO DON'T DO THIS HERE, IT GETS DONE IN
	# calc_loglike_for_optim()
	#b = BioGeoBEARS_model_object@params_table["b","est"]
	#phy$edge.length = phy$edge.length ^ b


	# Equal dispersal in all directions (unconstrained)
	# Equal extinction probability for all areas
#	dispersal_multipliers_matrix = matrix(1, nrow=length(areas), ncol=length(areas))

	# Multiply d by dispersal_multipliers_matrix (for relative distance)
# 	dmat_times_d = dispersal_multipliers_matrix * matrix(d, nrow=length(areas), ncol=length(areas))
# 	elist = rep(e, length(areas))
# 	amat = dispersal_multipliers_matrix * matrix(a, nrow=length(areas), ncol=length(areas))
	
	# Set up the instantaneous rate matrix (Q matrix)
#	Qmat = rcpp_states_list_to_DEmat(areas_list=areas_list, states_list=states_list, dmat=dmat_times_d, elist=elist, amat=amat, include_null_range=BioGeoBEARS_run_object$include_null_range, normalize_TF=TRUE, makeCOO_TF=force_sparse)


	#######################################################
	# Cladogenic model
	#######################################################
# 	j = BioGeoBEARS_model_object@params_table["j","est"]
# 	ysv = BioGeoBEARS_model_object@params_table["ysv","est"]
# 	v = BioGeoBEARS_model_object@params_table["v","est"]
# 	ys = BioGeoBEARS_model_object@params_table["ys","est"]
# 	y = BioGeoBEARS_model_object@params_table["y","est"]
# 	s = BioGeoBEARS_model_object@params_table["s","est"]
# 	sum_SPweights = y + s + j + v
# 
#	maxent_constraint_01 = BioGeoBEARS_model_object@params_table["mx01","est"]
	
	# Text version of speciation matrix	
#	maxent_constraint_01v = BioGeoBEARS_model_object@params_table["mx01v","est"]
	#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
		
	# Set the parameter controlling the size distribution of 
	# the smaller descendant species
# 	maxent01s_param = BioGeoBEARS_model_object@params_table["mx01s","est"]
# 	maxent01v_param = BioGeoBEARS_model_object@params_table["mx01v","est"]
# 	maxent01j_param = BioGeoBEARS_model_object@params_table["mx01j","est"]
# 	maxent01y_param = BioGeoBEARS_model_object@params_table["mx01y","est"]
# 
# 
# 	# Cladogenesis model inputs
# 	spPmat_inputs = NULL
	
	# This dmat is for dispersal multipliers, i.e. to apply to j events, 
	# NOT the dmat_times_d derived from the d parameter;
	# make sure there aren't others elsewhere!
#	dmat = dispersal_multipliers_matrix
# 	spPmat_inputs$dmat = dmat
# 
# 	states_indices = states_list
	
	# shorten the states_indices by 1 (cutting the 
	# null range state from the speciation matrix)
# 	if (BioGeoBEARS_run_object$include_null_range == TRUE)
# 		{
# 		states_indices[1] = NULL
# 		} # END if (include_null_range == TRUE)
# 	spPmat_inputs$l = states_indices
# 	spPmat_inputs$s = s
# 	spPmat_inputs$v = v
# 	spPmat_inputs$j = j
# 	spPmat_inputs$y = y
# 	spPmat_inputs$maxent01s_param = maxent01s_param
# 	spPmat_inputs$maxent01v_param = maxent01v_param
# 	spPmat_inputs$maxent01j_param = maxent01j_param
# 	spPmat_inputs$maxent01y_param = maxent01y_param

	outputs = BioGeoBEARS_model_object

	if ((is.numeric(BioGeoBEARS_run_object$timeperiods))) #&& (length(BioGeoBEARS_run_object$timeperiods) > 1))
		{
		# We need to put the params back into the inputs 
		# to get the reconstructed ancestors etc.
		# Note that fixlikes SHOULD be included here in the 
		# final results, if specified by the user at the beginning
		# (thanks to Julien for pointing out issue)
		return_condlikes_table = BioGeoBEARS_run_object$return_condlikes_table
		calc_ancprobs = BioGeoBEARS_run_object$calc_ancprobs
		
		fixnode = BioGeoBEARS_run_object$fixnode
		fixlikes = BioGeoBEARS_run_object$fixlikes
		
		# Need to store the model parameters in an inputs object to pass to calc_loglike_sp_stratified
		inputs2 = BioGeoBEARS_run_object
		inputs2$BioGeoBEARS_model_object = BioGeoBEARS_model_object
		
		calc_TTL_loglike_from_condlikes_table = BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table
		
		model_results = calc_loglike_sp_stratified(tip_condlikes_of_data_on_each_state, phy, Qmat=NULL, spPmat=NULL, min_branchlength=min_branchlength, return_what="all", probs_of_states_at_root=NULL, rootedge=TRUE, sparse=force_sparse, printlevel=BioGeoBEARS_run_object$printlevel, use_cpp=TRUE, input_is_COO=FALSE, spPmat_inputs=NULL, cppSpMethod=3, cluster_already_open=cluster_already_open, calc_ancprobs=calc_ancprobs, include_null_range=BioGeoBEARS_run_object$include_null_range, fixnode=fixnode, fixlikes=fixlikes, inputs=inputs2, allareas=allareas, all_states_list=all_states_list, return_condlikes_table=return_condlikes_table, calc_TTL_loglike_from_condlikes_table=calc_TTL_loglike_from_condlikes_table)
		} else {
		#print(params)
		params = BioGeoBEARS_model_object_to_est_params(BioGeoBEARS_model_object)
		#print(params)
		# Calculate the log-likelihood of the data, given the model parameters during this iteration	
		#model_results = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="all", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=1, calc_ancprobs=TRUE, include_null_range=BioGeoBEARS_run_object$include_null_range)

		# We need to put the params back into the inputs 
		# to get the reconstructed ancestors etc.
		# Note that fixlikes SHOULD be included here in the 
		# final results, if specified by the user at the beginning
		# (thanks to Julien for pointing out issue)
		calc_ancprobs = BioGeoBEARS_run_object$calc_ancprobs
		# (originally, to do local ancestral states, you would set
		#  calc_ancprobs to FALSE and use the subsequent optim_result
		#  $ fvalue to get the LnL optimal on that node state.
		#  I am now changing it to always use the fixlikes.)
		
		#print("Calculating final LnL...")
		model_results = calc_loglike_for_optim(params=params, BioGeoBEARS_run_object=BioGeoBEARS_run_object, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=BioGeoBEARS_run_object$print_optim, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, return_what="all", calc_ancprobs=calc_ancprobs)
		#print("model_results:")
		#print(model_results)
		}


	if (cluster_was_open == FALSE)
		{
		if (exists("cluster_open") && (cluster_open == TRUE))
			{
			cat("\n\nStopping cluster with ", num_cores_to_use, " cores.\n\n", sep="")
			stopCluster(cluster_already_open)
			}
		}
	
	
	
	
	#######################################################
	# Store results in a BEARS result
	#######################################################
	bears_output = model_results
	bears_output$inputs = BioGeoBEARS_run_object
	#bears_output$spPmat_inputs = spPmat_inputs
	bears_output$outputs = outputs
	bears_output$optim_result = optim_result2
	
	return(bears_output)
	} # END bears_optim_run <- function(BioGeoBEARS_run_object = define_BioGeoBEARS_run(), skip_optim=FALSE)


