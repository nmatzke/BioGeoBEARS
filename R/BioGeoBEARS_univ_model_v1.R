require("ape")
require("rexpokit")
require("cladoRcpp")



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
#' @param BioGeoBEARS_model_object Object containing the model.
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
#' @return \code{ttl_loglike} The log-likelihood of the data under the input model and parameters.
#' @export
#' @seealso \code{\link[stats]{convolve}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
calc_loglike_for_optim <- function(params, BioGeoBEARS_model_object, phy, tip_condlikes_of_data_on_each_state, print_optim=TRUE, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open)
	{
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
	distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

	dmat = matrix(d, nrow=length(areas), ncol=length(areas))
	elist = rep(e, length(areas))
	
	# Set up the instantaneous rate matrix (Q matrix)
	Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)

	# Cladogenic model
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
	spPmat_inputs$dmat = distances_mat
	spPmat_inputs$maxent01s_param = maxent01s_param
	spPmat_inputs$maxent01v_param = maxent01v_param
	spPmat_inputs$maxent01j_param = maxent01j_param
	spPmat_inputs$maxent01y_param = maxent01y_param


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
	ttl_loglike = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="loglike", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, cluster_already_open=cluster_already_open)
	ttl_loglike

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









# 6-parameter model, adding j, v (vicariance proportion), maxent_constraint_01 (for non-vicariant subsets), maxent_constraint_01v (weighting for size of smaller offspring)
#######################################################
# bears_optim_run
#######################################################
#' Run ML search from \code{BioGeoBEARS_run} object
#' 
#' Uses a BioGeoBEARS_run_object to simplify input.
#'
#' @param BioGeoBEARS_run_object Contains all inputs
#' @return \code{bears_output} A list of outputs.  bears_output$optim_result
#' @export
#' @seealso \code{\link{bears_2param_standard_fast}}, \code{\link[cladoRcpp]{numstates_from_numareas}}, \code{\link{getranges_from_LagrangePHYLIP}}, \code{\link[ape]{read.tree}}, \code{\link{calc_loglike_sp}}
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
#' extdata_dir = system.file("extdata/", package="BioGeoBEARS")
#' # tmp hard code: extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"
#'
#' # Set the filenames (Hawaiian Psychotria from Ree & Smith 2008)
#' trfn = paste(extdata_dir, "Psychotria_5.2.newick", sep="")
#' tr = read.tree(file=trfn)
#' 
#' geogfn = paste(extdata_dir, "Psychotria_geog.data", sep="")
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
bears_optim_run <- function(BioGeoBEARS_run_object = define_BioGeoBEARS_run())
	{
	defaults='	
	BioGeoBEARS_run_object = define_BioGeoBEARS_run()
	BioGeoBEARS_run_object
	'
	
	require(cladoRcpp)
	require(rexpokit)
	
	inputs = BioGeoBEARS_run_object

	
	# Load the model object
	BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object

	# Get geographic ranges at tips
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=BioGeoBEARS_run_object$geogfn)


	# Get the list of geographic areas
	areas = getareas_from_tipranges_object(tipranges)
	areas_list = seq(0, length(areas)-1, 1)		# 0-base indexes

	# Change the names to tipranges@df:
	# this doesn't make sense if areas_list is 0-based indexes
	names(tipranges@df) = areas_list

	#######################################################
	# Set the maximum range size (can be thought of as
	# a free parameter
	#######################################################
	if (is.na(BioGeoBEARS_run_object$max_range_size))
		{
		if (is.na(BioGeoBEARS_run_object$states_list))
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
	
	# Old/slow way of getting the list of states and speciation matrix (slow)
	# old_states_list = areas_list_to_states_list(areas, include_null_range=TRUE)
	# old_states_list
	# spmat = make_relprob_matrix_bi(old_states_list)
	# spmat
	
	# max_numareas = max(sapply(X=old_states_list, FUN=length), na.rm=TRUE)
	# max_numareas


	# Take the list of areas, and get list of possible states
	# (the user can manually input states if they like)
	if (is.na(BioGeoBEARS_run_object$states_list))
		{
		states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
		states_list
		} else {
		states_list = BioGeoBEARS_run_object$states_list
		}

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

	if (force_sparse == TRUE)
		{
		cat("\nNote: force_sparse is set to TRUE; length(states_list)=", length(states_list), "\n", sep="")
		}
	
	
	
	#######################################################
	# Load the phylogenetic tree
	#######################################################
	trfn = BioGeoBEARS_run_object$trfn
	phy = read.tree(file=trfn)

	# The likelihood of each state at the tips
	# Change this, if you have observations instead of presence/absence at the tips
	tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges, phy, states_list=states_list, maxareas=max_numareas)
	tip_condlikes_of_data_on_each_state

	

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
	if (.Platform$GUI != "AQUA" && ((is.na(num_cores_to_use) == TRUE) || (num_cores_to_use > 1)) )
		{
		# We are doing manual, optional processing on several cores;
		# this seems to have less overhead/hassle/incompatibility issues
		# than using mcmapply, mclapply, etc...
		require("parallel") #<- do this higher up

		num_cores_computer_has = detectCores()
		
		if (is.null(num_cores_to_use))
			{
			num_cores_to_use = num_cores_computer_has
			}

		# Don't do this, if the cluster is already open
		cat("\nYour computer has ", num_cores_computer_has, " cores. You have chosen to use:\nnum_cores_to_use = ", num_cores_to_use, " cores for the matrix exponentiations in the likelihood calculations.\n", sep="")

		cluster_already_open = makeCluster(rep("localhost",num_cores_to_use), type = "SOCK")
		cat("Started cluster with ", num_cores_to_use, " cores.\n\n", sep="")
		
		# Flag so that you remember to close cluster at the end
		cluster_open=TRUE
		} else {
		cluster_already_open = NULL
		}

	
	
	#######################################################
	# Check if there are multiple time periods
	#######################################################
	# i.e., timeperiods must exist (not be null and be numeric) and must be of length > 1
	if ((is.numeric(inputs$timeperiods))) #&& (length(inputs$timeperiods) > 1))
		{
		# Run optimization on a STRATIFIED tree
		allareas = areas_list
		all_states_list = states_list
		
		use_optimx = BioGeoBEARS_run_object$use_optimx
		if (use_optimx == FALSE)
			{
			optim_result2 = optim(par=params, fn=calc_loglike_for_optim_stratified, BioGeoBEARS_run_object=inputs, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
		
		#optim_result2 = nlminb(start=params, objective=calc_loglike_for_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, lower=lower, upper=upper, control=list(iter.max=50, trace=1, abs.tol=0.001))# method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
			} else {
			# Compare methods with optimx
			require(optimx)
			
			optim_result2 = optimx::optimx(par=params, fn=calc_loglike_for_optim_stratified, lower=lower, upper=upper, BioGeoBEARS_run_object=inputs, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))
	
			# Run with all methods, for testing:
			# optim_result2 = optimx::optimx(par=params, fn=calc_loglike_for_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))
	
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
			}
		} else {
		
		# Run optimization on a SINGLE tree
		use_optimx = BioGeoBEARS_run_object$use_optimx
		if (use_optimx == FALSE)
			{
			optim_result2 = optim(par=params, fn=calc_loglike_for_optim, BioGeoBEARS_model_object=BioGeoBEARS_model_object, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
		
		#optim_result2 = nlminb(start=params, objective=calc_loglike_for_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, lower=lower, upper=upper, control=list(iter.max=50, trace=1, abs.tol=0.001))# method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
			} else {
			# Compare methods with optimx
			require(optimx)
			
			optim_result2 = optimx::optimx(par=params, fn=calc_loglike_for_optim, BioGeoBEARS_model_object=BioGeoBEARS_model_object, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))
	
			# Run with all methods, for testing:
			# optim_result2 = optimx::optimx(par=params, fn=calc_loglike_for_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))
	
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
			}

		}
	

	
	#######################################################
	# Summarize results 
	#######################################################

	# Set the dispersal and extinction rate
	if (use_optimx == FALSE)
		{
		BioGeoBEARS_model_object = params_into_BioGeoBEARS_model_object(BioGeoBEARS_model_object, params=optim_result2$par)
		BioGeoBEARS_model_object = calc_linked_params_BioGeoBEARS_model_object(BioGeoBEARS_model_object)
		# # Set the dispersal and extinction rate
		# 			d = optim_result2$par[1]
		# 			e = optim_result2$par[2]
		# 			j = optim_result2$par[3]
		# 			v = optim_result2$par[4]
		# 			maxent01_param = optim_result2$par[5]
		# 			maxent01v_param = optim_result2$par[6]
		} else {
		BioGeoBEARS_model_object = params_into_BioGeoBEARS_model_object(BioGeoBEARS_model_object, params=optim_result2$par[[1]])
		BioGeoBEARS_model_object = calc_linked_params_BioGeoBEARS_model_object(BioGeoBEARS_model_object)
# 		d = optim_result2$par[[1]][1]
# 		e = optim_result2$par[[1]][2]
# 		j = optim_result2$par[[1]][3]
# 		v = optim_result2$par[[1]][4]
# 		maxent01_param = optim_result2$par[[1]][5]
# 		maxent01v_param = optim_result2$par[[1]][6]
		}
	
	
	# Update the output
	BioGeoBEARS_run_object$BioGeoBEARS_model_object = BioGeoBEARS_model_object
	
	# Set the dispersal and extinction rate
	d = BioGeoBEARS_model_object@params_table["d","est"]
	e = BioGeoBEARS_model_object@params_table["e","est"]

	# Equal dispersal in all directions (unconstrained)
	# Equal extinction probability for all areas
	dispersal_multipliers_matrix = matrix(1, nrow=length(areas), ncol=length(areas))

	dmat = matrix(d, nrow=length(areas), ncol=length(areas)) * dispersal_multipliers_matrix
	elist = rep(e, length(areas))
	
	# Set up the instantaneous rate matrix (Q matrix)
	Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)

	# Cladogenic model
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
	spPmat_inputs$dmat = dmat

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

	inputs$BioGeoBEARS_model_object = BioGeoBEARS_model_object

	if ((is.numeric(inputs$timeperiods))) #&& (length(inputs$timeperiods) > 1))
		{
		model_results = calc_loglike_sp_stratified(tip_condlikes_of_data_on_each_state, phy, Qmat=NULL, spPmat=NULL, min_branchlength=1e-21, return_what="loglike", probs_of_states_at_root=NULL, rootedge=TRUE, sparse=force_sparse, printlevel=0, use_cpp=TRUE, input_is_COO=FALSE, spPmat_inputs=NULL, cppSpMethod=3, cluster_already_open=cluster_already_open, calc_ancprobs=TRUE, null_range_allowed=TRUE, fixnode=NULL, fixlikes=NULL, inputs=inputs, allareas=allareas, all_states_list=all_states_list)
		} else {
		# Calculate the log-likelihood of the data, given the model parameters during this iteration	
		model_results = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="all", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, calc_ancprobs=TRUE)
		}



	if (exists("cluster_open") && (cluster_open == TRUE))
		{
		cat("\n\nStopping cluster with ", num_cores_to_use, " cores.\n\n", sep="")
		stopCluster(cluster_already_open)
		}
	
	
	
	
	#######################################################
	# Store results in a BEARS result
	#######################################################
	bears_output = model_results
	bears_output$inputs = inputs
	bears_output$spPmat_inputs = spPmat_inputs
	bears_output$outputs = BioGeoBEARS_run_object
	bears_output$optim_result = optim_result2
	
	return(bears_output)
	}


