#######################################################
# Functions in this file implement various standard 
# models of interest
#######################################################
require("ape")
require("rexpokit")
require("cladoRcpp")



#######################################################
# bears_2param_standard_fast
#######################################################
#' 2-parameter model, fixed cladogenesis model (as in LAGRANGE)
#' 
#' This function implements a biogeographical model with 2 free parameters (\emph{d}, rate of dispersal/range addition, and
#' \emph{e}, rate of extinction/range contraction), and a fixed cladogenesis model with equal probability of vicariance, sympatric-subset,
#' and sympatric-range-copying events, and with the smaller descendant always having a range size of 1 area.  Once the model is set up,
#' it is input into the optimization routine \code{\link[optimx]{optimx}} (the more common \code{\link[stats]{optim}} can also be used
#' by editing the function), and \code{\link{calc_loglike_sp}} is used to calculate the log-likelihood of 
#' each set of parameters.  Once the parameter values that give the data the maximum likelihood are found, they are reported back to the function
#' and returned to the user.
#'
#' This duplicates the model used in the standard LAGRANGE implementation (\cite{ReeSmith2008}, \cite{Ree2009configurator}, \cite{SmithRee2010_CPPversion},
#' with no constraints on dispersal or range size.
#' 
#' Here, all of the fastest processing options have been used.
#'
#' Model implementations are provided to show the user how a specific model can be set up and optimized.  This is preferable compared to the "black-box"
#' nature of most other inference packages.  Users are encouraged to experiment.  Useful models can be added to later versions of \code{BioGeoBEARS}.
#'
#' @param trfn The filename of the phylogenetic tree, in NEWICK format (\url{http://evolution.genetics.washington.edu/phylip/newicktree.html}).  
#' Tipnames should match the names in geogfn.  See \code{\link[ape]{read.tree}} in APE for reading in phylogenetic trees.
#' @param geogfn A PHYLIP-style file with geographic range data (see \code{\link{getranges_from_LagrangePHYLIP}}) for each tipname. This is the same format
#' used by C++ LAGRANGE (\cite{SmithRee2010_CPPversion}).
#' @param max_range_size The maximum rangesize, in number of areas.  Having a smaller maximum range size means that you can have more areas (the size of the
#' state space is greatly reduced; see \code{\link[cladoRcpp]{numstates_from_numareas}}.
#' @param num_cores_to_use If >1, parallel processing will be attempted. \bold{Note:} parallel processing via \code{library (parallel)} will work in Mac command-line
#' R, but not in Mac GUI \code{R.app}.
#' @return \code{bears_output} A list of outputs.  bears_output$optim_result
#' @export
#' @seealso \code{\link[cladoRcpp]{numstates_from_numareas}}, \code{\link{getranges_from_LagrangePHYLIP}}, \code{\link[ape]{read.tree}}, \code{\link{calc_loglike_sp}}
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
#' bears_output = bears_2param_standard_fast(trfn=trfn, geogfn=geogfn)
#' bears_output
#' }
#'
bears_2param_standard_fast <- function(trfn = "Psychotria_5.2.newick", geogfn = "Psychotria_geog.data", max_range_size=NULL,  num_cores_to_use=NULL)
	{
	defaults='
	trfn="Psychotria_5.2.newick"
	geogfn="Psychotria_geog.data"
	max_range_size=NULL
	num_cores_to_use=NULL
	'
	
	require(cladoRcpp)
	require(rexpokit)




	# Get geographic ranges at tips
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)

	# Get the list of geographic areas
	areas = getareas_from_tipranges_object(tipranges)
	areas_list = seq(0, length(areas)-1, 1)

	# Change the names to tipranges@df:
	names(tipranges@df) = areas_list


	#######################################################
	# Set the maximum range size (can be thought of as
	# a free parameter
	#######################################################
	if (is.null(max_range_size))
		{
		max_range_size = length(areas)
		need_to_fix_rangesize = TRUE
		} else {
		need_to_fix_rangesize = FALSE
		}

	
	# Old/slow way of getting the list of states and speciation matrix (slow)
	# old_states_list = areas_list_to_states_list(areas, include_null_range=TRUE)
	# old_states_list
	# spmat = make_relprob_matrix_bi(old_states_list)
	# spmat
	
	# max_numareas = max(sapply(X=old_states_list, FUN=length), na.rm=TRUE)
	# max_numareas


	# Take the list of areas, and get list of possible states
	states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
	states_list

	if (need_to_fix_rangesize == FALSE)
		{
		max_numareas = max(sapply(X=states_list, FUN=length), na.rm=TRUE)
		max_numareas

		# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
		if (max_numareas >= 8)
			{
			force_sparse = TRUE
			} else {
			force_sparse = FALSE
			}	

		} else {
		max_numareas = max_range_size
		}
	# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
	if (max_numareas >= 8)
		{
		force_sparse = TRUE
		} else {
		force_sparse = FALSE
		}	
	
	
	# Load the phylogenetic tree
	phy = read.tree(file=trfn)


	# The likelihood of each state at the tips
	tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges, phy, states_list=states_list, maxareas=max_numareas)
	tip_condlikes_of_data_on_each_state

	maxent_constraint_01 = 0.0001


	#######################################################
	# Set up the function for optimization
	#######################################################	
	function_to_optim <- function(params, phy, tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open)
		{
		junk='
		cluster_already_open=NULL
		'
		
		# Lagrange-like (all new species are within the geographic range of the ancestor)

		# Set the dispersal and extinction rate
		d = params[1]
		e = params[2]

		# Equal dispersal in all directions (unconstrained)
		# Equal extinction probability for all areas
		distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

		dmat = matrix(d, nrow=length(areas), ncol=length(areas))
		elist = rep(e, length(areas))
		
		# Set up the instantaneous rate matrix (Q matrix)
		Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)
	
		# Cladogenic model
		ys = 0.5
		j = 0
		v = 0.5
		sum_SPweights = ys + j + v

		# Text version of speciation matrix	
		maxent_constraint_01v = maxent_constraint_01
		#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
			
		# Set the parameter controlling the size distribution of 
		# the smaller descendant species
		maxent01s_param = maxent_constraint_01
		maxent01v_param = maxent_constraint_01
		maxent01j_param = maxent_constraint_01
		maxent01y_param = maxent_constraint_01


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
			# Before calculating the log likelihood, print it, in case there is e.g. a bug
			cat("d=", d, "; e=", e, "; j=", j, "; ys=", ys, "; v=", v, "; sum=", sum_SPweights, "; LnL=", sep="")
			}

		# Calculate the log-likelihood of the data, given the model parameters during this iteration	
		ttl_loglike = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="loglike", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, cluster_already_open=cluster_already_open)
		ttl_loglike
	
		if (print_optim == TRUE)
			{
			# If the log likelihood is successful, print it
			cat(ttl_loglike, "\n", sep="")
			}

		#return(-1*ttl_loglike)
		return(ttl_loglike)
		}

	
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
	params = c(0.01, 0.01)
	minj = 1e-05
	# start on Lagrange results
	#params = c(3.11882,  2.51741)
	lower = c(1e-15,1e-15)
	upper = c(1,1)
	
	# High performance computing
	# HPC using parallel package in R 2.15 or higher, which allows
	# mcmapply (multicore apply)
	# Don't use multicore if using R.app ('AQUA')
	if (.Platform$GUI != "AQUA" && ((is.null(num_cores_to_use) == TRUE) || (num_cores_to_use > 1)) )
		{
		# We are doing manual, optional processing on several cores;
		# this seems to have less overhead/hassle/incompatibility issues
		# than using mcmapply, mclapply, etc...
		#require("parallel") #<- do this higher up

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
		
	
	
	
	# Run optimization	
	use_optimx = TRUE
	if (use_optimx == FALSE)
		{
		optim_result2 = optim(par=params, fn=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
	
	#optim_result2 = nlminb(start=params, objective=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, lower=lower, upper=upper, control=list(iter.max=50, trace=1, abs.tol=0.001))# method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
		} else {
		# Compare methods with optimx
		#require(optimx)
		
		
		
		optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

		# Run with all methods, for testing:
		# optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

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
	

	if (exists("cluster_open") && (exists("cluster_open")))
		{
		cat("\n\nStopping cluster with ", num_cores_to_use, " cores.\n\n", sep="")
		stopCluster(cluster_already_open)
		}
	
	#######################################################
	# Summarize results 
	#######################################################

	if (use_optimx == FALSE)
		{
		# Set the dispersal and extinction rate
		d = optim_result2$par[1]
		e = optim_result2$par[2]
		} else {
		d = optim_result2$par[[1]][1]
		e = optim_result2$par[[1]][2]
		}

	# Equal dispersal in all directions (unconstrained)
	# Equal extinction probability for all areas
	distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

	dmat = matrix(d, nrow=length(areas), ncol=length(areas))
	elist = rep(e, length(areas))
	
	# Set up the instantaneous rate matrix (Q matrix)
	Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)

	# Cladogenic model
	ys = 0.5
	j = 0
	v = 0.5
	sum_SPweights = ys + j + v

	# Text version of speciation matrix	
	maxent_constraint_01v = maxent_constraint_01
		
	# Set the parameter controlling the size distribution of 
	# the smaller descendant species
	maxent01s_param = maxent_constraint_01
	maxent01v_param = maxent_constraint_01
	maxent01j_param = maxent_constraint_01
	maxent01y_param = maxent_constraint_01

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


	# Calculate the log-likelihood of the data, given the model parameters during this iteration	
	model_results = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="all", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, calc_ancprobs=TRUE)



	#######################################################
	# Store results in a BEARS result
	#######################################################
	bears_output = model_results
	bears_output$optim_result = optim_result2
	bears_output$spPmat_inputs = spPmat_inputs
	
	return(bears_output)
	}



#######################################################
# bears_2param_standard_fast
#######################################################
#' 2-parameter model, fixed cladogenesis model (as in LAGRANGE)
#' 
#' This function implements a biogeographical model with 2 free parameters (\emph{d}, rate of dispersal/range addition, and
#' \emph{e}, rate of extinction/range contraction), and a fixed cladogenesis model with equal probability of vicariance, sympatric-subset,
#' and sympatric-range-copying events, and with the smaller descendant always having a range size of 1 area.  Once the model is set up,
#' it is input into the optimization routine \code{\link[optimx]{optimx}} (the more common \code{\link[stats]{optim}} can also be used
#' by editing the function), and \code{\link{calc_loglike_sp}} is used to calculate the log-likelihood of 
#' each set of parameters.  Once the parameter values that give the data the maximum likelihood are found, they are reported back to the function
#' and returned to the user.
#'
#' This duplicates the model used in the standard LAGRANGE implementation (\cite{ReeSmith2008}, \cite{Ree2009configurator}, \cite{SmithRee2010_CPPversion},
#' with no constraints on dispersal or range size.
#' 
#' Here, all of the fastest processing options have been used.
#'
#' Model implementations are provided to show the user how a specific model can be set up and optimized.  This is preferable compared to the "black-box"
#' nature of most other inference packages.  Users are encouraged to experiment.  Useful models can be added to later versions of \code{BioGeoBEARS}.
#'
#' @param trfn The filename of the phylogenetic tree, in NEWICK format (\url{http://evolution.genetics.washington.edu/phylip/newicktree.html}).  
#' Tipnames should match the names in geogfn.  See \code{\link[ape]{read.tree}} in APE for reading in phylogenetic trees.
#' @param geogfn A PHYLIP-style file with geographic range data (see \code{\link{getranges_from_LagrangePHYLIP}}) for each tipname. This is the same format
#' used by C++ LAGRANGE (\cite{SmithRee2010_CPPversion}).
#' @param max_range_size The maximum rangesize, in number of areas.  Having a smaller maximum range size means that you can have more areas (the size of the
#' state space is greatly reduced; see \code{\link[cladoRcpp]{numstates_from_numareas}}.
#' @param num_cores_to_use If >1, parallel processing will be attempted. \bold{Note:} parallel processing via \code{library (parallel)} will work in Mac command-line
#' R, but not in Mac GUI \code{R.app}.
#' @param fixnode If the state at a particular node is going to be fixed (e.g. for ML marginal ancestral states), give the node number.
#' @param fixlikes The state likelihoods to be used at the fixed node.  I.e. 1 for the fixed state, and 0 for the others.
#' @return \code{bears_output} A list of outputs.  bears_output$optim_result
#' @export
#' @seealso \code{\link[cladoRcpp]{numstates_from_numareas}}, \code{\link{getranges_from_LagrangePHYLIP}}, \code{\link[ape]{read.tree}}, \code{\link{calc_loglike_sp}}
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
#' @examples
#' test=1
#' 
#' # Get the example files directory
#' extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
#' # tmp hard code: extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"
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
#' bears_output = bears_2param_standard_fast(trfn=trfn, geogfn=geogfn)
#' bears_output
#' }
#'
bears_2param_standard_fast_fixnode <- function(trfn = "Psychotria_5.2.newick", geogfn = "Psychotria_geog.data", max_range_size=NULL,  num_cores_to_use=NULL, fixnode=NULL, fixlikes=NULL)
	{
	defaults='
	trfn="Psychotria_5.2.newick"
	geogfn="Psychotria_geog.data"
	max_range_size=NULL
	num_cores_to_use=NULL
	'
	
	require(cladoRcpp)
	require(rexpokit)




	# Get geographic ranges at tips
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)

	# Get the list of geographic areas
	areas = getareas_from_tipranges_object(tipranges)
	areas_list = seq(0, length(areas)-1, 1)

	# Change the names to tipranges@df:
	names(tipranges@df) = areas_list


	#######################################################
	# Set the maximum range size (can be thought of as
	# a free parameter
	#######################################################
	if (is.null(max_range_size))
		{
		max_range_size = length(areas)
		need_to_fix_rangesize = TRUE
		} else {
		need_to_fix_rangesize = FALSE
		}

	
	# Old/slow way of getting the list of states and speciation matrix (slow)
	# old_states_list = areas_list_to_states_list(areas, include_null_range=TRUE)
	# old_states_list
	# spmat = make_relprob_matrix_bi(old_states_list)
	# spmat
	
	# max_numareas = max(sapply(X=old_states_list, FUN=length), na.rm=TRUE)
	# max_numareas


	# Take the list of areas, and get list of possible states
	states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
	states_list

	if (need_to_fix_rangesize == FALSE)
		{
		max_numareas = max(sapply(X=states_list, FUN=length), na.rm=TRUE)
		max_numareas

		# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
		if (max_numareas >= 8)
			{
			force_sparse = TRUE
			} else {
			force_sparse = FALSE
			}	

		} else {
		max_numareas = max_range_size
		}
	# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
	if (max_numareas >= 8)
		{
		force_sparse = TRUE
		} else {
		force_sparse = FALSE
		}	
	
	
	# Load the phylogenetic tree
	phy = read.tree(file=trfn)


	# The likelihood of each state at the tips
	tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges, phy, states_list=states_list, maxareas=max_numareas)
	tip_condlikes_of_data_on_each_state

	maxent_constraint_01 = 0.0001


	#######################################################
	# Set up the function for optimization
	#######################################################	
	function_to_optim <- function(params, phy, tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, fixnode=fixnode, fixlikes=fixlikes)
		{
		junk='
		cluster_already_open=NULL
		'
		
		# Lagrange-like (all new species are within the geographic range of the ancestor)

		# Set the dispersal and extinction rate
		d = params[1]
		e = params[2]

		# Equal dispersal in all directions (unconstrained)
		# Equal extinction probability for all areas
		distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

		dmat = matrix(d, nrow=length(areas), ncol=length(areas))
		elist = rep(e, length(areas))
		
		# Set up the instantaneous rate matrix (Q matrix)
		Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)
	
		# Cladogenic model
		ys = 0.5
		j = 0
		v = 0.5
		sum_SPweights = ys + j + v

		# Text version of speciation matrix	
		maxent_constraint_01v = maxent_constraint_01
		#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
			
		# Set the parameter controlling the size distribution of 
		# the smaller descendant species
		maxent01s_param = maxent_constraint_01
		maxent01v_param = maxent_constraint_01
		maxent01j_param = maxent_constraint_01
		maxent01y_param = maxent_constraint_01


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
			# Before calculating the log likelihood, print it, in case there is e.g. a bug
			cat("d=", d, "; e=", e, "; j=", j, "; ys=", ys, "; v=", v, "; sum=", sum_SPweights, "; LnL=", sep="")
			}

		# Calculate the log-likelihood of the data, given the model parameters during this iteration	
		ttl_loglike = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="loglike", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, cluster_already_open=cluster_already_open, fixnode=fixnode, fixlikes=fixlikes)
		ttl_loglike
	
		if (print_optim == TRUE)
			{
			# If the log likelihood is successful, print it
			cat(ttl_loglike, "\n", sep="")
			}

		#return(-1*ttl_loglike)
		return(ttl_loglike)
		}

	
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
	params = c(0.01, 0.01)
	minj = 1e-05
	# start on Lagrange results
	#params = c(3.11882,  2.51741)
	lower = c(1e-15,1e-15)
	upper = c(1,1)
	
	# High performance computing
	# HPC using parallel package in R 2.15 or higher, which allows
	# mcmapply (multicore apply)
	# Don't use multicore if using R.app ('AQUA')
	if (.Platform$GUI != "AQUA" && ((is.null(num_cores_to_use) == TRUE) || (num_cores_to_use > 1)) )
		{
		# We are doing manual, optional processing on several cores;
		# this seems to have less overhead/hassle/incompatibility issues
		# than using mcmapply, mclapply, etc...
		#require("parallel") #<- do this higher up

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
		
	
	
	
	# Run optimization	
	use_optimx = TRUE
	if (use_optimx == FALSE)
		{
		optim_result2 = optim(par=params, fn=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, fixnode=fixnode, fixlikes=fixlikes, method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
	
	#optim_result2 = nlminb(start=params, objective=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, lower=lower, upper=upper, control=list(iter.max=50, trace=1, abs.tol=0.001))# method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
		} else {
		# Compare methods with optimx
		#require(optimx)
		
		
		
		optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, fixnode=fixnode, fixlikes=fixlikes, itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

		# Run with all methods, for testing:
		# optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

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
	

	if (exists("cluster_open") && (exists("cluster_open")))
		{
		cat("\n\nStopping cluster with ", num_cores_to_use, " cores.\n\n", sep="")
		stopCluster(cluster_already_open)
		}
	


	#######################################################
	# Summarize results 
	#######################################################

	# We DON'T need this for marginal optimization over a particular fixed node
# 
# 	if (use_optimx == FALSE)
# 		{
# 		# Set the dispersal and extinction rate
# 		d = optim_result2$par[1]
# 		e = optim_result2$par[2]
# 		} else {
# 		d = optim_result2$par[[1]][1]
# 		e = optim_result2$par[[1]][2]
# 		}
# 
# 	# Equal dispersal in all directions (unconstrained)
# 	# Equal extinction probability for all areas
# 	distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))
# 
# 	dmat = matrix(d, nrow=length(areas), ncol=length(areas))
# 	elist = rep(e, length(areas))
# 	
# 	# Set up the instantaneous rate matrix (Q matrix)
# 	Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)
# 
# 	# Cladogenic model
# 	ys = 0.5
# 	j = 0
# 	v = 0.5
# 	sum_SPweights = ys + j + v
# 
# 	# Text version of speciation matrix	
# 	maxent_constraint_01v = maxent_constraint_01
# 		
# 	# Set the parameter controlling the size distribution of 
# 	# the smaller descendant species
# 	maxent01s_param = maxent_constraint_01
# 	maxent01v_param = maxent_constraint_01
# 	maxent01j_param = maxent_constraint_01
# 	maxent01y_param = maxent_constraint_01
# 
# 	# Cladogenesis model inputs
# 	spPmat_inputs = NULL
# 	states_indices = states_list
# 	states_indices[1] = NULL	# shorten the states_indices by 1 (cutting the null range state from the speciation matrix)
# 	spPmat_inputs$l = states_indices
# 	spPmat_inputs$s = ys
# 	spPmat_inputs$v = v
# 	spPmat_inputs$j = j
# 	spPmat_inputs$y = ys
# 	spPmat_inputs$dmat = distances_mat
# 	spPmat_inputs$maxent01s_param = maxent01s_param
# 	spPmat_inputs$maxent01v_param = maxent01v_param
# 	spPmat_inputs$maxent01j_param = maxent01j_param
# 	spPmat_inputs$maxent01y_param = maxent01y_param
# 
# 
# 	# Calculate the log-likelihood of the data, given the model parameters during this iteration	
# 	model_results = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="all", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, calc_ancprobs=TRUE)
# 
# 
# 
# 	#######################################################
# 	# Store results in a BEARS result
# 	#######################################################
# 	bears_output = model_results
# 	bears_output$optim_result = optim_result2
# 	bears_output$spPmat_inputs = spPmat_inputs
# 	
	return(optim_result2)
	}










#######################################################
# bears_2param_DIVA_fast
#######################################################
#' 2-parameter model, fixed cladogenesis model (as in LAGRANGE)
#' 
#' This function implements a biogeographical model with 2 free parameters (\emph{d}, rate of dispersal/range addition, and
#' \emph{e}, rate of extinction/range contraction), and a fixed cladogenesis model copying the DIVA model (\cite{Ronquist_1997_DIVA}.  This
#' has: equal probability of vicariance at all range sizes, but NO sympatric-subset speciation, no jump/founder-event speciation,
#' and sympatric-range-copying events are limited to the smaller descendant always having a range size of
#' 1 area(\cite{Ronquist_Sanmartin_2011}).
#'
#' Once the model is set up,
#' it is input into the optimization routine \code{\link[optimx]{optimx}} (the more common \code{\link[stats]{optim}} can also be used
#' by editing the function), and \code{\link{calc_loglike_sp}} is used to calculate the log-likelihood of 
#' each set of parameters.  Once the parameter values that give the data the maximum likelihood are found, they are reported back
#' to the function and returned to the user.
#'
#' This duplicates the model used in the standard DIVA implementation (\cite{ReeSmith2008}, \cite{Ree2009configurator}, \cite{SmithRee2010_CPPversion},
#' with no constraints on dispersal or range size.
#' 
#' Here, all of the fastest processing options have been used.
#'
#' Model implementations are provided to show the user how a specific model can be set up and optimized.  This is preferable compared to the "black-box"
#' nature of most other inference packages.  Users are encouraged to experiment.  Useful models can be added to later versions of \code{BioGeoBEARS}.
#'
#' @param trfn The filename of the phylogenetic tree, in NEWICK format (\url{http://evolution.genetics.washington.edu/phylip/newicktree.html}).  
#' Tipnames should match the names in geogfn.  See \code{\link[ape]{read.tree}} in APE for reading in phylogenetic trees.
#' @param geogfn A PHYLIP-style file with geographic range data (see \code{\link{getranges_from_LagrangePHYLIP}}) for each tipname. This is the same format
#' used by C++ LAGRANGE (\cite{SmithRee2010_CPPversion}).
#' @param max_range_size The maximum rangesize, in number of areas.  Having a smaller maximum range size means that you can have more areas (the size of the
#' state space is greatly reduced; see \code{\link[cladoRcpp]{numstates_from_numareas}}.
#' @param num_cores_to_use If >1, parallel processing will be attempted. \bold{Note:} parallel processing via \code{library (parallel)} will work in Mac command-line
#' R, but not in Mac GUI \code{R.app}.
#' @return \code{bears_output} A list of outputs.  bears_output$optim_result
#' @export
#' @seealso \code{\link[cladoRcpp]{numstates_from_numareas}}, \code{\link{getranges_from_LagrangePHYLIP}}, \code{\link[ape]{read.tree}}, \code{\link{calc_loglike_sp}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' Felsenstein, Joe.  The Newick tree format.  \url{http://evolution.genetics.washington.edu/phylip/newicktree.html}
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Ronquist_1997_DIVA
#'	 @cite Ronquist1996_DIVA
#'	 @cite Ronquist_Sanmartin_2011
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
#' bears_output = bears_2param_standard_fast(trfn=trfn, geogfn=geogfn)
#' bears_output
#' }
#'
bears_2param_DIVA_fast <- function(trfn = "Psychotria_5.2.newick", geogfn = "Psychotria_geog.data", max_range_size=NULL,  num_cores_to_use=NULL)
	{
	defaults='
	trfn="Psychotria_5.2.newick"
	geogfn="Psychotria_geog.data"
	max_range_size=NULL
	num_cores_to_use=NULL
	'
	
	require(cladoRcpp)
	require(rexpokit)




	# Get geographic ranges at tips
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)

	# Get the list of geographic areas
	areas = getareas_from_tipranges_object(tipranges)
	areas_list = seq(0, length(areas)-1, 1)

	# Change the names to tipranges@df:
	names(tipranges@df) = areas_list


	#######################################################
	# Set the maximum range size (can be thought of as
	# a free parameter
	#######################################################
	if (is.null(max_range_size))
		{
		max_range_size = length(areas)
		need_to_fix_rangesize = TRUE
		} else {
		need_to_fix_rangesize = FALSE
		}

	
	# Old/slow way of getting the list of states and speciation matrix (slow)
	# old_states_list = areas_list_to_states_list(areas, include_null_range=TRUE)
	# old_states_list
	# spmat = make_relprob_matrix_bi(old_states_list)
	# spmat
	
	# max_numareas = max(sapply(X=old_states_list, FUN=length), na.rm=TRUE)
	# max_numareas


	# Take the list of areas, and get list of possible states
	states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
	states_list

	if (need_to_fix_rangesize == FALSE)
		{
		max_numareas = max(sapply(X=states_list, FUN=length), na.rm=TRUE)
		max_numareas

		# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
		if (max_numareas >= 8)
			{
			force_sparse = TRUE
			} else {
			force_sparse = FALSE
			}	

		} else {
		max_numareas = max_range_size
		}
	# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
	if (max_numareas >= 8)
		{
		force_sparse = TRUE
		} else {
		force_sparse = FALSE
		}	
	
	
	# Load the phylogenetic tree
	phy = read.tree(file=trfn)


	# The likelihood of each state at the tips
	tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges, phy, states_list=states_list, maxareas=max_numareas)
	tip_condlikes_of_data_on_each_state

	maxent_constraint_01 = 0.0001


	#######################################################
	# Set up the function for optimization
	#######################################################	
	function_to_optim <- function(params, phy, tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open)
		{
		junk='
		cluster_already_open=NULL
		'
		
		# Lagrange-like (all new species are within the geographic range of the ancestor)

		# Set the dispersal and extinction rate
		d = params[1]
		e = params[2]

		# Equal dispersal in all directions (unconstrained)
		# Equal extinction probability for all areas
		distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

		dmat = matrix(d, nrow=length(areas), ncol=length(areas))
		elist = rep(e, length(areas))
		
		# Set up the instantaneous rate matrix (Q matrix)
		Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)
	
		# Cladogenic model
		ys = 0.5
		j = 0
		v = 0.5
		sum_SPweights = ys + j + v

		# Text version of speciation matrix	
		maxent_constraint_01v = 0.5
		#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
			
		# Set the parameter controlling the size distribution of 
		# the smaller descendant species
		maxent01s_param = maxent_constraint_01
		maxent01v_param = maxent_constraint_01v
		maxent01j_param = maxent_constraint_01
		maxent01y_param = maxent_constraint_01


		# Cladogenesis model inputs
		spPmat_inputs = NULL
		states_indices = states_list
		states_indices[1] = NULL	# shorten the states_indices by 1 (cutting the null range state from the speciation matrix)
		spPmat_inputs$l = states_indices
		spPmat_inputs$s = 0
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
			# Before calculating the log likelihood, print it, in case there is e.g. a bug
			cat("d=", d, "; e=", e, "; j=", j, "; y=", spPmat_inputs$y, "; s=", spPmat_inputs$s, "; v=", v, "; maxent01=", maxent_constraint_01, "; maxent01v=", maxent_constraint_01v, "; sum=", sum_SPweights, "; LnL=", sep="")
			}

		# Calculate the log-likelihood of the data, given the model parameters during this iteration	
		ttl_loglike = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="loglike", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, cluster_already_open=cluster_already_open)
		ttl_loglike
	
		if (print_optim == TRUE)
			{
			# If the log likelihood is successful, print it
			cat(ttl_loglike, "\n", sep="")
			}

		#return(-1*ttl_loglike)
		return(ttl_loglike)
		}

	
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
	params = c(0.01, 0.01)
	minj = 1e-05
	# start on Lagrange results
	#params = c(3.11882,  2.51741)
	lower = c(1e-15,1e-15)
	upper = c(1,1)
	
	# High performance computing
	# HPC using parallel package in R 2.15 or higher, which allows
	# mcmapply (multicore apply)
	# Don't use multicore if using R.app ('AQUA')
	if (.Platform$GUI != "AQUA" && ((is.null(num_cores_to_use) == TRUE) || (num_cores_to_use > 1)) )
		{
		# We are doing manual, optional processing on several cores;
		# this seems to have less overhead/hassle/incompatibility issues
		# than using mcmapply, mclapply, etc...
		#require("parallel") #<- do this higher up

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
		
	
	
	
	# Run optimization	
	use_optimx = TRUE
	if (use_optimx == FALSE)
		{
		optim_result2 = optim(par=params, fn=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
	
	#optim_result2 = nlminb(start=params, objective=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, lower=lower, upper=upper, control=list(iter.max=50, trace=1, abs.tol=0.001))# method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
		} else {
		# Compare methods with optimx
		#require(optimx)
		
		
		
		optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

		# Run with all methods, for testing:
		# optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

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
	

	if (exists("cluster_open") && (exists("cluster_open")))
		{
		cat("\n\nStopping cluster with ", num_cores_to_use, " cores.\n\n", sep="")
		stopCluster(cluster_already_open)
		}
	
	#######################################################
	# Summarize results 
	#######################################################

	if (use_optimx == FALSE)
		{
		# Set the dispersal and extinction rate
		d = optim_result2$par[1]
		e = optim_result2$par[2]
		} else {
		d = optim_result2$par[[1]][1]
		e = optim_result2$par[[1]][2]
		}

	# Equal dispersal in all directions (unconstrained)
	# Equal extinction probability for all areas
	distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

	dmat = matrix(d, nrow=length(areas), ncol=length(areas))
	elist = rep(e, length(areas))
	
	# Set up the instantaneous rate matrix (Q matrix)
	Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)

	# Cladogenic model
	ys = 0.5
	j = 0
	v = 0.5
	sum_SPweights = ys + j + v

	# Text version of speciation matrix	
	maxent_constraint_01v = maxent_constraint_01
		
	# Set the parameter controlling the size distribution of 
	# the smaller descendant species
	maxent01s_param = maxent_constraint_01
	maxent01v_param = maxent_constraint_01
	maxent01j_param = maxent_constraint_01
	maxent01y_param = maxent_constraint_01

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


	# Calculate the log-likelihood of the data, given the model parameters during this iteration	
	model_results = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="all", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, calc_ancprobs=TRUE)



	#######################################################
	# Store results in a BEARS result
	#######################################################
	bears_output = model_results
	bears_output$optim_result = optim_result2
	bears_output$spPmat_inputs = spPmat_inputs
	
	return(bears_output)
	}










#######################################################
# bears_2param_standard_fast_fortest
#######################################################
#' 2-parameter model, fixed cladogenesis model (as in LAGRANGE) -- older test version
#' 
#' This is an older, test version of \code{\link{bears_2param_standard_fast}}.
#'
#' @param trfn The filename of the phylogenetic tree, in NEWICK format (\url{http://evolution.genetics.washington.edu/phylip/newicktree.html}).  
#' Tipnames should match the names in geogfn.  See \code{\link[ape]{read.tree}} in APE for reading in phylogenetic trees.
#' @param geogfn A PHYLIP-style file with geographic range data (see \code{\link{getranges_from_LagrangePHYLIP}}) for each tipname. This is the same format
#' used by C++ LAGRANGE (\cite{SmithRee2010_CPPversion}).
#' @return \code{bears_output} A list of outputs.  bears_output$optim_result
#' @export
#' @seealso \code{\link[cladoRcpp]{numstates_from_numareas}}, \code{\link{getranges_from_LagrangePHYLIP}}, \code{\link[ape]{read.tree}}, \code{\link{calc_loglike_sp}}
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
#' bears_output = bears_2param_standard_fast_fortest(trfn=trfn, geogfn=geogfn)
#' bears_output
#' }
#'
bears_2param_standard_fast_fortest <- function(trfn = "test.newick", geogfn = "test.data")
	{
	defaults='
	trfn="Psychotria_5.2.newick"
	geogfn="Psychotria_geog.data"
	'
	
	require(cladoRcpp)
	require(rexpokit)

	# Get geographic ranges at tips
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)

	# Get the list of geographic areas
	areas = getareas_from_tipranges_object(tipranges)
	areas_list = seq(0, length(areas)-1, 1)

	# Change the names to tipranges@df:
	names(tipranges@df) = areas_list


	#######################################################
	# Set the maximum range size (can be thought of as
	# a free parameter
	#######################################################
	if (is.null(max_range_size))
		{
		max_range_size = length(areas)
		}
	
	# Old/slow way of getting the list of states and speciation matrix (slow)
	# old_states_list = areas_list_to_states_list(areas, include_null_range=TRUE)
	# old_states_list
	# spmat = make_relprob_matrix_bi(old_states_list)
	# spmat
	
	# max_numareas = max(sapply(X=old_states_list, FUN=length), na.rm=TRUE)
	# max_numareas

	#######################################################
	# Set the maximum range size (can be thought of as
	# a free parameter
	#######################################################
	if (is.null(max_range_size))
		{
		max_range_size = length(areas)
		need_to_fix_rangesize = TRUE
		} else {
		need_to_fix_rangesize = FALSE
		}


	# Take the list of areas, and get list of possible states
	states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
	states_list

	if (need_to_fix_rangesize == FALSE)
		{
		max_numareas = max(sapply(X=states_list, FUN=length), na.rm=TRUE)
		max_numareas

		# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
		if (max_numareas >= 8)
			{
			force_sparse = TRUE
			} else {
			force_sparse = FALSE
			}	

		} else {
		max_numareas = max_range_size
		}
	# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
	if (max_numareas >= 8)
		{
		force_sparse = TRUE
		} else {
		force_sparse = FALSE
		}	
	
	
	# Load the phylogenetic tree
	phy = read.tree(file=trfn)


	# The likelihood of each state at the tips
	tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges, phy, states_list=states_list, maxareas=max_numareas)
	tip_condlikes_of_data_on_each_state

	maxent_constraint_01 = 0.0001


	#######################################################
	# Set up the function for optimization
	#######################################################	
	function_to_optim <- function(params, phy, tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open)
		{
		# Lagrange-like (all new species are within the geographic range of the ancestor)

		# Set the dispersal and extinction rate
		d = params[1]
		e = params[2]

		# Equal dispersal in all directions (unconstrained)
		# Equal extinction probability for all areas
		distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

		dmat = matrix(d, nrow=length(areas), ncol=length(areas))
		elist = rep(e, length(areas))
		
		# Set up the instantaneous rate matrix (Q matrix)
		Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)
	
		# Cladogenic model
		ys = 0.5
		j = 0
		v = 0.5
		sum_SPweights = ys + j + v

		# Text version of speciation matrix	
		maxent_constraint_01v = maxent_constraint_01
		#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
			
		# Set the parameter controlling the size distribution of 
		# the smaller descendant species
		maxent01s_param = maxent_constraint_01
		maxent01v_param = maxent_constraint_01
		maxent01j_param = maxent_constraint_01
		maxent01y_param = maxent_constraint_01


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
			# Before calculating the log likelihood, print it, in case there is e.g. a bug
			cat("d=", d, "; e=", e, "; j=", j, "; ys=", ys, "; v=", v, "; sum=", sum_SPweights, "; LnL=", sep="")
			}

		# Calculate the log-likelihood of the data, given the model parameters during this iteration	
		ttl_loglike = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="loglike", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, cluster_already_open=cluster_already_open)
		ttl_loglike
	
		if (print_optim == TRUE)
			{
			# If the log likelihood is successful, print it
			cat(ttl_loglike, "\n", sep="")
			}
	
		return(ttl_loglike)
		}

	
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
	params = c(0.01, 0.01)
	minj = 1e-05
	# start on Lagrange results
	#params = c(3.11882,  2.51741)
	lower = c(1e-15,1e-15)
	upper = c(10,10)
	
	# High performance computing
	# HPC using parallel package in R 2.15 or higher, which allows
	# mcmapply (multicore apply)
	# Don't use multicore if using R.app ('AQUA')
	if (.Platform$GUI != "AQUA" && ((is.null(num_cores_to_use) == TRUE) || (num_cores_to_use > 1)) )
		{
		# We are doing manual, optional processing on several cores;
		# this seems to have less overhead/hassle/incompatibility issues
		# than using mcmapply, mclapply, etc...
		#require("parallel") #<- do this higher up

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

	
	
	
	# Run optimization	
	use_optimx = TRUE
	if (use_optimx == FALSE)
		{
		optim_result2 = optim(par=params, fn=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
	
	#optim_result2 = nlminb(start=params, objective=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, lower=lower, upper=upper, control=list(iter.max=50, trace=1, abs.tol=0.001))# method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
		} else {
		# Compare methods with optimx
		#require(optimx)
		
		
		
		optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

		# Run with all methods, for testing:
		# optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

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
	

	if (exists("cluster_open") && (cluster_open == TRUE))
		{
		cat("\n\nStopping cluster with ", num_cores_to_use, " cores.\n\n", sep="")
		stopCluster(cluster_already_open)
		}
	
	
	
	
	
	#######################################################
	# Summarize results 
	#######################################################

	# Set the dispersal and extinction rate
	if (use_optimx == FALSE)
		{
		# Set the dispersal and extinction rate
		d = optim_result2$par[1]
		e = optim_result2$par[2]
		} else {
		d = optim_result2$par[[1]][1]
		e = optim_result2$par[[1]][2]
		}

	# Equal dispersal in all directions (unconstrained)
	# Equal extinction probability for all areas
	distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

	dmat = matrix(d, nrow=length(areas), ncol=length(areas))
	elist = rep(e, length(areas))
	
	# Set up the instantaneous rate matrix (Q matrix)
	Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)

	# Cladogenic model
	ys = 0.5
	j = 0
	v = 0.5
	sum_SPweights = ys + j + v

	# Text version of speciation matrix	
	maxent_constraint_01v = maxent_constraint_01
		
	# Set the parameter controlling the size distribution of 
	# the smaller descendant species
	maxent01s_param = maxent_constraint_01
	maxent01v_param = maxent_constraint_01
	maxent01j_param = maxent_constraint_01
	maxent01y_param = maxent_constraint_01

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


	# Calculate the log-likelihood of the data, given the model parameters during this iteration	
	model_results = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="all", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, calc_ancprobs=TRUE)



	#######################################################
	# Store results in a BEARS result
	#######################################################
	bears_output = model_results
	bears_output$optim_result = optim_result2
	bears_output$spPmat_inputs = spPmat_inputs
	
	return(bears_output)
	}

















# Allow only range-copying (complete sympatry) during cladogenesis

#######################################################
# bears_2param_standard_fast_symOnly
#######################################################
#' 2-parameter model, no cladogenesis model (as in BayArea or other purely continuous-time model)
#' 
#' This implements a 2-parameter model, as in LAGRANGE or \code{\link{bears_2param_standard_fast}}, but omits the speciation/cladogenesis
#' model.  This means that the model is purely continuous-time, as when biogeographic range is treated as a discrete character in 
#' software designed for inference on morphological () or molecular data ().  This model is that implemented in \code{BayArea}, if no 
#' distance-dependent effect on dispersal probability is assumed.  Such distance-dependence could easily be added with a third parameter, 
#' however.
#' 
#' \code{BayArea} is a new program by Landis, Matzke, Moore, and Huelsenbeck; see \cite{Landis_Matzke_etal_2013_BayArea}. 
#' However, BayArea does not currently implement cladogenesis models; it only has continuous-time model for evolutionary change along branches.  In effect,
#' this means that the cladogenesis model is sympatric speciation with complete range copying with probability 1.
#' 
#' @param trfn The filename of the phylogenetic tree, in NEWICK format (\url{http://evolution.genetics.washington.edu/phylip/newicktree.html}).  
#' Tipnames should match the names in geogfn.  See \code{\link[ape]{read.tree}} in APE for reading in phylogenetic trees.
#' @param geogfn A PHYLIP-style file with geographic range data (see \code{\link{getranges_from_LagrangePHYLIP}}) for each tipname. This is the same format
#' used by C++ LAGRANGE (\cite{SmithRee2010_CPPversion}).
#' @param max_range_size The maximum rangesize, in number of areas.  Having a smaller maximum range size means that you can have more areas (the size of the
#' state space is greatly reduced; see \code{\link[cladoRcpp]{numstates_from_numareas}}.
#' @param num_cores_to_use If >1, parallel processing will be attempted. \bold{Note:} parallel processing via \code{library (parallel)} will work in Mac command-line
#' R, but not in Mac GUI \code{R.app}.
#' @return \code{bears_output} A list of outputs.  bears_output$optim_result
#' @export
#' @seealso \code{\link[cladoRcpp]{numstates_from_numareas}}, \code{\link{getranges_from_LagrangePHYLIP}}, \code{\link[ape]{read.tree}}, \code{\link{calc_loglike_sp}}
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
#' bears_output = bears_2param_standard_fast_symOnly(trfn=trfn, geogfn=geogfn)
#' bears_output
#' }
#'
bears_2param_standard_fast_symOnly <- function(trfn = "Psychotria_5.2.newick", geogfn = "Psychotria_geog.data", max_range_size=NULL, num_cores_to_use=NULL)
	{
	defaults='
	trfn = "Psychotria_5.2.newick"
	geogfn = "Psychotria_geog.data"
	'
	
	require(cladoRcpp)
	require(rexpokit)

	# Get geographic ranges at tips
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)


	# Get the list of geographic areas
	areas = getareas_from_tipranges_object(tipranges)
	areas_list = seq(0, length(areas)-1, 1)

	# Change the names to tipranges@df:
	names(tipranges@df) = areas_list

	
	# Old/slow way of getting the list of states and speciation matrix (slow)
	# old_states_list = areas_list_to_states_list(areas, include_null_range=TRUE)
	# old_states_list
	# spmat = make_relprob_matrix_bi(old_states_list)
	# spmat
	
	# max_numareas = max(sapply(X=old_states_list, FUN=length), na.rm=TRUE)
	# max_numareas

	#######################################################
	# Set the maximum range size (can be thought of as
	# a free parameter
	#######################################################
	if (is.null(max_range_size))
		{
		max_range_size = length(areas)
		need_to_fix_rangesize = TRUE
		} else {
		need_to_fix_rangesize = FALSE
		}


	# Take the list of areas, and get list of possible states
	states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
	states_list

	if (need_to_fix_rangesize == FALSE)
		{
		max_numareas = max(sapply(X=states_list, FUN=length), na.rm=TRUE)
		max_numareas

		# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
		if (max_numareas >= 8)
			{
			force_sparse = TRUE
			} else {
			force_sparse = FALSE
			}	

		} else {
		max_numareas = max_range_size
		}
	# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
	if (max_numareas >= 8)
		{
		force_sparse = TRUE
		} else {
		force_sparse = FALSE
		}	
	
	
	# Load the phylogenetic tree
	phy = read.tree(file=trfn)


	# The likelihood of each state at the tips
	tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges, phy, states_list=states_list, maxareas=max_numareas)
	tip_condlikes_of_data_on_each_state

	maxent_constraint_01 = 0.9999


	#######################################################
	# Set up the function for optimization
	#######################################################	
	function_to_optim <- function(params, phy, tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open)
		{
		# Lagrange-like (all new species are within the geographic range of the ancestor)

		# Set the dispersal and extinction rate
		d = params[1]
		e = params[2]

		# Equal dispersal in all directions (unconstrained)
		# Equal extinction probability for all areas
		distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

		dmat = matrix(d, nrow=length(areas), ncol=length(areas))
		elist = rep(e, length(areas))
		
		# Set up the instantaneous rate matrix (Q matrix)
		Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)
	
		# Cladogenic model
		ys = 1
		j = 0
		v = 0
		sum_SPweights = ys + j + v

		# Text version of speciation matrix	
		maxent_constraint_01v = maxent_constraint_01
		#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
			
		# Set the parameter controlling the size distribution of 
		# the smaller descendant species
		maxent01s_param = 0.0001
		maxent01v_param = 0.0001
		maxent01j_param = 0.0001
		maxent01y_param = maxent_constraint_01


		# Cladogenesis model inputs
		spPmat_inputs = NULL
		states_indices = states_list
		states_indices[1] = NULL	# shorten the states_indices by 1 (cutting the null range state from the speciation matrix)
		spPmat_inputs$l = states_indices
		spPmat_inputs$s = 0.0001
		spPmat_inputs$v = v
		spPmat_inputs$j = j
		spPmat_inputs$y = ys - spPmat_inputs$s
		spPmat_inputs$dmat = distances_mat
		spPmat_inputs$maxent01s_param = maxent01s_param
		spPmat_inputs$maxent01v_param = maxent01v_param
		spPmat_inputs$maxent01j_param = maxent01j_param
		spPmat_inputs$maxent01y_param = maxent01y_param


		if (print_optim == TRUE)
			{
			# Before calculating the log likelihood, print it, in case there is e.g. a bug
			cat("d=", d, "; e=", e, "; j=", j, "; ys=", ys, "; v=", v, "; sum=", sum_SPweights, "; LnL=", sep="")
			}

		# Calculate the log-likelihood of the data, given the model parameters during this iteration	
		ttl_loglike = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="loglike", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, cluster_already_open=cluster_already_open)
		ttl_loglike
	
		if (print_optim == TRUE)
			{
			# If the log likelihood is successful, print it
			cat(ttl_loglike, "\n", sep="")
			}

		return(ttl_loglike)
		}

	
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
	params = c(0.01, 0.01)
	minj = 1e-05
	# start on Lagrange results
	#params = c(3.11882,  2.51741)
	lower = c(1e-15,1e-15)
	upper = c(1,1)
	
	# High performance computing
	# HPC using parallel package in R 2.15 or higher, which allows
	# mcmapply (multicore apply)
	# Don't use multicore if using R.app ('AQUA')
	if (.Platform$GUI != "AQUA" && ((is.null(num_cores_to_use) == TRUE) || (num_cores_to_use > 1)) )
		{
		# We are doing manual, optional processing on several cores;
		# this seems to have less overhead/hassle/incompatibility issues
		# than using mcmapply, mclapply, etc...
		#require("parallel") #<- do this higher up

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

	
	
	
	# Run optimization	
	use_optimx = TRUE
	if (use_optimx == FALSE)
		{
		optim_result2 = optim(par=params, fn=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
	
	#optim_result2 = nlminb(start=params, objective=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, lower=lower, upper=upper, control=list(iter.max=50, trace=1, abs.tol=0.001))# method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
		} else {
		# Compare methods with optimx
		#require(optimx)
		
		
		
		optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

		# Run with all methods, for testing:
		# optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

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
	

	if (exists("cluster_open") && (cluster_open == TRUE))
		{
		cat("\n\nStopping cluster with ", num_cores_to_use, " cores.\n\n", sep="")
		stopCluster(cluster_already_open)
		}
	
	
	
	
	#######################################################
	# Summarize results 
	#######################################################

	# Set the dispersal and extinction rate
	if (use_optimx == FALSE)
		{
		# Set the dispersal and extinction rate
		d = optim_result2$par[1]
		e = optim_result2$par[2]
		} else {
		d = optim_result2$par[[1]][1]
		e = optim_result2$par[[1]][2]
		}

	# Equal dispersal in all directions (unconstrained)
	# Equal extinction probability for all areas
	distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

	dmat = matrix(d, nrow=length(areas), ncol=length(areas))
	elist = rep(e, length(areas))
	
	# Set up the instantaneous rate matrix (Q matrix)
	Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)

	# Cladogenic model
	ys = 1
	j = 0
	v = 0
	sum_SPweights = ys + j + v

	# Text version of speciation matrix	
	maxent_constraint_01v = maxent_constraint_01
	#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
		
	# Set the parameter controlling the size distribution of 
	# the smaller descendant species
	maxent01s_param = 0.0001
	maxent01v_param = 0.0001
	maxent01j_param = 0.0001
	maxent01y_param = maxent_constraint_01


	# Cladogenesis model inputs
	spPmat_inputs = NULL
	states_indices = states_list
	states_indices[1] = NULL	# shorten the states_indices by 1 (cutting the null range state from the speciation matrix)
	spPmat_inputs$l = states_indices
	spPmat_inputs$s = 0.0001
	spPmat_inputs$v = v
	spPmat_inputs$j = j
	spPmat_inputs$y = ys - spPmat_inputs$s
	spPmat_inputs$dmat = distances_mat
	spPmat_inputs$maxent01s_param = maxent01s_param
	spPmat_inputs$maxent01v_param = maxent01v_param
	spPmat_inputs$maxent01j_param = maxent01j_param
	spPmat_inputs$maxent01y_param = maxent01y_param


	# Calculate the log-likelihood of the data, given the model parameters during this iteration	
	model_results = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="all", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, calc_ancprobs=TRUE)




	#######################################################
	# Store results in a BEARS result
	#######################################################
	bears_output = model_results
	bears_output$optim_result = optim_result2
	bears_output$spPmat_inputs = spPmat_inputs
	
	return(bears_output)
	}








#######################################################
# bears_2param_standard_fast_symOnly_simp
#######################################################
#' 2-parameter model, no cladogenesis model (as in BayArea or other purely continuous-time model)
#' 
#' (Forcing no speciation model.) This implements a 2-parameter model, as in LAGRANGE or \code{\link{bears_2param_standard_fast}}, 
#' but omits the speciation/cladogenesis
#' model.  This means that the model is purely continuous-time, as when biogeographic range is treated as a discrete character in 
#' software designed for inference on morphological () or molecular data ().  This model is that implemented in \code{BayArea}, if no 
#' distance-dependent effect on dispersal probability is assumed.  Such distance-dependence could easily be added with a third parameter, 
#' however.
#' 
#' \code{BayArea} is a new program by Landis, Matzke, Moore, and Huelsenbeck; see \cite{Landis_Matzke_etal_2013_BayArea}. 
#' However, BayArea does not currently implement cladogenesis models; it only has continuous-time model for evolutionary change along branches.  In effect,
#' this means that the cladogenesis model is sympatric speciation with complete range copying with probability 1.
#' 
#' @param trfn The filename of the phylogenetic tree, in NEWICK format (\url{http://evolution.genetics.washington.edu/phylip/newicktree.html}).  
#' Tipnames should match the names in geogfn.  See \code{\link[ape]{read.tree}} in APE for reading in phylogenetic trees.
#' @param geogfn A PHYLIP-style file with geographic range data (see \code{\link{getranges_from_LagrangePHYLIP}}) for each tipname. This is the same format
#' used by C++ LAGRANGE (\cite{SmithRee2010_CPPversion}).
#' @param max_range_size The maximum rangesize, in number of areas.  Having a smaller maximum range size means that you can have more areas (the size of the
#' state space is greatly reduced; see \code{\link[cladoRcpp]{numstates_from_numareas}}.
#' @param num_cores_to_use If >1, parallel processing will be attempted. \bold{Note:} parallel processing via \code{library (parallel)} will work in Mac command-line
#' R, but not in Mac GUI \code{R.app}.
#' @return \code{bears_output} A list of outputs.  bears_output$optim_result
#' @export
#' @seealso \code{\link[cladoRcpp]{numstates_from_numareas}}, \code{\link{getranges_from_LagrangePHYLIP}}, \code{\link[ape]{read.tree}}, \code{\link{calc_loglike_sp}}
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
#' bears_output = bears_2param_standard_fast_symOnly(trfn=trfn, geogfn=geogfn)
#' bears_output
#' }
#'
bears_2param_standard_fast_symOnly_simp <- function(trfn = "Psychotria_5.2.newick", geogfn = "Psychotria_geog.data", max_range_size=NULL, num_cores_to_use=NULL)
	{
	defaults='
	trfn = "Psychotria_5.2.newick"
	geogfn = "Psychotria_geog.data"
	'
	
	require(cladoRcpp)
	require(rexpokit)

	# Get geographic ranges at tips
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)


	# Get the list of geographic areas
	areas = getareas_from_tipranges_object(tipranges)
	areas_list = seq(0, length(areas)-1, 1)

	# Change the names to tipranges@df:
	names(tipranges@df) = areas_list

	
	# Old/slow way of getting the list of states and speciation matrix (slow)
	# old_states_list = areas_list_to_states_list(areas, include_null_range=TRUE)
	# old_states_list
	# spmat = make_relprob_matrix_bi(old_states_list)
	# spmat
	
	# max_numareas = max(sapply(X=old_states_list, FUN=length), na.rm=TRUE)
	# max_numareas

	#######################################################
	# Set the maximum range size (can be thought of as
	# a free parameter
	#######################################################
	if (is.null(max_range_size))
		{
		max_range_size = length(areas)
		need_to_fix_rangesize = TRUE
		} else {
		need_to_fix_rangesize = FALSE
		}


	# Take the list of areas, and get list of possible states
	states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
	states_list

	if (need_to_fix_rangesize == FALSE)
		{
		max_numareas = max(sapply(X=states_list, FUN=length), na.rm=TRUE)
		max_numareas

		# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
		if (max_numareas >= 8)
			{
			force_sparse = TRUE
			} else {
			force_sparse = FALSE
			}	

		} else {
		max_numareas = max_range_size
		}
	# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
	if (max_numareas >= 8)
		{
		force_sparse = TRUE
		} else {
		force_sparse = FALSE
		}	
	
	
	# Load the phylogenetic tree
	phy = read.tree(file=trfn)


	# The likelihood of each state at the tips
	tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges, phy, states_list=states_list, maxareas=max_numareas)
	tip_condlikes_of_data_on_each_state

	maxent_constraint_01 = 0.9999


	#######################################################
	# Set up the function for optimization
	#######################################################	
	function_to_optim <- function(params, phy, tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open)
		{
		# Lagrange-like (all new species are within the geographic range of the ancestor)

		# Set the dispersal and extinction rate
		d = params[1]
		e = params[2]

		# Equal dispersal in all directions (unconstrained)
		# Equal extinction probability for all areas
		distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

		dmat = matrix(d, nrow=length(areas), ncol=length(areas))
		elist = rep(e, length(areas))
		
		# Set up the instantaneous rate matrix (Q matrix)
		Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)
	
		# Cladogenic model
		ys = 1
		j = 0
		v = 0
		sum_SPweights = ys + j + v

		# Text version of speciation matrix	
		maxent_constraint_01v = maxent_constraint_01
		#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
			
		# Set the parameter controlling the size distribution of 
		# the smaller descendant species
		maxent01s_param = 0.0001
		maxent01v_param = 0.0001
		maxent01j_param = 0.0001
		maxent01y_param = maxent_constraint_01


		# Cladogenesis model inputs
		spPmat_inputs = NULL
		states_indices = states_list
		states_indices[1] = NULL	# shorten the states_indices by 1 (cutting the null range state from the speciation matrix)
		spPmat_inputs$l = states_indices
		spPmat_inputs$s = 0.0001
		spPmat_inputs$v = v
		spPmat_inputs$j = j
		spPmat_inputs$y = ys - spPmat_inputs$s
		spPmat_inputs$dmat = distances_mat
		spPmat_inputs$maxent01s_param = maxent01s_param
		spPmat_inputs$maxent01v_param = maxent01v_param
		spPmat_inputs$maxent01j_param = maxent01j_param
		spPmat_inputs$maxent01y_param = maxent01y_param


		if (print_optim == TRUE)
			{
			# Before calculating the log likelihood, print it, in case there is e.g. a bug
			cat("d=", d, "; e=", e, "; j=", j, "; ys=", ys, "; v=", v, "; sum=", sum_SPweights, "; LnL=", sep="")
			}

		# Calculate the log-likelihood of the data, given the model parameters during this iteration	
		ttl_loglike = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="loglike", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=NULL, printlevel=0, cluster_already_open=cluster_already_open)
		ttl_loglike
	
		if (print_optim == TRUE)
			{
			# If the log likelihood is successful, print it
			cat(ttl_loglike, "\n", sep="")
			}

		return(ttl_loglike)
		}

	
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
	params = c(0.01, 0.01)
	minj = 1e-05
	# start on Lagrange results
	#params = c(3.11882,  2.51741)
	lower = c(1e-15,1e-15)
	upper = c(1,1)
	
	# High performance computing
	# HPC using parallel package in R 2.15 or higher, which allows
	# mcmapply (multicore apply)
	# Don't use multicore if using R.app ('AQUA')
	if (.Platform$GUI != "AQUA" && ((is.null(num_cores_to_use) == TRUE) || (num_cores_to_use > 1)) )
		{
		# We are doing manual, optional processing on several cores;
		# this seems to have less overhead/hassle/incompatibility issues
		# than using mcmapply, mclapply, etc...
		#require("parallel") #<- do this higher up

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

	
	
	
	# Run optimization	
	use_optimx = TRUE
	if (use_optimx == FALSE)
		{
		optim_result2 = optim(par=params, fn=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
	
	#optim_result2 = nlminb(start=params, objective=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, lower=lower, upper=upper, control=list(iter.max=50, trace=1, abs.tol=0.001))# method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
		} else {
		# Compare methods with optimx
		#require(optimx)
		
		
		
		optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

		# Run with all methods, for testing:
		# optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

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
	

	if (exists("cluster_open") && (cluster_open == TRUE))
		{
		cat("\n\nStopping cluster with ", num_cores_to_use, " cores.\n\n", sep="")
		stopCluster(cluster_already_open)
		}
	
	
	
	
	#######################################################
	# Summarize results 
	#######################################################

	# Set the dispersal and extinction rate
	if (use_optimx == FALSE)
		{
		# Set the dispersal and extinction rate
		d = optim_result2$par[1]
		e = optim_result2$par[2]
		} else {
		d = optim_result2$par[[1]][1]
		e = optim_result2$par[[1]][2]
		}

	# Equal dispersal in all directions (unconstrained)
	# Equal extinction probability for all areas
	distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

	dmat = matrix(d, nrow=length(areas), ncol=length(areas))
	elist = rep(e, length(areas))
	
	# Set up the instantaneous rate matrix (Q matrix)
	Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)

	# Cladogenic model
	ys = 1
	j = 0
	v = 0
	sum_SPweights = ys + j + v

	# Text version of speciation matrix	
	maxent_constraint_01v = maxent_constraint_01
	#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
		
	# Set the parameter controlling the size distribution of 
	# the smaller descendant species
	maxent01s_param = 0.0001
	maxent01v_param = 0.0001
	maxent01j_param = 0.0001
	maxent01y_param = maxent_constraint_01


	# Cladogenesis model inputs
	spPmat_inputs = NULL
	states_indices = states_list
	states_indices[1] = NULL	# shorten the states_indices by 1 (cutting the null range state from the speciation matrix)
	spPmat_inputs$l = states_indices
	spPmat_inputs$s = 0.0001
	spPmat_inputs$v = v
	spPmat_inputs$j = j
	spPmat_inputs$y = ys - spPmat_inputs$s
	spPmat_inputs$dmat = distances_mat
	spPmat_inputs$maxent01s_param = maxent01s_param
	spPmat_inputs$maxent01v_param = maxent01v_param
	spPmat_inputs$maxent01j_param = maxent01j_param
	spPmat_inputs$maxent01y_param = maxent01y_param


	# Calculate the log-likelihood of the data, given the model parameters during this iteration	
	model_results = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="all", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, calc_ancprobs=FALSE)




	#######################################################
	# Store results in a BEARS result
	#######################################################
	bears_output = model_results
	bears_output$optim_result = optim_result2
	bears_output$spPmat_inputs = spPmat_inputs
	
	return(bears_output)
	}











# 3-parameter model, adding j (founder-event speciation)
#######################################################
# bears_3param_standard_fast_fixnode
#######################################################
#' 3-parameter model, adding j (founder-event speciation)
#' 
#' This implements a 3-parameter model, basically \code{LAGRANGE} or \code{\link{bears_2param_standard_fast}}, but with a parameter \emph{j}
#' controlling the relative weight of "founder-event speciation" (\cite{Matzke_2012_IBS}) versus vicariance+sympatric speciation (which are mandated in
#' \code{LAGRANGE} and \code{\link{bears_2param_standard_fast}}.
#' 
#' @param trfn The filename of the phylogenetic tree, in NEWICK format (\url{http://evolution.genetics.washington.edu/phylip/newicktree.html}).  
#' Tipnames should match the names in geogfn.  See \code{\link[ape]{read.tree}} in APE for reading in phylogenetic trees.
#' @param geogfn A PHYLIP-style file with geographic range data (see \code{\link{getranges_from_LagrangePHYLIP}}) for each tipname. This is the same format
#' used by C++ LAGRANGE (\cite{SmithRee2010_CPPversion}).
#' @param max_range_size The maximum rangesize, in number of areas.  Having a smaller maximum range size means that you can have more areas (the size of the
#' state space is greatly reduced; see \code{\link[cladoRcpp]{numstates_from_numareas}}.
#' @param num_cores_to_use If >1, parallel processing will be attempted. \bold{Note:} parallel processing via \code{library (parallel)} will work in Mac command-line
#' R, but not in Mac GUI \code{R.app}.
#' @param fixnode If the state at a particular node is going to be fixed (e.g. for ML marginal ancestral states), give the node number.
#' @param fixlikes The state likelihoods to be used at the fixed node.  I.e. 1 for the fixed state, and 0 for the others.
#' @return \code{bears_output} A list of outputs.  bears_output$optim_result
#' @export
#' @seealso \code{\link[cladoRcpp]{numstates_from_numareas}}, \code{\link{getranges_from_LagrangePHYLIP}}, \code{\link[ape]{read.tree}}, \code{\link{calc_loglike_sp}}
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
#' # tmp hard code: extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"
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
#' bears_output = bears_3param_standard_fast(trfn=trfn, geogfn=geogfn)
#' bears_output
#' }
#'
bears_3param_standard_fast_fixnode <- function(trfn = "Psychotria_5.2.newick", geogfn = "Psychotria_geog.data", max_range_size=NULL, num_cores_to_use=NULL, fixnode=fixnode, fixlikes=fixlikes)
	{
	defaults='
	trfn = "Psychotria_5.2.newick"
	geogfn = "Psychotria_geog.data"
	'
	
	require(cladoRcpp)
	require(rexpokit)

	# Get geographic ranges at tips
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)


	# Get the list of geographic areas
	areas = getareas_from_tipranges_object(tipranges)
	areas_list = seq(0, length(areas)-1, 1)

	# Change the names to tipranges@df:
	names(tipranges@df) = areas_list

	#######################################################
	# Set the maximum range size (can be thought of as
	# a free parameter
	#######################################################
	if (is.null(max_range_size))
		{
		max_range_size = length(areas)
		need_to_fix_rangesize = TRUE
		} else {
		need_to_fix_rangesize = FALSE
		}
	
	# Old/slow way of getting the list of states and speciation matrix (slow)
	# old_states_list = areas_list_to_states_list(areas, include_null_range=TRUE)
	# old_states_list
	# spmat = make_relprob_matrix_bi(old_states_list)
	# spmat
	
	# max_numareas = max(sapply(X=old_states_list, FUN=length), na.rm=TRUE)
	# max_numareas


	# Take the list of areas, and get list of possible states
	states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
	states_list

	if (need_to_fix_rangesize == FALSE)
		{
		max_numareas = max(sapply(X=states_list, FUN=length), na.rm=TRUE)
		max_numareas

		# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
		if (max_numareas >= 8)
			{
			force_sparse = TRUE
			} else {
			force_sparse = FALSE
			}	

		} else {
		max_numareas = max_range_size
		}
	# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
	if (max_numareas >= 8)
		{
		force_sparse = TRUE
		} else {
		force_sparse = FALSE
		}	
	
	
	# Load the phylogenetic tree
	phy = read.tree(file=trfn)


	# The likelihood of each state at the tips
	tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges, phy, states_list=states_list, maxareas=max_numareas)
	tip_condlikes_of_data_on_each_state

	maxent_constraint_01 = 0.0001


	#######################################################
	# Set up the function for optimization
	#######################################################	
	function_to_optim <- function(params, phy, tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, fixnode=fixnode, fixlikes=fixlikes)
		{
		# Lagrange-like (all new species are within the geographic range of the ancestor)

		# Set the dispersal and extinction rate
		d = params[1]
		e = params[2]

		# Equal dispersal in all directions (unconstrained)
		# Equal extinction probability for all areas
		distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

		dmat = matrix(d, nrow=length(areas), ncol=length(areas))
		elist = rep(e, length(areas))
		
		# Set up the instantaneous rate matrix (Q matrix)
		Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)
	
		# Cladogenic model
		j = params[3]
		ysv = 1-j
		ys = 0.5 * ysv
		v = 0.5 * ysv
		sum_SPweights = ys + j + v

		# Text version of speciation matrix	
		maxent_constraint_01v = maxent_constraint_01
		#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
			
		# Set the parameter controlling the size distribution of 
		# the smaller descendant species
		maxent01s_param = maxent_constraint_01
		maxent01v_param = maxent_constraint_01
		maxent01j_param = maxent_constraint_01
		maxent01y_param = maxent_constraint_01


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
			# Before calculating the log likelihood, print it, in case there is e.g. a bug
			cat("d=", d, "; e=", e, "; j=", j, "; ys=", ys, "; v=", v, "; sum=", sum_SPweights, "; LnL=", sep="")
			}

		# Calculate the log-likelihood of the data, given the model parameters during this iteration	
		ttl_loglike = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="loglike", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, cluster_already_open=cluster_already_open, fixnode=fixnode, fixlikes=fixlikes)
		ttl_loglike
	
		if (print_optim == TRUE)
			{
			# If the log likelihood is successful, print it
			cat(ttl_loglike, "\n", sep="")
			}
	
		return(ttl_loglike)
		}

	
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
	params = c(0.01, 0.01, 0.1)
	minj = 1e-04
	# start on Lagrange results
	#params = c(3.11882,  2.51741)
	lower = c(1e-15,1e-15, minj)
	upper = c(1,1, 1-minj)
	
	# High performance computing
	# HPC using parallel package in R 2.15 or higher, which allows
	# mcmapply (multicore apply)
	# Don't use multicore if using R.app ('AQUA')
	if (.Platform$GUI != "AQUA" && ((is.null(num_cores_to_use) == TRUE) || (num_cores_to_use > 1)) )
		{
		# We are doing manual, optional processing on several cores;
		# this seems to have less overhead/hassle/incompatibility issues
		# than using mcmapply, mclapply, etc...
		#require("parallel") #<- do this higher up

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

	
	
	
	# Run optimization	
	use_optimx = TRUE
	if (use_optimx == FALSE)
		{
		optim_result2 = optim(par=params, fn=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, fixnode=fixnode, fixlikes=fixlikes, method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
	
	#optim_result2 = nlminb(start=params, objective=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, lower=lower, upper=upper, control=list(iter.max=50, trace=1, abs.tol=0.001))# method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
		} else {
		# Compare methods with optimx
		#require(optimx)
		
		
		
		optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, fixnode=fixnode, fixlikes=fixlikes, itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

		# Run with all methods, for testing:
		# optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

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
	

	if (exists("cluster_open") && (cluster_open == TRUE))
		{
		cat("\n\nStopping cluster with ", num_cores_to_use, " cores.\n\n", sep="")
		stopCluster(cluster_already_open)
		}
	
	
	
	
	
	#######################################################
	# Summarize results 
	#######################################################
# 
# 	# Set the dispersal and extinction rate
# 	if (use_optimx == FALSE)
# 		{
# 		# Set the dispersal and extinction rate
# 		d = optim_result2$par[1]
# 		e = optim_result2$par[2]
# 		j = optim_result2$par[3]
# 		} else {
# 		d = optim_result2$par[[1]][1]
# 		e = optim_result2$par[[1]][2]
# 		j = optim_result2$par[[1]][3]
# 		}
# 
# 	# Equal dispersal in all directions (unconstrained)
# 	# Equal extinction probability for all areas
# 	distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))
# 
# 	dmat = matrix(d, nrow=length(areas), ncol=length(areas))
# 	elist = rep(e, length(areas))
# 	
# 	# Set up the instantaneous rate matrix (Q matrix)
# 	Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)
# 
# 	# Cladogenic model
# 	ysv = 1-j
# 	ys = 0.5 * ysv
# 	v = 0.5 * ysv
# 	sum_SPweights = ys + j + v
# 
# 	# Text version of speciation matrix	
# 	maxent_constraint_01v = maxent_constraint_01
# 	#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
# 		
# 	# Set the parameter controlling the size distribution of 
# 	# the smaller descendant species
# 	maxent01s_param = maxent_constraint_01
# 	maxent01v_param = maxent_constraint_01
# 	maxent01j_param = maxent_constraint_01
# 	maxent01y_param = maxent_constraint_01
# 
# 	# Cladogenesis model inputs
# 	spPmat_inputs = NULL
# 	states_indices = states_list
# 	states_indices[1] = NULL	# shorten the states_indices by 1 (cutting the null range state from the speciation matrix)
# 	spPmat_inputs$l = states_indices
# 	spPmat_inputs$s = ys
# 	spPmat_inputs$v = v
# 	spPmat_inputs$j = j
# 	spPmat_inputs$y = ys
# 	spPmat_inputs$dmat = distances_mat
# 	spPmat_inputs$maxent01s_param = maxent01s_param
# 	spPmat_inputs$maxent01v_param = maxent01v_param
# 	spPmat_inputs$maxent01j_param = maxent01j_param
# 	spPmat_inputs$maxent01y_param = maxent01y_param
# 
# 	# Calculate the log-likelihood of the data, given the model parameters during this iteration	
# 	model_results = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="all", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, calc_ancprobs=TRUE)
# 
# 
# 
# 	#######################################################
# 	# Store results in a BEARS result
# 	#######################################################
# 	bears_output = model_results
# 	bears_output$optim_result = optim_result2
# 	bears_output$spPmat_inputs = spPmat_inputs
	
	return(optim_result2)
	}





# 3-parameter model, adding j (founder-event speciation)
#######################################################
# bears_3param_standard_fast
#######################################################
#' 3-parameter model, adding j (founder-event speciation)
#' 
#' This implements a 3-parameter model, basically \code{LAGRANGE} or \code{\link{bears_2param_standard_fast}}, but with a parameter \emph{j}
#' controlling the relative weight of "founder-event speciation" (\cite{Matzke_2012_IBS}) versus vicariance+sympatric speciation (which are mandated in
#' \code{LAGRANGE} and \code{\link{bears_2param_standard_fast}}.
#' 
#' @param trfn The filename of the phylogenetic tree, in NEWICK format (\url{http://evolution.genetics.washington.edu/phylip/newicktree.html}).  
#' Tipnames should match the names in geogfn.  See \code{\link[ape]{read.tree}} in APE for reading in phylogenetic trees.
#' @param geogfn A PHYLIP-style file with geographic range data (see \code{\link{getranges_from_LagrangePHYLIP}}) for each tipname. This is the same format
#' used by C++ LAGRANGE (\cite{SmithRee2010_CPPversion}).
#' @param max_range_size The maximum rangesize, in number of areas.  Having a smaller maximum range size means that you can have more areas (the size of the
#' state space is greatly reduced; see \code{\link[cladoRcpp]{numstates_from_numareas}}.
#' @param num_cores_to_use If >1, parallel processing will be attempted. \bold{Note:} parallel processing via \code{library (parallel)} will work in Mac command-line
#' R, but not in Mac GUI \code{R.app}.
#' @return \code{bears_output} A list of outputs.  bears_output$optim_result
#' @export
#' @seealso \code{\link[cladoRcpp]{numstates_from_numareas}}, \code{\link{getranges_from_LagrangePHYLIP}}, \code{\link[ape]{read.tree}}, \code{\link{calc_loglike_sp}}
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
#' # tmp hard code: extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"
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
#' bears_output = bears_3param_standard_fast(trfn=trfn, geogfn=geogfn)
#' bears_output
#' }
#'
bears_3param_standard_fast <- function(trfn = "Psychotria_5.2.newick", geogfn = "Psychotria_geog.data", max_range_size=NULL, num_cores_to_use=NULL)
	{
	defaults='
	trfn = "Psychotria_5.2.newick"
	geogfn = "Psychotria_geog.data"
	'
	
	require(cladoRcpp)
	require(rexpokit)

	# Get geographic ranges at tips
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)


	# Get the list of geographic areas
	areas = getareas_from_tipranges_object(tipranges)
	areas_list = seq(0, length(areas)-1, 1)

	# Change the names to tipranges@df:
	names(tipranges@df) = areas_list

	#######################################################
	# Set the maximum range size (can be thought of as
	# a free parameter
	#######################################################
	if (is.null(max_range_size))
		{
		max_range_size = length(areas)
		need_to_fix_rangesize = TRUE
		} else {
		need_to_fix_rangesize = FALSE
		}
	
	# Old/slow way of getting the list of states and speciation matrix (slow)
	# old_states_list = areas_list_to_states_list(areas, include_null_range=TRUE)
	# old_states_list
	# spmat = make_relprob_matrix_bi(old_states_list)
	# spmat
	
	# max_numareas = max(sapply(X=old_states_list, FUN=length), na.rm=TRUE)
	# max_numareas


	# Take the list of areas, and get list of possible states
	states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
	states_list

	if (need_to_fix_rangesize == FALSE)
		{
		max_numareas = max(sapply(X=states_list, FUN=length), na.rm=TRUE)
		max_numareas

		# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
		if (max_numareas >= 8)
			{
			force_sparse = TRUE
			} else {
			force_sparse = FALSE
			}	

		} else {
		max_numareas = max_range_size
		}
	# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
	if (max_numareas >= 8)
		{
		force_sparse = TRUE
		} else {
		force_sparse = FALSE
		}	
	
	
	# Load the phylogenetic tree
	phy = read.tree(file=trfn)


	# The likelihood of each state at the tips
	tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges, phy, states_list=states_list, maxareas=max_numareas)
	tip_condlikes_of_data_on_each_state

	maxent_constraint_01 = 0.0001


	#######################################################
	# Set up the function for optimization
	#######################################################	
	function_to_optim <- function(params, phy, tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open)
		{
		# Lagrange-like (all new species are within the geographic range of the ancestor)

		# Set the dispersal and extinction rate
		d = params[1]
		e = params[2]

		# Equal dispersal in all directions (unconstrained)
		# Equal extinction probability for all areas
		distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

		dmat = matrix(d, nrow=length(areas), ncol=length(areas))
		elist = rep(e, length(areas))
		
		# Set up the instantaneous rate matrix (Q matrix)
		Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)
	
		# Cladogenic model
		j = params[3]
		ysv = 1-j
		ys = 0.5 * ysv
		v = 0.5 * ysv
		sum_SPweights = ys + j + v

		# Text version of speciation matrix	
		maxent_constraint_01v = maxent_constraint_01
		#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
			
		# Set the parameter controlling the size distribution of 
		# the smaller descendant species
		maxent01s_param = maxent_constraint_01
		maxent01v_param = maxent_constraint_01
		maxent01j_param = maxent_constraint_01
		maxent01y_param = maxent_constraint_01


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
			# Before calculating the log likelihood, print it, in case there is e.g. a bug
			cat("d=", d, "; e=", e, "; j=", j, "; ys=", ys, "; v=", v, "; sum=", sum_SPweights, "; LnL=", sep="")
			}

		# Calculate the log-likelihood of the data, given the model parameters during this iteration	
		ttl_loglike = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="loglike", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, cluster_already_open=cluster_already_open)
		ttl_loglike
	
		if (print_optim == TRUE)
			{
			# If the log likelihood is successful, print it
			cat(ttl_loglike, "\n", sep="")
			}
	
		return(ttl_loglike)
		}

	
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
	params = c(0.01, 0.01, 0.1)
	minj = 1e-04
	# start on Lagrange results
	#params = c(3.11882,  2.51741)
	lower = c(1e-15,1e-15, minj)
	upper = c(1,1, 1-minj)
	
	# High performance computing
	# HPC using parallel package in R 2.15 or higher, which allows
	# mcmapply (multicore apply)
	# Don't use multicore if using R.app ('AQUA')
	if (.Platform$GUI != "AQUA" && ((is.null(num_cores_to_use) == TRUE) || (num_cores_to_use > 1)) )
		{
		# We are doing manual, optional processing on several cores;
		# this seems to have less overhead/hassle/incompatibility issues
		# than using mcmapply, mclapply, etc...
		#require("parallel") #<- do this higher up

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

	
	
	
	# Run optimization	
	use_optimx = TRUE
	if (use_optimx == FALSE)
		{
		optim_result2 = optim(par=params, fn=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
	
	#optim_result2 = nlminb(start=params, objective=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, lower=lower, upper=upper, control=list(iter.max=50, trace=1, abs.tol=0.001))# method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
		} else {
		# Compare methods with optimx
		#require(optimx)
		
		
		
		optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

		# Run with all methods, for testing:
		# optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

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
	

	if (exists("cluster_open") && (cluster_open == TRUE))
		{
		cat("\n\nStopping cluster with ", num_cores_to_use, " cores.\n\n", sep="")
		stopCluster(cluster_already_open)
		}
	
	
	
	
	
	#######################################################
	# Summarize results 
	#######################################################

	# Set the dispersal and extinction rate
	if (use_optimx == FALSE)
		{
		# Set the dispersal and extinction rate
		d = optim_result2$par[1]
		e = optim_result2$par[2]
		j = optim_result2$par[3]
		} else {
		d = optim_result2$par[[1]][1]
		e = optim_result2$par[[1]][2]
		j = optim_result2$par[[1]][3]
		}

	# Equal dispersal in all directions (unconstrained)
	# Equal extinction probability for all areas
	distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

	dmat = matrix(d, nrow=length(areas), ncol=length(areas))
	elist = rep(e, length(areas))
	
	# Set up the instantaneous rate matrix (Q matrix)
	Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)

	# Cladogenic model
	ysv = 1-j
	ys = 0.5 * ysv
	v = 0.5 * ysv
	sum_SPweights = ys + j + v

	# Text version of speciation matrix	
	maxent_constraint_01v = maxent_constraint_01
	#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
		
	# Set the parameter controlling the size distribution of 
	# the smaller descendant species
	maxent01s_param = maxent_constraint_01
	maxent01v_param = maxent_constraint_01
	maxent01j_param = maxent_constraint_01
	maxent01y_param = maxent_constraint_01

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

	# Calculate the log-likelihood of the data, given the model parameters during this iteration	
	model_results = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="all", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, calc_ancprobs=TRUE)



	#######################################################
	# Store results in a BEARS result
	#######################################################
	bears_output = model_results
	bears_output$optim_result = optim_result2
	bears_output$spPmat_inputs = spPmat_inputs
	
	return(bears_output)
	}












# 3-parameter model, adding v (vicariance proportion), but no j (founder-event speciation)
#######################################################
# bears_3param_standard_fast_noJ
#######################################################
#' 3-parameter model, adding v (vicariance proportion), but no j (founder-event speciation)
#' 
#' This implements a 3-parameter model, basically \code{LAGRANGE} or \code{\link{bears_2param_standard_fast}}, but with a parameter \emph{v}
#' controlling the relative weight of vicariance versus the range-copying and range-subset forms of sympatric speciation utilized in \code{LAGRANGE}
#' and \code{\link{bears_2param_standard_fast}}.
#' 
#' @param trfn The filename of the phylogenetic tree, in NEWICK format (\url{http://evolution.genetics.washington.edu/phylip/newicktree.html}).  
#' Tipnames should match the names in geogfn.  See \code{\link[ape]{read.tree}} in APE for reading in phylogenetic trees.
#' @param geogfn A PHYLIP-style file with geographic range data (see \code{\link{getranges_from_LagrangePHYLIP}}) for each tipname. This is the same format
#' used by C++ LAGRANGE (\cite{SmithRee2010_CPPversion}).
#' @param max_range_size The maximum rangesize, in number of areas.  Having a smaller maximum range size means that you can have more areas (the size of the
#' state space is greatly reduced; see \code{\link[cladoRcpp]{numstates_from_numareas}}.
#' @param num_cores_to_use If >1, parallel processing will be attempted. \bold{Note:} parallel processing via \code{library (parallel)} will work in Mac command-line
#' R, but not in Mac GUI \code{R.app}.
#' @return \code{bears_output} A list of outputs.  bears_output$optim_result
#' @export
#' @seealso \code{\link[cladoRcpp]{numstates_from_numareas}}, \code{\link{getranges_from_LagrangePHYLIP}}, \code{\link[ape]{read.tree}}, \code{\link{calc_loglike_sp}}
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
#' # tmp hard code: extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"
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
#' bears_output = bears_3param_standard_fast_noJ(trfn=trfn, geogfn=geogfn)
#' bears_output#' 
#' }
#'
bears_3param_standard_fast_noJ <- function(trfn = "Psychotria_5.2.newick", geogfn = "Psychotria_geog.data", max_range_size=NULL, num_cores_to_use=NULL)
	{
	defaults='
	trfn = "Psychotria_5.2.newick"
	geogfn = "Psychotria_geog.data"
	'
	
	require(cladoRcpp)
	require(rexpokit)

	# Get geographic ranges at tips
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)


	# Get the list of geographic areas
	areas = getareas_from_tipranges_object(tipranges)
	areas_list = seq(0, length(areas)-1, 1)

	# Change the names to tipranges@df:
	names(tipranges@df) = areas_list

	#######################################################
	# Set the maximum range size (can be thought of as
	# a free parameter
	#######################################################
	if (is.null(max_range_size))
		{
		max_range_size = length(areas)
		need_to_fix_rangesize = TRUE
		} else {
		need_to_fix_rangesize = FALSE
		}
	
	# Old/slow way of getting the list of states and speciation matrix (slow)
	# old_states_list = areas_list_to_states_list(areas, include_null_range=TRUE)
	# old_states_list
	# spmat = make_relprob_matrix_bi(old_states_list)
	# spmat
	
	# max_numareas = max(sapply(X=old_states_list, FUN=length), na.rm=TRUE)
	# max_numareas


	# Take the list of areas, and get list of possible states
	states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
	states_list

	if (need_to_fix_rangesize == FALSE)
		{
		max_numareas = max(sapply(X=states_list, FUN=length), na.rm=TRUE)
		max_numareas

		# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
		if (max_numareas >= 8)
			{
			force_sparse = TRUE
			} else {
			force_sparse = FALSE
			}	

		} else {
		max_numareas = max_range_size
		}
	# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
	if (max_numareas >= 8)
		{
		force_sparse = TRUE
		} else {
		force_sparse = FALSE
		}	
	
	
	# Load the phylogenetic tree
	phy = read.tree(file=trfn)


	# The likelihood of each state at the tips
	tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges, phy, states_list=states_list, maxareas=max_numareas)
	tip_condlikes_of_data_on_each_state

	maxent_constraint_01 = 0.0001


	#######################################################
	# Set up the function for optimization
	#######################################################	
	function_to_optim <- function(params, phy, tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open)
		{
		# Lagrange-like (all new species are within the geographic range of the ancestor)

		# Set the dispersal and extinction rate
		d = params[1]
		e = params[2]

		# Equal dispersal in all directions (unconstrained)
		# Equal extinction probability for all areas
		distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

		dmat = matrix(d, nrow=length(areas), ncol=length(areas))
		elist = rep(e, length(areas))
		
		# Set up the instantaneous rate matrix (Q matrix)
		Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)
	
		# Cladogenic model
		j = 0
		ysv = 1-j
		v = ysv * params[3]
		ys = ysv - v
		sum_SPweights = ys + j + v

		# Text version of speciation matrix	
		maxent_constraint_01v = maxent_constraint_01
		#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
			
		# Set the parameter controlling the size distribution of 
		# the smaller descendant species
		maxent01s_param = maxent_constraint_01
		maxent01v_param = maxent_constraint_01
		maxent01j_param = maxent_constraint_01
		maxent01y_param = maxent_constraint_01


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
			# Before calculating the log likelihood, print it, in case there is e.g. a bug
			cat("d=", d, "; e=", e, "; j=", j, "; ys=", ys, "; v=", v, "; sum=", sum_SPweights, "; LnL=", sep="")
			}

		# Calculate the log-likelihood of the data, given the model parameters during this iteration	
		ttl_loglike = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="loglike", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, cluster_already_open=cluster_already_open)
		ttl_loglike
	
		if (print_optim == TRUE)
			{
			# If the log likelihood is successful, print it
			cat(ttl_loglike, "\n", sep="")
			}
	
		return(ttl_loglike)
		}

	
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
	params = c(0.01, 0.01, 0.1)
	minj = 1e-05
	# start on Lagrange results
	#params = c(3.11882,  2.51741)
	lower = c(1e-15, 1e-15, minj)
	upper = c(1, 1, 1-minj)
	
	# High performance computing
	# HPC using parallel package in R 2.15 or higher, which allows
	# mcmapply (multicore apply)
	# Don't use multicore if using R.app ('AQUA')
	if (.Platform$GUI != "AQUA" && ((is.null(num_cores_to_use) == TRUE) || (num_cores_to_use > 1)) )
		{
		# We are doing manual, optional processing on several cores;
		# this seems to have less overhead/hassle/incompatibility issues
		# than using mcmapply, mclapply, etc...
		#require("parallel") #<- do this higher up

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

	
	
	
	# Run optimization	
	use_optimx = TRUE
	if (use_optimx == FALSE)
		{
		optim_result2 = optim(par=params, fn=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
	
	#optim_result2 = nlminb(start=params, objective=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, lower=lower, upper=upper, control=list(iter.max=50, trace=1, abs.tol=0.001))# method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
		} else {
		# Compare methods with optimx
		#require(optimx)
		
		
		
		optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

		# Run with all methods, for testing:
		# optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

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
	

	if (exists("cluster_open") && (cluster_open == TRUE))
		{
		cat("\n\nStopping cluster with ", num_cores_to_use, " cores.\n\n", sep="")
		stopCluster(cluster_already_open)
		}
	
	
	
	
	#######################################################
	# Summarize results 
	#######################################################

	# Set the dispersal and extinction rate
	if (use_optimx == FALSE)
		{
		# Set the dispersal and extinction rate
		d = optim_result2$par[1]
		e = optim_result2$par[2]
		v = optim_result2$par[3]
		} else {
		d = optim_result2$par[[1]][1]
		e = optim_result2$par[[1]][2]
		v = optim_result2$par[[1]][3]
		}
	j = 0
	# Equal dispersal in all directions (unconstrained)
	# Equal extinction probability for all areas
	distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

	dmat = matrix(d, nrow=length(areas), ncol=length(areas))
	elist = rep(e, length(areas))
	
	# Set up the instantaneous rate matrix (Q matrix)
	Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)

	# Cladogenic model
	ysv = 1-j
	v = ysv * v
	ys = ysv - v
	sum_SPweights = ys + j + v

	# Text version of speciation matrix	
	maxent_constraint_01v = maxent_constraint_01
	#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)

	# Set the parameter controlling the size distribution of 
	# the smaller descendant species
	maxent01s_param = maxent_constraint_01
	maxent01v_param = maxent_constraint_01
	maxent01j_param = maxent_constraint_01
	maxent01y_param = maxent_constraint_01

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

	# Calculate the log-likelihood of the data, given the model parameters during this iteration	
	model_results = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="all", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, calc_ancprobs=TRUE)



	#######################################################
	# Store results in a BEARS result
	#######################################################
	bears_output = model_results
	bears_output$optim_result = optim_result2
	bears_output$spPmat_inputs = spPmat_inputs
	
	return(bears_output)
	}
















# 4-parameter model, adding j (founder-event speciation) and v (vicariance proportion)
#######################################################
# bears_4param_standard_fast
#######################################################
#' 4-parameter model, adding j (founder-event speciation) and v (vicariance proportion)
#' 
#' This implements a 4-parameter model, basically \code{LAGRANGE} or \code{\link{bears_2param_standard_fast}}, but with a parameter \emph{j}
#' controlling the relative weight of "founder-event speciation" (\cite{Matzke_2012_IBS}) and another parameter \emph{v}
#' controlling the relative weight of vicariance.  The remainder of the weight (weights must sum to 1) is taken up by the range-copying and 
#' range-subset forms of sympatric speciation utilized in \code{LAGRANGE} and \code{\link{bears_2param_standard_fast}}.
#' 
#' @param trfn The filename of the phylogenetic tree, in NEWICK format (\url{http://evolution.genetics.washington.edu/phylip/newicktree.html}).  
#' Tipnames should match the names in geogfn.  See \code{\link[ape]{read.tree}} in APE for reading in phylogenetic trees.
#' @param geogfn A PHYLIP-style file with geographic range data (see \code{\link{getranges_from_LagrangePHYLIP}}) for each tipname. This is the same format
#' used by C++ LAGRANGE (\cite{SmithRee2010_CPPversion}).
#' @param max_range_size The maximum rangesize, in number of areas.  Having a smaller maximum range size means that you can have more areas (the size of the
#' state space is greatly reduced; see \code{\link[cladoRcpp]{numstates_from_numareas}}.
#' @param num_cores_to_use If >1, parallel processing will be attempted. \bold{Note:} parallel processing via \code{library (parallel)} will work in Mac command-line
#' R, but not in Mac GUI \code{R.app}.
#' @return \code{bears_output} A list of outputs.  bears_output$optim_result
#' @export
#' @seealso \code{\link[cladoRcpp]{numstates_from_numareas}}, \code{\link{getranges_from_LagrangePHYLIP}}, \code{\link[ape]{read.tree}}, \code{\link{calc_loglike_sp}}
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
#' # tmp hard code: extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"
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
#' bears_output = bears_4param_standard_fast(trfn=trfn, geogfn=geogfn)
#' bears_output
#' }
#'
bears_4param_standard_fast <- function(trfn = "Psychotria_5.2.newick", geogfn = "Psychotria_geog.data", max_range_size=NULL, num_cores_to_use=NULL)
	{
	defaults='
	trfn = "Psychotria_5.2.newick"
	geogfn = "Psychotria_geog.data"
	'
	
	require(cladoRcpp)
	require(rexpokit)

	# Get geographic ranges at tips
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)


	# Get the list of geographic areas
	areas = getareas_from_tipranges_object(tipranges)
	areas_list = seq(0, length(areas)-1, 1)

	# Change the names to tipranges@df:
	names(tipranges@df) = areas_list

	#######################################################
	# Set the maximum range size (can be thought of as
	# a free parameter
	#######################################################
	if (is.null(max_range_size))
		{
		max_range_size = length(areas)
		need_to_fix_rangesize = TRUE
		} else {
		need_to_fix_rangesize = FALSE
		}
	
	# Old/slow way of getting the list of states and speciation matrix (slow)
	# old_states_list = areas_list_to_states_list(areas, include_null_range=TRUE)
	# old_states_list
	# spmat = make_relprob_matrix_bi(old_states_list)
	# spmat
	
	# max_numareas = max(sapply(X=old_states_list, FUN=length), na.rm=TRUE)
	# max_numareas


	# Take the list of areas, and get list of possible states
	states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
	states_list

	if (need_to_fix_rangesize == FALSE)
		{
		max_numareas = max(sapply(X=states_list, FUN=length), na.rm=TRUE)
		max_numareas

		# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
		if (max_numareas >= 8)
			{
			force_sparse = TRUE
			} else {
			force_sparse = FALSE
			}	

		} else {
		max_numareas = max_range_size
		}
	# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
	if (max_numareas >= 8)
		{
		force_sparse = TRUE
		} else {
		force_sparse = FALSE
		}	
	
	
	# Load the phylogenetic tree
	phy = read.tree(file=trfn)


	# The likelihood of each state at the tips
	tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges, phy, states_list=states_list, maxareas=max_numareas)
	tip_condlikes_of_data_on_each_state

	maxent_constraint_01 = 0.0001


	#######################################################
	# Set up the function for optimization
	#######################################################	
	function_to_optim <- function(params, phy, tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open)
		{
		# Lagrange-like (all new species are within the geographic range of the ancestor)

		# Set the dispersal and extinction rate
		d = params[1]
		e = params[2]

		# Equal dispersal in all directions (unconstrained)
		# Equal extinction probability for all areas
		distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

		dmat = matrix(d, nrow=length(areas), ncol=length(areas))
		elist = rep(e, length(areas))
		
		# Set up the instantaneous rate matrix (Q matrix)
		Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)
	
		# Cladogenic model
		j = params[3]
		ysv = 1-j
		v = ysv * params[4]
		ys = ysv - v
		sum_SPweights = ys + j + v

		# Text version of speciation matrix	
		maxent_constraint_01v = maxent_constraint_01
		#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
			
		# Set the parameter controlling the size distribution of 
		# the smaller descendant species
		maxent01s_param = maxent_constraint_01
		maxent01v_param = maxent_constraint_01
		maxent01j_param = maxent_constraint_01
		maxent01y_param = maxent_constraint_01


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
			# Before calculating the log likelihood, print it, in case there is e.g. a bug
			cat("d=", d, "; e=", e, "; j=", j, "; ys=", ys, "; v=", v, "; maxent_constraint_01=", maxent_constraint_01, "; sum=", sum_SPweights, "; LnL=", sep="")
			}

		# Calculate the log-likelihood of the data, given the model parameters during this iteration	
		ttl_loglike = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="loglike", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, cluster_already_open=cluster_already_open)
		ttl_loglike
	
		if (print_optim == TRUE)
			{
			# If the log likelihood is successful, print it
			cat(ttl_loglike, "\n", sep="")
			}
	
	
		return(ttl_loglike)
		#return(-1*ttl_loglike)
		}

	
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
	params = c(0.01, 0.01, 0.1, 0.1)
	minj = 1e-05
	# start on Lagrange results
	#params = c(3.11882,  2.51741)
	lower = c(1e-15, 1e-15, minj, minj)
	upper = c(1, 1, 1-minj, 1-minj)
	
	
	# High performance computing
	# HPC using parallel package in R 2.15 or higher, which allows
	# mcmapply (multicore apply)
	# Don't use multicore if using R.app ('AQUA')
	if (.Platform$GUI != "AQUA" && ((is.null(num_cores_to_use) == TRUE) || (num_cores_to_use > 1)) )
		{
		# We are doing manual, optional processing on several cores;
		# this seems to have less overhead/hassle/incompatibility issues
		# than using mcmapply, mclapply, etc...
		#require("parallel") #<- do this higher up

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

	
	
	
	# Run optimization	
	use_optimx = TRUE
	if (use_optimx == FALSE)
		{
		optim_result2 = optim(par=params, fn=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
	
	#optim_result2 = nlminb(start=params, objective=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, lower=lower, upper=upper, control=list(iter.max=50, trace=1, abs.tol=0.001))# method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
		} else {
		# Compare methods with optimx
		#require(optimx)
		
		
		
		optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

		# Run with all methods, for testing:
		# optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

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
	

	if (exists("cluster_open") && (cluster_open == TRUE))
		{
		cat("\n\nStopping cluster with ", num_cores_to_use, " cores.\n\n", sep="")
		stopCluster(cluster_already_open)
		}
	
	# Compare methods with optimx
	#optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, control=list(all.methods=TRUE, maximize=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))
	
	# This might work
	# doesn't seem to find optimum
	# optim_result2 = nlminb(start=params, objective=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, lower=lower, upper=upper, control=list(iter.max=50, trace=1))# method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
	
	
	# This doesn't hit optimum
	#optim_result2 = nlminb(start=params, objective=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, lower=lower, upper=upper, control=list(iter.max=50, trace=1, abs.tol=0.001))# method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
	
	# Neither does this
	#optim_result2 = nlminb(start=params, objective=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, lower=lower, upper=upper, control=list(iter.max=50, trace=1, abs.tol=0.001, step.min=0.01, step.max=0.01))# method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
	
	
	#######################################################
	# Summarize results 
	#######################################################

	if (use_optimx == FALSE)
		{
		# Set the dispersal and extinction rate
		d = optim_result2$par[1]
		e = optim_result2$par[2]
		j = optim_result2$par[3]
		v = optim_result2$par[4]
		} else {
		d = optim_result2$par[[1]][1]
		e = optim_result2$par[[1]][2]
		j = optim_result2$par[[1]][3]
		v = optim_result2$par[[1]][4]
		}
	# Equal dispersal in all directions (unconstrained)
	# Equal extinction probability for all areas
	distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

	dmat = matrix(d, nrow=length(areas), ncol=length(areas))
	elist = rep(e, length(areas))
	
	# Set up the instantaneous rate matrix (Q matrix)
	Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)

	# Cladogenic model
	ysv = 1-j
	v = ysv * v
	ys = ysv - v
	sum_SPweights = ys + j + v

	# Text version of speciation matrix	
	maxent_constraint_01v = maxent_constraint_01
	#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
		
	# Set the parameter controlling the size distribution of 
	# the smaller descendant species
	maxent01s_param = maxent_constraint_01
	maxent01v_param = maxent_constraint_01
	maxent01j_param = maxent_constraint_01
	maxent01y_param = maxent_constraint_01

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

	# Calculate the log-likelihood of the data, given the model parameters during this iteration	
	model_results = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="all", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, calc_ancprobs=TRUE)



	#######################################################
	# Store results in a BEARS result
	#######################################################
	bears_output = model_results
	bears_output$optim_result = optim_result2
	bears_output$spPmat_inputs = spPmat_inputs
	
	return(bears_output)
	}









# 5-parameter model, adding j, v (vicariance proportion), maxent_constraint_01 (weighting for size of smaller offspring)
#######################################################
# bears_5param_standard_fast
#######################################################
#' 5-parameter model, adding j (founder-event speciation), v (vicariance proportion), and maxent_constraint_01 (weighting for size of smaller-ranged descendant lineage)
#' 
#' This implements a 5-parameter model, basically \code{LAGRANGE} or \code{\link{bears_2param_standard_fast}}, but with a parameter \emph{j}
#' controlling the relative weight of "founder-event speciation" (\cite{Matzke_2012_IBS}), and another parameter \emph{v}
#' controlling the relative weight of vicariance.  The remainder of the weight (weights must sum to 1) is taken up by the range-copying and 
#' range-subset forms of sympatric speciation utilized in \code{LAGRANGE} and \code{\link{bears_2param_standard_fast}}.  A fifth parameter, 
#' \emph{maxent_constraint_01}, controls the relative probability of 
#' daughter lineages of different rangesizes.  If maxent_constraint_01=0.0001, the smaller-ranged daughter lineage will have size 1 area, with probability 1.  If 
#' maxent_constraint_01=0.5, all different rangesizes will have equal probability, and if maxent_constraint_01=0.9999, the largest possible range will have
#' probability 1.  
#' 
#' @param trfn The filename of the phylogenetic tree, in NEWICK format (\url{http://evolution.genetics.washington.edu/phylip/newicktree.html}).  
#' Tipnames should match the names in geogfn.  See \code{\link[ape]{read.tree}} in APE for reading in phylogenetic trees.
#' @param geogfn A PHYLIP-style file with geographic range data (see \code{\link{getranges_from_LagrangePHYLIP}}) for each tipname. This is the same format
#' used by C++ LAGRANGE (\cite{SmithRee2010_CPPversion}).
#' @param max_range_size The maximum rangesize, in number of areas.  Having a smaller maximum range size means that you can have more areas (the size of the
#' state space is greatly reduced; see \code{\link[cladoRcpp]{numstates_from_numareas}}.
#' @param num_cores_to_use If >1, parallel processing will be attempted. \bold{Note:} parallel processing via \code{library (parallel)} will work in Mac command-line
#' R, but not in Mac GUI \code{R.app}.
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
#' extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
#' # tmp hard code: extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"
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
#' bears_output = bears_5param_standard_fast(trfn=trfn, geogfn=geogfn)
#' bears_output
#' }
#'
bears_5param_standard_fast <- function(trfn = "Psychotria_5.2.newick", geogfn = "Psychotria_geog.data", max_range_size=NULL, num_cores_to_use=NULL)
	{
	defaults='
	trfn = "Psychotria_5.2.newick"
	geogfn = "Psychotria_geog.data"
	'
	
	require(cladoRcpp)
	require(rexpokit)

	# Get geographic ranges at tips
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)


	# Get the list of geographic areas
	areas = getareas_from_tipranges_object(tipranges)
	areas_list = seq(0, length(areas)-1, 1)

	# Change the names to tipranges@df:
	names(tipranges@df) = areas_list

	#######################################################
	# Set the maximum range size (can be thought of as
	# a free parameter
	#######################################################
	if (is.null(max_range_size))
		{
		max_range_size = length(areas)
		need_to_fix_rangesize = TRUE
		} else {
		need_to_fix_rangesize = FALSE
		}
	
	# Old/slow way of getting the list of states and speciation matrix (slow)
	# old_states_list = areas_list_to_states_list(areas, include_null_range=TRUE)
	# old_states_list
	# spmat = make_relprob_matrix_bi(old_states_list)
	# spmat
	
	# max_numareas = max(sapply(X=old_states_list, FUN=length), na.rm=TRUE)
	# max_numareas


	# Take the list of areas, and get list of possible states
	states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
	states_list

	if (need_to_fix_rangesize == FALSE)
		{
		max_numareas = max(sapply(X=states_list, FUN=length), na.rm=TRUE)
		max_numareas

		# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
		if (max_numareas >= 8)
			{
			force_sparse = TRUE
			} else {
			force_sparse = FALSE
			}	

		} else {
		max_numareas = max_range_size
		}
	# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
	if (max_numareas >= 8)
		{
		force_sparse = TRUE
		} else {
		force_sparse = FALSE
		}	
	
	
	# Load the phylogenetic tree
	phy = read.tree(file=trfn)


	# The likelihood of each state at the tips
	tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges, phy, states_list=states_list, maxareas=max_numareas)
	tip_condlikes_of_data_on_each_state

	maxent_constraint_01 = 0.0001


	#######################################################
	# Set up the function for optimization
	#######################################################	
	function_to_optim <- function(params, phy, tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open)
		{
		# Lagrange-like (all new species are within the geographic range of the ancestor)

		# Set the dispersal and extinction rate
		d = params[1]
		e = params[2]

		# Equal dispersal in all directions (unconstrained)
		# Equal extinction probability for all areas
		distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

		dmat = matrix(d, nrow=length(areas), ncol=length(areas))
		elist = rep(e, length(areas))
		
		# Set up the instantaneous rate matrix (Q matrix)
		Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)
	
		# Cladogenic model
		j = params[3]
		ysv = 1-j
		v = ysv * params[4]
		ys = ysv - v
		sum_SPweights = ys + j + v

		maxent_constraint_01 = params[5]
		
		# Text version of speciation matrix	
		maxent_constraint_01v = maxent_constraint_01
		#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
			
		# Set the parameter controlling the size distribution of 
		# the smaller descendant species
		maxent01s_param = maxent_constraint_01
		maxent01v_param = maxent_constraint_01
		maxent01j_param = maxent_constraint_01
		maxent01y_param = maxent_constraint_01


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
			# Before calculating the log likelihood, print it, in case there is e.g. a bug
			cat("d=", d, "; e=", e, "; j=", j, "; ys=", ys, "; v=", v, "; maxent_constraint_01=", maxent_constraint_01, "; sum=", sum_SPweights, "; LnL=", sep="")
			}

		# Calculate the log-likelihood of the data, given the model parameters during this iteration	
		ttl_loglike = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="loglike", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, cluster_already_open=cluster_already_open)
		ttl_loglike
	
		if (print_optim == TRUE)
			{
			# If the log likelihood is successful, print it
			cat(ttl_loglike, "\n", sep="")
			}
	
	
		return(ttl_loglike)
		}

	
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
	params = c(0.01, 0.01, 0.1, 0.1, 0.0001)
	minj = 1e-05
	# start on Lagrange results
	#params = c(3.11882,  2.51741)
	lower = c(1e-15, 1e-15, minj, minj, minj)
	upper = c(1, 1, 1-minj, 1-minj, 1-minj)
	
	# High performance computing
	# HPC using parallel package in R 2.15 or higher, which allows
	# mcmapply (multicore apply)
	# Don't use multicore if using R.app ('AQUA')
	if (.Platform$GUI != "AQUA" && ((is.null(num_cores_to_use) == TRUE) || (num_cores_to_use > 1)) )
		{
		# We are doing manual, optional processing on several cores;
		# this seems to have less overhead/hassle/incompatibility issues
		# than using mcmapply, mclapply, etc...
		#require("parallel") #<- do this higher up

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

	
	
	
	
	# Run optimization	
	use_optimx = TRUE
	if (use_optimx == FALSE)
		{
		optim_result2 = optim(par=params, fn=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
	
	#optim_result2 = nlminb(start=params, objective=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, lower=lower, upper=upper, control=list(iter.max=50, trace=1, abs.tol=0.001))# method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
		} else {
		# Compare methods with optimx
		#require(optimx)
		
		
		
		optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

		# Run with all methods, for testing:
		# optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

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
	

	if (exists("cluster_open") && (cluster_open == TRUE))
		{
		cat("\n\nStopping cluster with ", num_cores_to_use, " cores.\n\n", sep="")
		stopCluster(cluster_already_open)
		}
	

	
	#######################################################
	# Summarize results 
	#######################################################

	# Set the dispersal and extinction rate
	if (use_optimx == FALSE)
		{
		# Set the dispersal and extinction rate
		d = optim_result2$par[1]
		e = optim_result2$par[2]
		j = optim_result2$par[3]
		v = optim_result2$par[4]
		maxent_constraint_01 = optim_result2$par[5]
		} else {
		d = optim_result2$par[[1]][1]
		e = optim_result2$par[[1]][2]
		j = optim_result2$par[[1]][3]
		v = optim_result2$par[[1]][4]
		maxent_constraint_01 = optim_result2$par[[1]][5]
		}
	# Equal dispersal in all directions (unconstrained)
	# Equal extinction probability for all areas
	distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

	dmat = matrix(d, nrow=length(areas), ncol=length(areas))
	elist = rep(e, length(areas))
	
	# Set up the instantaneous rate matrix (Q matrix)
	Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)

	# Cladogenic model
	ysv = 1-j
	v = ysv * v
	ys = ysv - v
	sum_SPweights = ys + j + v

	# Text version of speciation matrix	
	maxent_constraint_01v = maxent_constraint_01
	#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
		
	# Set the parameter controlling the size distribution of 
	# the smaller descendant species
	maxent01s_param = maxent_constraint_01
	maxent01v_param = maxent_constraint_01
	maxent01j_param = maxent_constraint_01
	maxent01y_param = maxent_constraint_01

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

	# Calculate the log-likelihood of the data, given the model parameters during this iteration	
	model_results = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="all", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, calc_ancprobs=TRUE)



	#######################################################
	# Store results in a BEARS result
	#######################################################
	bears_output = model_results
	bears_output$optim_result = optim_result2
	bears_output$spPmat_inputs = spPmat_inputs
	
	return(bears_output)
	}



# 5-parameter model, adding j, v (vicariance proportion), maxent_constraint_01 (weighting for size of smaller offspring)
#######################################################
# bears_5param_standard_fast_diffstart
#######################################################
#' 5-parameter model, with different starting points for optimization
#' 
#' This implements the same model as \code{\link{bears_5param_standard_fast}}, but uses different starting points and slightly different constraints.
#'
#' As the number of parameters increases, the importance of starting ML optimization runs from different places increases.  Several starting points should
#' be tried, especially if the likelihood surface seems flat.
#' 
#' @param trfn The filename of the phylogenetic tree, in NEWICK format (\url{http://evolution.genetics.washington.edu/phylip/newicktree.html}).  
#' Tipnames should match the names in geogfn.  See \code{\link[ape]{read.tree}} in APE for reading in phylogenetic trees.
#' @param geogfn A PHYLIP-style file with geographic range data (see \code{\link{getranges_from_LagrangePHYLIP}}) for each tipname. This is the same format
#' used by C++ LAGRANGE (\cite{SmithRee2010_CPPversion}).
#' @param max_range_size The maximum rangesize, in number of areas.  Having a smaller maximum range size means that you can have more areas (the size of the
#' state space is greatly reduced; see \code{\link[cladoRcpp]{numstates_from_numareas}}.
#' @param num_cores_to_use If >1, parallel processing will be attempted. \bold{Note:} parallel processing via \code{library (parallel)} will work in Mac command-line
#' R, but not in Mac GUI \code{R.app}.
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
#' extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
#' # tmp hard code: extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"
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
#' bears_output = bears_5param_standard_fast_diffstart(trfn=trfn, geogfn=geogfn)
#' bears_output
#' }
#'
bears_5param_standard_fast_diffstart <- function(trfn = "Psychotria_5.2.newick", geogfn = "Psychotria_geog.data", max_range_size=NULL, num_cores_to_use=NULL)
	{
	defaults='
	trfn = "Psychotria_5.2.newick"
	geogfn = "Psychotria_geog.data"
	max_range_size=NULL
	num_cores_to_use=NULL
	'
	
	require(cladoRcpp)
	require(rexpokit)

	# Get geographic ranges at tips
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)


	# Get the list of geographic areas
	areas = getareas_from_tipranges_object(tipranges)
	areas_list = seq(0, length(areas)-1, 1)

	# Change the names to tipranges@df:
	names(tipranges@df) = areas_list

	#######################################################
	# Set the maximum range size (can be thought of as
	# a free parameter
	#######################################################
	if (is.null(max_range_size))
		{
		max_range_size = length(areas)
		need_to_fix_rangesize = TRUE
		} else {
		need_to_fix_rangesize = FALSE
		}
	
	# Old/slow way of getting the list of states and speciation matrix (slow)
	# old_states_list = areas_list_to_states_list(areas, include_null_range=TRUE)
	# old_states_list
	# spmat = make_relprob_matrix_bi(old_states_list)
	# spmat
	
	# max_numareas = max(sapply(X=old_states_list, FUN=length), na.rm=TRUE)
	# max_numareas


	# Take the list of areas, and get list of possible states
	states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
	states_list

	if (need_to_fix_rangesize == FALSE)
		{
		max_numareas = max(sapply(X=states_list, FUN=length), na.rm=TRUE)
		max_numareas

		# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
		if (max_numareas >= 8)
			{
			force_sparse = TRUE
			} else {
			force_sparse = FALSE
			}	

		} else {
		max_numareas = max_range_size
		}
	# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
	if (max_numareas >= 8)
		{
		force_sparse = TRUE
		} else {
		force_sparse = FALSE
		}	
	
	
	# Load the phylogenetic tree
	phy = read.tree(file=trfn)


	# The likelihood of each state at the tips
	tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges, phy, states_list=states_list, maxareas=max_numareas)
	tip_condlikes_of_data_on_each_state

	maxent_constraint_01 = 0.0001


	#######################################################
	# Set up the function for optimization
	#######################################################	
	function_to_optim <- function(params, phy, tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open)
		{
		# Lagrange-like (all new species are within the geographic range of the ancestor)

		# Set the dispersal and extinction rate
		d = params[1]
		e = params[2]

		# Equal dispersal in all directions (unconstrained)
		# Equal extinction probability for all areas
		distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

		dmat = matrix(d, nrow=length(areas), ncol=length(areas))
		elist = rep(e, length(areas))
		
		# Set up the instantaneous rate matrix (Q matrix)
		Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)
	
		# Cladogenic model
		j = params[3]
		ysv = 1-j
		v = ysv * params[4]
		ys = ysv - v
		sum_SPweights = ys + j + v

		# test d=0.01588153; e=0.006352611; j=0.07886359; ys=0.7831373; v=0.1379991;

		maxent_constraint_01 = params[5]
		# maxent_constraint_01=0.464636; 		
		
		# Text version of speciation matrix	
		maxent_constraint_01v = maxent_constraint_01
		#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
			
		# Set the parameter controlling the size distribution of 
		# the smaller descendant species
		maxent01s_param = maxent_constraint_01
		maxent01v_param = maxent_constraint_01
		maxent01j_param = maxent_constraint_01
		maxent01y_param = maxent_constraint_01


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
			# Before calculating the log likelihood, print it, in case there is e.g. a bug
			cat("d=", d, "; e=", e, "; j=", j, "; ys=", ys, "; v=", v, "; maxent_constraint_01=", maxent_constraint_01, "; sum=", sum_SPweights, "; LnL=", sep="")
			}

		# Calculate the log-likelihood of the data, given the model parameters during this iteration	
		ttl_loglike = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="loglike", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, cluster_already_open=cluster_already_open)
		ttl_loglike
	
		if (print_optim == TRUE)
			{
			# If the log likelihood is successful, print it
			cat(ttl_loglike, "\n", sep="")
			}
	
	
		return(ttl_loglike)
		}

	
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
	params = c(0.05, 0.02, 0.2, 0.3, 0.5)
	minj = 1e-05
	# start on Lagrange results
	#params = c(3.11882,  2.51741)
	lower = c(1e-15, 1e-15, minj, minj, 0.0002)
	upper = c(1, 1, 1-minj, 1-minj, 1-0.0002)
	
	# High performance computing
	# HPC using parallel package in R 2.15 or higher, which allows
	# mcmapply (multicore apply)
	# Don't use multicore if using R.app ('AQUA')
	if (.Platform$GUI != "AQUA" && ((is.null(num_cores_to_use) == TRUE) || (num_cores_to_use > 1)) )
		{
		# We are doing manual, optional processing on several cores;
		# this seems to have less overhead/hassle/incompatibility issues
		# than using mcmapply, mclapply, etc...
		#require("parallel") #<- do this higher up

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

	
	
	
	
	# Run optimization	
	use_optimx = TRUE
	if (use_optimx == FALSE)
		{
		optim_result2 = optim(par=params, fn=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
	
	#optim_result2 = nlminb(start=params, objective=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, lower=lower, upper=upper, control=list(iter.max=50, trace=1, abs.tol=0.001))# method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
		} else {
		# Compare methods with optimx
		#require(optimx)
		
		
		
		optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

		# Run with all methods, for testing:
		# optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

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
	

	if (exists("cluster_open") && (cluster_open == TRUE))
		{
		cat("\n\nStopping cluster with ", num_cores_to_use, " cores.\n\n", sep="")
		stopCluster(cluster_already_open)
		}
	

	
	#######################################################
	# Summarize results 
	#######################################################

	# Set the dispersal and extinction rate
	if (use_optimx == FALSE)
		{
		# Set the dispersal and extinction rate
		d = optim_result2$par[1]
		e = optim_result2$par[2]
		j = optim_result2$par[3]
		v = optim_result2$par[4]
		maxent_constraint_01 = optim_result2$par[5]
		} else {
		d = optim_result2$par[[1]][1]
		e = optim_result2$par[[1]][2]
		j = optim_result2$par[[1]][3]
		v = optim_result2$par[[1]][4]
		maxent_constraint_01 = optim_result2$par[[1]][5]
		}
	# Equal dispersal in all directions (unconstrained)
	# Equal extinction probability for all areas
	distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

	dmat = matrix(d, nrow=length(areas), ncol=length(areas))
	elist = rep(e, length(areas))
	
	# Set up the instantaneous rate matrix (Q matrix)
	Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)

	# Cladogenic model
	ysv = 1-j
	v = ysv * v
	ys = ysv - v
	sum_SPweights = ys + j + v

	# Text version of speciation matrix	
	maxent_constraint_01v = maxent_constraint_01
	#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
		
	# Set the parameter controlling the size distribution of 
	# the smaller descendant species
	maxent01s_param = maxent_constraint_01
	maxent01v_param = maxent_constraint_01
	maxent01j_param = maxent_constraint_01
	maxent01y_param = maxent_constraint_01

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

	# Calculate the log-likelihood of the data, given the model parameters during this iteration	
	model_results = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="all", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, calc_ancprobs=TRUE)



	#######################################################
	# Store results in a BEARS result
	#######################################################
	bears_output = model_results
	bears_output$optim_result = optim_result2
	bears_output$spPmat_inputs = spPmat_inputs
	
	return(bears_output)
	}









# 5-parameter model, adding j, v (vicariance proportion), maxent_constraint_01v (weighting for size of smaller offspring)
#######################################################
# bears_5param_standard_fast_v
#######################################################
#' 5-parameter model, adding j (founder-event speciation), v (vicariance proportion), and maxent_constraint_01v (vicariance daughter sizes)
#' 
#' This implements a 5-parameter model, basically \code{LAGRANGE} or \code{\link{bears_2param_standard_fast}}, but with a parameter \emph{j}
#' controlling the relative weight of "founder-event speciation" (\cite{Matzke_2012_IBS}), and another parameter \emph{v}
#' controlling the relative weight of vicariance.  The remainder of the weight (weights must sum to 1) is taken up by the range-copying and 
#' range-subset forms of sympatric speciation utilized in \code{LAGRANGE} and \code{\link{bears_2param_standard_fast}}.  A fifth parameter, 
#' \emph{maxent_constraint_01v}, controls the relative probability of 
#' daughter lineages of different rangesizes, but only for the vicariance events, which are rather different from other types of speciation events.  
#' If maxent_constraint_01v=0.0001, the smaller-ranged daughter lineage will have size 1 area, with probability 1.  If 
#' maxent_constraint_01v=0.5, all different rangesizes will have equal probability, and if maxent_constraint_01v=0.9999, the largest possible range will have
#' probability 1 -- but note that in a vicariance context, this would mean at maximum rangesize of 50% of the areas.
#'
#' Non-vicariance events have hard-coded \code{maxent_constraint_01=0.0001}
#' 
#' @param trfn The filename of the phylogenetic tree, in NEWICK format (\url{http://evolution.genetics.washington.edu/phylip/newicktree.html}).  
#' Tipnames should match the names in geogfn.  See \code{\link[ape]{read.tree}} in APE for reading in phylogenetic trees.
#' @param geogfn A PHYLIP-style file with geographic range data (see \code{\link{getranges_from_LagrangePHYLIP}}) for each tipname. This is the same format
#' used by C++ LAGRANGE (\cite{SmithRee2010_CPPversion}).
#' @param max_range_size The maximum rangesize, in number of areas.  Having a smaller maximum range size means that you can have more areas (the size of the
#' state space is greatly reduced; see \code{\link[cladoRcpp]{numstates_from_numareas}}.
#' @param num_cores_to_use If >1, parallel processing will be attempted. \bold{Note:} parallel processing via \code{library (parallel)} will work in Mac command-line
#' R, but not in Mac GUI \code{R.app}.
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
#' extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
#' # tmp hard code: extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"
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
#' bears_output = bears_5param_standard_fast_v(trfn=trfn, geogfn=geogfn)
#' bears_output
#' }
#'
bears_5param_standard_fast_v <- function(trfn = "Psychotria_5.2.newick", geogfn = "Psychotria_geog.data", max_range_size=NULL, num_cores_to_use=NULL)
	{
	defaults='
	trfn = "Psychotria_5.2.newick"
	geogfn = "Psychotria_geog.data"
	'
	
	require(cladoRcpp)
	require(rexpokit)

	# Get geographic ranges at tips
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)


	# Get the list of geographic areas
	areas = getareas_from_tipranges_object(tipranges)
	areas_list = seq(0, length(areas)-1, 1)

	# Change the names to tipranges@df:
	names(tipranges@df) = areas_list

	#######################################################
	# Set the maximum range size (can be thought of as
	# a free parameter
	#######################################################
	if (is.null(max_range_size))
		{
		max_range_size = length(areas)
		need_to_fix_rangesize = TRUE
		} else {
		need_to_fix_rangesize = FALSE
		}
	
	# Old/slow way of getting the list of states and speciation matrix (slow)
	# old_states_list = areas_list_to_states_list(areas, include_null_range=TRUE)
	# old_states_list
	# spmat = make_relprob_matrix_bi(old_states_list)
	# spmat
	
	# max_numareas = max(sapply(X=old_states_list, FUN=length), na.rm=TRUE)
	# max_numareas


	# Take the list of areas, and get list of possible states
	states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
	states_list

	if (need_to_fix_rangesize == FALSE)
		{
		max_numareas = max(sapply(X=states_list, FUN=length), na.rm=TRUE)
		max_numareas

		# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
		if (max_numareas >= 8)
			{
			force_sparse = TRUE
			} else {
			force_sparse = FALSE
			}	

		} else {
		max_numareas = max_range_size
		}
	# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
	if (max_numareas >= 8)
		{
		force_sparse = TRUE
		} else {
		force_sparse = FALSE
		}	
	
	
	# Load the phylogenetic tree
	phy = read.tree(file=trfn)


	# The likelihood of each state at the tips
	tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges, phy, states_list=states_list, maxareas=max_numareas)
	tip_condlikes_of_data_on_each_state

	maxent_constraint_01 = 0.0001


	#######################################################
	# Set up the function for optimization
	#######################################################	
	function_to_optim <- function(params, phy, tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open)
		{
		# Lagrange-like (all new species are within the geographic range of the ancestor)

		# Set the dispersal and extinction rate
		d = params[1]
		e = params[2]

		# Equal dispersal in all directions (unconstrained)
		# Equal extinction probability for all areas
		distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

		dmat = matrix(d, nrow=length(areas), ncol=length(areas))
		elist = rep(e, length(areas))
		
		# Set up the instantaneous rate matrix (Q matrix)
		Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)
	
		# Cladogenic model
		j = params[3]
		ysv = 1-j
		v = ysv * params[4]
		ys = ysv - v
		sum_SPweights = ys + j + v

		maxent_constraint_01 = 0.0001
		
		# Text version of speciation matrix	
		maxent_constraint_01v = maxent_constraint_01
		#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
			
		# Set the parameter controlling the size distribution of 
		# the smaller descendant species
		maxent01s_param = maxent_constraint_01
		maxent01v_param = params[5]
		maxent01j_param = maxent_constraint_01
		maxent01y_param = maxent_constraint_01


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
			# Before calculating the log likelihood, print it, in case there is e.g. a bug
			cat("d=", d, "; e=", e, "; j=", j, "; ys=", ys, "; v=", v, "; maxent01v=", maxent_constraint_01v, "; sum=", sum_SPweights, "; LnL=", sep="")
			}

		# Calculate the log-likelihood of the data, given the model parameters during this iteration	
		ttl_loglike = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="loglike", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, cluster_already_open=cluster_already_open)
		ttl_loglike
	
		if (print_optim == TRUE)
			{
			# If the log likelihood is successful, print it
			cat(ttl_loglike, "\n", sep="")
			}
	
		
		return(ttl_loglike)
		}

	
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
	params = c(0.01, 0.01, 0.1, 0.1, 0.0001)
	minj = 1e-05
	# start on Lagrange results
	#params = c(3.11882,  2.51741)
	lower = c(1e-15, 1e-15, minj, minj, minj)
	upper = c(1, 1, 1-minj, 1-minj, 1-minj)
	
	# High performance computing
	# HPC using parallel package in R 2.15 or higher, which allows
	# mcmapply (multicore apply)
	# Don't use multicore if using R.app ('AQUA')
	if (.Platform$GUI != "AQUA" && ((is.null(num_cores_to_use) == TRUE) || (num_cores_to_use > 1)) )
		{
		# We are doing manual, optional processing on several cores;
		# this seems to have less overhead/hassle/incompatibility issues
		# than using mcmapply, mclapply, etc...
		#require("parallel") #<- do this higher up

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

	
	
	
	# Run optimization	
	use_optimx = TRUE
	if (use_optimx == FALSE)
		{
		optim_result2 = optim(par=params, fn=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
	
	#optim_result2 = nlminb(start=params, objective=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, lower=lower, upper=upper, control=list(iter.max=50, trace=1, abs.tol=0.001))# method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
		} else {
		# Compare methods with optimx
		#require(optimx)
		
		
		
		optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

		# Run with all methods, for testing:
		# optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

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
	

	if (exists("cluster_open") && (cluster_open == TRUE))
		{
		cat("\n\nStopping cluster with ", num_cores_to_use, " cores.\n\n", sep="")
		stopCluster(cluster_already_open)
		}
	
	
	
	
	
	#######################################################
	# Summarize results 
	#######################################################

	# Set the dispersal and extinction rate
	if (use_optimx == FALSE)
		{
		# Set the dispersal and extinction rate
		d = optim_result2$par[1]
		e = optim_result2$par[2]
		j = optim_result2$par[3]
		v = optim_result2$par[4]
		maxent01v_param = optim_result2$par[5]
		} else {
		d = optim_result2$par[[1]][1]
		e = optim_result2$par[[1]][2]
		j = optim_result2$par[[1]][3]
		v = optim_result2$par[[1]][4]
		maxent01v_param = optim_result2$par[[1]][5]
		}
	# Equal dispersal in all directions (unconstrained)
	# Equal extinction probability for all areas
	distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

	dmat = matrix(d, nrow=length(areas), ncol=length(areas))
	elist = rep(e, length(areas))
	
	# Set up the instantaneous rate matrix (Q matrix)
	Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)

	# Cladogenic model
	ysv = 1-j
	v = ysv * v
	ys = ysv - v
	sum_SPweights = ys + j + v

	maxent_constraint_01 = 0.0001

	# Text version of speciation matrix	
	maxent_constraint_01v = maxent01v_param
	#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
		
	# Set the parameter controlling the size distribution of 
	# the smaller descendant species
	maxent01s_param = maxent_constraint_01
	maxent01v_param = maxent_constraint_01v
	maxent01j_param = maxent_constraint_01
	maxent01y_param = maxent_constraint_01

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


	# Calculate the log-likelihood of the data, given the model parameters during this iteration	
	model_results = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="all", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, calc_ancprobs=TRUE)



	#######################################################
	# Store results in a BEARS result
	#######################################################
	bears_output = model_results
	bears_output$optim_result = optim_result2
	bears_output$spPmat_inputs = spPmat_inputs
	
	return(bears_output)
	}









# 6-parameter model, adding j, v (vicariance proportion), maxent_constraint_01 (for non-vicariant subsets), maxent_constraint_01v (weighting for size of smaller offspring)
#######################################################
# bears_6param_standard_fast_ys_v
#######################################################
#' 6-parameter model, adding j (founder-event speciation), v (vicariance proportion), and both maxent_constraint_01 and maxent_constraint_01v
#' 
#' This implements a 6-parameter model, basically \code{LAGRANGE} or \code{\link{bears_2param_standard_fast}}, but with a parameter \emph{j}
#' controlling the relative weight of "founder-event speciation" (\cite{Matzke_2012_IBS}), and another parameter \emph{v}
#' controlling the relative weight of vicariance.  The remainder of the weight (weights must sum to 1) is taken up by the range-copying and 
#' range-subset forms of sympatric speciation utilized in \code{LAGRANGE} and \code{\link{bears_2param_standard_fast}}.  A fifth parameter, 
#' \emph{maxent_constraint_01}, controls the relative probability of 
#' daughter lineages of different rangesizes.  A sixth parameter, 
#' \emph{maxent_constraint_01v}, controls the relative probability of 
#' daughter lineages of different rangesizes, but only for the vicariance events, which are rather different from other types of speciation events.  
#' If maxent_constraint_01v=0.0001, the smaller-ranged daughter lineage will have size 1 area, with probability 1.  If 
#' maxent_constraint_01v=0.5, all different rangesizes will have equal probability, and if maxent_constraint_01v=0.9999, the largest possible range will have
#' probability 1 -- but note that in a vicariance context, this would mean at maximum rangesize of 50% of the areas.
#'
#' @param trfn The filename of the phylogenetic tree, in NEWICK format (\url{http://evolution.genetics.washington.edu/phylip/newicktree.html}).  
#' Tipnames should match the names in geogfn.  See \code{\link[ape]{read.tree}} in APE for reading in phylogenetic trees.
#' @param geogfn A PHYLIP-style file with geographic range data (see \code{\link{getranges_from_LagrangePHYLIP}}) for each tipname. This is the same format
#' used by C++ LAGRANGE (\cite{SmithRee2010_CPPversion}).
#' @param max_range_size The maximum rangesize, in number of areas.  Having a smaller maximum range size means that you can have more areas (the size of the
#' state space is greatly reduced; see \code{\link[cladoRcpp]{numstates_from_numareas}}.
#' @param num_cores_to_use If >1, parallel processing will be attempted. \bold{Note:} parallel processing via \code{library (parallel)} will work in Mac command-line
#' R, but not in Mac GUI \code{R.app}.
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
#' extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
#' # tmp hard code: extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"
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
#' bears_output = bears_6param_standard_fast_ys_v(trfn=trfn, geogfn=geogfn)
#' bears_output
#' }
#'
bears_6param_standard_fast_ys_v <- function(trfn = "Psychotria_5.2.newick", geogfn = "Psychotria_geog.data", max_range_size=NULL, num_cores_to_use=NULL)
	{
	defaults='
	trfn = "Psychotria_5.2.newick"
	geogfn = "Psychotria_geog.data"
	'
	
	require(cladoRcpp)
	require(rexpokit)

	# Get geographic ranges at tips
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)


	# Get the list of geographic areas
	areas = getareas_from_tipranges_object(tipranges)
	areas_list = seq(0, length(areas)-1, 1)

	# Change the names to tipranges@df:
	names(tipranges@df) = areas_list

	#######################################################
	# Set the maximum range size (can be thought of as
	# a free parameter
	#######################################################
	if (is.null(max_range_size))
		{
		max_range_size = length(areas)
		need_to_fix_rangesize = TRUE
		} else {
		need_to_fix_rangesize = FALSE
		}
	
	# Old/slow way of getting the list of states and speciation matrix (slow)
	# old_states_list = areas_list_to_states_list(areas, include_null_range=TRUE)
	# old_states_list
	# spmat = make_relprob_matrix_bi(old_states_list)
	# spmat
	
	# max_numareas = max(sapply(X=old_states_list, FUN=length), na.rm=TRUE)
	# max_numareas


	# Take the list of areas, and get list of possible states
	states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
	states_list

	if (need_to_fix_rangesize == FALSE)
		{
		max_numareas = max(sapply(X=states_list, FUN=length), na.rm=TRUE)
		max_numareas

		# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
		if (max_numareas >= 8)
			{
			force_sparse = TRUE
			} else {
			force_sparse = FALSE
			}	

		} else {
		max_numareas = max_range_size
		}
	# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
	if (max_numareas >= 8)
		{
		force_sparse = TRUE
		} else {
		force_sparse = FALSE
		}	
	
	
	# Load the phylogenetic tree
	phy = read.tree(file=trfn)


	# The likelihood of each state at the tips
	tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges, phy, states_list=states_list, maxareas=max_numareas)
	tip_condlikes_of_data_on_each_state

	maxent_constraint_01 = 0.0001


	#######################################################
	# Set up the function for optimization
	#######################################################	
	function_to_optim <- function(params, phy, tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open)
		{
		# Lagrange-like (all new species are within the geographic range of the ancestor)

		# Set the dispersal and extinction rate
		d = params[1]
		e = params[2]

		# Equal dispersal in all directions (unconstrained)
		# Equal extinction probability for all areas
		distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

		dmat = matrix(d, nrow=length(areas), ncol=length(areas))
		elist = rep(e, length(areas))
		
		# Set up the instantaneous rate matrix (Q matrix)
		Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)
	
		# Cladogenic model
		j = params[3]
		ysv = 1-j
		v = ysv * params[4]
		ys = ysv - v
		sum_SPweights = ys + j + v

		maxent_constraint_01 = params[5]
		
		# Text version of speciation matrix	
		maxent_constraint_01v = params[6]
		#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
			
		# Set the parameter controlling the size distribution of 
		# the smaller descendant species
		maxent01s_param = maxent_constraint_01
		maxent01v_param = maxent_constraint_01v
		maxent01j_param = maxent_constraint_01
		maxent01y_param = maxent_constraint_01


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
			# Before calculating the log likelihood, print it, in case there is e.g. a bug
			cat("d=", d, "; e=", e, "; j=", j, "; ys=", ys, "; v=", v, "; maxent01=", maxent_constraint_01, "; maxent01v=", maxent_constraint_01v, "; sum=", sum_SPweights, "; LnL=", sep="")
			}

		# Calculate the log-likelihood of the data, given the model parameters during this iteration	
		ttl_loglike = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="loglike", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, cluster_already_open=cluster_already_open)
		ttl_loglike
	
		if (print_optim == TRUE)
			{
			# If the log likelihood is successful, print it
			cat(ttl_loglike, "\n", sep="")
			}
		
		return(ttl_loglike)
		}

	
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
	params = c(0.01, 0.01, 0.1, 0.1, 0.0001, 0.0001)
	minj = 1e-05
	# start on Lagrange results
	#params = c(3.11882,  2.51741)
	lower = c(1e-15, 1e-15, minj, minj, minj, minj)
	upper = c(1, 1, 1-minj, 1-minj, 1-minj, 1-minj)
	
	# High performance computing
	# HPC using parallel package in R 2.15 or higher, which allows
	# mcmapply (multicore apply)
	# Don't use multicore if using R.app ('AQUA')
	if (.Platform$GUI != "AQUA" && ((is.null(num_cores_to_use) == TRUE) || (num_cores_to_use > 1)) )
		{
		# We are doing manual, optional processing on several cores;
		# this seems to have less overhead/hassle/incompatibility issues
		# than using mcmapply, mclapply, etc...
		#require("parallel") #<- do this higher up

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

	
	
	
	# Run optimization	
	use_optimx = TRUE
	if (use_optimx == FALSE)
		{
		optim_result2 = optim(par=params, fn=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
	
	#optim_result2 = nlminb(start=params, objective=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, lower=lower, upper=upper, control=list(iter.max=50, trace=1, abs.tol=0.001))# method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
		} else {
		# Compare methods with optimx
		#require(optimx)
		
		
		
		optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

		# Run with all methods, for testing:
		# optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

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
	

	if (exists("cluster_open") && (cluster_open == TRUE))
		{
		cat("\n\nStopping cluster with ", num_cores_to_use, " cores.\n\n", sep="")
		stopCluster(cluster_already_open)
		}
	
	
	
	
	
	#######################################################
	# Summarize results 
	#######################################################

	# Set the dispersal and extinction rate
	if (use_optimx == FALSE)
		{
		# Set the dispersal and extinction rate
		d = optim_result2$par[1]
		e = optim_result2$par[2]
		j = optim_result2$par[3]
		v = optim_result2$par[4]
		maxent01_param = optim_result2$par[5]
		maxent01v_param = optim_result2$par[6]
		} else {
		d = optim_result2$par[[1]][1]
		e = optim_result2$par[[1]][2]
		j = optim_result2$par[[1]][3]
		v = optim_result2$par[[1]][4]
		maxent01_param = optim_result2$par[[1]][5]
		maxent01v_param = optim_result2$par[[1]][6]
		}
	# Equal dispersal in all directions (unconstrained)
	# Equal extinction probability for all areas
	distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

	dmat = matrix(d, nrow=length(areas), ncol=length(areas))
	elist = rep(e, length(areas))
	
	# Set up the instantaneous rate matrix (Q matrix)
	Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)

	# Cladogenic model
	ysv = 1-j
	v = ysv * v
	ys = ysv - v
	sum_SPweights = ys + j + v

	maxent_constraint_01 = maxent01_param

	# Text version of speciation matrix	
	maxent_constraint_01v = maxent01v_param
	#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
		
	# Set the parameter controlling the size distribution of 
	# the smaller descendant species
	maxent01s_param = maxent_constraint_01
	maxent01v_param = maxent_constraint_01v
	maxent01j_param = maxent_constraint_01
	maxent01y_param = maxent_constraint_01

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


	# Calculate the log-likelihood of the data, given the model parameters during this iteration	
	model_results = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="all", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, calc_ancprobs=TRUE)



	#######################################################
	# Store results in a BEARS result
	#######################################################
	bears_output = model_results
	bears_output$optim_result = optim_result2
	bears_output$spPmat_inputs = spPmat_inputs
	
	return(bears_output)
	}











# 6-parameter model, adding j, v (vicariance proportion), maxent_constraint_01 (for non-vicariant subsets), maxent_constraint_01v (weighting for size of smaller offspring)
#######################################################
# bears_5param_standard_fast_v
#######################################################
#' 6-parameter model, adding j (founder-event speciation), v (vicariance proportion), and both maxent_constraint_01 and maxent_constraint_01v
#' 
#' This implements a 6-parameter model, basically \code{LAGRANGE} or \code{\link{bears_2param_standard_fast}}, but with a parameter \emph{j}
#' controlling the relative weight of "founder-event speciation" (\cite{Matzke_2012_IBS}), and another parameter \emph{v}
#' controlling the relative weight of vicariance.  The remainder of the weight (weights must sum to 1) is taken up by the range-copying and 
#' range-subset forms of sympatric speciation utilized in \code{LAGRANGE} and \code{\link{bears_2param_standard_fast}}.  A fifth parameter, 
#' \emph{maxent_constraint_01}, controls the relative probability of 
#' daughter lineages of different rangesizes.  A sixth parameter, 
#' \emph{maxent_constraint_01v}, controls the relative probability of 
#' daughter lineages of different rangesizes, but only for the vicariance events, which are rather different from other types of speciation events.  
#' If maxent_constraint_01v=0.0001, the smaller-ranged daughter lineage will have size 1 area, with probability 1.  If 
#' maxent_constraint_01v=0.5, all different rangesizes will have equal probability, and if maxent_constraint_01v=0.9999, the largest possible range will have
#' probability 1 -- but note that in a vicariance context, this would mean at maximum rangesize of 50% of the areas.
#'
#' @param trfn The filename of the phylogenetic tree, in NEWICK format (\url{http://evolution.genetics.washington.edu/phylip/newicktree.html}).  
#' Tipnames should match the names in geogfn.  See \code{\link[ape]{read.tree}} in APE for reading in phylogenetic trees.
#' @param geogfn A PHYLIP-style file with geographic range data (see \code{\link{getranges_from_LagrangePHYLIP}}) for each tipname. This is the same format
#' used by C++ LAGRANGE (\cite{SmithRee2010_CPPversion}).
#' @param max_range_size The maximum rangesize, in number of areas.  Having a smaller maximum range size means that you can have more areas (the size of the
#' state space is greatly reduced; see \code{\link[cladoRcpp]{numstates_from_numareas}}.
#' @param num_cores_to_use If >1, parallel processing will be attempted. \bold{Note:} parallel processing via \code{library (parallel)} will work in Mac command-line
#' R, but not in Mac GUI \code{R.app}.
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
#' extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
#' # tmp hard code: extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"
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
#' bears_output = bears_6param_standard_fast_ys_v(trfn=trfn, geogfn=geogfn)
#' bears_output
#' }
#' # Add c (direct changes), b (branch length exponent), x (distance exponent)
bears_9param_standard_fast_ys_v_cb <- function(trfn = "Psychotria_5.2.newick", geogfn = "Psychotria_geog.data", max_range_size=NULL, num_cores_to_use=NULL)
	{
	defaults='
	trfn = "Psychotria_5.2.newick"
	geogfn = "Psychotria_geog.data"
	'
	
	require(cladoRcpp)
	require(rexpokit)

	# Get geographic ranges at tips
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)


	# Get the list of geographic areas
	areas = getareas_from_tipranges_object(tipranges)
	areas_list = seq(0, length(areas)-1, 1)

	# Change the names to tipranges@df:
	names(tipranges@df) = areas_list

	#######################################################
	# Set the maximum range size (can be thought of as
	# a free parameter
	#######################################################
	if (is.null(max_range_size))
		{
		max_range_size = length(areas)
		need_to_fix_rangesize = TRUE
		} else {
		need_to_fix_rangesize = FALSE
		}
	
	# Old/slow way of getting the list of states and speciation matrix (slow)
	# old_states_list = areas_list_to_states_list(areas, include_null_range=TRUE)
	# old_states_list
	# spmat = make_relprob_matrix_bi(old_states_list)
	# spmat
	
	# max_numareas = max(sapply(X=old_states_list, FUN=length), na.rm=TRUE)
	# max_numareas


	# Take the list of areas, and get list of possible states
	states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
	states_list

	if (need_to_fix_rangesize == FALSE)
		{
		max_numareas = max(sapply(X=states_list, FUN=length), na.rm=TRUE)
		max_numareas

		# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
		if (max_numareas >= 8)
			{
			force_sparse = TRUE
			} else {
			force_sparse = FALSE
			}	

		} else {
		max_numareas = max_range_size
		}
	# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
	if (max_numareas >= 8)
		{
		force_sparse = TRUE
		} else {
		force_sparse = FALSE
		}	
	
	
	# Load the phylogenetic tree
	phy = read.tree(file=trfn)


	# The likelihood of each state at the tips
	tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges, phy, states_list=states_list, maxareas=max_numareas)
	tip_condlikes_of_data_on_each_state

	maxent_constraint_01 = 0.0001


	#######################################################
	# Set up the function for optimization
	#######################################################	
	function_to_optim <- function(params, phy, tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open)
		{
		# Lagrange-like (all new species are within the geographic range of the ancestor)

		# Set the dispersal and extinction rate
		d = params[1]
		e = params[2]

		# Equal dispersal in all directions (unconstrained)
		# Equal extinction probability for all areas
		distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

		dmat = matrix(d, nrow=length(areas), ncol=length(areas))
		elist = rep(e, length(areas))
		
		# Set up the instantaneous rate matrix (Q matrix)
		Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)
	
		# Cladogenic model
		j = params[3]
		ysv = 1-j
		v = ysv * params[4]
		ys = ysv - v
		sum_SPweights = ys + j + v

		maxent_constraint_01 = params[5]
		
		# Text version of speciation matrix	
		maxent_constraint_01v = params[6]
		#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
			
		# Set the parameter controlling the size distribution of 
		# the smaller descendant species
		maxent01s_param = maxent_constraint_01
		maxent01v_param = maxent_constraint_01v
		maxent01j_param = maxent_constraint_01
		maxent01y_param = maxent_constraint_01


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
			# Before calculating the log likelihood, print it, in case there is e.g. a bug
			cat("d=", d, "; e=", e, "; j=", j, "; ys=", ys, "; v=", v, "; maxent01=", maxent_constraint_01, "; maxent01v=", maxent_constraint_01v, "; sum=", sum_SPweights, "; LnL=", sep="")
			}

		# Calculate the log-likelihood of the data, given the model parameters during this iteration	
		ttl_loglike = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="loglike", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, cluster_already_open=cluster_already_open)
		ttl_loglike
	
		if (print_optim == TRUE)
			{
			# If the log likelihood is successful, print it
			cat(ttl_loglike, "\n", sep="")
			}
		
		return(ttl_loglike)
		}

	
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
	params = c(0.01, 0.01, 0.1, 0.1, 0.0001, 0.0001)
	minj = 1e-05
	# start on Lagrange results
	#params = c(3.11882,  2.51741)
	lower = c(1e-15, 1e-15, minj, minj, minj, minj)
	upper = c(1, 1, 1-minj, 1-minj, 1-minj, 1-minj)
	
	# High performance computing
	# HPC using parallel package in R 2.15 or higher, which allows
	# mcmapply (multicore apply)
	# Don't use multicore if using R.app ('AQUA')
	if (.Platform$GUI != "AQUA" && ((is.null(num_cores_to_use) == TRUE) || (num_cores_to_use > 1)) )
		{
		# We are doing manual, optional processing on several cores;
		# this seems to have less overhead/hassle/incompatibility issues
		# than using mcmapply, mclapply, etc...
		#require("parallel") #<- do this higher up

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

	
	
	
	# Run optimization	
	use_optimx = TRUE
	if (use_optimx == FALSE)
		{
		optim_result2 = optim(par=params, fn=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open, method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
	
	#optim_result2 = nlminb(start=params, objective=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, lower=lower, upper=upper, control=list(iter.max=50, trace=1, abs.tol=0.001))# method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
		} else {
		# Compare methods with optimx
		#require(optimx)
		
		
		
		optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

		# Run with all methods, for testing:
		# optim_result2 = optimx(par=params, fn=function_to_optim, lower=lower, upper=upper, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, areas_list=areas_list, states_list=states_list, force_sparse=force_sparse, cluster_already_open=cluster_already_open,itnmax=250, method=c("bobyqa"), control=list(all.methods=FALSE, maximize=TRUE, save.failures=TRUE))# method="L-BFGS-B", control=list(fnscale=-1, trace=2, maxit=500))

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
	

	if (exists("cluster_open") && (cluster_open == TRUE))
		{
		cat("\n\nStopping cluster with ", num_cores_to_use, " cores.\n\n", sep="")
		stopCluster(cluster_already_open)
		}
	
	
	
	
	
	#######################################################
	# Summarize results 
	#######################################################

	# Set the dispersal and extinction rate
	if (use_optimx == FALSE)
		{
		# Set the dispersal and extinction rate
		d = optim_result2$par[1]
		e = optim_result2$par[2]
		j = optim_result2$par[3]
		v = optim_result2$par[4]
		maxent01_param = optim_result2$par[5]
		maxent01v_param = optim_result2$par[6]
		} else {
		d = optim_result2$par[[1]][1]
		e = optim_result2$par[[1]][2]
		j = optim_result2$par[[1]][3]
		v = optim_result2$par[[1]][4]
		maxent01_param = optim_result2$par[[1]][5]
		maxent01v_param = optim_result2$par[[1]][6]
		}
	# Equal dispersal in all directions (unconstrained)
	# Equal extinction probability for all areas
	distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

	dmat = matrix(d, nrow=length(areas), ncol=length(areas))
	elist = rep(e, length(areas))
	
	# Set up the instantaneous rate matrix (Q matrix)
	Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)

	# Cladogenic model
	ysv = 1-j
	v = ysv * v
	ys = ysv - v
	sum_SPweights = ys + j + v

	maxent_constraint_01 = maxent01_param

	# Text version of speciation matrix	
	maxent_constraint_01v = maxent01v_param
	#spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
		
	# Set the parameter controlling the size distribution of 
	# the smaller descendant species
	maxent01s_param = maxent_constraint_01
	maxent01v_param = maxent_constraint_01v
	maxent01j_param = maxent_constraint_01
	maxent01y_param = maxent_constraint_01

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


	# Calculate the log-likelihood of the data, given the model parameters during this iteration	
	model_results = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=NULL, return_what="all", sparse=force_sparse, use_cpp=TRUE, input_is_COO=force_sparse, spPmat_inputs=spPmat_inputs, printlevel=0, calc_ancprobs=TRUE)



	#######################################################
	# Store results in a BEARS result
	#######################################################
	bears_output = model_results
	bears_output$optim_result = optim_result2
	bears_output$spPmat_inputs = spPmat_inputs
	
	return(bears_output)
	}




















#######################################################
# Some basic types of runs
#######################################################

#######################################################
# bears_2param_standard_slowQ_slowSP
#######################################################
#' 2-parameter model, fixed cladogenesis model -- slow version
#' 
#' This implements the same 2-parameter model found in \code{LAGRANGE} or \code{\link{bears_2param_standard_fast}}, but using the original
#' slower options for matrix exponentiation and cladogenesis events.
#' 
#' @param trfn The filename of the phylogenetic tree, in NEWICK format (\url{http://evolution.genetics.washington.edu/phylip/newicktree.html}).  
#' Tipnames should match the names in geogfn.  See \code{\link[ape]{read.tree}} in APE for reading in phylogenetic trees.
#' @param geogfn A PHYLIP-style file with geographic range data (see \code{\link{getranges_from_LagrangePHYLIP}}) for each tipname. This is the same format
#' used by C++ LAGRANGE (\cite{SmithRee2010_CPPversion}).
#' @return \code{bears_output} A list of outputs.  bears_output$optim_result
#' @param max_range_size The maximum rangesize, in number of areas.  Having a smaller maximum range size means that you can have more areas (the size of the
#' state space is greatly reduced; see \code{\link[cladoRcpp]{numstates_from_numareas}}.
#' @export
#' @seealso \code{\link[cladoRcpp]{numstates_from_numareas}}, \code{\link{getranges_from_LagrangePHYLIP}}, \code{\link[ape]{read.tree}}, \code{\link{calc_loglike_sp}}
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
#' # tmp hard code: extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"
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
#' bears_output = bears_2param_standard_slowQ_slowSP(trfn=trfn, geogfn=geogfn)
#' bears_output
#' }
#'
bears_2param_standard_slowQ_slowSP <- function(trfn = "Psychotria_5.2.newick", geogfn = "Psychotria_geog.data", max_range_size=NULL)
	{
	defaults='
	trfn = "Psychotria_5.2.newick"
	geogfn = "Psychotria_geog.data"
	'
	
	require(cladoRcpp)
	require(rexpokit)

	# Get geographic ranges at tips
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)


	# Get the list of geographic areas
	areas = getareas_from_tipranges_object(tipranges)
	areas_list = seq(0, length(areas)-1, 1)

	# Change the names to tipranges@df:
	names(tipranges@df) = areas_list

	
	# Old/slow way of getting the list of states and speciation matrix (slow)
	old_states_list = cladoRcpp:::areas_list_to_states_list_old(areas, include_null_range=TRUE)
	old_states_list
	spmat = make_relprob_matrix_bi(old_states_list)
	spmat
	
	# max_numareas = max(sapply(X=old_states_list, FUN=length), na.rm=TRUE)
	# max_numareas

	#######################################################
	# Set the maximum range size (can be thought of as
	# a free parameter
	#######################################################
	if (is.null(max_range_size))
		{
		max_range_size = length(areas)
		need_to_fix_rangesize = TRUE
		} else {
		need_to_fix_rangesize = FALSE
		}



	if (need_to_fix_rangesize == FALSE)
		{
		max_numareas = max(sapply(X=states_list, FUN=length), na.rm=TRUE)
		max_numareas

		# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
		if (max_numareas >= 8)
			{
			force_sparse = TRUE
			} else {
			force_sparse = FALSE
			}	

		} else {
		max_numareas = max_range_size
		}
	# Check for large rate matrix (>= 7 areas; here, sparse matrices will work better)
	if (max_numareas >= 8)
		{
		force_sparse = TRUE
		} else {
		force_sparse = FALSE
		}	

	# Take the list of areas, and get list of possible states
	states_list = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
	states_list
	
	
	# Load the phylogenetic tree
	phy = read.tree(file=trfn)


	# The likelihood of each state at the tips
	tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges, phy, states_list=states_list, maxareas=max_numareas)
	tip_condlikes_of_data_on_each_state

	maxent_constraint_01 = 0.0001


	#######################################################
	# Set up the function for optimization
	#######################################################	
	function_to_optim <- function(params, phy, tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, spmat=spmat, areas_list=areas_list, states_list=states_list)
		{
		# Lagrange-like (all new species are within the geographic range of the ancestor)

		# Set the dispersal and extinction rate
		d = params[1]
		e = params[2]

		# Equal dispersal in all directions (unconstrained)
		# Equal extinction probability for all areas
		distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

		dmat = matrix(d, nrow=length(areas), ncol=length(areas))
		elist = rep(e, length(areas))
		
		# Set up the instantaneous rate matrix (Q matrix)
		Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=FALSE)
	
		# Cladogenic model
		ys = 0.5
		j = 0
		v = 0.5
		sum_SPweights = ys + j + v

		# Text version of speciation matrix	
		maxent_constraint_01v = maxent_constraint_01
		spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
			
		# Set the parameter controlling the size distribution of 
		# the smaller descendant species
		maxent01s_param = maxent_constraint_01
		maxent01v_param = maxent_constraint_01
		maxent01j_param = maxent_constraint_01
		maxent01y_param = maxent_constraint_01
	
		# Calculate the log-likelihood of the data, given the model parameters during this iteration	
		ttl_loglike = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=spPmat, return_what="loglike", sparse=FALSE, use_cpp=FALSE, input_is_COO=FALSE, spPmat_inputs=NULL, printlevel=0)
		ttl_loglike
	
		if (print_optim == TRUE)
			{
			cat("d=", d, "; e=", e, "; j=", j, "; ys=", ys, "; v=", v, ";\tsum=", sum_SPweights, ";\tLnL=", ttl_loglike, "\n", sep="")
			}
	
		return(ttl_loglike)
		}

	
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
	params = c(0.01, 0.01)
	minj = 1e-05
	# start on Lagrange results
	#params = c(3.11882,  2.51741)
	lower = c(1e-15,1e-15)
	upper = c(1,1)
	
	# Run optimization	
	optim_result2 = optim(par=params, fn=function_to_optim, phy=phy, tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, print_optim=TRUE, maxent_constraint_01=maxent_constraint_01, spmat=spmat, areas_list=areas_list, states_list=states_list, method="L-BFGS-B", lower=lower, upper=upper, control=list(fnscale=-1, trace=2, maxit=500))
	
	
	
	
	
	#######################################################
	# Summarize results 
	#######################################################

	# Set the dispersal and extinction rate
	d = optim_result2$par[1]
	e = optim_result2$par[2]

	# Equal dispersal in all directions (unconstrained)
	# Equal extinction probability for all areas
	distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))

	dmat = matrix(d, nrow=length(areas), ncol=length(areas))
	elist = rep(e, length(areas))
	
	# Set up the instantaneous rate matrix (Q matrix)
	Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=FALSE)

	# Cladogenic model
	ys = 0.5
	j = 0
	v = 0.5
	sum_SPweights = ys + j + v

	# Text version of speciation matrix	
	maxent_constraint_01v = maxent_constraint_01
	spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", mergesym="*", ys=ys, j=j, v=v, maxent_constraint_01=maxent_constraint_01, maxent_constraint_01v=maxent_constraint_01v, max_numareas=max_numareas)
		
	# Set the parameter controlling the size distribution of 
	# the smaller descendant species
	maxent01s_param = maxent_constraint_01
	maxent01v_param = maxent_constraint_01
	maxent01j_param = maxent_constraint_01
	maxent01y_param = maxent_constraint_01

	# Calculate the log-likelihood of the data, given the model parameters during this iteration	
	model_results = calc_loglike_sp(tip_condlikes_of_data_on_each_state=tip_condlikes_of_data_on_each_state, phy=phy, Qmat=Qmat, spPmat=spPmat, return_what="all", sparse=FALSE, use_cpp=FALSE, input_is_COO=FALSE, spPmat_inputs=NULL, printlevel=0)


	#######################################################
	# Store results in a BEARS result
	#######################################################
	bears_output = model_results
	bears_output$optim_result = optim_result2
	
	return(bears_output)
	}









