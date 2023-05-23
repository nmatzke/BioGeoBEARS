
#require("ape")
#require("rexpokit")
#require("cladoRcpp")


#######################################################
# Defining Classes
#######################################################
# Based on advice here:
# http://stackoverflow.com/questions/7368262/how-to-properly-document-s4-class-slots-using-roxygen2




# Set up a default BioGeoBEARS model
# continuous parameters
# implements LAGRANGE model


#######################################################
# BioGeoBEARS_model_defaults
#######################################################
#' Set up a default BioGeoBEARS model object
#' 
#' The default starting model is the DEC model (DEC=Dispersal-Extinction-Cladogenesis)
#' of Ree & Smith (2008). The starting parameters are d=0.01 and e=0.01. All of the 
#' other BioGeoBEARS models can be created by modifying the setup of the DEC model 
#' (adding free parameters, changing the weight of different cladogenesis events, etc.).
#' This does \emph{not} mean that all of these models are special cases of DEC. The 
#' phrase "X is a special case of Y" means that model X can be produced by taking a 
#' free parameter in Y and setting that parameter to a fixed value equal to the fixed
#' value assumed in model X. Model X is said to be "nested" in Y.
#' 
#' Thus, DEC is a special case of DEC+J (when j=0), but not the reverse.
#' DEC, DIVALIKE and BAYAREALIKE are not special cases of each other, but they are
#' each, respectively, special cases of DEC+x, DIVALIKE+x, and BAYAREALIKE+x, for 
#' instance, when x is fixed to x=0. 
#'
#' In the general way that all of these models involve dispersal, 
#' extinction/extirpation, and cladogenesis, they can reasonably be termed "DEC-type"
#' models, but that does not mean they are the same model as the DEC model originally
#' implemented in the program \code{lagrange}. Similarly, the idea that DEC is a 
#' "natural" or "default" model for biogeography must be questioned -- there are 
#' many other models that could have been (and still could be) used in biogeography;
#' the ubiquity of the DEC model is basically a historial accident.
#' 
#' @param minval_anagenesis Minimum value above zero for d, e, a, b parameters.
#' @param minval_cladogenesis Minimum value above zero for j, v, etc.
#' @param maxval Maximum value for d, e, a
#' @return \code{param_table} Return the parameter table object
#' @export
#' @seealso \code{\link[base]{rbind}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @examples
#' test=1
BioGeoBEARS_model_defaults <- function(minval_anagenesis=1e-12, minval_cladogenesis=1e-5, maxval=5)
	{
	defaults='
	minval_anagenesis = 1e-12
	minval_cladogenesis = 1e-5
	maxval = 50
	'
	
	param_data_starter = as.data.frame(matrix(data=0, nrow=1, ncol=7), stringsAsFactors=FALSE)
	names(param_data_starter) = c("type", "init", "min", "max", "est", "note", "desc")
	param_table = NULL
	param_names = NULL
	
	# For a parameter on a branch (anagenesis)
	
	# dispersal rate
	param_name = "d"
	param_data = param_data_starter
	param_data$type = "free"
	param_data$init = 0.01
	param_data$min = 0 + minval_anagenesis
	param_data$max = maxval - minval_anagenesis
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "anagenesis: rate of 'dispersal' (range expansion)"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table
	
	# extinction rate
	param_name = "e"
	param_data = param_data_starter
	param_data$type = "free"
	param_data$init = 0.01
	param_data$min = 0 + minval_anagenesis
	param_data$max = maxval - minval_anagenesis
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "anagenesis: rate of 'extinction' (range contraction)"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table

	# between-ranges transition rate (anagenic; i.e., with a standard morphology model)
	param_name = "a"
	param_data = param_data_starter
	param_data$type = "fixed"
	param_data$init = 0.0
	param_data$min = 0 + minval_anagenesis
	param_data$max = maxval - minval_anagenesis
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "anagenesis: rate of range-switching (i.e. for a standard char.)"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table

	# branch-length exponent (0: all branches=1; 1: use branchlengths; 2: use branchlengths^2, not meaningful)
	param_name = "b"
	param_data = param_data_starter
	param_data$type = "fixed"
	param_data$init = 1.0
	param_data$min = 0 + minval_anagenesis
	param_data$max = 1 - minval_anagenesis
	param_data$est = param_data$init
	param_data$note = "non-stratified only"
	param_data$desc = "anagenesis: exponent on branch lengths"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table

	# Exponent on distance (dist^-1 = inverse distance weighting; dist^-0 = 1 (equal probs); dist^-2 = exponential dispersal decay)
	param_name = "x"
	param_data = param_data_starter
	param_data$type = "fixed"
	param_data$init = 0.0
	param_data$min = -2.5
	param_data$max = 2.5
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "exponent on distance (modifies d, j, a)"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table

	# Exponent on environmental distance (dist^-1 = inverse distance weighting; 
	# dist^-0 = 1 (equal probs); dist^-2 = exponential dispersal decay)
	param_name = "n"
	param_data = param_data_starter
	param_data$type = "fixed"
	param_data$init = 0.0
	param_data$min = -10
	param_data$max = 10
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "exponent on environmental distance (modifies d, j, a)"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table

	# Exponent on manual dispersal multipliers (mult^<0 = reverse multipliers; 
	# mult^0 = 1 (equal probs); mult^>0 = strengthen multipliers)
	param_name = "w"
	param_data = param_data_starter
	param_data$type = "fixed"
	param_data$init = 1.0
	param_data$min = -10
	param_data$max = 10
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "exponent on manual dispersal multipliers (modifies d, j, a)"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table


	# Exponent on extinction risk with area (area^-0, equiprobable extinction; area^-1, prob(ext) is inverse to area; area^-2, prob declines exponentially)
	param_name = "u"
	param_data = param_data_starter
	param_data$type = "fixed"
	param_data$init = 0.0
	param_data$min = -10
	param_data$max = 10
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "anagenesis: exponent on extinction risk with area (modifies e)"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table

	
	# For a parameter at a branch point (cladogenesis)
	# (these are slightly different as I sometimes get errors running optim()/optimx() when 0 is 
	#  included as a possibility)
	param_name = "j"
	param_data = param_data_starter
	param_data$type = "fixed"
	param_data$init = 0.0
	param_data$min = 0 + minval_cladogenesis
	param_data$max = 3 - minval_cladogenesis
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "cladogenesis: relative per-event weight of jump dispersal"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table

	param_name = "ysv"
	param_data = param_data_starter
	param_data$type = "3-j"
	param_data$init = 3 - minval_cladogenesis
	param_data$min = 0 + minval_cladogenesis
	param_data$max = 3
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "cladogenesis: y+s+v"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table


	param_name = "ys"
	param_data = param_data_starter
	param_data$type = "ysv*2/3"
	param_data$init = 2 - minval_cladogenesis
	param_data$min = 0 + minval_cladogenesis
	param_data$max = 2
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "cladogenesis: y+s"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table

	param_name = "y"
	param_data = param_data_starter
	param_data$type = "ysv*1/3"
	param_data$init = 1
	param_data$min = 0 + minval_cladogenesis
	param_data$max = 1
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "cladogenesis: relative per-event weight of sympatry (range-copying)"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table

	param_name = "s"
	param_data = param_data_starter
	param_data$type = "ysv*1/3"
	param_data$init = 1
	param_data$min = 0 + minval_cladogenesis
	param_data$max = 1
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "cladogenesis: relative per-event weight of subset speciation"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table

	param_name = "v"
	param_data = param_data_starter
	param_data$type = "ysv*1/3"
	param_data$init = 1
	param_data$min = 0 + minval_cladogenesis
	param_data$max = 1
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "cladogenesis: relative per-event weight of vicariant speciation"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table



	# Maxent parameters (prob of smaller offspring | size)
	param_name = "mx01"
	param_data = param_data_starter
	param_data$type = "fixed"
	param_data$init = 0.0001
	param_data$min = 0.0001
	param_data$max = 0.9999
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "cladogenesis: controls range size of smaller daughter"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table

	param_name = "mx01j"
	param_data = param_data_starter
	param_data$type = "mx01"
	param_data$init = param_table[param_data$type,"init"]
	param_data$min = 0.0001
	param_data$max = 0.9999
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "cladogenesis: controls range size of smaller daughter"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table

	param_name = "mx01y"
	param_data = param_data_starter
	param_data$type = "mx01"
	param_data$init = param_table[param_data$type,"init"]
	param_data$min = 0.0001
	param_data$max = 0.9999
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "cladogenesis: controls range size of smaller daughter"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table

	param_name = "mx01s"
	param_data = param_data_starter
	param_data$type = "mx01"
	param_data$init = param_table[param_data$type,"init"]
	param_data$min = 0.0001
	param_data$max = 0.9999
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "cladogenesis: controls range size of smaller daughter"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table

	param_name = "mx01v"
	param_data = param_data_starter
	param_data$type = "mx01"
	param_data$init = param_table[param_data$type,"init"]
	param_data$min = 0.0001
	param_data$max = 0.9999
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "cladogenesis: controls range size of smaller daughter"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table

	# Prob distribution on range sizes at root
	param_name = "mx01r"
	param_data = param_data_starter
	param_data$type = "fixed"
	param_data$init = 0.5
	param_data$min = 0.0001
	param_data$max = 0.9999
	param_data$est = param_data$init
	param_data$note = "no"
	param_data$desc = "root: controls range size probabilities of root"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table


	# Mean frequency of truly sampling OTU of interest
	param_name = "mf"
	param_data = param_data_starter
	param_data$type = "fixed"
	param_data$init = 0.1
	param_data$min = 0.005
	param_data$max = 1-param_data$min
	param_data$est = param_data$init
	param_data$note = "yes"
	param_data$desc = "mean frequency of truly sampling OTU of interest"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table


	# Detection probability (per true sample of OTU of interest)
	param_name = "dp"
	param_data = param_data_starter
	param_data$type = "fixed"
	param_data$init = 1
	param_data$min = 0.005
	param_data$max = 1-param_data$min
	param_data$est = param_data$init
	param_data$note = "yes"
	param_data$desc = "detection probability per true sample of OTU of interest"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table

	# False detection probability (per true sample of non-OTU of interest)
	param_name = "fdp"
	param_data = param_data_starter
	param_data$type = "fixed"
	param_data$init = 0
	param_data$min = 0.005
	param_data$max = 1-param_data$min
	param_data$est = param_data$init
	param_data$note = "yes"
	param_data$desc = "false detection of OTU probability per true taphonomic control sample"
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table


	
	return(param_table)
	}




#######################################################
# get_clado_perEvent_weights
#######################################################
#' Get the per-event weights at cladogenesis
#' 
#' This function calculates the per-event weight as a proportion of some total
#' weight, e.g. default 1.  This is mostly useful for displaying the relative 
#' weight of different processes on a per-event basis.  This is NOT the same 
#' thing as the actual calculation of P(event|ancestral range), the probability
#' of a particular cladogenetic range-changing event (e.g., AB->A,B) conditional 
#' on an ancestral range (e.g., AB).
#'
#' For example, if model has per-event weights of j=0, s=1, y=1, v=1, the 
#' \code{get_clado_perEvent_weights()} result would be 0, 0.333, 0.333, 0.333.
#'
#' Background on how per-event weights are used to calculate conditional probabilities:
#' 
#' At a cladogenesis event, a large number of events are possible, conditional on a 
#' particular ancestral range.  (The possible events, and thus the weights and 
#' probabilities, will be different, depending on whether the ancestor is "A", "AB", 
#' "BCD", or whatever.) The simplest way to
#' compute these is just to assign some weight to each event, then sum all the events
#' and divide by the sum to get the probabilities.  More complex schemes can be 
#' imagined, but these are fairly pointless as they would all break down once
#' e.g. distance-dependence, user-specified connectivities, etc., are imposed.
#' 
#' In addition, one could imagine trying to assign total probabilities to each category of
#' event, but each row of the cladogenesis matrix may have a different count of the different
#' types of events (one row may have 1 y event and 2 j events; another row may have 4 j, 2 v, and
#' 2 s, and 0 y events; etc.).
#' 
#' One thing that IS meaningful is the per-event weight, i.e. the values that
#' the program is using for j, v, y, and s.  These ARE meaningful, as long as
#' they are forced to sum to some value (default 4).  This ensures that they are
#' identifiable (otherwise, j,v,y,s=1 and j,v,y,s=2 would be the same model).
#' 
#' @param params_table The \code{params_table} from a \code{BioGeoBEARS_model_object}, e.g. 
#' from \code{BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table}
#' @param sumval Default=1.
#' @param plotwhat Default "est", use "init" to get the initial starting values instead.
#' @return \code{wts} Return the per-event weights
#' @export
#' @seealso \code{\link[base]{rbind}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @examples
#' 
#' # default DEC+J model
#' BioGeoBEARS_run_object = define_BioGeoBEARS_run()
#' BioGeoBEARS_run_object$BioGeoBEARS_model_object@@params_table
#' params_table = BioGeoBEARS_run_object$BioGeoBEARS_model_object@@params_table
#' params_table
#' 
#' get_clado_perEvent_weights(params_table)
#' 
#' 
#' # DEC+J model
#' BioGeoBEARS_run_object = define_BioGeoBEARS_run()
#' BioGeoBEARS_run_object$BioGeoBEARS_model_object@@params_table["j","type"] = "free"
#' BioGeoBEARS_run_object$BioGeoBEARS_model_object@@params_table["j","init"] = 1
#' BioGeoBEARS_run_object$BioGeoBEARS_model_object@@params_table["j","est"] = 1
#' 
#' BioGeoBEARS_run_object$BioGeoBEARS_model_object@@params_table
#' 
#' BioGeoBEARS_run_object$BioGeoBEARS_model_object = 
#' calc_linked_params_BioGeoBEARS_model_object(
#' BioGeoBEARS_model_object=BioGeoBEARS_run_object$BioGeoBEARS_model_object, 
#' update_init=TRUE)
#' 
#' 
#' BioGeoBEARS_run_object$BioGeoBEARS_model_object@@params_table
#' params_table = BioGeoBEARS_run_object$BioGeoBEARS_model_object@@params_table
#' 
#' get_clado_perEvent_weights(params_table)
#' 
#'
get_clado_perEvent_weights <- function(params_table, sumval=1, plotwhat="est")
	{
	j = params_table["j", plotwhat]
	y = params_table["y", plotwhat]
	s = params_table["s", plotwhat]
	v = params_table["v", plotwhat]
	
	jysv_val = j+y+s+v
	
	j = j / jysv_val * sumval
	y = y / jysv_val * sumval
	s = s / jysv_val * sumval
	v = v / jysv_val * sumval
	
	vals = c(j, y, s, v)
	tmpmat = matrix(vals, nrow=1)
	wts = adf2(tmpmat)
	names(wts) = c("j", "y", "s", "v")
	
	return(wts)
	} # END get_clado_perEvent_weights()




#######################################################
# Take a BioGeoBEARS results_object, and get the LnL,
# numparams, and ML parameter estimates either as 
# a string or a table
#######################################################


#' Extract ML search results from a BioGeoBEARS_results_object
#'
#' Takes a \code{results_object}, and gets the LnL,
#' numparams, and ML parameter estimates either as 
#' a string or a table.
#'
#' @param results_object A \code{BioGeoBEARS_results_object} that is the result of a 
#' \code{\link{bears_optim_run}}. Typically named \code{res}, 
#' \code{resDEC}, \code{resDECj}, etc.)
#' @param returnwhat One of "table" (default), "string", or "param_names". 
#' @param addl_params The names of any additional (non-free) parameters that you 
#' want to extract, in addition to the free parameters.
#' @param paramsstr_digits The number of digits to return for the lnL and parameters, in string output.
#' @return Output is a \code{param_ests} table, a \code{paramstr} string, 
#' or \code{param_names} parameter names, depending on input \code{return_what}.
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
extract_params_from_BioGeoBEARS_results_object <- function(results_object, returnwhat="table", addl_params=NULL, paramsstr_digits=4)
	{
	defaults='
	results_object = res
	returnwhat="table"
	addl_params=NULL
	paramsstr_digits=4
	'
	
	# Error check on inputs
	if ( (returnwhat != "table")  && (returnwhat != "string") && (returnwhat != "param_names") )
		{
		errortxt = paste("\n\nERROR IN :\n\nInput 'returnwhat' must equal 'table' or 'string',\nbut you have returnwhat='", returnwhat, "'.\n\n", sep="")
		cat(errortxt)
		
		stop("STOPPING on error in extract_params_from_BioGeoBEARS_results_object().")
		}
	
	# Extract the log-likelihood
	# Fixing a problem where use of optim instead of optimx causes a problem
	# This will work with optim, optimx2012, or optimx2013
	LnL = get_LnL_from_optim_result(optimx_result=results_object$optim_result, use_optimx=results_object$inputs$use_optimx)

	params_table = results_object$outputs@params_table
	params_table

	#get_clado_perEvent_weights(params_table)


	# PLOT TITLE
	# What should be on the plot title
	
	params_free_TF = params_table$type == "free"
	params_free = (rownames(params_table))[params_free_TF]
	numparams = sum(params_free_TF)
	
	# Add additional user-specified parameters, if desired
	if (length(addl_params) > 0)
		{
		params_free = c(params_free, unlist(addl_params))
		params_free = unique(params_free)
		}
	
	# Write the string of FREE parameters
	paramstrs = rep("", length(params_free))
	param_names = NULL
	param_ests = NULL
	for (i in 1:length(params_free))
		{
		param_name = params_free[i]
		param_est = params_table[param_name,"est"]
		param_print = round(param_est, digits=paramsstr_digits)
		paramstrs[i] = paste(param_name, "=", param_print, "; ", sep="")
		
		# Store for output
		param_names = c(param_names, param_name)
		param_ests = c(param_ests, param_est)
		}

	# Return just the string
	if (returnwhat == "param_names")
		{
		return(param_names)
		}

	paramstrs = c(paramstrs, "LnL=", sprintf(fmt="%.2f", LnL) )
	paramstr = paste0(paramstrs, collapse="")
	
	# Return just the string
	if (returnwhat == "string")
		{
		return(paramstr)
		}
	
	# Store for output
	param_names = c("LnL", "numparams", param_names)
	param_ests = c(LnL, numparams, param_ests)

	# Handy summary outputs
	param_ests = matrix(data=param_ests, nrow=1)
	param_ests = adf2(param_ests)
	names(param_ests) = param_names
	
	param_ests = dfnums_to_numeric(param_ests)
	
	return(param_ests)
	} # END extract_params_from_BioGeoBEARS_results_object()
	








#######################################################
# BioGeoBEARS_model
#######################################################
#' An object of class BioGeoBEARS_model holding the model inputs
#'
#'@section Slots: 
#'  \describe{
#'    \item{\code{df}:}{Data.frame of class \code{"numeric"}, containing data from df}
#'  }
#'
#' @name BioGeoBEARS_model 
#' @rdname BioGeoBEARS_model-class
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @seealso \code{\link{define_tipranges_object}}, \code{\link{getareas_from_tipranges_object}},
#' \code{\link[cladoRcpp]{areas_list_to_states_list_old}}, \code{\link{areas_list_to_states_list_new}},
#' \code{\link{tipranges_to_tip_condlikes_of_data_on_each_state}}
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @examples
#' tipranges_object = define_tipranges_object()
#' tipranges_object
#'
setClass(Class="BioGeoBEARS_model", representation=representation(params_table="data.frame"))

#setClass(Class="BioGeoBEARS_model", representation=representation(params_table="data.frame"), contains="numeric")


#setClass(Class="BioGeoBEARS_model", representation=representation(params_table="data.frame"), prototype=BioGeoBEARS_model_defaults())


#setClass("tipranges", representation(df="data.frame"),
#    contains = "numeric"
#)


#######################################################
# define_BioGeoBEARS_model_object
#######################################################
#' Define a BioGeoBEARS_model class and object
#' 
#' Class \code{BioGeoBEARS_model} is an extension of the \code{\link{data.frame}} class.  It is used for holding
#' discrete geographic range data for the tips on a phylogeny. Geographic ranges are represented with bit 
#' encoding (0/1) indicating absence or presence in each possible area.
#' 
#' This is just a data.frame with:
#' rows = taxanames\cr
#' columns = area names\cr
#' cells = 0/1 representing empty/occupied\cr
#' 
#' @param minval_anagenesis Minimum value above zero for d, e, a, b parameters. Default is 1e-12.
#' @param minval_cladogenesis Minimum value above zero for j, v, etc. Default is 1e-5.
#' @param maxval Maximum value for d, e, a. Default is 5.
#' @return \code{BioGeoBEARS_model_object} The BioGeoBEARS_model object, of class \code{BioGeoBEARS_model}
#' @export
#' @seealso \code{\link[cladoRcpp]{areas_list_to_states_list_old}}, \code{\link{areas_list_to_states_list_new}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @examples
#' testval=1
#' BioGeoBEARS_model_object = define_BioGeoBEARS_model_object()
#' BioGeoBEARS_model_object
#' define_BioGeoBEARS_model_object()
define_BioGeoBEARS_model_object <- function(minval_anagenesis=1e-12, minval_cladogenesis=1e-5, maxval=5)
	{
	# Define the BioGeoBEARS_model class;

	#BioGeoBEARS_model_object = new("BioGeoBEARS_model", df=tmpdf)
	BioGeoBEARS_model_object = new("BioGeoBEARS_model", params_table=BioGeoBEARS_model_defaults(minval_anagenesis, minval_cladogenesis, maxval))
	BioGeoBEARS_model_object
	
	# you can get the dataframe with
	# BioGeoBEARS_model_object@df

	return(BioGeoBEARS_model_object)
	}



#######################################################
# BioGeoBEARS_model_object_to_init_params
#######################################################
#' Produce initial parameters from a BioGeoBEARS model object
#' 
#' This function returns the initial values of the (free) parameters from
#' a \code{BioGeoBEARS_model_object}.
#' 
#' @param BioGeoBEARS_model_object The BioGeoBEARS_model object, of class 
#' \code{BioGeoBEARS_model}. E.g. from \code{BioGeoBEARS_run_object$BioGeoBEARS_model_object}
#' @return \code{params} parameter vector
#' @export
#' @seealso \code{\link[BioGeoBEARS]{define_BioGeoBEARS_model_object}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @examples
#' test=1
#' BioGeoBEARS_run_object = define_BioGeoBEARS_run()
#' BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
#' 
#' BioGeoBEARS_model_object_to_init_params(BioGeoBEARS_model_object)
#' 
#' # [1] 0.01 0.01
BioGeoBEARS_model_object_to_init_params <- function(BioGeoBEARS_model_object)
	{
	ex='
	BioGeoBEARS_run_object = define_BioGeoBEARS_run()
	BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
	BioGeoBEARS_model_object_to_init_params(BioGeoBEARS_model_object)
	'
	
	free_params_TF = BioGeoBEARS_model_object@params_table$type == "free"
	num_free_params = sum(free_params_TF, na.rm=TRUE)
	
	params = rep(NA, num_free_params)
	
	params = BioGeoBEARS_model_object@params_table$init[free_params_TF]
	
	return(params)	
	}



#######################################################
# BioGeoBEARS_model_object_to_est_params
#######################################################
#' Extract estimated parameters from a BioGeoBEARS model object
#' 
#' This function returns the estimated values of the (free) parameters from
#' a \code{BioGeoBEARS_model_object}.
#' 
#' @param BioGeoBEARS_model_object The BioGeoBEARS_model object, of class 
#' \code{BioGeoBEARS_model}. E.g. from \code{BioGeoBEARS_run_object$BioGeoBEARS_model_object}
#' @return \code{params} parameter vector
#' @export
#' @seealso \code{\link[BioGeoBEARS]{define_BioGeoBEARS_model_object}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @examples
#' test=1
#' BioGeoBEARS_run_object = define_BioGeoBEARS_run()
#' BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
#' 
#' BioGeoBEARS_model_object_to_est_params(BioGeoBEARS_model_object)
#' 
#' # [1] 0.01 0.01
BioGeoBEARS_model_object_to_est_params <- function(BioGeoBEARS_model_object)
	{
	ex='
	BioGeoBEARS_run_object = define_BioGeoBEARS_run()
	BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
	BioGeoBEARS_model_object_to_est_params(BioGeoBEARS_model_object)
	'
	free_params_TF = BioGeoBEARS_model_object@params_table$type == "free"
	num_free_params = sum(free_params_TF, na.rm=TRUE)
	
	params = rep(NA, num_free_params)
	
	params = BioGeoBEARS_model_object@params_table$est[free_params_TF]
	
	return(params)	
	}




#######################################################
# BioGeoBEARS_model_object_to_params_lower
#######################################################
#' Extract the lower limit on the parameters from a BioGeoBEARS model object
#' 
#' This function returns the lower limits of the (free) parameters from
#' a \code{BioGeoBEARS_model_object}.
#' 
#' @param BioGeoBEARS_model_object The BioGeoBEARS_model object, of class 
#' \code{BioGeoBEARS_model}. E.g. from \code{BioGeoBEARS_run_object$BioGeoBEARS_model_object}
#' @return \code{params} parameter vector
#' @export
#' @seealso \code{\link[BioGeoBEARS]{define_BioGeoBEARS_model_object}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @examples
#' test=1
#' BioGeoBEARS_run_object = define_BioGeoBEARS_run()
#' BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
#' 
#' BioGeoBEARS_model_object_to_est_params(BioGeoBEARS_model_object)
#' 
#' # [1] 1e-12 1e-12
BioGeoBEARS_model_object_to_params_lower <- function(BioGeoBEARS_model_object)
	{
	ex='
	BioGeoBEARS_run_object = define_BioGeoBEARS_run()
	BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
	BioGeoBEARS_model_object_to_params_lower(BioGeoBEARS_model_object)
	'
	free_params_TF = BioGeoBEARS_model_object@params_table$type == "free"
	num_free_params = sum(free_params_TF, na.rm=TRUE)
	
	params = rep(NA, num_free_params)
	
	params = BioGeoBEARS_model_object@params_table$min[free_params_TF]
	
	return(params)	
	}


#######################################################
# BioGeoBEARS_model_object_to_params_upper
#######################################################
#' Extract the upper limit on the parameters from a BioGeoBEARS model object
#' 
#' This function returns the upper limits of the (free) parameters from
#' a \code{BioGeoBEARS_model_object}.
#' 
#' @param BioGeoBEARS_model_object The BioGeoBEARS_model object, of class 
#' \code{BioGeoBEARS_model}. E.g. from \code{BioGeoBEARS_run_object$BioGeoBEARS_model_object}
#' @return \code{params} parameter vector
#' @export
#' @seealso \code{\link[BioGeoBEARS]{define_BioGeoBEARS_model_object}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @examples
#' test=1
#' BioGeoBEARS_run_object = define_BioGeoBEARS_run()
#' BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
#' 
#' BioGeoBEARS_model_object_to_params_upper(BioGeoBEARS_model_object)
#' 
#' # [1] 5 5
BioGeoBEARS_model_object_to_params_upper <- function(BioGeoBEARS_model_object)
	{
	ex='
	BioGeoBEARS_run_object = define_BioGeoBEARS_run()
	BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
	BioGeoBEARS_model_object_to_params_upper(BioGeoBEARS_model_object)
	'
	free_params_TF = BioGeoBEARS_model_object@params_table$type == "free"
	num_free_params = sum(free_params_TF, na.rm=TRUE)
	
	params = rep(NA, num_free_params)
	
	params = BioGeoBEARS_model_object@params_table$max[free_params_TF]
	
	return(params)	
	}


#######################################################
# params_into_BioGeoBEARS_model_object
#######################################################
#' Feed modified parameters back into a BioGeoBEARS model object
#' 
#' This function takes a list of parameter values (typically, the 
#' current iteration of parameters in a Maximum Likelihood search)
#' and inputs them, in order, into the free parameters specified
#' in the \code{BioGeoBEARS_model_object@params_table}.
#' 
#' @param BioGeoBEARS_model_object The BioGeoBEARS_model object, of class 
#' \code{BioGeoBEARS_model}. E.g. from \code{BioGeoBEARS_run_object$BioGeoBEARS_model_object}
#' @param params parameter vector
#' @param initTF If \code{TRUE} (default), update the "init" (initial values) column with these parameters.
#' @param estTF If \code{TRUE} (default), update the "est" (estimates) column with these parameters.
#' @return \code{BioGeoBEARS_model_object} The BioGeoBEARS_model object, of class \code{BioGeoBEARS_model}
#' @export
#' @seealso \code{\link[BioGeoBEARS]{define_BioGeoBEARS_model_object}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @examples
#' test=1
#' BioGeoBEARS_run_object = define_BioGeoBEARS_run()
#' BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
#' BioGeoBEARS_model_object
#' BioGeoBEARS_model_object2 = params_into_BioGeoBEARS_model_object(BioGeoBEARS_model_object, params=c(0.34,0.21))
#' BioGeoBEARS_model_object2
#' 
params_into_BioGeoBEARS_model_object <- function(BioGeoBEARS_model_object, params, initTF=TRUE, estTF=TRUE)
	{
	ex='
	BioGeoBEARS_run_object = define_BioGeoBEARS_run()
	BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
	BioGeoBEARS_model_object
	BioGeoBEARS_model_object2 = params_into_BioGeoBEARS_model_object(BioGeoBEARS_model_object, params=c(0.34,0.21))
	BioGeoBEARS_model_object2
	'

	free_params_TF = BioGeoBEARS_model_object@params_table$type == "free"
	num_free_params = sum(free_params_TF, na.rm=TRUE)
	
	if (initTF == TRUE)
		{
		BioGeoBEARS_model_object@params_table$init[free_params_TF] = params	
		}
		
	if (estTF == TRUE)
		{
		BioGeoBEARS_model_object@params_table$est[free_params_TF] = params
		}
	
	if ( (initTF == FALSE) && (estTF == FALSE))
		{
		warning("WARNING in params_into_BioGeoBEARS_model_object -- neither column 'init' nor 'est' were updated, check inputs 'initTF' and 'estTF'.")
		}
	
	return(BioGeoBEARS_model_object)
	}

#######################################################
# merge_words_nonwords
#######################################################
#' Merge lists of words and nonwords (numbers) that may be of different length
#' 
#' Utility function for printing BioGeoBEARS parameters to screen or
#' text files. The list of \code{nonwords} might be 1 longer than
#' the list of \code{words}. (Longer produces an error.)
#' 
#' @param words A list of words
#' @param nonwords A list of nonwords
#' @return \code{sentence} A text string.
#' @export
#' @seealso \code{\link[BioGeoBEARS]{define_BioGeoBEARS_model_object}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @examples
#' test=1
#' words = c("A", "B", "C")
#' nonwords = c(1,2,3,4)
#' merge_words_nonwords(words, nonwords)
#' 
#' \dontrun{
#' # Produces an error
#' words = c("A", "B", "C")
#' nonwords = c(1,2,3,4, 5)
#' merge_words_nonwords(words, nonwords)
#' }
merge_words_nonwords <- function(words, nonwords)
	{
	ex='
	words = c("A", "B", "C")
	nonwords = c(1,2,3,4)
	merge_words_nonwords(words, nonwords)
	
	# Produces an error
	words = c("A", "B", "C")
	nonwords = c(1,2,3,4, 5)
	merge_words_nonwords(words, nonwords)
	'
	
	if (length(nonwords) == ((length(words) + 1)))
		{
		words = c(words, "")
		}
	
	wordsmat = cbind(nonwords, words)
	paste1 = apply(X=wordsmat, MARGIN=1, FUN=paste, sep="", collapse="")
	paste1
	
	sentence = paste(paste1, sep="", collapse="")
	return(sentence)
	}








# Get the parameter results from an optim search

#' Get the parameter results from an optim search
#'
#' Extracts the parameters from the results object 
#' produced by 
#' the \code{\link[stats]{optim}} function for 
#' Maximum Likelihood optimization.
#'
#' This is necessary, because, annoyingly, \code{optim}, \code{optimx} from pre-2013,
#' \code{optimx} from post-2013, and \code{GenSA}, all return their
#' parameter estimates in slightly different objects.
#'
#' @param optim_result The result produced by \code{\link[stats]{optim}}()
#' @return optim_param_results A vector of the ML parameters
#' @export
#' @seealso \code{\link[stats]{optim}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
get_params_from_optim <- function(optim_result)
	{
	optim_param_results = as.numeric(optim_result$par)
	
	return(optim_param_results)
	}


# Get the parameter results from an optim search
#' Get the parameter results from a GenSA search
#'
#' Extracts the parameters from the results object 
#' produced by 
#' the \code{\link[GenSA]{GenSA}} function for 
#' Maximum Likelihood optimization (GenSA = Generalized
#' Simulated Annealing, which seems to do a more thorough
#' search when the number of parameters increases to 5+).
#'
#' This is necessary, because, annoyingly, \code{optim}, \code{optimx} from pre-2013,
#' \code{optimx} from post-2013, and \code{GenSA}, all return their
#' parameter estimates in slightly different objects.
#'
#' @param GenSA_result The result produced by \code{\link[stats]{optim}}()
#' @return GenSA_param_results A vector of the ML parameters
#' @export
#' @seealso \code{\link[GenSA]{GenSA}} 
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
get_params_from_GenSA <- function(GenSA_result)
	{
	GenSA_param_results = as.numeric(GenSA_result$par)
	
	return(GenSA_param_results)
	}



# Get the parameter results from an optimx 2012 search


#' Get the parameter results from an optimx 2012 search
#'
#' Extracts the parameters from the results object 
#' produced by 
#' the \code{\link[optimx]{optimx}} function for 
#' Maximum Likelihood optimization (2012 and earlier).
#'
#' This is necessary, because, annoyingly, \code{optim}, \code{optimx} from pre-2013,
#' \code{optimx} from post-2013, and \code{GenSA}, all return their
#' parameter estimates in slightly different objects.
#'
#' @param optimx_result The result produced by \code{\link[optimx]{optimx}}()
#' @return optimx_param_results A vector of the ML parameters
#' @export
#' @seealso \code{\link[optimx]{optimx}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
get_params_from_optimx2012 <- function(optimx_result)
	{
	optimx_param_results = as.numeric(optimx_result$par[[1]])
	
	return(optimx_param_results)
	}


#' Get the parameter results from an optimx 2013 search
#'
#' Extracts the parameters from the results object 
#' produced by 
#' the \code{\link[optimx]{optimx}} function for 
#' Maximum Likelihood optimization (2013 and later).
#'
#' This is necessary because, annoyingly, \code{optim}, \code{optimx} from pre-2013,
#' \code{optimx} from post-2013, and \code{GenSA}, all return their
#' parameter estimates in slightly different objects.
#'
#' @param optimx_result The result produced by \code{\link[optimx]{optimx}}()
#' @return optimx_param_results A vector of the ML parameters
#' @export
#' @seealso \code{\link[optimx]{optimx}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
get_params_from_optimx2013 <- function(optimx_result)
	{
	nparams = nparams_from_optimx2013(optimx_result)
	optimx_param_results = as.numeric(optimx_result[1:nparams])
	
	return(optimx_param_results)
	}


#' Get the number of free parameters in an ML search result
#'
#' Utility function for \code{get_params_from_optimx2013}, which figures
#' out the number of free parameters in the results object.
#'
#' This is necessary because, annoyingly, \code{optim}, \code{optimx} from pre-2013,
#' \code{optimx} from post-2013, and \code{GenSA}, all return their
#' parameter estimates in slightly different objects.
#'
#' @param optimx_result The result produced by \code{\link[optimx]{optimx}}()
#' @return nparams The number of free parameters
#' @export
#' @seealso \code{\link[optimx]{optimx}}, \code{\link{get_params_from_optimx2013}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
nparams_from_optimx2013 <- function(optimx_result)
	{
	# Find the colnums and names of the free parameters
	value_TF = names(optimx_result) == "value"
	tmp_colnums = 1:length(value_TF)
	
	# The log-likelihood (LnL) is in this column
	LnL_colnum = tmp_colnums[value_TF]
	
	# This is the number of parameters
	nparams = LnL_colnum - 1
	return(nparams)
	}


#######################################################
# Put optim or optimx result into a BioGeoBEARS_model_object
#######################################################

#' Get the parameter estimate from a optim, optimx, or GenSA result
#'
#' This function runs the appropriate \code{get_} function to 
#' extract the parameter inferences from ML searches with
#' \code{optim}, \code{optimx}, or \code{GenSA}.
#'
#' @param optimx_result The result object returned by \code{optim}, 
#' \code{optimx}, or \code{GenSA}
#' @param use_optimx If \code{TRUE} (default) or \code{"optimx"}, use \code{\link[optimx]{optimx}} rather than \code{\link[stats]{optim}}. 
#' If \code{FALSE} or \code{"optim"}, use \code{\link[stats]{optim}}. If \code{"GenSA"}, use \code{\link[GenSA]{GenSA}}.
#' @return params_from_ML A vector of parameters from the ML search
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
get_params_from_optim_or_optimx_result <- function(optimx_result, use_optimx=TRUE)
	{
	# Get the params from ML search to	
	# set the dispersal and extinction rate (and j, etc)
	# in the BioGeoBEARS_model_object
	# (includes updating linked params)
	
	if ( (use_optimx == FALSE) || (use_optimx == "optim") )
		{
		params_from_ML = get_params_from_optim(optimx_result)
		} 
		
	if ( (use_optimx == TRUE) || (use_optimx == "optimx") )
		{
		# USING OPTIMX()
		# optimx has 2012 and 2013 version
		# optimx 2013 has different output; params are in p1, p2, etc.
		
		# Check for optimx 2012 or 2013, and extract parameters accordingly
		if (packageVersion("optimx") < 2013)
			{
			# optimx 2012
			params_from_ML = get_params_from_optimx2012(optimx_result)
			} else {
			# optimx 2013
			params_from_ML = get_params_from_optimx2013(optimx_result)
			} # END optimx2012 versus optimx2013
		} # END optim() versus optimx()

	if (use_optimx == "GenSA")
		{
		params_from_ML = get_params_from_GenSA(optimx_result)
		} 

	
	return(params_from_ML)
	}




#######################################################
# Put optim or optimx result into a BioGeoBEARS_model_object
#######################################################

#' Put the parameter estimate from a optim, optimx, or GenSA result into a BioGeoBEARS_model_object
#'
#' Updates the current estimates (column \code{est}) in a \code{BioGeoBEARS_model_object@params_table},
#' with the parameter inferences from ML searches with
#' \code{optim}, \code{optimx}, or \code{GenSA}.
#'
#' @param BioGeoBEARS_model_object The BioGeoBEARS_model object, of class 
#' \code{BioGeoBEARS_model}. E.g. from \code{BioGeoBEARS_run_object$BioGeoBEARS_model_object}
#' @param optimx_result The result object returned by \code{optim}, 
#' \code{optimx}, or \code{GenSA}
#' @param use_optimx If \code{TRUE} (default) or \code{"optimx"}, use \code{\link[optimx]{optimx}} rather than \code{\link[stats]{optim}}. 
#' If \code{FALSE} or \code{"optim"}, use \code{\link[stats]{optim}}. If \code{"GenSA"}, use \code{\link[GenSA]{GenSA}}.
#' @return BioGeoBEARS_model_object A BioGeoBEARS_model object, of class 
#' \code{BioGeoBEARS_model}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
update_BioGeoBEARS_model_object_w_optimx_result <- function(BioGeoBEARS_model_object, optimx_result, use_optimx=TRUE)
	{
	# Get the params from ML search to	
	# set the dispersal and extinction rate (and j, etc)
	# in the BioGeoBEARS_model_object
	# (includes updating linked params)
	
	# Get ML params from each of the types of optim/optimx result:
	params_from_ML = get_params_from_optim_or_optimx_result(optimx_result, use_optimx=use_optimx)
	
	
	# Update the model object and linked parameters
	BioGeoBEARS_model_object = params_into_BioGeoBEARS_model_object(BioGeoBEARS_model_object, params=params_from_ML)
	BioGeoBEARS_model_object = calc_linked_params_BioGeoBEARS_model_object(BioGeoBEARS_model_object)
	
	return(BioGeoBEARS_model_object)
	}






#' Put the initial parameters into the optim_results format
#'
#' If you wish to calculate the data likelihood under a certain set of pre-specified
#' parameters, those initial parameters (in the init column of a \code{BioGeoBEARS_model_object@params_table},
#' need to be input into a results object from \code{optim}, \code{optimx}, or \code{GenSA}. This function
#' constructs the necessary object. 
#'
#' The function is called when, in bears_optim_run, skip_optim_option=="return_all". 
#' The function converts the 'init' parameters into the optim_results format, for 
#' further use in \code{bears_optim_run}.
#' 
#' @param BioGeoBEARS_model_object The BioGeoBEARS_model object, of class 
#' \code{BioGeoBEARS_model}. E.g. from \code{BioGeoBEARS_run_object$BioGeoBEARS_model_object}
#' @param total_loglikelihood The total_loglikelihood (one of the parts of an ML search result)
#' @param use_optimx If \code{TRUE} (default) or \code{"optimx"}, use \code{\link[optimx]{optimx}} rather than \code{\link[stats]{optim}}. 
#' If \code{FALSE} or \code{"optim"}, use \code{\link[stats]{optim}}. If \code{"GenSA"}, use \code{\link[GenSA]{GenSA}}.
#' @return BioGeoBEARS_model_object A BioGeoBEARS_model object, of class 
#' \code{BioGeoBEARS_model}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
put_params_into_optim_or_optimx_result <- function(BioGeoBEARS_model_object, total_loglikelihood, use_optimx)
	{

	# Get the free parameters from the model object
	free_params_TF = BioGeoBEARS_model_object@params_table$type == "free"
	num_free_params = sum(free_params_TF, na.rm=TRUE)
	params = BioGeoBEARS_model_object@params_table$init[free_params_TF]
	
	
	# Get the params from ML search to	
	# set the dispersal and extinction rate (and j, etc)
	# in the BioGeoBEARS_model_object
	# (includes updating linked params)
	
	# BASIC OPTIM
	if ( (use_optimx == FALSE) || (use_optimx == "optim") )
		{
		optim_result = list()
		optim_result$par = params
		optim_result$value = total_loglikelihood
		optim_result$counts = rep("user", 2)
		names(optim_result$counts) = c("function", "gradient")
		optim_result$convergence = total_loglikelihood
		optim_result$message = "These parameter values were manually input by the user (so, not directly output by optim() )."
		} # END if ( (use_optimx == FALSE) || (use_optimx == "optim") )
		

	# USING OPTIMX
	if ( (use_optimx == TRUE) || (use_optimx == "optimx") )
		{
		# optimx has 2012 and 2013 version
		# optimx 2013 has different output; params are in p1, p2, etc.
		
		# Check for optimx 2012 or 2013, and extract parameters accordingly
		if (packageVersion("optimx") < 2013)
			{
			###########################
			# optimx 2012
			###########################
			#params_from_ML = get_params_from_optimx2012(optimx_result)
			stoptxt = paste0("STOP ERROR in put_params_into_optim_or_optimx_result(): The function 'put_params_into_optim_or_optimx_result' will not work for versions of 'optimx' from 2012 and before. Upgrade to a newer version of the 'optimx' package to use.")
			cat("\n\n")
			cat(stoptxt)
			cat("\n\n")
			stop(stoptxt)
			} else {
			###########################
			# optimx 2013+
			###########################
			#params_from_ML = get_params_from_optim(optimx_result)
			num_elements = num_free_params + 8
			param_names = paste0("p", 1:num_free_params)
		
			optim_result = as.data.frame(matrix(data=rep(0, times=num_elements), nrow=1, ncol=num_elements), stringsAsFactors=FALSE)
			names(optim_result) = c(param_names, "value", "fevals", "gevals", "niter", "convcode", "kkt1", "kkt2", "xtimes")
			optim_result$value = total_loglikelihood
			optim_result$fevals = "user"
			optim_result$gevals = "user"
			optim_result$niter = "user"
			optim_result$convcode = "user"
			optim_result$kkt1 = "user"
			optim_result$kkt2 = "user"
			optim_result$xtimes = "user"
			row.names(optim_result) == "bobyqa"
			class(optim_result) = c("optimx", "data.frame")
		
			# Input the free parameters that you estimated
			for (i in 1:length(param_names))
				{
				cmdtxt = paste0("optim_result$", param_names[i], " = ", params[i])
				eval(parse(text=cmdtxt))
				}
			return(optim_result)
			#params_from_ML = get_params_from_optimx2013(optimx_result)
			} # END optimx2012 versus optimx2013
		} # END optim() versus optimx()


	# USING GenSA
	if (use_optimx == "GenSA")
		{
		#params_from_ML = get_params_from_GenSA(optimx_result)
		optim_result = list()
		optim_result$value = -1 * total_loglikelihood
		optim_result$par = params
		optim_result$trace.mat = "These parameter values were manually input by the user; so, no $trace.mat available from GenSA."
		optim_result$counts = "user"
		} # END if (use_optimx == "GenSA")

	return(optim_result)
	} # END put_params_into_optim_or_optimx_result <- function(res=NA, BioGeoBEARS_model_object=NA, use_optimx=NA, total_loglikelihood=NA)







# 
#' Get the log-likelihood from an optim search
#'
#' Extracts the log-likelihood value from an ML search with \code{\link[stats]{optim}}.
#'
#' @param optimx_result The result from an \code{\link[stats]{optim}} search
#' @return LnL The log-likelihood
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
get_LnL_from_optim <- function(optimx_result)
	{
	LnL = as.numeric(optimx_result$value)
	return(LnL)
	}

# Get the log-likelihood from a GenSA search
#' Get the log-likelihood from a GenSA search
#'
#' Extracts the log-likelihood value from an ML search with \code{\link[GenSA]{GenSA}}.
#'
#' @param optim_result The result from an \code{\link[GenSA]{GenSA}} search
#' @return LnL The log-likelihood
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
get_LnL_from_GenSA <- function(optimx_result)
	{
	LnL = -1 * as.numeric(optimx_result$value)
	return(LnL)
	}

# Get the log-likelihood result from an optimx pre-2012 search

#' Get the log-likelihood from an optimx pre-2012 search
#'
#' Extracts the log-likelihood value from an ML search with the pre-2013 \code{\link[optimx]{optimx}}.
#'
#' @param optimx_result The result from an \code{\link[optimx]{optimx}} search
#' @return LnL The log-likelihood
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
get_LnL_from_optimx2012 <- function(optimx_result)
	{
	LnL = as.numeric(optimx_result$fvalues)
	return(LnL)
	}

# Get the parameter results from an optimx 2013 search
#' Get the log-likelihood from an optimx post-2013 search
#'
#' Extracts the log-likelihood value from an ML search with the post-2013 \code{\link[optimx]{optimx}}.
#'
#' @param optimx_result The result from an \code{\link[optimx]{optimx}} search
#' @return LnL The log-likelihood
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
get_LnL_from_optimx2013 <- function(optimx_result)
	{
	LnL = as.numeric(optimx_result$value)
	return(LnL)
	}


#' Get the log-likelihood from a optim, optimx, or GenSA result
#'
#' This function runs the appropriate \code{get_LnL_} function to 
#' extract the log-likelihood from ML searches with
#' \code{optim}, \code{optimx}, or \code{GenSA}.
#'
#' @param optimx_result The result object returned by \code{optim}, 
#' \code{optimx}, or \code{GenSA}
#' @param use_optimx If \code{TRUE} (default) or \code{"optimx"}, use \code{\link[optimx]{optimx}} rather than \code{\link[stats]{optim}}. 
#' If \code{FALSE} or \code{"optim"}, use \code{\link[stats]{optim}}. If \code{"GenSA"}, use \code{\link[GenSA]{GenSA}}.
#' @return LnL The log-likelihood
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
get_LnL_from_optim_result <- function(optimx_result, use_optimx=TRUE)
	{
	# Optimx
	if ( (use_optimx == TRUE) || (use_optimx == "optimx") )
		{
		# Using optimx() results
		if (packageVersion("optimx") < 2013)
			{
			# optimx 2012
			LnL = get_LnL_from_optimx2012(optimx_result)
			} else {
			# optimx 2013
			LnL = get_LnL_from_optimx2013(optimx_result)
			} # end optimx 2012 vs. 2013
		}
	
	# Optim
	if ( (use_optimx == FALSE) || (use_optimx == "optim") )
		{
		# Using optim() results
		LnL = get_LnL_from_optim(optimx_result)
		} # end optim vs. optimx

	if (use_optimx == "GenSA")
		{
		# Using optim() results
		LnL = get_LnL_from_GenSA(optimx_result)
		} # end optim vs. optimx

	return(LnL)
	}

#' Get the log-likelihood from a BioGeoBEARS_results_object
#'
#' This function takes a \code{BioGeoBEARS_results_object}, and runs 
#' the appropriate \code{get_LnL_} function to 
#' extract the log-likelihood from the embededded ML searche result with
#' \code{optim}, \code{optimx}, or \code{GenSA}.
#'
#' @param res A \code{BioGeoBEARS_results_object} that is the result of a 
#' \code{\link{bears_optim_run}}. (typically \code{res}, 
#' \code{resDEC}, \code{resDECj}, etc.)
#' @param use_optimx If \code{TRUE} (default) or \code{"optimx"}, use \code{\link[optimx]{optimx}} rather than \code{\link[stats]{optim}}. 
#' If \code{FALSE} or \code{"optim"}, use \code{\link[stats]{optim}}. If \code{"GenSA"}, use \code{\link[GenSA]{GenSA}}.
#' @return LnL The log-likelihood
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
get_LnL_from_BioGeoBEARS_results_object <- function(res)
	{
	BioGeoBEARS_run_object = res$inputs
	optimx_result = res$optim_result
	LnL = get_LnL_from_optim_result(optimx_result, use_optimx=BioGeoBEARS_run_object$use_optimx)
	return(LnL)
	}





#######################################################
# calc_linked_params_BioGeoBEARS_model_object
#######################################################
#' Update parameters that are deterministic functions of free parameters
#' 
#' This function updates the linked parameters (which are listed as neither "fixed" nor
#' "free" in \code{params_table$type}; i.e., they are equations which are calculated from #' the fixed and free parameters,
#' which should have already been set by other functions).
#'
#' \code{params_table$type} is typically stored in: \code{BioGeoBEARS_run_object$BioGeoBEARS_model_object@@params_table}.
#' 
#' @param BioGeoBEARS_model_object The BioGeoBEARS_model object, of class 
#' \code{BioGeoBEARS_model}. E.g. from \code{BioGeoBEARS_run_object$BioGeoBEARS_model_object}
#' @param update_init If \code{TRUE}, put the estimates into the initial values in the \code{params_table}.
#' Default: \code{FALSE}.
#' @return \code{BioGeoBEARS_model_object} Updated version of the BioGeoBEARS_model
#' object, of class \code{BioGeoBEARS_model}.
#' @export
#' @seealso \code{\link[BioGeoBEARS]{define_BioGeoBEARS_model_object}} \code{\link[BioGeoBEARS]{define_BioGeoBEARS_run}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @examples
#' # Define a BioGeoBEARS run object
#' BioGeoBEARS_run_object = define_BioGeoBEARS_run()
#' BioGeoBEARS_run_object$BioGeoBEARS_model_object@@params_table
#' 
#' # Set 'j' to be free, i.e. as in a DEC+J model (adding jump dispersal
#' # to the LAGRANGE DEC model)
#' BioGeoBEARS_run_object$BioGeoBEARS_model_object@@params_table["j","type"] = "free"
#' BioGeoBEARS_run_object$BioGeoBEARS_model_object@@params_table["j","init"] = 0.25
#' BioGeoBEARS_run_object$BioGeoBEARS_model_object@@params_table["j","est"] = 0.25
#' 
#' # Display result
#' BioGeoBEARS_run_object$BioGeoBEARS_model_object@@params_table
#' 
#' # Update the other parameters
#' BioGeoBEARS_run_object$BioGeoBEARS_model_object = 
#' calc_linked_params_BioGeoBEARS_model_object(
#' BioGeoBEARS_model_object=BioGeoBEARS_run_object$BioGeoBEARS_model_object)
#' 
#' # Display result
#' BioGeoBEARS_run_object$BioGeoBEARS_model_object@@params_table
#'
calc_linked_params_BioGeoBEARS_model_object <- function(BioGeoBEARS_model_object, update_init=FALSE)
	{
	# identify free parameters
	linked_TF1 = BioGeoBEARS_model_object@params_table$type != "free"

	# identify fixed parameters
	linked_TF2 = BioGeoBEARS_model_object@params_table$type != "fixed"

	# the other parameters (linked parameters) are calculated from the fixed and free parameters
	# which should have already been set
	linked_TF = (linked_TF1 + linked_TF2) == 2
	
	# If none are linked, skip
	if (sum(linked_TF) == 0)
		{
		return(BioGeoBEARS_model_object)
		}
	
	linked_types = unique(BioGeoBEARS_model_object@params_table$type[linked_TF])
	
	# Get the list of variables
	mstr = paste(rownames(BioGeoBEARS_model_object@params_table), sep="", collapse="|")
	
	for (ll in 1:length(linked_types))
		{
		# Match this type
		linked_type = linked_types[ll]
		this_type_TF = BioGeoBEARS_model_object@params_table$type == linked_type
		
		# Get the formula
		tmp_formula = linked_type
		
		# Extract the words
		words = regmatches(x=tmp_formula, m=gregexpr(mstr, tmp_formula), invert=FALSE)[[1]]
		nonwords = regmatches(x=tmp_formula, m=gregexpr(mstr, tmp_formula), invert=TRUE)[[1]]
		
		
		# Get the values for each word
		tmpfun <- function(word, params_table)
			{
			return(params_table[word, "est"])
			}
		vals = sapply(X=words, FUN=tmpfun, params_table=BioGeoBEARS_model_object@params_table)

		# Make the formula, with numbers put in
		new_formula = merge_words_nonwords(vals, nonwords)
		new_formula
		
		# Calculate the new value
		newval = NA
		cmdstr = paste("newval = ", new_formula, sep="")
		eval(parse(text=cmdstr))
		
		# Input it into the matched cells
		BioGeoBEARS_model_object@params_table$est[this_type_TF] = newval
		
		# Also update the initial values, if desired
		if (update_init == TRUE)
			{
			BioGeoBEARS_model_object@params_table$init = BioGeoBEARS_model_object@params_table$est
			BioGeoBEARS_model_object@params_table$init[this_type_TF] = newval
			}
		
		}
	return(BioGeoBEARS_model_object)
	}




# Define a MAXIMUM LIKELIHOOD search

#######################################################
# define_BioGeoBEARS_run
#######################################################
#' Define a maximum likelihood search, perhaps stratified
#' 
#' Set up the inputs object for an ML search.  See parameter descriptions for defaults.
#' 
#' @param abbr Text abbreviation of run, e.g. "default"
#' @param description Text description of run, e.g. "defaults"
#' @param BioGeoBEARS_model_object Default is \code{define_BioGeoBEARS_model_object()}
#' @param trfn The filename of the phylogenetic tree, in NEWICK format (\url{http://evolution.genetics.washington.edu/phylip/newicktree.html}).  
#' Tipnames should match the names in geogfn.  See \code{\link[ape]{read.tree}} in APE for reading in phylogenetic trees. Default is \code{NULL}, which sets the file to the default "Psychotria_5.2.newick" in the BioGeoBEARS extension data directory, find with \code{system.file("extdata", package="BioGeoBEARS")}.
#' @param geogfn A PHYLIP-style file with geographic range data (see \code{\link{getranges_from_LagrangePHYLIP}}) for each tipname. This is the same format
#' used by C++ LAGRANGE (\cite{SmithRee2010_CPPversion}). The default is \code{NULL}, which sets the file to the default "Psychotria_geog.data" in the BioGeoBEARS extension data directory, find with \code{system.file("extdata", package="BioGeoBEARS")}
#' @param timesfn Filename for the stratified times.
#' @param distsfn Filename for the changing distances.
#' @param envdistsfn Filename for the changing environmental distances.
#' @param dispersal_multipliers_fn Filename for the changing hard-coded dispersal multipliers
#' @param area_of_areas_fn Filename for the area of each area
#' @param areas_allowed_fn Filename for the allowed areas for single-species ranges (for turning off areas back in time, e.g. as you go back in time there are fewer Hawaiian islands).
#' @param areas_adjacency_fn Filename indicating area adjacency, *in the sense of allowing/disallowing certain ranges in the state space* (which is a different thing than modifying the dispersal rates between areas).
#' @param detects_fn Filename for the counts of detections of OTUs of interest. See \code{\link{calc_obs_like}}.
#' @param controls_fn Filename for the counts of taphonomic controls (which INCLUDE the OTUs of interest). See \code{\link{calc_obs_like}}.
#' @param max_range_size The maximum rangesize, in number of areas.  Having a smaller maximum range size means that you can have more areas (the size of the
#' state space is greatly reduced; see \code{\link[cladoRcpp]{numstates_from_numareas}}.
#' @param states_list A list of the possible states/geographic ranges, in 0-based index form.
#' @param force_sparse Should sparse matrix exponentiation be used?  Default \code{FALSE}, which means dense matrix exponentiation is always used. 
#' If \code{NA}, the program will
#' use sparse matrix exponentiation for transition matrices above rank 128 (size 128x128).  NOTE: Sparse matrix exponentiation seems to give correlated, but 
#' not exact, results, and these errors may accumulate.  Presumably the problems become less with larger matrices, but I have not explored this in detail.
#' @param use_detection_model If TRUE, use the detection model (with parameters mf, dp, and fdp) and counts of detections and counts of taphonomic controls to 
#' calculate the \code{tip_condlikes_of_data_on_each_state}.
#' @param print_optim If TRUE (default), print the optimization steps as ML estimation progresses.
#' @param num_cores_to_use If >1, parallel processing will be attempted. \bold{Note:} parallel processing via \code{library (parallel)} will work in Mac command-line
#' R, but not in Mac GUI \code{R.app}.
#' @param cluster_already_open If the user wants to distribute the matrix exponentiation calculations from all the branches across a number of processors/nodes on 
#' a cluster, specify the cluster here.  E.g. \code{cluster_already_open = makeCluster(rep("localhost",num_cores_to_use), type = "SOCK")}.  Note: this will work on 
#' most platforms, including Macs running R from command line, but will NOT work on Macs running the R GUI \code{R.app}, because parallel processing functions like
#' \code{MakeCluster} from e.g. \code{library(parallel)} for some reason crash R.app.  The program runs a check for R.app and will just run on 1 node if found. 
#' @param use_optimx If \code{TRUE} (default) or \code{"optimx"}, use \code{\link[optimx]{optimx}} rather than \code{\link[stats]{optim}}. 
#' If \code{FALSE} or \code{"optim"}, use \code{\link[stats]{optim}}. If \code{"GenSA"}, use \code{\link[GenSA]{GenSA}}.
#' @param rescale_params If FALSE (default), parameters are used as-is. If TRUE, they are re-scaled, 
#' with the observed min (the min of "init", "est", and "min") in BioGeoBEARS_model_object@params_table 
#' re-set to zero, and the observed max (the max of "init", "est", and "max") set to 1. See functions
#' \code{\link{scale_BGB_params}} and \code{\link{unscale_BGB_params}} for the algorithm used to
#' scale and unscale the parameters. Scaling \em{might} be helpful for the ML search when 
#' parameters have much different sizes. Google e.g. discussions by Ben Bolker on scaling in ML
#' searches.
#' @param return_condlikes_table If \code{TRUE}, return the table of ALL conditional likelihood results, including at branch subsections
#' (only some should be used in calculating the final log-likelihood of the geography range data on the tree!)
#' @param calc_TTL_loglike_from_condlikes_table If TRUE, force making of the condlikes table, and use it to calculate the log-likelihood
#' (default=TRUE; matches LAGRANGE).
#' @param calc_ancprobs If \code{TRUE} (default), calculate and return the necessary pieces (uppass and downpass probs) for ancestral states.
#' @param fixnode If the state at a particular node is going to be fixed (e.g. for ML marginal ancestral states), give the node number.
#' @param fixlikes The state likelihoods to be used at the fixed node.  I.e. 1 for the fixed state, and 0 for the others.
#' @param speedup If \code{TRUE} (default), set the maximum number of iterations to \code{itnmax=50*(number of free parameters)}, instead of the
#' \code{optimx} default, 250.  Also set \code{optimx} \code{reltol} parameter to 0.001 (instead of the default, ~1e-8).
#' @param include_null_range Should the null range (no areas) be included as a state?  Default \code{TRUE}.
#' @param min_branchlength Nodes with branches below this branchlength will not be treated as cladogenesis events; instead, they will be treated as 
#' if an OTU had been sampled from an anagenetic lineage, i.e. as if you had a direct ancestor.  This is useful for putting fossils into the biogeography analysis,
#' when you have fossil species that range through time. (Note: the proper way to obtain such trees, given that most phylogenetic methods force all OTUs to be tips 
#' rather than direct ancestors, is another question subject to active research.  However, one method might be to just set a branch-length cutoff, and treat any
#' branches sufficiently small as direct ancestors.)
#' @param tmpwd The working directory in which the input and output files will be placed. Default is \code{\link[base]{getwd}}. This is stored 
#' mostly for future reference; users are responsible for manually navigating to the appropriate directory ahead of time, using \code{\link[base]{setwd}}.
#' @param printlevel Default 1 prints basic alerts.  2 prints even more. 0 will remove most of them.
#' @return \code{inputs} Inputs for ML search.
#' @export
#' @seealso \code{\link{readfiles_BioGeoBEARS_run}}, \code{\link[BioGeoBEARS]{define_BioGeoBEARS_model_object}}, \code{\link[base]{setwd}}, \code{\link[base]{getwd}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @examples
#' test=1
#' 
define_BioGeoBEARS_run <- function(abbr="default", description="defaults", BioGeoBEARS_model_object=define_BioGeoBEARS_model_object(minval_anagenesis=1e-12, minval_cladogenesis=1e-5, maxval=5), trfn=NULL, geogfn=NULL, timesfn=NA, distsfn=NA, dispersal_multipliers_fn=NA, area_of_areas_fn=NA, areas_allowed_fn=NA, areas_adjacency_fn=NA, detects_fn=NA, controls_fn=NA, max_range_size=NA, states_list=NULL, force_sparse=FALSE, use_detection_model=FALSE, print_optim=TRUE, printlevel=0, on_NaN_error=-1e50, num_cores_to_use=NA, cluster_already_open=FALSE, use_optimx=TRUE, rescale_params=FALSE, return_condlikes_table=FALSE, calc_TTL_loglike_from_condlikes_table=TRUE, calc_ancprobs=TRUE, fixnode=NULL, fixlikes=NULL, speedup=TRUE, include_null_range=TRUE, useAmbiguities=FALSE, min_branchlength=0.000001, tmpwd=getwd())
	{
	inputs = list()
	
	# Load the default input files
	#extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
	extdata_dir = system.file("extdata", package="BioGeoBEARS")
	if (is.null(geogfn))
		{
		geogfn = np(paste(addslash(extdata_dir), "Psychotria_geog.data", sep=""))
		inputs$geogfn <- geogfn
		}

	if (is.null(trfn))
		{
		trfn = np(paste(addslash(extdata_dir), "Psychotria_5.2.newick", sep=""))
		inputs$trfn <- trfn
		}

	inputs$trfn = trfn
	inputs$geogfn = geogfn

	
	inputs$abbr = abbr
	inputs$description = description
	inputs$BioGeoBEARS_model_object = BioGeoBEARS_model_object
	inputs$timesfn = timesfn
	inputs$distsfn = distsfn										# distance between areas, for dispersal ~ dist^x
	inputs$dispersal_multipliers_fn = dispersal_multipliers_fn		# hard-coded dispersal multiplier (or 0s/1s for constraints)
	inputs$area_of_areas_fn = area_of_areas_fn						# area of each areas (for extinction ~ area^u)
	inputs$areas_allowed_fn = areas_allowed_fn						# if ONLY these areas are allowed at certain timepoints
	inputs$areas_adjacency_fn = areas_adjacency_fn					# if ONLY connected areas are allowed
	inputs$detects_fn = detects_fn									# used only if use_detection_model==TRUE
	inputs$controls_fn = controls_fn								# used only if use_detection_model==TRUE
	inputs$max_range_size = max_range_size
	inputs$states_list = states_list
	inputs$force_sparse = force_sparse
	inputs$use_detection_model = use_detection_model
	inputs$print_optim = print_optim
	inputs$printlevel = printlevel
	inputs$on_NaN_error = on_NaN_error
	inputs$wd = tmpwd												# Store the working directory you are in
	inputs$num_cores_to_use = num_cores_to_use
	inputs$cluster_already_open = cluster_already_open
	inputs$use_optimx = use_optimx
	inputs$rescale_params = rescale_params
	inputs$return_condlikes_table = return_condlikes_table
	inputs$calc_TTL_loglike_from_condlikes_table = calc_TTL_loglike_from_condlikes_table
	inputs$calc_ancprobs = calc_ancprobs
	inputs$fixnode = fixnode
	inputs$fixlikes = fixlikes
	inputs$speedup = speedup
	
	# New, 2014-02-16
	inputs$include_null_range = include_null_range
	
	# 2014-05-07
	inputs$useAmbiguities = useAmbiguities

	# 2016-02-24
	inputs$min_branchlength = min_branchlength
	
	return(inputs)
	}


#######################################################
# BioGeoBEARS_run
#######################################################
# An object of class BioGeoBEARS_run holding the model inputs
#
#
#@section Slots: 
#  \describe{
#    \item{\code{list}:}{List of class \code{"list"}, containing inputs list from define_BioGeoBEARS_run}
#  }
#
# @note Go BEARS!
# @name BioGeoBEARS_run 
# @rdname BioGeoBEARS_run-class
# @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
# @seealso \code{\link{define_tipranges_object}}, \code{\link{getareas_from_tipranges_object}},
# \code{\link[cladoRcpp]{areas_list_to_states_list_old}}, \code{\link{areas_list_to_states_list_new}},
# \code{\link{tipranges_to_tip_condlikes_of_data_on_each_state}}
# @references
# \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
# @examples
# test=1
#
# setClass(Class="BioGeoBEARS_run", representation=representation(inputs="list"), prototype=define_BioGeoBEARS_run())

# setClass(Class="BioGeoBEARS_run", contains="list", slots=list(abbr = "character", 
# description = "character", 
# trfn = "character", 
# geogfn = "character", 
# BioGeoBEARS_model_object = "BioGeoBEARS_model", 
# timesfn = "character", 
# distsfn = "character", 
# dispersal_multipliers_fn = "character", 
# area_of_areas_fn = "character", 
# areas_allowed_fn = "character", 
# areas_adjacency_fn = "character", 
# detects_fn = "character", 
# controls_fn = "character", 
# max_range_size = "numeric", 
# force_sparse = "logical", 
# use_detection_model = "logical", 
# print_optim = "logical", 
# printlevel = "numeric", 
# on_NaN_error = "numeric", 
# num_cores_to_use = "numeric", 
# cluster_already_open = "character", 
# speedup = "logical", 
# include_null_range = "logical", 
# useAmbiguities = "logical", 
# min_branchlength = "numeric",
# use_optimx = "character", 
# rescale_params = "logical", 
# return_condlikes_table = "logical", 
# calc_TTL_loglike_from_condlikes_table = "logical", 
# calc_ancprobs = "logical"))



# dispersal_multipliers file is just a list of distance matrices, separated by blank lines,
# from youngest to oldest
#######################################################
# read_dispersal_multipliers_fn
#######################################################
#' Read in the hard-coded dispersal multipliers from file
#' 
#' dispersal_multipliers file is just a list of distance matrices, separated by blank lines,
#' from youngest to oldest
#' 
#' @param inputs The inputs list
#' @param dispersal_multipliers_fn The dispersal multipliers filename.
#' @return \code{list_of_dispersal_multipliers_mats} A list object
#' @export
#' @seealso \code{\link[stats]{convolve}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @examples
#' test=1
read_dispersal_multipliers_fn <- function(inputs=NULL, dispersal_multipliers_fn=NULL)
	{
	defaults='
	dispersal_multipliers_fn = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/Psychotria_dists_stratified/Hawaii_KOMH_dispersal_multipliers.txt"
	'

	if (!is.null(inputs))
		{
		dispersal_multipliers_fn = inputs$dispersal_multipliers_fn
		}
	if (!is.null(dispersal_multipliers_fn))
		{
		dispersal_multipliers_fn = dispersal_multipliers_fn
		}
	
	
	tmplines = readLines(dispersal_multipliers_fn)
	list_of_dispersal_multipliers_mats = list()
	lnum = 1
	
	newmat = TRUE
	for (i in 1:length(tmplines))
		{
		# Using a quick "trim" operation to avoid the 
		# especially the annoying end-of-line issues
		tmplines[i] = quicktrim(tmplines[i])

		if (tmplines[i] == "END")
			{
			return(list_of_dispersal_multipliers_mats)
			}
		
		if (tmplines[i] == "")
			{
			# You've hit the end of a matrix,
			# increment list and add
			tmpmat = as.matrix(tmprows)
			tmpmat = adf2(tmpmat)
			names(tmpmat) = nameslist

			# Check that matrices are square
			dims = dim(tmpmat)
			if (dims[1] != dims[2])
				{
				stoptxt = paste0("\nFATAL ERROR: the matrix is not square! Instead, yours has ", dims[1], " rows x ", dims[2], " columns. Check your input file. Printing the matrix below.\n")
		
				cat(stoptxt)
				print(tmpmat)
				cat("\n\n")
				stop(stoptxt)
				} # END if (dims[1] != dims[2])
	


			list_of_dispersal_multipliers_mats[[lnum]] = tmpmat
			lnum = lnum + 1
			
			newmat = TRUE
			} else {
			if (newmat == TRUE)
				{
				nameslist = strsplit(tmplines[i], split="\t")[[1]]
				newmat = FALSE
				tmprows = NULL
				} else {
				tmprow = as.numeric(strsplit(tmplines[i], split="\t")[[1]])
				
				# Error check
				if (length(nameslist) == length(tmprow))
					{
					tmprows = rbind(tmprows, tmprow)
					} else {
					errortxt = paste("\n\nERROR in read_dispersal_multipliers_fn():\n\nIn line ", i, " of '", dispersal_multipliers_fn, "', you have ", length(nameslist), " area names, but ", length(tmprow), " numbers.\n\nLook at the file and fix this mismatch to proceed. Have a nice day!\n\n", sep="")
					cat(errortxt)
					
					stop("STOP ERROR in read_dispersal_multipliers_fn().\n\n")
					
					} # END if (length(nameslist) == length(tmprow))

				}
			}
		}
	list_of_dispersal_multipliers_mats
	
	

	
	return(list_of_dispersal_multipliers_mats)
	}




# Timeperiods file is just a list of times, 1 per line, from
# youngest to oldest
#######################################################
# read_times_fn
#######################################################
#' Read in the stratification time breakpoints
#' 
#' The timeperiods file is just a list of times, 1 per line, from
#' youngest to oldest. Do not include time=0.
#' 
#' @param inputs The inputs list
#' @param timesfn The times filename.
#' @return \code{timeperiods} A list object
#' @export
#' @seealso \code{\link[stats]{convolve}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @examples
#' test=1
read_times_fn <- function(inputs=NULL, timesfn=NULL)
	{
	defaults='
	timesfn = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/Psychotria_dists_stratified/Hawaii_timeperiods.txt"
	'
	
	if (!is.null(inputs))
		{
		timesfn = inputs$timesfn
		}
	if (!is.null(timesfn))
		{
		timesfn = timesfn
		timeperiods = as.numeric(readLines(timesfn))
		
		# Error check
		if ( all(timeperiods==sort(timeperiods)) == FALSE || (0 %in% timeperiods)==TRUE )
			{
			txt = "STOP ERROR in read_times_fn(): 'times' in a time-stratified analysis must be sorted from youngest to oldest in the input timeperiods file (the dispersal matrices etc. should be in the same order, also). Also, time '0' should *not* be included. Either of these will cause this error. So, check your timeperiods input file."
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			stop(txt)
			} # END if ( all( times == sort(times) ) )

		} else {
		timeperiods = NULL
		}
	return(timeperiods)
	}



#' library(ape)
#' library(BioGeoBEARS)
#' 
#' extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
#' 
#' # Load timeperiods from file
#' times_fn = np(paste(addslash(extdata_dir), "examples/Psychotria_M3strat/BGB/timeperiods.txt", sep=""))
#' inputs = NULL
#' inputs$times_fn = times_fn
#' timeperiods = read_times_fn(inputs$times_fn)
#' 
#' # Or set the timeperiods manually
#' timeperiods = c(0.5, 1.9, 3.7, 5.1, 10)
#' 
#' # Load tree from file
#' trfn = np(paste(addslash(extdata_dir), "Psychotria_5.2.newick", sep=""))
#' tr = read.tree(trfn)
#' 
#' # Run the function
#' counts_df = sum_nodes_branchlengths_by_timeperiod(tr, timeperiods)
#' counts_df
#' 
#' # Check that the sums add up
#' sum(tr$edge.length)
#' sum(counts_df$branchlength_sums)
#' 
#' tr$Nnode
#' sum(counts_df$node_counts)
#' 
#' #   time_tops time_bots node_counts branchlength_sums
#' # 1       0.0       0.5           0          9.500000
#' # 2       0.5       1.9           6         23.908071
#' # 3       1.9       3.7           9         13.601658
#' # 4       3.7       5.1           2          4.490857
#' # 5       5.1      10.0           1          0.200000
#' 
sum_nodes_branchlengths_by_timeperiod <- function(tr, timeperiods)
	{
	# Check if "tr" is useable
	if (("phylo" %in% class(tr)) == TRUE)
		{
		trfn = "temp_tree.newick"
		write.tree(tr, file=trfn)
		} else if (("character" %in% class(tr)) == TRUE) {
		trfn = tr
		} else {
		txt = "STOP ERROR in sum_nodes_branchlengths_by_timeperiod(): argument 'tr' must be either a 'phylo' tree object, or a 'character' filename."
		stop(txt)
		}
	
	
	
	inputs = NULL
	inputs$trfn = trfn
	inputs$timeperiods = timeperiods

	tree_pieces = section_the_tree(inputs, make_master_table=TRUE, plot_pieces=FALSE, cut_fossils=FALSE)
	
	# WATCH OUT: Because of this, when cut_fossils=FALSE, the sum of edge lengths will be wrong...
	# 
	# This will extend ALL tips up to time_bp=0 my.  Keep track of true tip age through orig_tr_table$fossils and orig_tr_table$time_bp
	# 	phy_as_it_is_chopped_down = extend_tips_to_ultrametricize(obj=phy_as_it_is_chopped_down, age_of_root=0, tips_end_at_this_date=NA)
	#
	# ...check each tree pice
	
	
	tree_pieces

	node_counts = NULL
	branchlength_sums = NULL

	for (i in 1:length(tree_pieces$tree_sections_list))
		{
		tree_piece = tree_pieces$tree_sections_list[[i]]
		branchlength_sum = 0
		nodes_sum = 0
	
		for (j in 1:length(tree_piece$return_pieces_list))
			{
			if ("numeric" %in% class(tree_piece$return_pieces_list[[j]]))
				{
				# Check if this branch segment actually is actually alive during or after this time bin
				stratum = i
				piecenum = j
				TF1 = tree_pieces$master_table$stratum == stratum
				TF2 = tree_pieces$master_table$piecenum == piecenum
				TF = (TF1 + TF2) == 2
				master_row = tree_pieces$master_table[TF,]
				if (master_row$time_bp > master_row$time_bot)
					{
					next() # skip; this fossil tip is older than the time bin
					}
				else if ( (master_row$time_bp < master_row$time_bot) && (master_row$time_bp > master_row$time_top) )
					{
					# This branch tip ends in this time bin; add it's segmental length
					length_of_segment = master_row$time_bot - master_row$time_bp
					branchlength_sum = branchlength_sum + length_of_segment
					}
				else
					{
					# This branch is continuing through the time bin; add it
					branchlength_sum = branchlength_sum + tree_piece$return_pieces_list[[j]]
					}
				} else {
				tmptr = tree_piece$return_pieces_list[[j]]
				nodes_sum = nodes_sum + tmptr$Nnode
				if (is.null(tmptr$root.edge))
					{
					branchlength_sum = branchlength_sum + sum(tmptr$edge.length) + 0.0
					} else {
					branchlength_sum = branchlength_sum + sum(tmptr$edge.length) + tmptr$root.edge
					}
				} # END if ("numeric" %in% class
			} # END for (j in 1:length(tree_piece$return_pieces_list))
		node_counts = c(node_counts, nodes_sum)
		branchlength_sums = c(branchlength_sums, branchlength_sum)	
		} # END for (i in 1:length(tree_pieces))
	
	time_tops = c(0.0, timeperiods[1:(length(timeperiods)-1)])
	time_bots = timeperiods
	tmpmat = cbind(time_tops, time_bots, node_counts, branchlength_sums)
	counts_df = as.data.frame(tmpmat, stringsAsFactors=FALSE)
	names(counts_df) = c("time_tops", "time_bots", "node_counts", "branchlength_sums")
	counts_df
	} # END sum_nodes_branchlengths_by_timeperiod <- function(tr, timeperiods)





subset_distmats <- function(distmats_list, rows_to_keep_TF, replace_NAs_with=0.0)
	{
	defaults='
	tmpfn = "geological_distances_v3_div100_stay_same.txt"
	distmats_list = read_distances_fn(inputs=NULL, distsfn=tmpfn)
	rows_to_keep_TF = c(TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE)

	new_distmats_list = subset_distmats(distmats_list, rows_to_keep_TF=rows_to_keep_TF, replace_NAs_with=0.0)

	'
	new_distmats_list = list()

	for (i in 1:length(distmats_list))
		{
		# subset cols
		tmpmat = distmats_list[[i]]
		tmpmat2 = tmpmat[,rows_to_keep_TF]
		tmpmat3 = tmpmat2[rows_to_keep_TF,]
	
		# Replace NAs with 0
		tmpmat3[is.na(tmpmat3)] = replace_NAs_with
	
		new_distmats_list[[i]] = tmpmat3
		}
	return(new_distmats_list)
	} # END subset_distmats <- function(distmats_list, rows_to_keep_TF, replace_NAs_with=0.0)


write_distances_to_fn <- function(new_distmats_list, outfn)
	{
	# Write distances matrix
	defaults='
	tmpfn = "geological_distances_v3_div100_stay_same.txt"
	outfn = gsub(pattern=".txt", replacement="_subset.txt", x=tmpfn)
	'
	for (i in 1:length(new_distmats_list))
		{
		if (i == 1)
			{
			write.table(x=new_distmats_list[i], file=outfn, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE, append=FALSE)
			} else {
			write.table(x=new_distmats_list[i], file=outfn, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE, append=TRUE)		
			}
	
		# Line break	
		if (i < length(new_distmats_list))
			{
			write(x="", file=outfn, append=TRUE)
			} else {
			write(x="", file=outfn, append=TRUE)
			write(x="END", file=outfn, append=TRUE)
			}
		}
	
	return()
	} # END write_distances_to_fn <- function(new_distmats_list, outfn)



# Distances file is just a list of distance matrices, separated by blank lines,
# from youngest to oldest
#######################################################
# read_distances_fn
#######################################################
#' Read in the distances by time
#' 
#' Distances file is just a list of distance matrices, separated by blank lines,
#' from youngest to oldest. Tab-delimited. 
#' 
#' Column headers: yes, row names: no. Ends with blank line, "END", blank line.
#' 
#' @param inputs The inputs list
#' @param distsfn The distances filename.
#' @return \code{list_of_distances_mats} A list object
#' @export
#' @seealso \code{\link[stats]{convolve}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @examples
#' test=1
read_distances_fn <- function(inputs=NULL, distsfn=NULL)
	{
	defaults='
	distsfn = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/Psychotria_dists_stratified/Hawaii_KOMH_distances.txt"
	'

	if (!is.null(inputs))
		{
		distsfn = inputs$distsfn
		}
	if (!is.null(distsfn))
		{
		distsfn = distsfn
		}
	
	
	tmplines = readLines(distsfn)
	list_of_distances_mats = list()
	lnum = 1
	
	newmat = TRUE
	for (i in 1:length(tmplines))
		{
		# Using a quick "trim" operation to avoid the 
		# especially the annoying end-of-line issues
		tmplines[i] = quicktrim(tmplines[i])
		
		print(tmplines[i])
		
		if (tmplines[i] == "END")
			{
			return(list_of_distances_mats)
			}
		
		if (tmplines[i] == "")
			{
			# You've hit the end of a matrix,
			# increment list and add
			tmpmat = as.matrix(tmprows)
			tmpmat = adf2(tmpmat)
			names(tmpmat) = nameslist

			# Correct distances if any are <= 0
			# (to ensure consistent behavior of the exponent, i.e. dist ^ (x_exponent)
			tmpmat = as.matrix(tmpmat)
			diag(tmpmat) = NA
			minval = min(tmpmat, na.rm=TRUE)
 			if (minval <= 0)
 				{
 				# Print
 				stoptxt = "STOP ERROR in read_distances_fn(): Minimum distance between regions must be > 0"
 				cat("\n\n")
 				cat(stoptxt)
 				cat("\n\n")
				cat("Printing locations of zero values in the matrix:\n\n") 				
 				# Find the error locations
				TF = tmpmat == 0
				rownums_x_colnums = expand.grid(1:nrow(tmpmat), 1:ncol(tmpmat))
				head(rownums_x_colnums)
				temp = rownums_x_colnums[c(TF),]
				
				rownums_x_colnums = temp[!is.na(temp[,"Var1"]),]
				namevals = colnames(tmpmat)
				
				
				for (j in 1:nrow(rownums_x_colnums))
					{
					val = tmpmat[rownums_x_colnums[j,1], rownums_x_colnums[j,2]]
					tmpstr = paste0(rownums_x_colnums[j,1], ":", namevals[rownums_x_colnums[j,1]], ", ", rownums_x_colnums[j,2], ":", namevals[rownums_x_colnums[j,2]], " = ", val)
					print(tmpstr)
					}
				
 				stop(stoptxt)
 				}

			# Check that matrices are square
			dims = dim(tmpmat)
			if (dims[1] != dims[2])
				{
				stoptxt = paste0("\nFATAL ERROR: the matrix is not square! Instead, yours has ", dims[1], " rows x ", dims[2], " columns. Check your input file. Printing the matrix below.\n")
		
				cat(stoptxt)
				print(tmpmat)
				cat("\n\n")
				stop(stoptxt)
				} # END if (dims[1] != dims[2])
				
			list_of_distances_mats[[lnum]] = tmpmat
			lnum = lnum + 1
			
			newmat = TRUE
			} else {
			if (newmat == TRUE)
				{
				nameslist = strsplit(tmplines[i], split="\\s+")[[1]]
				newmat = FALSE
				tmprows = NULL
				} else {
				tmprow = as.numeric(strsplit(tmplines[i], split="\t")[[1]])

				# Error check
				if (length(nameslist) == length(tmprow))
					{
					tmprows = rbind(tmprows, tmprow)
					} else {
					errortxt = paste("\n\nSTOP ERROR in read_distances_fn():\n\nIn line ", i, " of '", distsfn, "', you have ", length(nameslist), " area names, but ", length(tmprow), " numbers.\n\nLook at the file and fix this mismatch to proceed. Have a nice day!\n\n", sep="")
					cat(errortxt)
					
					stop("STOP ERROR in read_distances_fn().\n\n")
					
					} # END if (length(nameslist) == length(tmprow))
				}
			}
		} # END for (i in 1:length(tmplines))
	list_of_distances_mats
	return(list_of_distances_mats)
	}



# Envdistances file is just a list of envdistance matrices, separated by blank lines,
# from youngest to oldest
#######################################################
# read_envdistances_fn
#######################################################
#' Read in the envdistances by time
#' 
#' Distances file is just a list of envdistance matrices, separated by blank lines,
#' from youngest to oldest. Tab-delimited. 
#' 
#' Column headers: yes, row names: no. Ends with blank line, "END", blank line.
#' 
#' @param inputs The inputs list
#' @param envdistsfn The envdistances filename.
#' @return \code{list_of_envdistances_mats} A list object
#' @export
#' @seealso \code{\link[stats]{convolve}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @examples
#' test=1
read_envdistances_fn <- function(inputs=NULL, envdistsfn=NULL)
	{
	defaults='
	envdistsfn = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/Psychotria_envdists_stratified/Hawaii_KOMH_envdistances.txt"
	'

	if (!is.null(inputs))
		{
		envdistsfn = inputs$envdistsfn
		}
	if (!is.null(envdistsfn))
		{
		envdistsfn = envdistsfn
		}
	
	
	tmplines = readLines(envdistsfn)
	list_of_envdistances_mats = list()
	lnum = 1
	
	newmat = TRUE
	for (i in 1:length(tmplines))
		{
		# Using a quick "trim" operation to avoid the 
		# especially the annoying end-of-line issues
		tmplines[i] = quicktrim(tmplines[i])

		print(tmplines[i])
		
		if (tmplines[i] == "END")
			{
			return(list_of_envdistances_mats)
			}
		
		if (tmplines[i] == "")
			{
			# You've hit the end of a matrix,
			# increment list and add
			tmpmat = as.matrix(tmprows)
			tmpmat = adf2(tmpmat)
			names(tmpmat) = nameslist

			# Correct envdistances if any are <= 0
			# (to ensure consistent behavior of the exponent, i.e. envdist ^ (x_exponent)
			tmpmat = as.matrix(tmpmat)
			diag(tmpmat) = NA
			minval = min(tmpmat, na.rm=TRUE)
 			if (minval <= 0)
 				{
 				# Print
 				stoptxt = "STOP ERROR in read_envdistances_fn(): Minimum distance between regions must be > 0"
 				cat("\n\n")
 				cat(stoptxt)
 				cat("\n\n")
				cat("Printing locations of zero values in the matrix:\n\n") 				
 				# Find the error locations
				TF = tmpmat == 0
				rownums_x_colnums = expand.grid(1:nrow(tmpmat), 1:ncol(tmpmat))
				head(rownums_x_colnums)
				temp = rownums_x_colnums[c(TF),]
				
				rownums_x_colnums = temp[!is.na(temp[,"Var1"]),]
				namevals = colnames(tmpmat)
				
				
				for (j in 1:nrow(rownums_x_colnums))
					{
					val = tmpmat[rownums_x_colnums[j,1], rownums_x_colnums[j,2]]
					tmpstr = paste0(rownums_x_colnums[j,1], ":", namevals[rownums_x_colnums[j,1]], ", ", rownums_x_colnums[j,2], ":", namevals[rownums_x_colnums[j,2]], " = ", val)
					print(tmpstr)
					}
				
 				stop(stoptxt)
 				}

			# Check that matrices are square
			dims = dim(tmpmat)
			if (dims[1] != dims[2])
				{
				stoptxt = paste0("\nFATAL ERROR: the matrix is not square! Instead, yours has ", dims[1], " rows x ", dims[2], " columns. Check your input file. Printing the matrix below.\n")
		
				cat(stoptxt)
				print(tmpmat)
				cat("\n\n")
				stop(stoptxt)
				} # END if (dims[1] != dims[2])
				
			list_of_envdistances_mats[[lnum]] = tmpmat
			lnum = lnum + 1
			
			newmat = TRUE
			} else {
			if (newmat == TRUE)
				{
				nameslist = strsplit(tmplines[i], split="\t")[[1]]
				newmat = FALSE
				tmprows = NULL
				} else {
				tmprow = as.numeric(strsplit(tmplines[i], split="\t")[[1]])

				# Error check
				if (length(nameslist) == length(tmprow))
					{
					tmprows = rbind(tmprows, tmprow)
					} else {
					errortxt = paste("\n\nERROR in read_envdistances_fn():\n\nIn line ", i, " of '", envdistsfn, "', you have ", length(nameslist), " area names, but ", length(tmprow), " numbers.\n\nLook at the file and fix this mismatch to proceed. Have a nice day!\n\n", sep="")
					cat(errortxt)
					
					stop("STOP ERROR in read_envdistances_fn().\n\n")
					
					} # END if (length(nameslist) == length(tmprow))
				}
			}
		}
	list_of_envdistances_mats
	return(list_of_envdistances_mats)
	}


# area_areas file is just a list of distance matrices, separated by blank lines,
# from youngest to oldest

#######################################################
# read_area_of_areas_fn
#######################################################
#' Read in the area areas by time
#' 
#' area_areas file is just a list of vectors of areas, separated by blank lines,
#' from youngest to oldest. (e.g., like a distances file, but with
#' only 2 lines per time period: (1) the area labels, (2) the areas of each area.
#'
#' The area-of-areas functionality was included for experimental purposes, but
#' given the unreliability of estimating <i>e</i> in DEC-type models, it is 
#' probably of dubious utility.
#' 
#' @param inputs The inputs list
#' @param area_of_areas_fn The area-of-areas filename.
#' @return \code{list_of_area_areas_mats} A list object
#' @export
#' @seealso \code{\link[stats]{convolve}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @examples
#' test=1
read_area_of_areas_fn <- function(inputs=NULL, area_of_areas_fn=NULL)
	{
	defaults='
	area_areas_fn = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/Psychotria_dists_stratified/Hawaii_KOMH_area_of_areas.txt"
	'

	if (!is.null(inputs))
		{
		area_of_areas_fn = inputs$area_of_areas_fn
		}
	if (!is.null(area_of_areas_fn))
		{
		area_of_areas_fn = area_of_areas_fn
		}
	
	
	tmplines = readLines(area_of_areas_fn)
	list_of_area_areas_mats = list()
	lnum = 1
	
	newmat = TRUE
	for (i in 1:length(tmplines))
		{
		# Using a quick "trim" operation to avoid the 
		# especially the annoying end-of-line issues
		tmplines[i] = quicktrim(tmplines[i])

		if (tmplines[i] == "END")
			{
			return(list_of_area_areas_mats)
			}
		
		if (tmplines[i] == "")
			{
			# You've hit the end of a matrix,
			# increment list and add
			tmpmat = as.matrix(tmprows)
			tmpmat = adf2(tmpmat)
			names(tmpmat) = nameslist

			# Correct distances if any are <= 0
			# (to ensure consistent behavior of the exponent, i.e. dist ^ (x_exponent)
			tmpmat = as.matrix(tmpmat)
			minval = min(tmpmat, na.rm=TRUE)
 			if (minval <= 0)
 				{
 				stop("\n\nERROR: Minimum area of a region must be > 0")
 				}


			list_of_area_areas_mats[[lnum]] = tmpmat
			lnum = lnum + 1
			
			newmat = TRUE
			} else {
			if (newmat == TRUE)
				{
				nameslist = strsplit(tmplines[i], split="\t")[[1]]
				newmat = FALSE
				tmprows = NULL
				} else {
				tmprow = as.numeric(strsplit(tmplines[i], split="\t")[[1]])

				# Error check
				if (length(nameslist) == length(tmprow))
					{
					tmprows = rbind(tmprows, tmprow)
					} else {
					errortxt = paste("\n\nERROR in read_area_of_areas_fn():\n\nIn line ", i, " of '", area_of_areas_fn, "', you have ", length(nameslist), " area names, but ", length(tmprow), " numbers.\n\nLook at the file and fix this mismatch to proceed. Have a nice day!\n\n", sep="")
					cat(errortxt)
					
					stop("STOP ERROR in read_area_of_areas_fn().\n\n")
					
					} # END if (length(nameslist) == length(tmprow))
				}
			}
		}
	list_of_area_areas_mats
	return(list_of_area_areas_mats)
	}





# areas_allowed file is just a list of distance matrices, separated by blank lines,
# from youngest to oldest

#######################################################
# read_areas_allowed_fn
#######################################################
#' Read in the area areas by time
#' 
#' The areas_allowed file is just a list of 1/0 matrices, separated by blank lines,
#' from youngest to oldest. 1s represent allowed areas: ranges/states that 
#' include areas set to '0' will be disallowed from the state space.  Numbers 
#' within a matrix are tab-delimited.
#'
#' Technically, we could specify areas-allowed with just a single row, rather than a 
#' matrix, but it will be kept as a matrix for back-compatibility (just make sure the 
#' matrix is symmetric about the diagonal).
#'
#' Except for simple cases, a better way to specify different areas allowed is to
#' use a manual modification of the states_list; see instructions on PhyloWiki.
#' 
#' @param inputs The inputs list
#' @param areas_allowed_fn The areas-allowed filename.
#' @return \code{list_of_areas_allowed_mats} A list object
#' @export
#' @seealso \code{\link[stats]{convolve}} \code{\link{read_areas_adjacency_fn}} \code{\link[cladoRcpp]{rcpp_areas_list_to_states_list}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @examples
#' test=1
read_areas_allowed_fn <- function(inputs=NULL, areas_allowed_fn=NULL)
	{
	defaults='
	areas_allowed_fn = "/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/Psychotria_M3strat/areas_allowed.txt"
	'

	if (!is.null(inputs))
		{
		areas_allowed_fn = inputs$areas_allowed_fn
		}
	if (!is.null(areas_allowed_fn))
		{
		areas_allowed_fn = areas_allowed_fn
		}
	
	
	tmplines = readLines(areas_allowed_fn)
	list_of_areas_allowed_mats = list()
	lnum = 1
	
	newmat = TRUE
	for (i in 1:length(tmplines))
		{
		# Using a quick "trim" operation to avoid the 
		# especially the annoying end-of-line issues
		tmplines[i] = quicktrim(tmplines[i])

		if (tmplines[i] == "END")
			{
			return(list_of_areas_allowed_mats)
			}
			
		if (tmplines[i] == "")
			{
			# You've hit the end of a matrix,
			# increment list and add
			tmpmat = as.matrix(tmprows)
			tmpmat = adf2(tmpmat)
			names(tmpmat) = nameslist

			# Check that matrices are square
			dims = dim(tmpmat)
			if (dims[1] != dims[2])
				{
				stoptxt = paste0("\nFATAL ERROR: the matrix is not square! Instead, yours has ", dims[1], " rows x ", dims[2], " columns. Check your input file. Printing the matrix below.\n")
		
				cat(stoptxt)
				print(tmpmat)
				cat("\n\n")
				stop(stoptxt)
				} # END if (dims[1] != dims[2])
	
			list_of_areas_allowed_mats[[lnum]] = tmpmat
			lnum = lnum + 1
			
			newmat = TRUE
			} else {
			if (newmat == TRUE)
				{
				nameslist = strsplit(tmplines[i], split="\t")[[1]]
				newmat = FALSE
				tmprows = NULL
				} else {
				tmprow = as.numeric(strsplit(tmplines[i], split="\t")[[1]])

				# Error check
				if (length(nameslist) == length(tmprow))
					{
					tmprows = rbind(tmprows, tmprow)
					} else {
					errortxt = paste("\n\nERROR in read_areas_allowed_fn():\n\nIn line ", i, " of '", areas_allowed_fn, "', you have ", length(nameslist), " area names, but ", length(tmprow), " numbers.\n\nLook at the file and fix this mismatch to proceed. Have a nice day!\n\n", sep="")
					cat(errortxt)
					
					stop("STOP ERROR in read_areas_allowed_fn().\n\n")
					
					} # END if (length(nameslist) == length(tmprow))
				}
			}
		}
	list_of_areas_allowed_mats
	return(list_of_areas_allowed_mats)
	}



# areas_adjacency file is just a list of distance matrices, separated by blank lines,
# from youngest to oldest

#######################################################
# read_areas_adjacency_fn
#######################################################
#' Read in the area areas by time
#' 
#' The areas_adjacency file is just a list of 1/0 matrices, separated by blank lines,
#' from youngest to oldest. 1s represent allowed combinations of areas that are
#' allowed as ranges in the state space.  Numbers within a matrix are tab-delimited.
#'
#' This is meant to work like the areas_adjacency option in Python Lagrange. See 
#' discussion in the BioGeoBEARS Google Group, here: https://groups.google.com/forum/#!searchin/biogeobears/areas_allowed/biogeobears/7e7U9NtPTBU/4MyxckASMksJ
#'
#' It is possible that the areas_adjacency option might not cover every possible
#' constraint that users might want to put on the list of possible ranges/states. Users
#' can always edit the states_list manually or with a function before analysis.
#' 
#' The areas_adjacency matrix strategy can only encode very simple cases; for more complex cases, 
#' users should manually specify the states_list in each time stratum, excluding
#' area combinations they think are impossible. See PhyloWiki for instructions.
#'
#' @param inputs The inputs list
#' @param areas_adjacency_fn The areas-adjacency filename.
#' @return \code{list_of_areas_adjacency_mats} A list object
#' @export
#' @seealso \code{\link[stats]{convolve}} \code{\link{read_areas_allowed_fn}} \code{\link[cladoRcpp]{rcpp_areas_list_to_states_list}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://groups.google.com/forum/#!searchin/biogeobears/areas_allowed/biogeobears/7e7U9NtPTBU/4MyxckASMksJ}
#' @examples
#' test=1
read_areas_adjacency_fn <- function(inputs=NULL, areas_adjacency_fn=NULL)
	{
	defaults='
	areas_adjacency_fn = "/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/Psychotria_M3strat/BGB/areas_adjacency.txt"
	'

	if (!is.null(inputs))
		{
		areas_adjacency_fn = inputs$areas_adjacency_fn
		}
	if (!is.null(areas_adjacency_fn))
		{
		areas_adjacency_fn = areas_adjacency_fn
		}
	
	
	tmplines = readLines(areas_adjacency_fn)
	list_of_areas_adjacency_mats = list()
	lnum = 1
	
	newmat = TRUE
	for (i in 1:length(tmplines))
		{
		# Using a quick "trim" operation to avoid the 
		# especially the annoying end-of-line issues
		tmplines[i] = quicktrim(tmplines[i])

		if (tmplines[i] == "END")
			{
			return(list_of_areas_adjacency_mats)
			}
			
		if (tmplines[i] == "")
			{
			# You've hit the end of a matrix,
			# increment list and add
			tmpmat = as.matrix(tmprows)
			tmpmat = adf2(tmpmat)
			names(tmpmat) = nameslist

			# Check that matrices are square
			dims = dim(tmpmat)
			if (dims[1] != dims[2])
				{
				stoptxt = paste0("\nFATAL ERROR: the matrix is not square! Instead, yours has ", dims[1], " rows x ", dims[2], " columns. Check your input file. Printing the matrix below.\n")
		
				cat(stoptxt)
				print(tmpmat)
				cat("\n\n")
				stop(stoptxt)
				} # END if (dims[1] != dims[2])
	
			list_of_areas_adjacency_mats[[lnum]] = tmpmat
			lnum = lnum + 1
			
			newmat = TRUE
			} else {
			if (newmat == TRUE)
				{
				nameslist = strsplit(tmplines[i], split="\t")[[1]]
				newmat = FALSE
				tmprows = NULL
				} else {
				tmprow = as.numeric(strsplit(tmplines[i], split="\t")[[1]])

				# Error check
				if (length(nameslist) == length(tmprow))
					{
					tmprows = rbind(tmprows, tmprow)
					} else {
					errortxt = paste("\n\nERROR in read_areas_adjacency_fn():\n\nIn line ", i, " of '", areas_adjacency_fn, "', you have ", length(nameslist), " area names, but ", length(tmprow), " numbers.\n\nLook at the file and fix this mismatch to proceed. Have a nice day!\n\n", sep="")
					cat(errortxt)
					
					stop("STOP ERROR in read_areas_adjacency_fn().\n\n")
					
					} # END if (length(nameslist) == length(tmprow))
				}
			}
		}
	list_of_areas_adjacency_mats
	return(list_of_areas_adjacency_mats)
	} # END read_areas_adjacency_fn



#######################################################
# Cut down the states list according to areas_allowed_mat
#######################################################
#######################################################
# prune_states_list
#######################################################
#' Cut down the states list according to areas_allowed_mat
#' 
#' Go through a list of states.  Remove states that represent areas disallowed according
#' to areas_allowed_mat.  It is assumed (crucial!) that the areas in the \code{states_list}, 
#' and in the \code{areas_allowed_mat}, have the same order.
#' 
#' @param states_list_0based_index A states_list, 0-based, e.g. from \code{\link[cladoRcpp]{rcpp_areas_list_to_states_list}}
#' @param areas_allowed_mat The matrix of area combinations allowed (represented by 1s)
#' @return \code{states_list_0based_index_new} A 0-based list of allowed states/ranges
#' @export
#' @seealso \code{\link[cladoRcpp]{rcpp_areas_list_to_states_list}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @examples
#' test=1
prune_states_list <- function(states_list_0based_index, areas_allowed_mat)
	{
	defaults='
	areas_allowed_fn = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/Psychotria_M3strat/BGB/areas_allowed.txt"
	moref(areas_allowed_fn)
	
	tmpinputs = NULL
	tmpinputs$areas_allowed_fn = areas_allowed_fn
	list_of_areas_allowed_mats = read_areas_allowed_fn(inputs=tmpinputs)
	list_of_areas_allowed_mats
	areas_allowed_mat = list_of_areas_allowed_mats[[1]]
	areas_allowed_mat
	
	areas = c("K", "O", "M", "H")
	states_list_0based_index = rcpp_areas_list_to_states_list(areas)
	states_list_0based_index
	
	prune_states_list(states_list_0based_index, list_of_areas_allowed_mats[[1]])
	prune_states_list(states_list_0based_index, list_of_areas_allowed_mats[[2]])
	prune_states_list(states_list_0based_index, list_of_areas_allowed_mats[[3]])
	prune_states_list(states_list_0based_index, list_of_areas_allowed_mats[[4]])
	

	areas_allowed_Sara = matrix(data=c(1,1,0,0,1,1,1,0,0,1,1,1,0,0,1,1), ncol=4, byrow=TRUE)
	areas_allowed_Sara
	
	prune_states_list(states_list_0based_index, areas_allowed_Sara)
	

	
	'
	
	
	# First, eliminate the null range (represented by NA) if present
	nums = 1:length(states_list_0based_index)
	NA_TF = is.na(states_list_0based_index)
	nonNA_range_nums = nums[NA_TF == FALSE]
	
	mini_func <- function(areas_0based, areas_allowed_mat)
		{
		areas_1based = areas_0based + 1
		
		tmpmat = areas_allowed_mat[areas_1based, areas_1based]
		
		if (is.null(dim(tmpmat)))
			{
			if (tmpmat == 1)
				{
				return(TRUE)
				} else {
				return(FALSE)
				}
			} else {
			if (sum(tmpmat[1,]) == ncol(tmpmat))
				{
				return(TRUE)
				} else {
				return(FALSE)				
				}
			}
		}


	match_TF = sapply(X=states_list_0based_index[nonNA_range_nums], FUN=mini_func, areas_allowed_mat=areas_allowed_mat)
	match_TF
	
	if (sum(NA_TF) == 1)
		{
		match_TF = c(TRUE, match_TF)
		} else {
		pass = 1
		# match_TF = match_TF
		}
	
	states_list_0based_index_new = states_list_0based_index[match_TF]


	# Print to show it
	areas_allowed_mat
	states_list_0based_index_new

	return(states_list_0based_index_new)
	}





#######################################################
# Cut down the states list according to areas_adjacency_mat
#######################################################
#######################################################
# prune_states_list_by_adjacency
#######################################################
#' Cut down the states list according to areas_adjacency_mat
#' 
#' Go through a list of states.  Remove states that represent areas disadjacent 
#' according to areas_adjacency_mat.  It is assumed (crucial!) that the areas 
#' in the \code{states_list}, 
#' and in the \code{areas_adjacency_mat}, have the same order.
#' 
#' This is meant to work like the areas_adjacency option in Python Lagrange. See 
#' discussion in the BioGeoBEARS Google Group, here: https://groups.google.com/forum/#!searchin/biogeobears/areas_allowed/biogeobears/7e7U9NtPTBU/4MyxckASMksJ
#' 
#' @param states_list_0based_index A states_list, 0-based, e.g. from \code{\link[cladoRcpp]{rcpp_areas_list_to_states_list}}
#' @param areas_adjacency_mat The matrix of area combinations adjacency 
#' (represented by 1s)
#' @return \code{states_list_0based_index_new} A 0-based list of adjacency states/ranges
#' @export
#' @seealso \code{\link{read_areas_adjacency_fn}} \code{\link[cladoRcpp]{rcpp_areas_list_to_states_list}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://groups.google.com/forum/#!searchin/biogeobears/areas_allowed/biogeobears/7e7U9NtPTBU/4MyxckASMksJ}
#' @examples
#' test=1
prune_states_list_by_adjacency <- function(states_list_0based_index, areas_adjacency_mat)
	{
	defaults='
	areas_adjacency_fn = "/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/Psychotria_M3strat/BGB/areas_adjacency.txt"
	moref(areas_adjacency_fn)
	
	tmpinputs = NULL
	tmpinputs$areas_adjacency_fn = areas_adjacency_fn
	list_of_areas_adjacency_mats = read_areas_adjacency_fn(inputs=tmpinputs)
	list_of_areas_adjacency_mats
	areas_adjacency_mat = list_of_areas_adjacency_mats[[1]]
	areas_adjacency_mat
	
	areas = c("K", "O", "M", "H")
	states_list_0based_index = rcpp_areas_list_to_states_list(areas)
	states_list_0based_index
	
	prune_states_list_by_adjacency(states_list_0based_index, list_of_areas_adjacency_mats[[1]])
	prune_states_list_by_adjacency(states_list_0based_index, list_of_areas_adjacency_mats[[2]])
	prune_states_list_by_adjacency(states_list_0based_index, list_of_areas_adjacency_mats[[3]])
	prune_states_list_by_adjacency(states_list_0based_index, list_of_areas_adjacency_mats[[4]])
	

	areas_adjacency_Sara = matrix(data=c(1,1,0,0,1,1,1,0,0,1,1,1,0,0,1,1), ncol=4, byrow=TRUE)
	areas_adjacency_Sara
	
	prune_states_list(states_list_0based_index, areas_adjacency_Sara)
	prune_states_list_by_adjacency(states_list_0based_index, areas_adjacency_Sara)
	
	
	# Example from BioGeoBEARS Google Group
	# https://groups.google.com/forum/#!searchin/biogeobears/areas_allowed/biogeobears/7e7U9NtPTBU/4MyxckASMksJ
	areas = c("A", "B", "C", "D", "E", "F")
	states_list_0based_index = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
	states_list_0based_index
	areas_adjacency_mat = matrix(data=
		  c(1, 1, 1, 0, 1, 0,
			1, 1, 1, 0, 1, 1,
			1, 1, 1, 1, 1, 0,
			0, 0, 1, 1, 1, 0,
			1, 1, 1, 1, 1, 1,
			0, 1, 0, 0, 1, 1), nrow=6, ncol=6, byrow=TRUE)

	prune_states_list_by_adjacency(states_list_0based_index, areas_adjacency_mat)
	
	' # END defaults
	
	# Set up list of TRUEs (default is to keep all states)
	numstates = length(states_list_0based_index)
	states_to_keep = rep(TRUE, times=numstates)
	states_to_keep

	# Delete states that are not adjacent in the adjacency matrix
	# Skip the first state, if it's null
	if ( is.na(states_list_0based_index[[1]]) || states_list_0based_index[[1]] == "_" )
		{
		start_i = 2
		} else {
		start_i = 1
		}	

	for (i in start_i:length(states_list_0based_index))
		{	
		areas_in_this_state_range_0based = states_list_0based_index[[i]]
		areas_in_this_state_range_1based = areas_in_this_state_range_0based + 1
		adjacent_if_1s = areas_adjacency_mat[areas_in_this_state_range_1based, areas_in_this_state_range_1based]
		
		# Convert to a MATRIX so that length and sum can be EQUAL
		adjacent_if_1s = as.matrix(adjacent_if_1s)
		#print(adjacent_if_1s)
	
		if ( length(adjacent_if_1s) == sum(adjacent_if_1s) )
			{
			states_to_keep[i] = TRUE
			} else {
			states_to_keep[i] = FALSE
			}
		}

	# Modify the states_list_0based_index
	states_list_0based_index
	states_list_0based_index[states_to_keep]

	states_list_0based_index = states_list_0based_index[states_to_keep]
	
	# Put the new states_list_0based_index into the BioGeoBEARS_run_object
	states_list_0based_index_new = states_list_0based_index


	# Print to show it
	states_list_0based_index_new

	return(states_list_0based_index_new)
	}





#######################################################
# readfiles_BioGeoBEARS_run
#######################################################
#' Read in the extra input files, if any
#' 
#' This function reads input files for stratification, constraints, and detection, i.e.,
#' everything except the tree and geography files. E.g., \code{areas_allowed_fn} file is 
#' just a list of matrices, separated by blank lines,
#' from youngest to oldest. 
#' 
#' @param inputs The inputs list
#' @return \code{inputs} The modified inputs list
#' @export
#' @seealso \code{\link{define_BioGeoBEARS_run}}, \code{\link{read_times_fn}}, 
#' \code{\link{read_distances_fn}}, \code{\link{read_dispersal_multipliers_fn}}, 
#' \code{\link{read_area_of_areas_fn}}, \code{\link{read_areas_allowed_fn}}, \code{\link{read_areas_adjacency_fn}}, 
#' \code{\link{read_detections}}, \code{\link{read_controls}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @examples
#' test=1
readfiles_BioGeoBEARS_run <- function(inputs)
	{
	# Read in all the specified input files
	if (is.character(inputs$timesfn))
		{
		inputs$timeperiods = read_times_fn(inputs)
		}
	
	if (is.character(inputs$distsfn))
		{
		inputs$list_of_distances_mats = read_distances_fn(inputs)
		}

	if (is.character(inputs$envdistsfn))
		{
		inputs$list_of_envdistances_mats = read_envdistances_fn(inputs=NULL, envdistsfn=inputs$envdistsfn)
		}

	if (is.character(inputs$dispersal_multipliers_fn))
		{
		inputs$list_of_dispersal_multipliers_mats = read_dispersal_multipliers_fn(inputs)
		
		
# 		Error check
# 		for (i in 1:length(list_of_dispersal_multipliers_mats))
# 			{
# 		inputs$list_of_dispersal_multipliers_mats
# 			}
	
		}	
		
	if (is.character(inputs$area_of_areas_fn))
		{
		inputs$list_of_area_of_areas = read_area_of_areas_fn(inputs)
		}	
				
	if (is.character(inputs$areas_allowed_fn))
		{
		inputs$list_of_areas_allowed_mats = read_areas_allowed_fn(inputs)
		}

	if (is.character(inputs$areas_adjacency_fn))
		{
		inputs$list_of_areas_adjacency_mats = read_areas_adjacency_fn(inputs)
		}
	
	
	# If the detection model is turned on:
	if (inputs$use_detection_model == TRUE)
		{
		phy = check_trfn(trfn=inputs$trfn)
		
		if (is.character(inputs$detects_fn))
			{
			inputs$detects_df = read_detections(inputs$detects_fn, OTUnames=NULL, areanames=NULL, tmpskip=0, phy=phy)
			}
		if (is.character(inputs$detects_fn))
			{
			inputs$controls_df = read_controls(inputs$controls_fn, OTUnames=NULL, areanames=NULL, tmpskip=0, phy=phy)
			}
		}
		
	return(inputs)
	}






# Check a tree file --  uploaded 2017-11-02
check_trfn <- function(trfn)
	{
	
	# Load the tree
	# 2017-08-04: Added check for NEXUS input tree
	tmptr = try(expr=read.tree(trfn))
	
	if (class(tmptr) == "try-error")
		{
		stoptxt = paste0("STOP ERROR in check_BioGeoBEARS_run() or check_trfn(): There was an error in reading the tree file. You specified the trfn (TRee FileName) as '", inputs$trfn, "'.  Options to try:\n\n1. Read the error message to see if you can figure out what went wrong.\n\n.2. Try typing, in the R window, 'read.tree(trfn)', and see if this works.\n\n")
		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
		
		# Check for a NEXUS file
		error_txt = 'Error in if (tp[3] != "")'
		TF = grepl(pattern=error_txt, x=tmptr) == TRUE
		if (TF == TRUE)
			{
			txt = paste0("NOTE: Your error message begins with '", error_txt, "'.\n\nThis error typically is a result of trying to load a NEXUS phylogeny file, when what you need is a Newick phylogeny file. Some advice:\n\n(a) Open the file in a plain-text editor, and see what kind of file it is.\n\n(b) If it is a NEXUS-format phylogeny, read it in with FigTree or read.nexus() in the R package 'ape'.\n\n(c) Then save out to a Newick file (with FigTree) or with write.tree() (with the 'ape' package, in R).\n\n(d) Then use that Newick file as the input for BioGeoBEARS.")
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			} # END if (TF == TRUE)
		
		return(tmptr)
		} # END if (class(try_result) == "try-error")
	
	# Exit
	return(tmptr)
	} # END check_trfn <- function(trfn)
	





#######################################################
# check_BioGeoBEARS_run
#######################################################
#' Check the inputs for various problems
#' 
#' Numerous subtle mistakes in the input files for a BioGeoBEARS run can cause the run 
#' to crash.  As I come across these, I am putting in error checks for them.
#'
#' Some include: 
#' 
#' - Trees with negative branchlengths (as produced sometimes by e.g. BEAST MCC consensus trees
#' (MCC = majority clade consensus).  These trees are always fully resolved, but the median node
#' heights can sometimes be behind the node position in the tree.  Users should fix this manually,
#' pathological results or crashes will result otherwise.
#'
#' - Trees with polytomies.  \code{BioGeoBEARS} (and \code{LAGRANGE}, and \code{DIVA}) assume a model where lineages 
#' bifurcate, and never multifurcate.  Users can convert multifurcating trees to bifurcating trees
#' with \code{APE}'s \code{\link[ape]{multi2di}} (they will have to decide what branchlength to use for the new branches; it should be small, 
#' but bigger than the minimum branchlength used to identify fossils hooks (as hooks are considered to be
#' anagenetic members of a lineage, and thus are connected to the tree without a cladogenesis event invoked).  
#' Users can then run their analysis several times on differently-resolved trees.  
#' 
#' NOTE: After the above correction, users may wish to correct the tip branchlengths
#' (or make some other adjustment) so that all the tips are at age 0 my before present, as in an ultrametric tree.
#' (However, note that trees with fossil tips are not ultrametric according to APE's \code{\link[ape]{is.ultrametric}}, 
#' even though they are time-scaled.  To make living (nonfossil) tips line up to zero, see \code{\link{average_tr_tips}} or the 
#' (different!) .  
#' They should be used with care.  Alternately,
#' a small amount of error in tip heights will make very little difference in the likelihood calculations (e.g. if some
#' tips are 0.1 my too high, but the tree spans 200 my), which would be an argument for not requiring perfection after the
#' (crucial) corrections of negative branchlengths, zero-branchlengths, and polytomies have been made.
#'
#' - Check for an absurdly large number of states.  I've set the limit at 2500 (it starts getting slow around 200), users can
#' override with \code{allow_huge_ranges=TRUE}.
#' 
#' - Geography tipranges files should have same number of area labels as columns.
#' 
#' - Geography tipranges files should have same number of taxa as the tree, and with the (exact!!) same names.  This can be 
#' the source of many headaches, as different programs (\code{Mesquite}, etc.) treat spaces, periods, etc. in different ways, and 
#' re-write tipnames with/without quotes, underscores, etc.; and in my experience, my biologist colleagues find it very difficult
#' to guarantee that the tipnames in their tree and their data tables will match exactly.  The SAFEST approach is to NEVER use
#' these characters in tipnames or table names: space, comma, semicolon, dashes, parentheses, brackets, apostrophes or quote marks, or periods.  
#' Use ONLY letters, numbers, and underscores (_).  When plotting trees, \code{APE} automatically reads underscores as spaces, which
#' is nice for display.
#'
#' - There must be the same or more timeperiods than the other stratified items (distances matrices, etc.)
#' 
#' - The last time in the timeperiods file must be older than the root age of the tree
#'
#' @param inputs The inputs list (typically a BioGeoBEARS_run_object)
#' @param allow_huge_ranges Default FALSE, which will stop the run if there are more than 2500 states. If TRUE, this will just print a warning, and continue, at which point
#' you will wait for weeks or forever for the analysis to finish.  See \code{\link[cladoRcpp]{cladoRcpp}}'s \code{\link[cladoRcpp]{numstates_from_numareas}} function 
#' to calculate the size of the state space ahead of time, and links therein to see how the number of states scales with areas (2^number of areas, in an 
#' unconstrained analysis), how the size of the transition matrix you will be exponentiating scales (size = numstates * numstates), and the size of the 
#' ancestor/left-descendant/right-descendant cladogenesis matrix scales (numstates * numstates * numstates).  At 2500 states, this is 2500^3 = 1.5625e+10 
#' combinations of ancestor/left/right to check at every cladogenesis event, although \code{\link[cladoRcpp]{cladoRcpp}}'s tricks speed this up substantially.
#' @param allow_null_range_tips If FALSE (default), the function will call stop() if there is a null range in
#' the geography file (all zeros, e.g.: 0000).  However, since the null range is in the state space of DEC and similar
#' models, advanced users might wish to have the option (perhaps when a fossil lineage goes extinct), if so, set to TRUE.
#' @return \code{TRUE} if no errors found; otherwise a stop() is called.
#' @export
#' @seealso \code{\link{average_tr_tips}}, 
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @examples
#' test=1
check_BioGeoBEARS_run <- function(inputs, allow_huge_ranges=FALSE, allow_null_range_tips=NULL)
	{
	runjunk='
	inputs = BioGeoBEARS_run_object
	allow_huge_ranges=FALSE
	allow_null_range_tips=FALSE
	'
	
	if (is.null(inputs$allow_null_tips))
		{
		inputs$allow_null_tips = FALSE
		}
	
	if (is.null(allow_null_range_tips) == FALSE)
		{
		if (allow_null_range_tips == FALSE)
			{
			inputs$allow_null_tips = FALSE
			}
		if (allow_null_range_tips == TRUE)
			{
			inputs$allow_null_tips = TRUE
			}
		}
	
	# Check detection model setup
	if (inputs$use_detection_model == TRUE)
		{
		if ( (is.null(inputs$detects_fn) == TRUE) || (is.null(inputs$controls_fn) == TRUE) )
			{
			stoptxt = paste0("STOP ERROR in check_BioGeoBEARS_run(): You set inputs$use_detection_model == TRUE, but inputs$detects_fn and/or inputs$controls_fn is NULL. Using the detection model requires a detections counts file (inputs$detects_fn) and a taphonomic control counts file (inputs$controls_fn).\n\nBackground on finding files in R:\n\nCheck your inputs and your working directory (wd). Use commands like:\n - 'getwd()' to get your current R working directory (wd)\n - '?setwd' to see how to set your R working directory\n - 'list.files()' to list the files in your current R working directory\n - and use 'system('open .')' to open your file browsing program from the R command line (this works on macs, at least).")
			cat("\n\n")
			cat(stoptxt)
			cat("\n\n")
			stop(stoptxt)
			} # END if ( (is.null(inputs$detects_fn) == TRUE) || (is.null(inputs$controls_fn) == TRUE) )
		} # END if (inputs$use_detection_model == TRUE)
	
	
	
	# Load the tree
	# 2017-08-04: Added check for NEXUS input tree
	# tmptr = read.tree(inputs$trfn
	tmptr = check_trfn(trfn=inputs$trfn)
	
	
	# Make sure it exists
	if (exists("tmptr") == FALSE)
		{
		stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR in inputs: no readable Newick tree at:\n", inputs$trfn, "\n", sep="")
		cat(stoptxt)
		stop(stoptxt)
		}

	
	# Check that all tipnames are unique
	tipnames = tmptr$tip.label
	uniq_tipnames = unique(tipnames)
	if (length(uniq_tipnames) != length(tipnames))
		{
		stoptxt = paste("\n\ncheck_BioGeoBEARS_run() says: FATAL ERROR in inputs: your tree has non-unique tipnames. Make all tipnames unique in both your tree file and geography file.  Current tipnames are listed below:\n\n", sep="")
		
		for (i in 1:length(uniq_tipnames))
			{
			TF = uniq_tipnames[i] %in% tipnames
			if (sum(TF) > 1)
				{
				cat(uniq_tipnames[i])
				cat("\n")
				} # END if (sum(TF) > 1)
			} # END for (i in 1:length(uniq_tipnames))
		
		stop(stoptxt)
		} # END if (length(uniq_tipnames) != length(tipnames))


	# Check that all tipnames have no spaces
	TF = grepl(" ", tipnames)
	if (sum(TF) > 0)
		{
		stoptxt = paste("\n\ncheck_BioGeoBEARS_run() says: FATAL ERROR in inputs: your tree has tipnames with spaces. Take these out of your Newick file and re-run. Current tipnames are listed below.\n\n", sep="")
		
		print(tipnames[TF])
		stop(stoptxt)
		} # END if (sum(TF) > 0)



	# Check that all tipnames have no single-quotes (')
	TF = grepl("'", tipnames)
	if (sum(TF) > 0)
		{
		stoptxt = paste("\n\ncheck_BioGeoBEARS_run() says: FATAL ERROR in inputs: your tree has tipnames with apostrophes ('). Take these out of your Newick file and re-run. Current tipnames are listed below.\n\n", sep="")
		
		print(tipnames[TF])
		stop(stoptxt)
		} # END if (sum(TF) > 0)


	# Check for negative branchlengths
	brlen_equal_below_0_TF = tmptr$edge.length <= 0
	if (sum(brlen_equal_below_0_TF) > 0)
		{
		tmptr_table = prt(tmptr, printflag=FALSE)
		rows_w_BL0_TF = tmptr_table$edge.length <= 0
		rows_w_BL0_TF[is.na(rows_w_BL0_TF)] = FALSE
		
		nodenums = tmptr_table$node[rows_w_BL0_TF]
		branchlengths = tmptr_table$edge.length[rows_w_BL0_TF]
		edge_nums = tmptr_table$parent_br[rows_w_BL0_TF]
		tmptable = cbind(nodenums, branchlengths, edge_nums)
		tmptable = as.data.frame(tmptable, stringsAsFactors=FALSE)
		
		tmptxt = paste0(nodenums, collapse=",")
		tmptxt2 = paste0(edge_nums, collapse=",")
		
		stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR in check_BioGeoBEARS_run(): the input tree has branchlengths <= 0, at these nodes:\n\n", tmptxt, "\n\n...And, at these edge numbers:\n\n", tmptxt2, "\n\n")
		cat(stoptxt)
		
		print(tmptable)
		 
		cat("\nThis can sometimes happen in e.g. MCC (majority clade consensus) trees output by BEAST's TreeAnnotator.\nYou must fix the Newick file. See ?check_BioGeoBEARS_run, and PhyloWiki, for comments. One option is to try impose_min_brlen() and then write the tree to a file.\n", sep="")
		cat(stoptxt)
		stop(stoptxt)
		}

	# Check for polytomies
	if (is.binary(tmptr) == FALSE)
		{
		stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR in inputs: Your tree not bifurcating, i.e. is.binary(tmptr) returns FALSE.\n", 
		"\nYou must fix the Newick file. APE's multi2di() function is an option.  See ?check_BioGeoBEARS_run for comments.\n", sep="")
		cat(stoptxt)
		stop(stoptxt)
		}

	# Check for singletons
	if (has.singles(tmptr) == FALSE)
		{
		stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR in inputs: Your tree not bifurcating, because it as singletons (direct ancestor) nodes. I.e. has.singles(tmptr) returns FALSE.\n", 
		"\nYou must fix the Newick file. APE's collapse.singles() function is an option.  See ?check_BioGeoBEARS_run for comments.\n", sep="")
		cat(stoptxt)
		stop(stoptxt)
		}



	# Check for rooted tree
	if (is.rooted(tmptr) == FALSE)
		{
		stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR in inputs: Your tree is not rooted, i.e. is.rooted(tmptr) returns FALSE.\n", 
		"\nYou must fix the Newick file. APE's root() function is an option, or e.g. re-rooting by hand in FigTree.  However, you will want to make sure that all your tips still come up to the present (assuming you have a typical molecular tree, i.e. no fossils). \n", sep="")
		cat(stoptxt)
		stop(stoptxt)
		}
	
	
	# Load the tipranges file
	if (inputs$use_detection_model == FALSE)
		{
		tipranges = getranges_from_LagrangePHYLIP(inputs$geogfn)
		}
	if (inputs$use_detection_model == TRUE)
		{
		if (is.character(inputs$detects_fn))
			{
			inputs$detects_df = read_detections(inputs$detects_fn, OTUnames=NULL, areanames=NULL, tmpskip=0, phy=tmptr)
			tmp_blah = as.matrix(inputs$detects_df)
			tmp_blah[isblank_TF(tmp_blah)] = 0
			tmp_blah[tmp_blah > 0] = 1
			tmp_input = adf2(tmp_blah)
			tipranges_object = define_tipranges_object(tmpdf=tmp_input)
			tipranges_object@df = adf2(tipranges_object@df)
			rownames(tipranges_object@df) = rownames(tmp_blah)
			tipranges = tipranges_object
			} else {
			txt = paste0("STOP ERROR in bears_optim_run(): you set use_detection_model=TRUE, so an input file is required for input 'detects_fn'. This is required to set up the tipranges internally.")
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			stop(txt)
			}# END if (is.character(inputs$detects_fn))
		} # END if (inputs$use_detection_model == TRUE)



	# Make sure tipranges now exists
	if (exists("tipranges") == FALSE)
		{
		stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR in inputs: no readable tipranges text file at inputs$geogfn = '", inputs$geogfn, "', nor at inputs$detects_fn = '", inputs$detects_fn,  "'. Check your inputs and your working directory (wd). Use commands like:\n - 'getwd()' to get your current R working directory (wd)\n - '?setwd' to see how to set your R working directory\n - 'list.files()' to list the files in your current R working directory\n - and use 'system('open .')' to open your file browsing program from the R command line (this works on macs, at least).\n", sep="")
		cat(stoptxt)
		stop(stoptxt)
		}

	# Read in all the specified input files
	tipranges_colnames_TF = is.na(colnames(tipranges@df))
	if (sum(tipranges_colnames_TF) > 0)
		{
		catstr = paste(colnames(tipranges@df), collapse=" ", sep="")
		stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR in inputs: tipranges area names (columns) have NAs:\n", catstr, 
		"\nThis probably means your input tipranges file is missing ", sum(tipranges_colnames_TF), " areaname(s).\n", sep="")
		cat(stoptxt)
		moref(inputs$geogfn)
		stop(stoptxt)
		}
	
	
	###################################################################
	# Check that geography and tree files have identical taxa lists
	tipnames = sort(tmptr$tip.label)
	geogtaxa = sort(rownames(tipranges@df))
	###################################################################
	# Make a table of the matches and mis-matches	
	match_mismatch_table_df = compare_two_name_lists(names1=geogtaxa, names2=tipnames, listdesc1="geography file", listdesc2="Newick tree file", list_or_file_txt="files")
	
	# Check that tipranges taxa names match tree taxa names
	tipnames_in_geogfile_TF = tipnames %in% geogtaxa
	if ((sum(tipnames_in_geogfile_TF) == length(tipnames_in_geogfile_TF)) == FALSE)
		{
		tipnames_NOT_in_geogfile_TF = tipnames_in_geogfile_TF == FALSE
		numtips_not_in_geogfn = sum(tipnames_NOT_in_geogfile_TF)
		
		tips_not_in_geogfile = tipnames[tipnames_NOT_in_geogfile_TF]
		tips_not_in_geogfile_txt = paste(tips_not_in_geogfile, collapse=", ", sep="")
		
		if (length(geogtaxa) == length(tipnames))
			{
			match_TF = geogtaxa == tipnames
			match_TF_txt = paste(match_TF, collapse=" ", sep="")
			} else {
			match_TF_txt = paste("Cannot display TRUE/FALSE, as length(tipnames)=", length(tipnames), " & length(geogtaxa)=", length(geogtaxa), ".\n", sep="")
			}
		stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR in input files: Your geography file has a list of species, and your newick file has a list of species.  These lists have to match *exactly*.  This error message is saying that you have one or more species names that are missing, mis-spelled, differing due to underscores vs. spaces, contain inappropriate characters like apostrophes, etc. Error message: ", numtips_not_in_geogfn, " tree tip(s) are not in the geographic ranges file.  These are listed.\n",
		tips_not_in_geogfile_txt, "\n",
		"TRUE/FALSE between sort(tmptr$tip.label)==sort(rownames(tipranges@df)):\n",
		match_TF_txt, "\n", sep="")
		cat(stoptxt)
		
		cat("\n\nList of matches and mis-matches between the geography file and Newick tree file:\n\n")
		print(match_mismatch_table_df)
		
		cat("\n\nNOTE: Computers are very literal. Things like space vs. '_' vs. '.' matter, as do tabs versus spaces, and multiple spaces or tabs versus single ones. Also, you should NOT have ' or similar characters in either your Newick tree file or your geography file.\n\nTo edit these files, view them in a REAL plain-text editor, not in Word or whatever.\n\nInformation about plain-text editors: http://phylo.wikidot.com/biogeobears#texteditors .\n\nInformation about BioGeoBEARS file formats, with examples: http://phylo.wikidot.com/biogeobears#links_to_files .\n\n")
		
		stop(stoptxt)
		}
	

	# Check that tree taxa names match tipranges taxa names
	tipnames = sort(tmptr$tip.label)
	geogtaxa = sort(rownames(tipranges@df))
	
	geogtaxa_in_treetips_TF = geogtaxa %in% tipnames
	if ((sum(geogtaxa_in_treetips_TF) == length(geogtaxa_in_treetips_TF)) == FALSE)
		{
		geogtaxa_NOT_in_treetips_TF = geogtaxa_in_treetips_TF == FALSE
		num_geogtaxa_not_in_treetips = sum(geogtaxa_NOT_in_treetips_TF)
		
		geogtaxa_not_in_treetips = geogtaxa[geogtaxa_NOT_in_treetips_TF]
		geogtaxa_not_in_treetips_txt = paste(geogtaxa_not_in_treetips, collapse=", ", sep="")
		
		if (length(geogtaxa) == length(tipnames))
			{
			match_TF = geogtaxa == tipnames
			match_TF_txt = paste(match_TF, collapse=" ", sep="")
			} else {
			match_TF_txt = paste("Cannot display TRUE/FALSE, as length(tipnames)=", length(tipnames), " & length(geogtaxa)=", length(geogtaxa), ".\n", sep="")
			}
		stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR in inputs: ", num_geogtaxa_not_in_treetips, " taxa in the geography file are not in the tree tips.  These are listed.  This error message is saying you have one or more species names that are missing, mis-spelled, differing due to underscores vs. spaces, contain inappropriate characters like apostrophes, etc.\n",
		geogtaxa_not_in_treetips_txt, "\n",
		"TRUE/FALSE between sort(tmptr$tip.label)==sort(rownames(tipranges@df)):\n",
		match_TF_txt, "\n", sep="")
		
		cat("\n\nList of matches and mis-matches between the geography file and Newick tree file:\n\n")
		print(match_mismatch_table_df)
		
		cat(stoptxt)

		cat("\n\nNOTE: Computers are very literal. Things like space vs. '_' vs. '.' matter, as do tabs versus spaces, and multiple spaces or tabs versus single ones. Also, you should NOT have ' or similar characters in either your Newick tree file or your geography file.\n\nTo edit these files, view them in a REAL plain-text editor, not in Word or whatever.\n\nInformation about plain-text editors: http://phylo.wikidot.com/biogeobears#texteditors .\n\nInformation about BioGeoBEARS file formats, with examples: http://phylo.wikidot.com/biogeobears#links_to_files .\n\n")
	
		stop(stoptxt)
		}
	
	# Convert "?" to zero for purposes of checking for ranges too large
	tmp_tipranges = tipranges@df
	tmp_tipranges[tmp_tipranges == "?"] = 0
	tmp_tipranges = dfnums_to_numeric(tmp_tipranges)
	
	# Check that no tips have larger ranges than allowed by max, if there is a inputs$max_range_size specified
	if (!is.na(inputs$max_range_size))
		{

		
		max_tipsize = max(rowSums(tmp_tipranges))
		if (max_tipsize > inputs$max_range_size)
			{
			tips_too_big_TF = rowSums(tmp_tipranges) > inputs$max_range_size
			tipranges_too_big = tipranges@df[tips_too_big_TF, ]
			
			stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR in inputs: max_tipsize=", max_tipsize, " > inputs$max_range_size=", inputs$max_range_size, ". Examples:\n", sep="")
			cat(stoptxt)
			print(tipranges_too_big)
			cat("\n")
			
			stop(stoptxt)
			}
		}



	# Check for OTUs with ZERO (0) ranges coded
	rangesize_is_ZERO_TF = rowSums(tmp_tipranges) == 0
	ranges_have_NAs_TF = is.na(unlist(tipranges@df))
		
	if ((inputs$allow_null_tips==FALSE) && ( (sum(rangesize_is_ZERO_TF) > 0) || (sum(ranges_have_NAs_TF) > 0) ))
		{
		stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR in inputs: your tipranges file has NAs and/or tips with rangesize 0!\n See tipranges printed below:\n\n", sep="")
		cat(stoptxt)
		print(tipranges@df)
		
		stoptxt2 = paste("\nAnd your problematic tipranges rows are:\n", sep="")
		cat(stoptxt2)
		
		if (sum(rangesize_is_ZERO_TF) > 0)
			{
			TF = rangesize_is_ZERO_TF
			print(tipranges@df[TF, ])
			}
		
		if (sum(ranges_have_NAs_TF) > 0)
			{
			TF = rep(FALSE, nrow(tipranges@df))
			for (i in 1:nrow(tipranges@df))
				{
				tmprow = tipranges@df[i,]
				if (sum(is.na(tmprow)) > 0)
					{
					TF[i] = TRUE
					}
				}
			print(tipranges@df[TF,])
			}
		
		stop(stoptxt)

		}


	# Check for absurdly huge states_list
	tmp_numareas = ncol(tipranges@df)
	if (!is.na(inputs$max_range_size))
		{
		max_tipsize = max(rowSums(tmp_tipranges))
		} else {
		max_tipsize = tmp_numareas
		}
	tmp_numstates1 = length(inputs$states_list)
	if (length(tmp_numstates1) == 0)
		{
		tmp_numstates1 = numstates_from_numareas(numareas=tmp_numareas, maxareas=max_tipsize, include_null_range=TRUE)
		}
	if (tmp_numstates1 > 2500)
		{
		stoptxt = paste("\ncheck_BioGeoBEARS_run() says: Your setup has ", tmp_numstates1, " states (# states = # combinations of geographic ranges). This will be veerry slow. \n",
		"In check_BioGeoBEARS_run(), set allow_huge_ranges=TRUE to proceed, but you probably shouldn't bother.  See e.g. ?numstates_from_numareas.\n", sep="")
		cat(stoptxt)
		if (allow_huge_ranges == TRUE)
			{
			pass=1
			} else {
			stop(stoptxt)
			}
		}
	
	
	#######################################################
	# Check that all free parameters start inside the limits!
	#######################################################
	numrows = nrow(inputs$BioGeoBEARS_model_object@params_table)
	list_of_is = NULL
	for (i in 1:numrows)
		{
		if (inputs$BioGeoBEARS_model_object@params_table[i,"type"] == "free")
			{
			TF1 = inputs$BioGeoBEARS_model_object@params_table[i,"init"] >= inputs$BioGeoBEARS_model_object@params_table[i,"min"]
			TF2 = inputs$BioGeoBEARS_model_object@params_table[i,"init"] <= inputs$BioGeoBEARS_model_object@params_table[i,"max"]
			if ( (TF1 + TF2) == 2)
				{
				next()
				} else {
				list_of_is = c(list_of_is, i)
				} # END if ( (TF1 + TF2 + TF3 + TF4) == 4)
			} # END if (inputs$BioGeoBEARS_model_object@params_table[i,"type"] == "free")
		
		if (length(list_of_is) > 0)
			{
			stoptxt = paste0("check_BioGeoBEARS_run() says: STOP ERROR: ", length(list_of_is), " row(s) of the BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table are set to be 'free', but they have starting values ('init') outside of the specified min/max.\n\nFix manually, or run fix_BioGeoBEARS_params_minmax(). Printing these rows to screen....\n\n")
			cat("\n\n")
			cat(stoptxt)
			print(inputs$BioGeoBEARS_model_object@params_table[list_of_is,])
			
			stop(stoptxt)
			} # END if (length(list_of_is) > 0)
		} # END for (i in 1:numrows)
	
	
	#############################################################
	# If "x" is being used, there must be a distances file
	#############################################################
	TF1 = inputs$BioGeoBEARS_model_object@params_table["x","type"] == "free"
	TF2 = inputs$BioGeoBEARS_model_object@params_table["x","init"] != 0
	TF3 = inputs$BioGeoBEARS_model_object@params_table["x","est"] != 0

	if (TF1 || TF2 || TF3)
		{
		if (is.character(inputs$distsfn) == FALSE)
			{
			stoptxt = paste("\nFATAL ERROR: Your 'x' parameter is free or nonzero, but you have input\n",
		                "no distances file.\n\n", sep="")
		    cat(stoptxt)
		    stop(stoptxt)
		    }
		}



	# If "w" is being used (other than default 1), there must be a manual dispersal multipliers file
	TF1 = inputs$BioGeoBEARS_model_object@params_table["w","type"] == "free"
	TF2 = inputs$BioGeoBEARS_model_object@params_table["w","init"] != 1.0
	TF3 = inputs$BioGeoBEARS_model_object@params_table["w","est"] != 1.0

	if (TF1 || TF2 || TF3)
		{
		if (is.character(inputs$dispersal_multipliers_fn) == FALSE)
			{
			stoptxt = paste("\nFATAL ERROR: Your 'w' parameter is not set to '1', or is free, but you have input\n",
		                "no manual dispersal multipliers file.\n\n", sep="")
		    cat(stoptxt)
		    stop(stoptxt)
		    }
		}
		                
	# If "n" is being used, there must be an environmental distances file
	TF1 = inputs$BioGeoBEARS_model_object@params_table["n","type"] == "free"
	TF2 = inputs$BioGeoBEARS_model_object@params_table["n","init"] != 0
	TF3 = inputs$BioGeoBEARS_model_object@params_table["n","est"] != 0

	if (TF1 || TF2 || TF3)
		{
		if (is.character(inputs$envdistsfn) == FALSE)
			{
			stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: Your 'w' parameter is not set to '1', or is free, but you have input\n",
		                "no environmental distances file.\n\n", sep="")
		    cat(stoptxt)
		    stop(stoptxt)
		    }
		}
	
	

	####################################################################	
	# Check that areas for distances are in the same order
	####################################################################	
	
	# Check that matrices are square
	# Check that areas for distances are in the same order
	if (is.character(inputs$distsfn))
		{
		# Check that matrices are square
		dims = dim(inputs$list_of_distances_mats[[1]])
		if (dims[1] != dims[2])
			{
			stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the distance matrix is not square! Instead, yours has ", dims[1], " rows x ", dims[2], " columns. Check your input file. Printing the matrix below.\n")
			
			cat(stoptxt)
			print(inputs$list_of_distances_mats[[1]])
			cat("\n\n")
			stop(stoptxt)
			} # END if (dims[1] != dims[2])
		
		
		areanames_from_distsfn = colnames(inputs$list_of_distances_mats[[1]])
		areanames = names(tipranges@df)
		
		if (length(areanames) != length(areanames_from_distsfn))
			{
			stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the areas list in your geography file and in your distances file are not the same length!\n")
			
			cat(stoptxt)
			cat("\nPrinting the two lists below.\n\n")
			cat("Areas in geography file (names(tipranges@df)):\n")
			print(areanames)
			cat("Areas in distances file (colnames(inputs$list_of_distances_mats[[1]]) ):\n")
			print(areanames_from_distsfn)
			cat("\n\n")
			stop(stoptxt)
			} # END if (length(areanames) != length(areanames_from_distsfn))
		
		# If they are the same length, check that they match!
		TFs = areanames == areanames_from_distsfn
		if (all(TFs) == FALSE)
			{
			stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the area names in your geography file and in your distances file do not match! Either they are in a different order, or they are different. They must be the same, and be in the same order. In these matrices, the areas must be in the same order in BOTH rows and columns, but only the columns get area labels. For example files, see: http://phylo.wikidot.com/biogeobears#links_to_files . For advice on using a *PLAIN-TEXT* editor (not e.g. Word) to edit text files, see: http://phylo.wikidot.com/biogeobears#texteditors \n")
			
			cat(stoptxt)
			cat("\nPrinting the two lists below.\n\n")
			cat("Areas in geography file (names(tipranges@df)):\n")
			print(areanames)
			cat("Areas in distances file (colnames(inputs$list_of_distances_mats[[1]]) ):\n")
			print(areanames_from_distsfn)
			cat("\n\n")
			stop(stoptxt)
			} # END if (all(TFs) == FALSE)
		} # END if (is.character(inputs$distsfn))


	# Check that matrices are square
	# Check that areas for environmental distances are in the same order
	if (is.character(inputs$envdistsfn))
		{
		# Check that matrices are square
		dims = dim(inputs$list_of_envdistances_mats[[1]])
		if (dims[1] != dims[2])
			{
			stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the envdistances matrix is not square! Instead, yours has ", dims[1], " rows x ", dims[2], " columns. Check your input file. Printing the matrix below.\n")
			
			cat(stoptxt)
			print(inputs$list_of_envdistances_mats[[1]])
			cat("\n\n")
			stop(stoptxt)
			} # END if (dims[1] != dims[2])
		

		areanames_from_envdistsfn = colnames(inputs$list_of_envdistances_mats[[1]])
		areanames = names(tipranges@df)
		
		if (length(areanames) != length(areanames_from_envdistsfn))
			{
			stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the areas list in your geography file and in your environmental distances file are not the same length!\n")
			
			cat(stoptxt)
			cat("\nPrinting the two lists below.\n\n")
			cat("Areas in geography file (names(tipranges@df)):\n")
			print(areanames)
			cat("Areas in environmental distances file (colnames(inputs$list_of_envdistances_mats[[1]]) ):\n")
			print(areanames_from_envdistsfn)
			cat("\n\n")
			stop(stoptxt)
			} # END if (length(areanames) != length(areanames_from_envdistsfn))
		
		# If they are the same length, check that they match!
		TFs = areanames == areanames_from_envdistsfn
		if (all(TFs) == FALSE)
			{
			stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the area names in your geography file and in your environmental distances file do not match! Either they are in a different order, or they are different. They must be the same, and be in the same order. In these matrices, the areas must be in the same order in BOTH rows and columns, but only the columns get area labels. For example files, see: http://phylo.wikidot.com/biogeobears#links_to_files . For advice on using a *PLAIN-TEXT* editor (not e.g. Word) to edit text files, see: http://phylo.wikidot.com/biogeobears#texteditors \n")
			
			cat(stoptxt)
			cat("\nPrinting the two lists below.\n\n")
			cat("Areas in geography file (names(tipranges@df)):\n")
			print(areanames)
			cat("Areas in environmental distances file (colnames(inputs$list_of_envdistances_mats[[1]]) ):\n")
			print(areanames_from_envdistsfn)
			cat("\n\n")
			stop(stoptxt)
			} # END if (all(TFs) == FALSE)
		} # END if (is.character(inputs$envdistsfn))


	# Check that matrices are square
	# Check that areas for manual dispersal multipliers are in the same order
	if (is.character(inputs$dispersal_multipliers_fn))
		{
		# Check that matrices are square
		dims = dim(inputs$list_of_dispersal_multipliers_mats[[1]])
		if (dims[1] != dims[2])
			{
			stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the dispersal multipliers matrix is not square! Instead, yours has ", dims[1], " rows x ", dims[2], " columns. Check your input file. Printing the matrix below.\n")
			
			cat(stoptxt)
			print(inputs$list_of_dispersal_multipliers_mats[[1]])
			cat("\n\n")
			stop(stoptxt)
			} # END if (dims[1] != dims[2])


		areanames_from_dispersal_multipliers_fn = colnames(inputs$list_of_dispersal_multipliers_mats[[1]])
		areanames = names(tipranges@df)
		
		if (length(areanames) != length(areanames_from_dispersal_multipliers_fn))
			{
			stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the areas list in your geography file and in your manual dispersal multipliers file are not the same length!\n")
			
			cat(stoptxt)
			cat("\nPrinting the two lists below.\n\n")
			cat("Areas in geography file (names(tipranges@df)):\n")
			print(areanames)
			cat("Areas in manual dispersal multipliers file (colnames(inputs$list_of_dispersal_multipliers_mats[[1]]) ):\n")
			print(areanames_from_dispersal_multipliers_fn)
			cat("\n\n")
			stop(stoptxt)
			} # END if (length(areanames) != length(areanames_from_dispersal_multipliers_fn))
		
		# If they are the same length, check that they match!
		TFs = areanames == areanames_from_dispersal_multipliers_fn
		if (all(TFs) == FALSE)
			{
			stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the area names in your geography file and in your manual dispersal multipliers file do not match! Either they are in a different order, or they are different. They must be the same, and be in the same order. In these matrices, the areas must be in the same order in BOTH rows and columns, but only the columns get area labels. For example files, see: http://phylo.wikidot.com/biogeobears#links_to_files . For advice on using a *PLAIN-TEXT* editor (not e.g. Word) to edit text files, see: http://phylo.wikidot.com/biogeobears#texteditors \n")
			
			cat(stoptxt)
			cat("\nPrinting the two lists below.\n\n")
			cat("Areas in geography file (names(tipranges@df)):\n")
			print(areanames)
			cat("Areas in manual dispersal multipliers file (colnames(inputs$list_of_dispersal_multipliers_mats[[1]]) ):\n")
			print(areanames_from_dispersal_multipliers_fn)
			cat("\n\n")
			stop(stoptxt)
			} # END if (all(TFs) == FALSE)
		} # END if (is.character(inputs$dispersal_multipliers_fn))


	# Check that matrices are square
	# Check that areas for area of areas are in the same order
	if (is.character(inputs$area_of_areas_fn))
		{
		areanames_from_area_of_areas_fn = colnames(inputs$list_of_area_of_areas[[1]])
		areanames = names(tipranges@df)
		
		if (length(areanames) != length(areanames_from_area_of_areas_fn))
			{
			stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the areas list in your geography file and in your area of areas file are not the same length!\n")
			
			cat(stoptxt)
			cat("\nPrinting the two lists below.\n\n")
			cat("Areas in geography file (names(tipranges@df)):\n")
			print(areanames)
			cat("Areas in area of areas file (colnames(inputs$list_of_area_of_areas[[1]]) ):\n")
			print(areanames_from_area_of_areas_fn)
			cat("\n\n")
			stop(stoptxt)
			} # END if (length(areanames) != length(areanames_from_area_of_areas_fn))
		
		# If they are the same length, check that they match!
		TFs = areanames == areanames_from_area_of_areas_fn
		if (all(TFs) == FALSE)
			{
			stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the area names in your geography file and in your area of areas file do not match! Either they are in a different order, or they are different. They must be the same, and be in the same order. In these matrices, the areas must be in the same order in BOTH rows and columns, but only the columns get area labels. For example files, see: http://phylo.wikidot.com/biogeobears#links_to_files . For advice on using a *PLAIN-TEXT* editor (not e.g. Word) to edit text files, see: http://phylo.wikidot.com/biogeobears#texteditors \n")
			
			cat(stoptxt)
			cat("\nPrinting the two lists below.\n\n")
			cat("Areas in geography file (names(tipranges@df)):\n")
			print(areanames)
			cat("Areas in area of areas file (colnames(inputs$list_of_area_of_areas[[1]]) ):\n")
			print(areanames_from_area_of_areas_fn)
			cat("\n\n")
			stop(stoptxt)
			} # END if (all(TFs) == FALSE)
		} # END if (is.character(inputs$area_of_areas_fn))


	# Check that matrices are square
	# Check that areas for areas allowed are in the same order
	if (is.character(inputs$areas_allowed_fn))
		{
		# Check that matrices are square
		dims = dim(inputs$list_of_areas_allowed_mats[[1]])
		if (dims[1] != dims[2])
			{
			stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the areas allowed matrix is not square! Instead, yours has ", dims[1], " rows x ", dims[2], " columns. Check your input file. Printing the matrix below.\n")
			
			cat(stoptxt)
			print(inputs$list_of_areas_allowed_mats[[1]])
			cat("\n\n")
			stop(stoptxt)
			} # END if (dims[1] != dims[2])


		areanames_from_areas_allowed_fn = colnames(inputs$list_of_areas_allowed_mats[[1]])
		areanames = names(tipranges@df)
		
		if (length(areanames) != length(areanames_from_areas_allowed_fn))
			{
			stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the areas list in your geography file and in your areas allowed file are not the same length!\n")
			
			cat(stoptxt)
			cat("\nPrinting the two lists below.\n\n")
			cat("Areas in geography file (names(tipranges@df)):\n")
			print(areanames)
			cat("Areas in areas allowed file (colnames(inputs$list_of_areas_allowed_mats[[1]]) ):\n")
			print(areanames_from_areas_allowed_fn)
			cat("\n\n")
			stop(stoptxt)
			} # END if (length(areanames) != length(areanames_from_areas_allowed_fn))
		
		# If they are the same length, check that they match!
		TFs = areanames == areanames_from_areas_allowed_fn
		if (all(TFs) == FALSE)
			{
			stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the area names in your geography file and in your areas allowed file do not match! Either they are in a different order, or they are different. They must be the same, and be in the same order. In these matrices, the areas must be in the same order in BOTH rows and columns, but only the columns get area labels. For example files, see: http://phylo.wikidot.com/biogeobears#links_to_files . For advice on using a *PLAIN-TEXT* editor (not e.g. Word) to edit text files, see: http://phylo.wikidot.com/biogeobears#texteditors \n")
			
			cat(stoptxt)
			cat("\nPrinting the two lists below.\n\n")
			cat("Areas in geography file (names(tipranges@df)):\n")
			print(areanames)
			cat("Areas in areas allowed file (colnames(inputs$list_of_areas_allowed_mats[[1]]) ):\n")
			print(areanames_from_areas_allowed_fn)
			cat("\n\n")
			stop(stoptxt)
			} # END if (all(TFs) == FALSE)
		} # END if (is.character(inputs$areas_allowed_fn))

	# Check that matrices are square
	# Check that areas for areas adjacency are in the same order
	if (is.character(inputs$areas_adjacency_fn))
		{
		# Check that matrices are square
		dims = dim(inputs$list_of_areas_adjacency_mats[[1]])
		if (dims[1] != dims[2])
			{
			stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the areas adjacency matrix is not square! Instead, yours has ", dims[1], " rows x ", dims[2], " columns. Check your input file. Printing the matrix below.\n")
			
			cat(stoptxt)
			print(inputs$list_of_areas_adjacency_mats[[1]])
			cat("\n\n")
			stop(stoptxt)
			} # END if (dims[1] != dims[2])

		areanames_from_areas_adjacency_fn = colnames(inputs$list_of_areas_adjacency_mats[[1]])
		areanames = names(tipranges@df)
		
		if (length(areanames) != length(areanames_from_areas_adjacency_fn))
			{
			stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the areas list in your geography file and in your areas adjacency file are not the same length!\n")
			
			cat(stoptxt)
			cat("\nPrinting the two lists below.\n\n")
			cat("Areas in geography file (names(tipranges@df)):\n")
			print(areanames)
			cat("Areas in areas adjacency file (colnames(inputs$list_of_areas_adjacency_mats[[1]]) ):\n")
			print(areanames_from_areas_adjacency_fn)
			cat("\n\n")
			stop(stoptxt)
			} # END if (length(areanames) != length(areanames_from_areas_adjacency_fn))
		
		# If they are the same length, check that they match!
		TFs = areanames == areanames_from_areas_adjacency_fn
		if (all(TFs) == FALSE)
			{
			stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the area names in your geography file and in your areas adjacency file do not match! Either they are in a different order, or they are different. They must be the same, and be in the same order. In these matrices, the areas must be in the same order in BOTH rows and columns, but only the columns get area labels. For example files, see: http://phylo.wikidot.com/biogeobears#links_to_files . For advice on using a *PLAIN-TEXT* editor (not e.g. Word) to edit text files, see: http://phylo.wikidot.com/biogeobears#texteditors \n")
			
			cat(stoptxt)
			cat("\nPrinting the two lists below.\n\n")
			cat("Areas in geography file (names(tipranges@df)):\n")
			print(areanames)
			cat("Areas in areas adjacency file (colnames(inputs$list_of_areas_adjacency_mats[[1]]) ):\n")
			print(areanames_from_areas_adjacency_fn)
			cat("\n\n")
			stop(stoptxt)
			} # END if (all(TFs) == FALSE)
		} # END if (is.character(inputs$areas_adjacency_fn))







	####################################################################	
	# Check that the number of times, and the number of matrices, lines up
	####################################################################	
	
	# If there is a times file, it should have already been read in 
	# (via readfiles_BioGeoBEARS_run())
	if (is.character(inputs$timesfn))
		{
		# Extract the times
		timevals = inputs$timeperiods
		
		# Are times in order?
		if (identical(timevals, sort(timevals)) == FALSE)
			{
			stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: Your timeperiods are not in order from youngest to oldest. They need to be.\n")
			cat(stoptxt)
			stop(stoptxt)
			}
		
		
		# Check that the oldest time is older than the bottom of the tree
		trtable = prt(tmptr, printflag=FALSE)
		root_age = max(trtable$time_bp)
		oldest_time = timevals[length(timevals)]

		if (root_age >= oldest_time)
			{
			stoptxt = paste0("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: The root of your tree is age=", root_age, ". But your oldest time in your times file is only time=", oldest_time, ". The oldest time in the times file must be older than the root age.\n\n(Times in the times file represent time-bin bottoms. For example, if your first time is 10 Ma, this means the first time bin stretches from 0-10 Ma. (Ma = mega-annum = millions of years ago.)\n")
			cat(stoptxt)
			stop(stoptxt)
			}


		
		# Read distance matrices
		if (is.character(inputs$distsfn))
			{
			if (length(inputs$list_of_distances_mats) < length(inputs$timeperiods))
				{
				stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: fewer distances matrices than timeperiods.\n",
				"length(inputs$timeperiods)=", length(inputs$timeperiods), "\n",
				"length(inputs$list_of_distances_mats)=", length(inputs$list_of_distances_mats), "\n",
				sep="")
				cat(stoptxt)
				stop(stoptxt)
				}
			}

		# Read environmental distance matrices
		if (is.character(inputs$envdistsfn))
			{
			if (length(inputs$list_of_envdistances_mats) < length(inputs$timeperiods))
				{
				stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: fewer environmental distances matrices than timeperiods.\n",
				"length(inputs$timeperiods)=", length(inputs$timeperiods), "\n",
				"length(inputs$list_of_envdistances_mats)=", length(inputs$list_of_envdistances_mats), "\n",
				sep="")
				cat(stoptxt)
				stop(stoptxt)
				}
			}


		# Read dispersal multiplier matrices
		if (is.character(inputs$dispersal_multipliers_fn))
			{
			if (length(inputs$list_of_dispersal_multipliers_mats) < length(inputs$timeperiods))
				{
				stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: fewer manual dispersal multipliers matrices than timeperiods.\n",
				"length(inputs$timeperiods)=", length(inputs$timeperiods), "\n",
				"length(inputs$list_of_dispersal_multipliers_mats)=", length(inputs$list_of_dispersal_multipliers_mats), "\n",
				sep="")
				cat(stoptxt)
				stop(stoptxt)
				}
			}	
		
		# Read area-of-areas rows
		if (is.character(inputs$area_of_areas_fn))
			{
			if (length(inputs$list_of_area_of_areas) < length(inputs$timeperiods))
				{
				stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: fewer area-of-areas vectors than timeperiods.\n",
				"length(inputs$timeperiods)=", length(inputs$timeperiods), "\n",
				"length(inputs$list_of_area_of_areas)=", length(inputs$list_of_area_of_areas), "\n",
				sep="")
				cat(stoptxt)
				stop(stoptxt)
				}
			}	
				
		# Read areas-allowed matrices
		if (is.character(inputs$areas_allowed_fn))
			{
			if (length(inputs$list_of_areas_allowed_mats) != length(inputs$timeperiods))
				{
				stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: different number of area-allowed matrices than timeperiods.\n",
				"length(inputs$timeperiods)=", length(inputs$timeperiods), "\n",
				"length(inputs$list_of_areas_allowed_mats)=", length(inputs$list_of_areas_allowed_mats), "\n",
				sep="")
				cat(stoptxt)
				stop(stoptxt)
				}
			
			# Go through the areas-allowed matrices, check that they are square
			for (i in 1:length(inputs$list_of_areas_allowed_mats))
				{
				tmp_areas_allow_mat = inputs$list_of_areas_allowed_mats[[i]]
				if (nrow(tmp_areas_allow_mat) != ncol(tmp_areas_allow_mat))
					{
					errortxt = paste0("\n\ncheck_BioGeoBEARS_run() says: STOP ERROR:\n\nAreas-allowed matrices should be square, but your areas-allowed matrix #", i, " is not square:\n\nncol=", ncol(tmp_areas_allow_mat), "\nnrow=", nrow(tmp_areas_allow_mat), "\n\n")
					cat(errortxt)
				
					cat("Printing the offending areas-allowed matrix:\n\n")
					print(tmp_areas_allow_mat)
					cat("\n\n")

					
					stop(errortxt)
					} # END if (nrow(tmp_areas_allow_mat) != ncol())
				} # END for (i in 1:length(inputs$list_of_areas_allowed_mats))
			} # END if (is.character(inputs$areas_allowed_fn))


		# Read areas-adjacency matrices
		if (is.character(inputs$areas_adjacency_fn))
			{
			if (length(inputs$list_of_areas_adjacency_mats) != length(inputs$timeperiods))
				{
				stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: different number of areas-adjacency matrices than timeperiods.\n",
				"length(inputs$timeperiods)=", length(inputs$timeperiods), "\n",
				"length(inputs$list_of_areas_adjacency_mats)=", length(inputs$list_of_areas_adjacency_mats), "\n",
				sep="")
				cat(stoptxt)
				stop(stoptxt)
				}
			
			# Go through the areas-adjacency matrices, check that they are square
			for (i in 1:length(inputs$list_of_areas_adjacency_mats))
				{
				tmp_areas_adjacency_mat = inputs$list_of_areas_adjacency_mats[[i]]
				if (nrow(tmp_areas_adjacency_mat) != ncol(tmp_areas_adjacency_mat))
					{
					errortxt = paste0("\n\ncheck_BioGeoBEARS_run() says: STOP ERROR:\n\nAreas-adjacency matrices should be square, but your areas-adjacency matrix #", i, " is not square:\n\nncol=", ncol(tmp_areas_adjacency_mat), "\nnrow=", nrow(tmp_areas_adjacency_mat), "\n\n")
					cat(errortxt)
				
					cat("Printing the offending areas-adjacency matrix:\n\n")
					print(tmp_areas_adjacency_mat)
					cat("\n\n")

					
					stop(errortxt)
					} # END if (nrow(tmp_areas_adjacency_mat) != ncol())
				} # END for (i in 1:length(inputs$list_of_areas_adjacency_mats))
			} # END if (is.character(inputs$areas_adjacency_fn))




		####################################################################	
		# Check for tree sections, in a time-stratified analysis
		####################################################################	
		
		# Check for section_the_tree()
		if ("tree_sections_list" %in% names(inputs) == FALSE)
			{
			stoptxt = paste("\ncheck_BioGeoBEARS_run() says: FATAL ERROR: You have time slices, but you do not have 'inputs$tree_sections_list'.\n",
			"Run 'section_the_tree()' to add tree sections to your input BioGeoBEARS_run_object.\n\n(e.g., uncomment the line(s) starting '#section_the_tree()' in the example script!!) (You will have to do this for *each* model, e.g. 6 times for the 6 models in the PhyloWiki example script.)\n\n", sep="")
			cat(stoptxt)
			stop(stoptxt)
			}
		
		
		} # End time check	



	#######################################################
	# Check detections file
	#######################################################
	if (inputs$use_detection_model == TRUE)
		{
		# Check that the files exist and have been loaded
		if (is.character(inputs$detects_fn) == FALSE)
			{
			stoptxt = paste("\n\ncheck_BioGeoBEARS_run() says: FATAL ERROR: You have $use_detection_model set to TRUE, but you have no \n",
			"detections text file given in '$detects_fn'!\n\n", sep="")
			cat(stoptxt)
			stop(stoptxt)
			} else {
			# You have a file, so check if it is loaded
			if (class(inputs$detects_df) != "data.frame")
				{
				stoptxt = paste("\n\ncheck_BioGeoBEARS_run() says: FATAL ERROR: You have referenced a detections text file given in '$detects_fn' at:\n\n",
				inputs$detects_fn, "\n",
				"\n...but 'inputs$detects_df' is not a data.frame, perhaps empty!\n\n",
				"You should use 'readfiles_BioGeoBEARS_run()' to load these files.\n\n", 
				"Printing inputs$detects_df below:\n\n",
				"inputs$detects_df = \n\n", sep="")
				
				cat(stoptxt)
				print(inputs$detects_df)
				stop(stoptxt)				
				}
			
			# Check the order of the table rownames
			tr = check_trfn(trfn=inputs$trfn)
			tipnames = tr$tip.label
			
			table_rownames = row.names(inputs$detects_df)
			
			TF = table_rownames == tipnames
			if (sum(TF) != length(TF))
				{
				stoptxt = paste("\n\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the rownames in inputs$detects_df do not match the tip.labels in tr$tip.labels!!\n\n",
				"Printing both below:\n\n", sep="")
				
				cat(stoptxt)
				print("tipnames:")
				print(tipnames)

				print("table_rownames:")
				print(table_rownames)
				
				print("match TF:")
				print(TF)
				
				stop(stoptxt)
				}
			}
		}

	#######################################################
	# Check controls file
	#######################################################
	if (inputs$use_detection_model == TRUE)
		{
		# Check that the files exist and have been loaded
		if (is.character(inputs$controls_fn) == FALSE)
			{
			stoptxt = paste("\n\ncheck_BioGeoBEARS_run() says: FATAL ERROR: You have $use_detection_model set to TRUE, but you have no \n",
			"taphonomic controls text file given in '$controls_fn'!\n\n", sep="")
			cat(stoptxt)
			stop(stoptxt)
			} else {
			# You have a file, so check if it is loaded
			if (class(inputs$controls_df) != "data.frame")
				{
				stoptxt = paste("\n\ncheck_BioGeoBEARS_run() says: FATAL ERROR: You have referenced a taphonomic controls text file given in '$controls_fn' at:\n\n",
				inputs$controls_fn, "\n",
				"\n...but 'inputs$controls_df' is not a data.frame, perhaps empty!\n\n",
				"You should use 'readfiles_BioGeoBEARS_run()' to load these files.\n\n", 
				"Printing inputs$controls_df below:\n\n",
				"inputs$controls_df = \n\n", sep="")
				
				cat(stoptxt)
				print(inputs$controls_df)
				stop(stoptxt)				
				}
			
			# Check the order of the table rownames
			tr = check_trfn(trfn=inputs$trfn)
			tipnames = tr$tip.label
			
			table_rownames = row.names(inputs$controls_df)
			
			TF = table_rownames == tipnames
			if (sum(TF) != length(TF))
				{
				stoptxt = paste("\n\ncheck_BioGeoBEARS_run() says: FATAL ERROR: the rownames in inputs$controls_df do not match the tip.labels in tr$tip.labels!!\n\n",
				"Printing both below:\n\n", sep="")
				
				cat(stoptxt)
				print("tipnames:")
				print(tipnames)

				print("table_rownames:")
				print(table_rownames)
				
				print("match TF:")
				print(TF)
				
				stop(stoptxt)
				}
			}
		}



	return(TRUE)
	}



#' Ensure all parameters are inside the min/max limits
#'
#' This function checks for, and fixes, cases where the user has
#' input parameter values that are outside of the minimum/maximum
#' units.
#'
#' The function checks the \code{BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table},
#' comparing the \code{init} (initial values) and \code{est} (estimates) columns against the
#' columns \code{min} and \code{max}. Any parameters below the min are reset to the min, any
#' above the max are reset to the max.
#'
#' The function can take a \code{BioGeoBEARS_run_object}, \code{BioGeoBEARS_model_object}, 
#' or \code{params_table}. The first non-NULL input in that ordering is used to 
#' generate the output.
#'
#' @param BioGeoBEARS_run_object The inputs list (typically a BioGeoBEARS_run_object), 
#' derived from e.g. \code{define_BioGeoBEARS_run()}.
#' @param BioGeoBEARS_model_object The BioGeoBEARS_model object, of class 
#' \code{BioGeoBEARS_model}. E.g. from \code{BioGeoBEARS_run_object$BioGeoBEARS_model_object}
#' @param params_table The \code{params_table} from a \code{BioGeoBEARS_model_object}, e.g. 
#' from \code{BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table}
#' @return \code{BioGeoBEARS_run_object}, \code{BioGeoBEARS_model_object}, or \code{params_table}, depending
#' on the first non-\code{NULL} input.
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
fix_BioGeoBEARS_params_minmax <- function(BioGeoBEARS_run_object=NULL, BioGeoBEARS_model_object=NULL, params_table=NULL)
	{
	defaults='
	tmp_inputs = define_BioGeoBEARS_run()
	params_table = tmp_inputs$BioGeoBEARS_model_object@params_table
	fix_BioGeoBEARS_params_minmax(params_table=params_table)
	params_table["d", "init"] = 1e-15
	fix_BioGeoBEARS_params_minmax(params_table=params_table)
	'
	
	# Determine what to return -- the first non-null has priority
	if (is.null(BioGeoBEARS_run_object))
		{
		if (is.null(BioGeoBEARS_model_object))
			{
			if (is.null(params_table))
				{
				stoptxt = paste0("STOP ERROR in fix_BioGeoBEARS_params_minmax(): One of the inputs must be non-NULL.")
				cat("\n\n")
				cat(stoptxt)
				cat("\n\n")
				stop(stoptxt)
				} else {
				params_table = params_table
				returnwhat = "params_table"
				} # END if (is.null(params_table))
			} else {
			params_table = BioGeoBEARS_model_object@params_table
			returnwhat = "BioGeoBEARS_model_object"
			} # END if (is.null(BioGeoBEARS_model_object))
		} else {
		params_table = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
		returnwhat = "BioGeoBEARS_run_object"
		} # END if (is.null(BioGeoBEARS_run_object))
	
	
	#######################################################
	# Check that all free parameters start inside the limits!
	#######################################################
	numrows = nrow(params_table)
	list_of_is = NULL
	for (i in 1:numrows)
		{
		if (params_table[i,"type"] == "free")
			{
			TF1 = params_table[i,"init"] >= params_table[i,"min"]
			TF2 = params_table[i,"init"] <= params_table[i,"max"]
			TF3 = params_table[i,"est"] >= params_table[i,"min"]
			TF4 = params_table[i,"est"] <= params_table[i,"max"]
			if ( (TF1 + TF2 + TF3 + TF4) == 4)
				{
				next()
				} else {
				list_of_is = c(list_of_is, i)
				} # END if ( (TF1 + TF2 + TF3 + TF4) == 4)
			} # END if (params_table[i,"type"] == "free")
		} # END for (i in 1:numrows)
		
	if (length(list_of_is) > 0)
		{
		warningtxt = paste0("fix_BioGeoBEARS_params_minmax() is fixing limits issues for some rows of the params_table. Printing these rows to screen....\n\n")
		cat("\n\n")
		cat(warningtxt)
		print(params_table[list_of_is,])
		cat("\n\n")
		
		# Fix violations of min/max
		for (j in 1:length(list_of_is))
			{
			if (params_table[list_of_is[j],"init"] < params_table[list_of_is[j],"min"])
				{
				params_table[list_of_is[j],"init"] = params_table[list_of_is[j],"min"]
				params_table[list_of_is[j],"est"] = params_table[list_of_is[j],"min"]
				}
			if (params_table[list_of_is[j],"est"] < params_table[list_of_is[j],"min"])
				{
				params_table[list_of_is[j],"est"] = params_table[list_of_is[j],"min"]
				}
			if (params_table[list_of_is[j],"init"] > params_table[list_of_is[j],"max"])
				{
				params_table[list_of_is[j],"init"] = params_table[list_of_is[j],"max"]
				params_table[list_of_is[j],"est"] = params_table[list_of_is[j],"max"]
				}
			if (params_table[list_of_is[j],"est"] > params_table[list_of_is[j],"max"])
				{
				params_table[list_of_is[j],"est"] = params_table[list_of_is[j],"max"]
				}
			} # END for (j in 1:length(list_of_is))
		} # END if (length(list_of_is) > 0)
	
	
	# Return as appropriate
	if (returnwhat == "BioGeoBEARS_run_object")
		{
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table = params_table
		return(BioGeoBEARS_run_object)
		} # END if (returnwhat == "BioGeoBEARS_run_object")

	if (returnwhat == "BioGeoBEARS_model_object")
		{
		BioGeoBEARS_model_object@params_table = params_table
		return(BioGeoBEARS_model_object)
		} # END if (returnwhat == "BioGeoBEARS_model_object")

	if (returnwhat == "params_table")
		{
		params_table = params_table
		return(params_table)
		} # END if (returnwhat == "params_table")

	return(stop("fix_BioGeoBEARS_params_minmax(): shouldn't get here"))
	} # END fix_BioGeoBEARS_params_minmax <- function(BioGeoBEARS_run_object)





#######################################################
# tipranges
#######################################################
#' The tipranges class
#'
#' This class holds geographic range data for each tip in a phylogeny.
#'
#' Geographic range data can be read into a tipranges class object with BioGeoBEARS functions, e.g. \code{define_tipranges_object} or
#' \code{getareas_from_tipranges_object}.
#'
#' Class \code{tipranges} is an extension of the \code{\link{data.frame}} class.  It is used for holding
#' discrete geographic range data for the tips on a phylogeny. Geographic ranges are represented with bit 
#' encoding (0/1) indicating absence or presence in each possible area.
#' 
#' This is just a data.frame with:
#' rows = taxanames\cr
#' columns = area names\cr
#' cells = 0/1 representing empty/occupied\cr
#'
#'@section Slots: 
#'  \describe{
#'    \item{\code{df}:}{Data.frame of class \code{"numeric"}, containing data from df}
#'  }
#'
#' @name tipranges 
#' @rdname tipranges-class
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @seealso \code{\link{define_tipranges_object}}, \code{\link{getareas_from_tipranges_object}},
#' \code{\link[cladoRcpp]{areas_list_to_states_list_old}}, \code{\link{areas_list_to_states_list_new}},
#' \code{\link{tipranges_to_tip_condlikes_of_data_on_each_state}}
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @examples
#' tipranges_object = define_tipranges_object()
#' tipranges_object
#'
setClass("tipranges", representation(df="data.frame"),
    contains = "numeric"
)





#######################################################
# define_tipranges_object
#######################################################
#' Define a tipranges class and object
#' 
#' Class \code{tipranges} is an extension of the \code{\link{data.frame}} class.  It is used for holding
#' discrete geographic range data for the tips on a phylogeny. Geographic ranges are represented with bit 
#' encoding (0/1) indicating absence or presence in each possible area.
#' 
#' This is just a data.frame with:
#' rows = taxanames\cr
#' columns = area names\cr
#' cells = 0/1 representing empty/occupied\cr
#' 
#' @param tmpdf The user may input a \code{data.frame} holding the range data, if they like. Default is \code{NULL}, which
#' means the function will produce a temporary \code{data.frame} as an example.
#' @return \code{tipranges_object} The tipranges object, of class \code{tipranges}
#' @export
#' @seealso \code{\link{getareas_from_tipranges_object}}, \code{\link[cladoRcpp]{areas_list_to_states_list_old}},
#' \code{\link{areas_list_to_states_list_new}}, \code{\link{tipranges_to_tip_condlikes_of_data_on_each_state}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @examples
#' testval=1
#' tipranges_object = define_tipranges_object()
#' tipranges_object
#'
define_tipranges_object <- function(tmpdf=NULL)
	{
	# Define the tipranges class;
	# Now done separately
	#tipranges_class = setClass("tipranges", representation(df="data.frame"))

	# Create an instance of the class
	
	# junk:
	## A class that extends the built-in data type "numeric"
	# setClass("numWithId", representation(id = "character"),
    #     contains = "numeric")
	

	
	# Make a simple default example df if needed
#	if (is.null(tmpdf))
#		{
		# Default tip ranges data
		areanames = c("A", "B", "C")
		tipnames = c("tip1", "tip2", "tip3")
		tip1_geog = c(1, 0, 0)
		tip2_geog = c(0, 1, 0)
		tip3_geog = c(0, 0, 1)
	
		# Make the temporary data frame
		tmpdf2 = cbind(tip1_geog, tip2_geog, tip3_geog)
		tmpdf2 = adf2(data.matrix(tmpdf2))
		names(tmpdf2) = areanames
		row.names(tmpdf2) = tipnames
#		}
	
	if (is.null(tmpdf))
		{
		tipranges_object = new("tipranges", df=tmpdf2)
		} else {
		tipranges_object = new("tipranges", df=tmpdf2)
		tipranges_object@df = tmpdf
		}

	# you can get the dataframe with
	# tipranges_object@df


	return(tipranges_object)
	}


