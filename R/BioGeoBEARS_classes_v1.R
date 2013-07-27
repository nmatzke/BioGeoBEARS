
require("ape")
require("rexpokit")
require("cladoRcpp")


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
#' What it says.
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
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
BioGeoBEARS_model_defaults <- function(minval_anagenesis=1e-15, minval_cladogenesis=1e-5, maxval=5)
	{
	defaults='
	minval_anagenesis = 1e-15
	minval_cladogenesis = 1e-5
	maxval = 5
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
	param_data$min = -10
	param_data$max = 10
	param_data$est = param_data$init
	param_data$note = "works"
	param_data$desc = "anagenesis: exponent on distance (modifies d, j, a)"
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
	param_data$max = 1
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
	param_data$min = 0 + minval_anagenesis
	param_data$max = 1 - minval_anagenesis
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
	param_data$min = 0 + minval_anagenesis
	param_data$max = 1 - minval_anagenesis
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
	param_data$min = 0 + minval_anagenesis
	param_data$max = 1 - minval_anagenesis
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
# get_perEvent_probs
#######################################################
#' Get the per-event probabilities at cladogenesis
#' 
#' At a cladogenesis event, a large number of events are possible. The simplest way to
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
#' This function calculates the per-event weight as a proportion of some total
#' weight, e.g. default 1.  If the optim result was j=0, s=1, y=1, v=1, the \code{get_perEventprobs()}
#' result would be 0, 0.333, 0.333, 0.333.
#' 
#' @param params_table The \code{params_table} from a \code{BioGeoBEARS_model_object}.
#' @param sumval Default=1.
#' @param plotwhat Default "est", use "init" to get the initial starting values instead.
#' @return \code{wts} Return the per-event weights
#' @export
#' @seealso \code{\link[base]{rbind}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' 
#' # default DEC+J model
#' BioGeoBEARS_run_object = define_BioGeoBEARS_run()
#' BioGeoBEARS_run_object$BioGeoBEARS_model_object@@params_table
#' params_table = BioGeoBEARS_run_object$BioGeoBEARS_model_object@@params_table
#' params_table
#' 
#' get_perEvent_probs(params_table)
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
#' get_perEvent_probs(params_table)
#' 
#'
get_perEvent_probs <- function(params_table, sumval=1, plotwhat="est")
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
	}






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
#' @note Go BEARS!
#' @name BioGeoBEARS_model 
#' @rdname BioGeoBEARS_model-class
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @seealso \code{\link{define_tipranges_object}}, \code{\link{getareas_from_tipranges_object}},
#' \code{\link[cladoRcpp]{areas_list_to_states_list_old}}, \code{\link{areas_list_to_states_list_new}},
#' \code{\link{tipranges_to_tip_condlikes_of_data_on_each_state}}
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
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
#' @param minval_anagenesis Minimum value above zero for d, e, a, b parameters.
#' @param minval_cladogenesis Minimum value above zero for j, v, etc.
#' @param maxval Maximum value for d, e, a
#' @return \code{BioGeoBEARS_model_object} The BioGeoBEARS_model object, of class \code{BioGeoBEARS_model}
#' @export
#' @seealso \code{\link[cladoRcpp]{areas_list_to_states_list_old}}, \code{\link{areas_list_to_states_list_new}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' BioGeoBEARS_model_object = define_BioGeoBEARS_model_object()
#' BioGeoBEARS_model_object
#' define_BioGeoBEARS_model_object()
define_BioGeoBEARS_model_object <- function(minval_anagenesis=1e-15, minval_cladogenesis=1e-5, maxval=5)
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
#' What it says.
#' 
#' @param BioGeoBEARS_model_object The BioGeoBEARS_model object, of class \code{BioGeoBEARS_model}
#' @return \code{params} parameter vector
#' @export
#' @seealso \code{\link[BioGeoBEARS]{define_BioGeoBEARS_model_object}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
BioGeoBEARS_model_object_to_init_params <- function(BioGeoBEARS_model_object)
	{
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
#' What it says.
#' 
#' @param BioGeoBEARS_model_object The BioGeoBEARS_model object, of class \code{BioGeoBEARS_model}
#' @return \code{params} parameter vector
#' @export
#' @seealso \code{\link[BioGeoBEARS]{define_BioGeoBEARS_model_object}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
BioGeoBEARS_model_object_to_est_params <- function(BioGeoBEARS_model_object)
	{
	free_params_TF = BioGeoBEARS_model_object@params_table$type == "free"
	num_free_params = sum(free_params_TF, na.rm=TRUE)
	
	params = rep(NA, num_free_params)
	
	params = BioGeoBEARS_model_object@params_table$est[free_params_TF]
	
	return(params)	
	}




#######################################################
# BioGeoBEARS_model_object_to_params_lower
#######################################################
#' Produce the lower limit on the parameters from a BioGeoBEARS model object
#' 
#' What it says.
#' 
#' @param BioGeoBEARS_model_object The BioGeoBEARS_model object, of class \code{BioGeoBEARS_model}
#' @return \code{params} parameter vector
#' @export
#' @seealso \code{\link[BioGeoBEARS]{define_BioGeoBEARS_model_object}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
BioGeoBEARS_model_object_to_params_lower <- function(BioGeoBEARS_model_object)
	{
	free_params_TF = BioGeoBEARS_model_object@params_table$type == "free"
	num_free_params = sum(free_params_TF, na.rm=TRUE)
	
	params = rep(NA, num_free_params)
	
	params = BioGeoBEARS_model_object@params_table$min[free_params_TF]
	
	return(params)	
	}


#######################################################
# BioGeoBEARS_model_object_to_params_upper
#######################################################
#' Produce the upper limit on the parameters from a BioGeoBEARS model object
#' 
#' What it says.
#' 
#' @param BioGeoBEARS_model_object The BioGeoBEARS_model object, of class \code{BioGeoBEARS_model}
#' @return \code{params} parameter vector
#' @export
#' @seealso \code{\link[BioGeoBEARS]{define_BioGeoBEARS_model_object}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
BioGeoBEARS_model_object_to_params_upper <- function(BioGeoBEARS_model_object)
	{
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
#' What it says.
#' 
#' @param BioGeoBEARS_model_object The BioGeoBEARS_model object, of class \code{BioGeoBEARS_model}
#' @param params parameter vector
#' @return \code{BioGeoBEARS_model_object} The BioGeoBEARS_model object, of class \code{BioGeoBEARS_model}
#' @export
#' @seealso \code{\link[BioGeoBEARS]{define_BioGeoBEARS_model_object}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
params_into_BioGeoBEARS_model_object <- function(BioGeoBEARS_model_object, params)
	{
	free_params_TF = BioGeoBEARS_model_object@params_table$type == "free"
	num_free_params = sum(free_params_TF, na.rm=TRUE)
	
	BioGeoBEARS_model_object@params_table$est[free_params_TF] = params
	
	return(BioGeoBEARS_model_object)
	}

#######################################################
# merge_words_nonwords
#######################################################
#' Merge lists of words and nonwords (numbers) that may be of different length
#' 
#' Utility function.
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
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
merge_words_nonwords <- function(words, nonwords)
	{
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
#' @param BioGeoBEARS_model_object The BioGeoBEARS_model object, of class \code{BioGeoBEARS_model}
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
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
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
#' Tipnames should match the names in geogfn.  See \code{\link[ape]{read.tree}} in APE for reading in phylogenetic trees. Default "Psychotria_5.2.newick"
#' @param geogfn A PHYLIP-style file with geographic range data (see \code{\link{getranges_from_LagrangePHYLIP}}) for each tipname. This is the same format
#' used by C++ LAGRANGE (\cite{SmithRee2010_CPPversion}). Default "Psychotria_geog.data"
#' @param timesfn Filename for the stratified times.
#' @param distsfn Filename for the changing distances.
#' @param dispersal_multipliers_fn Filename for the changing hard-coded dispersal multipliers
#' @param area_of_areas_fn Filename for the area of each area
#' @param areas_allowed_fn Filename for the allowed connections between areas for single-species ranges.
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
#' @param tmpwd The working directory in which the input and output files will be placed. Default is \code{\link[base]{getwd}}. This is stored 
#' mostly for future reference; users are responsible for manually navigating to the appropriate directory ahead of time, using \code{\link[base]{setwd}}.
#' @param num_cores_to_use If >1, parallel processing will be attempted. \bold{Note:} parallel processing via \code{library (parallel)} will work in Mac command-line
#' R, but not in Mac GUI \code{R.app}.
#' @param cluster_already_open If the user wants to distribute the matrix exponentiation calculations from all the branches across a number of processors/nodes on 
#' a cluster, specify the cluster here.  E.g. \code{cluster_already_open = makeCluster(rep("localhost",num_cores_to_use), type = "SOCK")}.  Note: this will work on 
#' most platforms, including Macs running R from command line, but will NOT work on Macs running the R GUI \code{R.app}, because parallel processing functions like
#' \code{MakeCluster} from e.g. \code{library(parallel)} for some reason crash R.app.  The program runs a check for R.app and will just run on 1 node if found. 
#' @param use_optimx If TRUE, use \code{optimx} rather that \code{optim}.
#' @param return_condlikes_table If \code{TRUE}, return the table of ALL conditional likelihood results, including at branch subsections
#' (only some should be used in calculating the final log-likelihood of the geography range data on the tree!)
#' @param calc_TTL_loglike_from_condlikes_table If TRUE, force making of the condlikes table, and use it to calculate the log-likelihood
#' (default=TRUE; matches LAGRANGE).
#' @param calc_ancprobs If \code{TRUE} (default), calculate and return the necessary pieces (uppass and downpass probs) for ancestral states.
#' @param fixnode If the state at a particular node is going to be fixed (e.g. for ML marginal ancestral states), give the node number.
#' @param fixlikes The state likelihoods to be used at the fixed node.  I.e. 1 for the fixed state, and 0 for the others.
#' @param speedup If \code{TRUE} (default), set the maximum number of iterations to \code{itnmax=50*(number of free parameters)}, instead of the
#' \code{optimx} default, 250.  Also set \code{optimx} \code{reltol} parameter to 0.001 (instead of the default, ~1e-8).
#' @return \code{inputs} Inputs for ML search.
#' @export
#' @seealso \code{\link{readfiles_BioGeoBEARS_run}}, \code{\link[BioGeoBEARS]{define_BioGeoBEARS_model_object}}, \code{\link[base]{setwd}}, \code{\link[base]{getwd}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
#' 
define_BioGeoBEARS_run <- function(abbr="default", description="defaults", BioGeoBEARS_model_object=define_BioGeoBEARS_model_object(), trfn="Psychotria_5.2.newick", geogfn="Psychotria_geog.data", timesfn=NA, distsfn=NA, dispersal_multipliers_fn=NA, area_of_areas_fn=NA, areas_allowed_fn=NA, detects_fn=NA, controls_fn=NA, max_range_size=NA, states_list=NULL, force_sparse=FALSE, use_detection_model=FALSE, print_optim=TRUE, num_cores_to_use=NA, cluster_already_open=FALSE, use_optimx=TRUE, return_condlikes_table=FALSE, calc_TTL_loglike_from_condlikes_table=TRUE, calc_ancprobs=TRUE, fixnode=NULL, fixlikes=NULL, speedup=TRUE, tmpwd=getwd())
	{
	inputs = list()
	
	inputs$abbr = abbr
	inputs$description = description
	inputs$BioGeoBEARS_model_object = BioGeoBEARS_model_object
	inputs$trfn = trfn
	inputs$geogfn = geogfn
	inputs$timesfn = timesfn
	inputs$distsfn = distsfn										# distance between areas, for dispersal ~ dist^x
	inputs$dispersal_multipliers_fn = dispersal_multipliers_fn		# hard-coded dispersal multiplier (or 0s/1s for constraints)
	inputs$area_of_areas_fn = area_of_areas_fn						# area of each areas (for extinction ~ area^u)
	inputs$areas_allowed_fn = areas_allowed_fn						# if ONLY connected areas are allowed
	inputs$detects_fn = detects_fn									# used only if use_detection_model==TRUE
	inputs$controls_fn = controls_fn								# used only if use_detection_model==TRUE
	inputs$max_range_size = max_range_size
	inputs$states_list = states_list
	inputs$force_sparse = force_sparse
	inputs$use_detection_model = use_detection_model
	inputs$print_optim = print_optim
	inputs$wd = tmpwd												# Store the working directory you are in
	inputs$num_cores_to_use = num_cores_to_use
	inputs$cluster_already_open = cluster_already_open
	inputs$use_optimx = use_optimx
	inputs$return_condlikes_table = return_condlikes_table
	inputs$calc_TTL_loglike_from_condlikes_table = calc_TTL_loglike_from_condlikes_table
	inputs$calc_ancprobs = calc_ancprobs
	inputs$fixnode = fixnode
	inputs$fixlikes = fixlikes
	inputs$speedup = speedup

	
	return(inputs)
	}


#######################################################
# BioGeoBEARS_run
#######################################################
#' An object of class BioGeoBEARS_run holding the model inputs
#'
#'
#'@section Slots: 
#'  \describe{
#'    \item{\code{list}:}{List of class \code{"list"}, containing inputs list from define_BioGeoBEARS_run}
#'  }
#'
#' @note Go BEARS!
#' @name BioGeoBEARS_run 
#' @rdname BioGeoBEARS_run-class
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @seealso \code{\link{define_tipranges_object}}, \code{\link{getareas_from_tipranges_object}},
#' \code{\link[cladoRcpp]{areas_list_to_states_list_old}}, \code{\link{areas_list_to_states_list_new}},
#' \code{\link{tipranges_to_tip_condlikes_of_data_on_each_state}}
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
#'
setClass(Class="BioGeoBEARS_run", representation=representation(inputs="list"), prototype=define_BioGeoBEARS_run())



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
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
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
				tmprows = rbind(tmprows, tmprow)
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
#' youngest to oldest.
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
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
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
		}
	
	timeperiods = as.numeric(readLines(timesfn))
	
	return(timeperiods)
	}




# Distances file is just a list of distance matrices, separated by blank lines,
# from youngest to oldest
#######################################################
# read_distances_fn
#######################################################
#' Read in the distances by time
#' 
#' Distances file is just a list of distance matrices, separated by blank lines,
#' from youngest to oldest.
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
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
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
			# (to ensure consistent behavior of the exponent, i.e. dist ^ -1*(x_exponent)
			tmpmat = as.matrix(tmpmat)
			diag(tmpmat) = NA
			minval = min(tmpmat, na.rm=TRUE)
 			if (minval <= 0)
 				{
 				stop("\n\nERROR: Minimum distance between regions must be > 0")
 				}
			
			list_of_distances_mats[[lnum]] = tmpmat
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
				tmprows = rbind(tmprows, tmprow)
				}
			}
		}
	list_of_distances_mats
	return(list_of_distances_mats)
	}


# area_areas file is just a list of distance matrices, separated by blank lines,
# from youngest to oldest

#######################################################
# read_area_of_areas_fn
#######################################################
#' Read in the area areas by time
#' 
#' area_areas file is just a list of distance matrices, separated by blank lines,
#' from youngest to oldest.
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
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
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
			# (to ensure consistent behavior of the exponent, i.e. dist ^ -1*(x_exponent)
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
				tmprows = rbind(tmprows, tmprow)
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
#' areas_allowed file is just a list of 1/0 matrices, separated by blank lines,
#' from youngest to oldest. 1s represent allowed combinations of areas
#' 
#' @param inputs The inputs list
#' @param areas_allowed_fn The areas-allowed filename.
#' @return \code{list_of_areas_allowed_mats} A list object
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
read_areas_allowed_fn <- function(inputs=NULL, areas_allowed_fn=NULL)
	{
	defaults='
	areas_allowed_fn = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/Psychotria_dists_stratified/Hawaii_KOMH_areas_allowed.txt"
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
				tmprows = rbind(tmprows, tmprow)
				}
			}
		}
	list_of_areas_allowed_mats
	return(list_of_areas_allowed_mats)
	}


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
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
prune_states_list <- function(states_list_0based_index, areas_allowed_mat)
	{
	defaults='
	areas_allowed_fn = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/Psychotria_dists_stratified/Hawaii_KOMH_areas_allowed.txt"
	tmpinputs = NULL
	tmpinputs$areas_allowed_fn = areas_allowed_fn
	list_of_areas_allowed_mats = read_areas_allowed_fn(inputs=tmpinputs)
	list_of_areas_allowed_mats
	areas_allowed_mat = list_of_areas_allowed_mats[[1]]
	areas_allowed_mat
	
	areas = c("K", "O", "M", "H")
	states_list_0based_index = rcpp_areas_list_to_states_list(areas)
	states_list_0based_index
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
# readfiles_BioGeoBEARS_run
#######################################################
#' Read in the extra input files, if any
#' 
#' This function reads input files for stratification, constraints, and detection, i.e.,
#' everything except the tree and geography files. E.g., \code{areas_allowed_fn} file is 
#' just a list of distance matrices, separated by blank lines,
#' from youngest to oldest.
#' 
#' @param inputs The inputs list
#' @return \code{inputs} The modified inputs list
#' @export
#' @seealso \code{\link{define_BioGeoBEARS_run}}, \code{\link{read_times_fn}}, 
#' \code{\link{read_distances_fn}}, \code{\link{read_dispersal_multipliers_fn}}, 
#' \code{\link{read_area_of_areas_fn}}, \code{\link{read_areas_allowed_fn}}, 
#' \code{\link{read_detections}}, \code{\link{read_controls}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
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

	if (is.character(inputs$dispersal_multipliers_fn))
		{
		inputs$list_of_dispersal_multipliers_mats = read_dispersal_multipliers_fn(inputs)
		}	
		
	if (is.character(inputs$area_of_areas_fn))
		{
		inputs$list_of_area_of_areas = read_area_of_areas_fn(inputs)
		}	
				
	if (is.character(inputs$areas_allowed_fn))
		{
		inputs$list_of_areas_allowed_mats = read_areas_allowed_fn(inputs)
		}
	
	
	# If the detection model is turned on:
	if (inputs$use_detection_model == TRUE)
		{
		phy = read.tree(inputs$trfn)
		
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
#' - Check for an absurdly large number of states.  I've set the limit at 500 (it starts getting slow around 200), users can
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
#' @param inputs The inputs list
#' @param allow_huge_ranges Default FALSE, which will stop the run if there are more than 500 states. If TRUE, this will just print a warning, and continue, at which point
#' you will wait for weeks or forever for the analysis to finish.  See \code{\link[cladoRcpp]{cladoRcpp}}'s \code{\link[cladoRcpp]{numstates_from_numareas}} function 
#' to calculate the size of the state space ahead of time, and links therein to see how the number of states scales with areas (2^number of areas, in an 
#' unconstrained analysis), how the size of the transition matrix you will be exponentiating scales (size = numstates * numstates), and the size of the 
#' ancestor/left-descendant/right-descendant cladogenesis matrix scales (numstates * numstates * numstates).  At 500 states, this is 500^3 = 125,000,000 
#' combinations of ancestor/left/right to check at every cladogenesis event, although \code{\link[cladoRcpp]{cladoRcpp}}'s tricks speed this up substantially.
#' @return \code{TRUE} if no errors found; otherwise a stop() is called.
#' @export
#' @seealso \code{\link{average_tr_tips}}, 
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
check_BioGeoBEARS_run <- function(inputs, allow_huge_ranges=FALSE)
	{
	# Load the tree
	tmptr = read.tree(inputs$trfn)
	
	# Make sure it exists
	if (exists("tmptr") == FALSE)
		{
		stoptxt = paste("\nFATAL ERROR in inputs: no readable Newick tree at:\n", inputs$trfn, "\n", sep="")
		cat(stoptxt)
		stop(stoptxt)
		}

	# Check for negative branchlengths
	blren_equal_below_0_TF = tmptr$edge.length <= 0
	if (sum(blren_equal_below_0_TF) > 0)
		{
		tmptxt = paste(tmptr$edge.length[blren_equal_below_0_TF], collapse=", ", sep="")
		stoptxt = paste("\nFATAL ERROR in average_tr_tips(): the output tree has branchlengths <= 0:\n", tmptxt, 
		"\nThis can sometimes happen in e.g. MCC (majority clade consensus) trees output by BEAST's TreeAnnotator.\nYou must fix the Newick file. See ?check_BioGeoBEARS_run for comments.\n", sep="")
		cat(stoptxt)
		stop(stoptxt)
		}

	# Check for polytomies
	if (is.binary.tree(tmptr) == FALSE)
		{
		stoptxt = paste("\nFATAL ERROR in inputs: tree not bifurcating, i.e. is.binary.tree(tmptr) returns FALSE.\n", 
		"\nYou must fix the Newick file. APE's multi2di() function is an option.  See ?check_BioGeoBEARS_run for comments.\n", sep="")
		cat(stoptxt)
		stop(stoptxt)
		}
	
	
	# Load the tipranges file
	tipranges = getranges_from_LagrangePHYLIP(inputs$geogfn)
	# Make sure it exists
	if (exists("tipranges") == FALSE)
		{
		stoptxt = paste("\nFATAL ERROR in inputs: no readable tipranges text file at:\n", inputs$geogfn, "\n", sep="")
		cat(stoptxt)
		stop(stoptxt)
		}

	# Read in all the specified input files
	tipranges_colnames_TF = is.na(colnames(tipranges@df))
	if (sum(tipranges_colnames_TF) > 0)
		{
		catstr = paste(colnames(tipranges@df), collapse=" ", sep="")
		stoptxt = paste("\nFATAL ERROR in inputs: tipranges area names (columns) have NAs:\n", catstr, 
		"\nThis probably means your input tipranges file is missing ", sum(tipranges_colnames_TF), " areaname(s).\n", sep="")
		cat(stoptxt)
		moref(inputs$geogfn)
		stop(stoptxt)
		}
	
	# Check that tipranges taxa names match tree taxa names
	tipnames = sort(tmptr$tip.label)
	geogtaxa = sort(rownames(tipranges@df))
	
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
		stoptxt = paste("\nFATAL ERROR in inputs: ", numtips_not_in_geogfn, " tree tips are not in the geographic ranges file.  These are:\n",
		tips_not_in_geogfile_txt, "\n",
		"TRUE/FALSE between sort(tmptr$tip.label)==sort(rownames(tipranges@df)):\n",
		match_TF_txt, "\n", sep="")
		cat(stoptxt)
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
		
		geogtaxa_not_in_treetips = geogtaxa(geogtaxa_NOT_in_treetips_TF)
		geogtaxa_not_in_treetips_txt = paste(geogtaxa_not_in_treetips, collapse=", ", sep="")
		
		if (length(geogtaxa) == length(tipnames))
			{
			match_TF = geogtaxa == tipnames
			match_TF_txt = paste(match_TF, collapse=" ", sep="")
			} else {
			match_TF_txt = paste("Cannot display TRUE/FALSE, as length(tipnames)=", length(tipnames), " & length(geogtaxa)=", length(geogtaxa), ".\n", sep="")
			}
		stoptxt = paste("\nFATAL ERROR in inputs: ", num_geogtaxa_not_in_treetips, " taxa in the geography file are not in the tree tips.  These are:\n",
		geogtaxa_not_in_treetips_txt, "\n",
		"TRUE/FALSE between sort(tmptr$tip.label)==sort(rownames(tipranges@df)):\n",
		match_TF_txt, "\n", sep="")
		cat(stoptxt)
		stop(stoptxt)
		}
	

	# Check that no tips have larger ranges than allowed by max, if there is a inputs$max_range_size specified
	if (!is.na(inputs$max_range_size))
		{
		max_tipsize = max(rowSums(tipranges@df))
		if (max_tipsize > inputs$max_range_size)
			{
			tips_too_big_TF = rowSums(tipranges@df) > inputs$max_range_size
			tipranges_too_big = tipranges@df[tips_too_big_TF, ]
			
			stoptxt = paste("\nFATAL ERROR in inputs: max_tipsize=", max_tipsize, " > inputs$max_range_size=", inputs$max_range_size, ". Examples:\n", sep="")
			for (i in 1:sum(tips_too_big_TF))
				{
				cat(tipranges_too_big[i,])
				cat("\n")
				}
			stop(stoptxt)
			}
		}



	# Check for OTUs with ZERO (0) ranges coded
	rangesize_is_ZERO_TF = (rowSums(tipranges@df)) == 0
	ranges_have_NAs_TF = is.na(unlist(tipranges@df))
	if ( (sum(rangesize_is_ZERO_TF) > 0) || (sum(ranges_have_NAs_TF) > 0) )
		{
		stoptxt = paste("\nFATAL ERROR in inputs: your tipranges file has NAs and/or tips with rangesize 0!\n See tipranges printed below:\n\n", sep="")
		cat(stoptxt)
		print(tipranges@df)
		stop(stoptxt)

		}


	# Check for absurdly huge states_list
	tmp_numareas = ncol(tipranges@df)
	if (!is.na(inputs$max_range_size))
		{
		max_tipsize = max(rowSums(tipranges@df))
		} else {
		max_tipsize = tmp_numareas
		}
	tmp_numstates1 = length(inputs$states_list)
	if (length(tmp_numstates1) == 0)
		{
		tmp_numstates1 = numstates_from_numareas(numareas=tmp_numareas, maxareas=max_tipsize, include_null_range=TRUE)
		}
	if (tmp_numstates1 > 500)
		{
		stoptxt = paste("\nYour setup has ", tmp_numstates1, " states (# states = # combinations of geographic ranges). This will be veerry slow. \n",
		"In check_BioGeoBEARS_run(), set allow_huge_ranges=TRUE to proceed, but you probably shouldn't bother.  See e.g. ?numstates_from_numareas.\n", sep="")
		cat(stoptxt)
		if (allow_huge_ranges == TRUE)
			{
			pass=1
			} else {
			stop(stoptxt)
			}
		}
	
	
	# If there is a times file, it should have already been read in 
	# (via readfiles_BioGeoBEARS_run())
	if (is.character(inputs$timesfn))
		{
		# Extract the times
		timevals = inputs$timeperiods
	
		if (is.character(inputs$distsfn))
			{
			if (length(inputs$list_of_distances_mats) < length(inputs$timeperiods))
				{
				stoptxt = paste("\nFATAL ERROR: fewer distances matrices than timeperiods.\n",
				"length(inputs$timeperiods)=", length(inputs$timeperiods), "\n",
				"length(inputs$list_of_distances_mats)=", length(inputs$list_of_distances_mats), "\n",
				sep="")
				cat(stoptxt)
				stop(stoptxt)
				}
			}

		if (is.character(inputs$dispersal_multipliers_fn))
			{
			if (length(inputs$list_of_dispersal_multipliers_mats) < length(inputs$timeperiods))
				{
				stoptxt = paste("\nFATAL ERROR: fewer manual dispersal multipliers matrices than timeperiods.\n",
				"length(inputs$timeperiods)=", length(inputs$timeperiods), "\n",
				"length(inputs$list_of_dispersal_multipliers_mats)=", length(inputs$list_of_dispersal_multipliers_mats), "\n",
				sep="")
				cat(stoptxt)
				stop(stoptxt)
				}
			}	
		
		if (is.character(inputs$area_of_areas_fn))
			{
			if (length(inputs$list_of_area_of_areas) < length(inputs$timeperiods))
				{
				stoptxt = paste("\nFATAL ERROR: fewer area-of-areas vectors than timeperiods.\n",
				"length(inputs$timeperiods)=", length(inputs$timeperiods), "\n",
				"length(inputs$list_of_area_of_areas)=", length(inputs$list_of_area_of_areas), "\n",
				sep="")
				cat(stoptxt)
				stop(stoptxt)
				}
			}	
				
		if (is.character(inputs$areas_allowed_fn))
			{
			if (length(inputs$list_of_areas_allowed_mats) < length(inputs$timeperiods))
				{
				stoptxt = paste("\nFATAL ERROR: fewer area-allowed matrices than timeperiods.\n",
				"length(inputs$timeperiods)=", length(inputs$timeperiods), "\n",
				"length(inputs$list_of_areas_allowed_mats)=", length(inputs$list_of_areas_allowed_mats), "\n",
				sep="")
				cat(stoptxt)
				stop(stoptxt)
				}
			}
		}		

	return(TRUE)
	}








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
#' @note Go BEARS!
#' @name tipranges 
#' @rdname tipranges-class
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @seealso \code{\link{define_tipranges_object}}, \code{\link{getareas_from_tipranges_object}},
#' \code{\link[cladoRcpp]{areas_list_to_states_list_old}}, \code{\link{areas_list_to_states_list_new}},
#' \code{\link{tipranges_to_tip_condlikes_of_data_on_each_state}}
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
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
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
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


