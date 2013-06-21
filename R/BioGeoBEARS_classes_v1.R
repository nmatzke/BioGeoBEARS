
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
	
	param_data_starter = as.data.frame(matrix(data=0, nrow=1, ncol=5), stringsAsFactors=FALSE)
	names(param_data_starter) = c("type", "init", "min", "max", "est")
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
	param_data$max = 1 - minval_cladogenesis
	param_data$est = param_data$init
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table

	param_name = "ysv"
	param_data = param_data_starter
	param_data$type = "1-j"
	param_data$init = 0.5
	param_data$min = 0 + minval_cladogenesis
	param_data$max = 1 - minval_cladogenesis
	param_data$est = param_data$init
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table


	param_name = "ys"
	param_data = param_data_starter
	param_data$type = "ysv"
	param_data$init = 0.5
	param_data$min = 0 + minval_cladogenesis
	param_data$max = 1 - minval_cladogenesis
	param_data$est = param_data$init
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table

	param_name = "y"
	param_data = param_data_starter
	param_data$type = "ys"
	param_data$init = param_table[param_data$type,"init"]
	param_data$min = 0 + minval_cladogenesis
	param_data$max = 1 - minval_cladogenesis
	param_data$est = param_data$init
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table

	param_name = "s"
	param_data = param_data_starter
	param_data$type = "ys"
	param_data$init = param_table[param_data$type,"init"]
	param_data$min = 0 + minval_cladogenesis
	param_data$max = 1 - minval_cladogenesis
	param_data$est = param_data$init
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table

	param_name = "v"
	param_data = param_data_starter
	param_data$type = "ysv"
	param_data$init = param_table[param_data$type,"init"]
	param_data$min = 0 + minval_cladogenesis
	param_data$max = 1 - minval_cladogenesis
	param_data$est = param_data$init
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
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table


	# Detection probability (per observation)
	param_name = "dp"
	param_data = param_data_starter
	param_data$type = "fixed"
	param_data$init = 1
	param_data$min = 0 + minval_anagenesis
	param_data$max = 1 - minval_anagenesis
	param_data$est = param_data$init
	names(param_data) = names(param_data_starter)
	param_table = rbind(param_table, param_data)
	rownames(param_table) = (param_names = c(param_names, param_name))
	param_table


	
	return(param_table)
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
setClass(Class="BioGeoBEARS_model", representation=representation(params_table="data.frame"), contains="numeric")
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






#' new_BioGeoBEARS_model_object = calc_linked_params_BioGeoBEARS_model_object(BioGeoBEARS_model_object)
#' new_BioGeoBEARS_model_object
#' 

#######################################################
# calc_linked_params_BioGeoBEARS_model_object
#######################################################
#' Update parameters that are deterministic functions of free parameters
#' 
#' What it says.
#' 
#' @param BioGeoBEARS_model_object The BioGeoBEARS_model object, of class \code{BioGeoBEARS_model}
#' @return \code{BioGeoBEARS_model_object} (Updated version of) the BioGeoBEARS_model object, of class \code{BioGeoBEARS_model}
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
calc_linked_params_BioGeoBEARS_model_object <- function(BioGeoBEARS_model_object)
	{
	linked_TF1 = BioGeoBEARS_model_object@params_table$type != "free"
	linked_TF2 = BioGeoBEARS_model_object@params_table$type != "fixed"
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
		}
	return(BioGeoBEARS_model_object)
	}




# Define a MAXIMUM LIKELIHOOD search

#######################################################
# define_BioGeoBEARS_run
#######################################################
#' Define a MAXIMUM LIKELIHOOD search, perhaps stratified
#' 
#' Set up the inputs object for an ML search.
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
#' @param max_range_size The maximum rangesize, in number of areas.  Having a smaller maximum range size means that you can have more areas (the size of the
#' state space is greatly reduced; see \code{\link[cladoRcpp]{numstates_from_numareas}}.
#' @param states_list A list of the possible states/geographic ranges, in 0-based index form.
#' @param force_sparse Should sparse matrix exponentiation be used?
#' @param print_optim If TRUE (default), print the optimization steps as ML estimation progresses.
#' @param num_cores_to_use If >1, parallel processing will be attempted. \bold{Note:} parallel processing via \code{library (parallel)} will work in Mac command-line
#' R, but not in Mac GUI \code{R.app}.
#' @param cluster_already_open If the user wants to distribute the matrix exponentiation calculations from all the branches across a number of processors/nodes on 
#' a cluster, specify the cluster here.  E.g. \code{cluster_already_open = makeCluster(rep("localhost",num_cores_to_use), type = "SOCK")}.  Note: this will work on 
#' most platforms, including Macs running R from command line, but will NOT work on Macs running the R GUI \code{R.app}, because parallel processing functions like
#' \code{MakeCluster} from e.g. \code{library(parallel)} for some reason crash R.app.  The program runs a check for R.app and will just run on 1 node if found. 
#' @param use_optimx If TRUE, use \code{optimx} rather that \code{optim}.
#' @return \code{inputs} Inputs for ML search.
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
#' 
define_BioGeoBEARS_run <- function(abbr="default", description="defaults", BioGeoBEARS_model_object=define_BioGeoBEARS_model_object(), trfn="Psychotria_5.2.newick", geogfn="Psychotria_geog.data", timesfn=NA, distsfn=NA, dispersal_multipliers_fn=NA, area_of_areas_fn=NA, areas_allowed_fn=NA, max_range_size=NA, states_list=NA, force_sparse=NA, print_optim=TRUE, num_cores_to_use=NA, cluster_already_open=FALSE, use_optimx=TRUE)
	{
	inputs = list()
	
	inputs$abbr = abbr
	inputs$description = description
	inputs$BioGeoBEARS_model_object = BioGeoBEARS_model_object
	inputs$trfn = trfn
	inputs$geogfn = geogfn
	inputs$timesfn = timesfn
	inputs$distsfn = distsfn									# distance between areas, for dispersal ~ dist^x
	inputs$dispersal_multipliers_fn = dispersal_multipliers_fn	# hard-coded dispersal multiplier (or 0s/1s for constraints)
	inputs$area_of_areas_fn = area_of_areas_fn						# area of each areas (for extinction ~ area^u)
	inputs$areas_allowed_fn = areas_allowed_fn				# if ONLY connected areas are allowed
	inputs$max_range_size = max_range_size
	inputs$states_list = states_list
	inputs$force_sparse = force_sparse
	inputs$print_optim = print_optim
	inputs$num_cores_to_use = num_cores_to_use
	inputs$cluster_already_open = cluster_already_open
	inputs$use_optimx = use_optimx
	
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
#' areas_allowed file is just a list of distance matrices, separated by blank lines,
#' from youngest to oldest.
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
# readfiles_BioGeoBEARS_run
#######################################################
#' Read in the stratification input files, if any
#' 
#' areas_allowed file is just a list of distance matrices, separated by blank lines,
#' from youngest to oldest.
#' 
#' @param inputs The inputs list
#' @return \code{inputs} The modified inputs list
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
	inputs	
		
	return(inputs)
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
	
	# Default tip ranges data
	areanames = c("A", "B", "C")
	tipnames = c("tip1", "tip2", "tip3")
	tip1_geog = c(1, 0, 0)
	tip2_geog = c(0, 1, 0)
	tip3_geog = c(0, 0, 1)
	
	# Make a simple default example df if needed
	if (is.null(tmpdf))
		{
		# Make the temporary data frame
		tmpdf = cbind(tip1_geog, tip2_geog, tip3_geog)
		tmpdf = adf(tmpdf)
		names(tmpdf) = areanames
		row.names(tmpdf) = tipnames
		}
	
	tipranges_object = new("tipranges", df=tmpdf)

	# you can get the dataframe with
	# tipranges_object@df


	return(tipranges_object)
	}


