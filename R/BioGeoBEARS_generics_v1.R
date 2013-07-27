require("ape")
require("rexpokit")
require("cladoRcpp")


# source all .R files in a directory, except "compile" and "package" files

#######################################################
# sourceall
#######################################################
#' Source all .R files in a directory, except "compile" and "package" files
#' 
#' Utility function.
#' 
#' @param path The path to source
#' @param pattern Default is .R
#' @param ... Additional arguments to source
#' @return \code{path} The path that was sourced.
#' @export
#' @seealso \code{\link[base]{source}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
sourceall <- function(path=path, pattern="\\.R", ...)
	{
	tmppath = np(addslash(path))
	Rfiles = list.files(path=tmppath, pattern="\\.R", ...)
	
	# Files to remove
	Rfiles_remove_TF1 = grepl("compile", Rfiles)
	Rfiles_remove_TF2 = grepl("package", Rfiles)
	Rfiles_remove_TF = (Rfiles_remove_TF1 + Rfiles_remove_TF2) >= 1
	
	Rfiles = Rfiles[Rfiles_remove_TF == FALSE]

	cat("\nSourcing Rfiles in ", path, "...\n", sep="")

	
	for (Rfile in Rfiles)
		{
		cat("Sourcing Rfile: ", Rfile, "\n", sep="")
		fullfn = np(slashslash(paste(addslash(path), Rfile, sep="")))
		source(fullfn, chdir=TRUE, ...)
		}

	cat("\nDone sourcing Rfiles in ", path, "...\n", sep="")
	return(path)
	}




#######################################################
# printall
#######################################################
#' Print an entire table to screen
#' 
#' Utility function.  This prints a table to screen in chunks of \code{chunksize_toprint} 
#' (default=40).  This avoids the annoying situation of not being able to see the bottom 
#' of a table. Note that if you print something huge, you will be waiting for awhile (try
#' ESC or CTRL-C to cancel such an operation).
#'
#' Another option is to reset options to something like: \code{options(max.print=99999)}, but this
#' is hard to remember.  Your current setting is \code{getOption("max.print")}.
#' 
#' @param dtf The \code{\link[base]{data.frame}} to \code{\link[base]{print}}.
#' @param chunksize_toprint Number of lines to print. Default 50.
#' @param printflag For optional printing. Passed to \code{\link{prflag}}.
#' @return NULL
#' @export
#' @seealso \code{\link[base]{print}}, \code{\link{prflag}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
printall <- function(dtf, chunksize_toprint = 40, printflag=TRUE)
	{
	# Print everything in a data frame, in chunks of e.g. 50 rows
	if (nrow(dtf) <= chunksize_toprint)
		{
		prflag(dtf, printflag=printflag)
		return(dtf)
		}
	rows_toprint = seq(1, nrow(dtf), chunksize_toprint)
	
	if (printflag == TRUE)
		{
		for (i in 1 : (length(rows_toprint)-1) )
			{
			tmp_rows_toprint_start = rows_toprint[i]
			tmp_rows_toprint_end = rows_toprint[i+1]
			prflag(dtf[tmp_rows_toprint_start:tmp_rows_toprint_end, ])
			}
		
		# Then print the end
		tmp_rows_toprint_start = rows_toprint[length(rows_toprint)]
		tmp_rows_toprint_end = nrow(dtf)
		prflag(dtf[tmp_rows_toprint_start:tmp_rows_toprint_end, ])
		}	
	}





#######################################################
# prflag
#######################################################
#' Utility function to conditionally print intermediate results
#'
#' Just a handy shortcut function, allowing other functions to optionally 
#' print, depending on the value of \code{printflag}.
#' 
#' @param x What to print.
#' @param printflag If TRUE, do the printing
#' @return nothing
#' @export
#' @seealso \code{\link{get_daughters}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
prflag <- function(x, printflag=TRUE)
	{
	# A standard function to print (or not) certain variables,
	#   based on a master printflag
	# This avoids having to comment in/out various code chunks
	#   while debugging.
	if (printflag == TRUE)
		{
		# CAT instead of PRINT if it's a string or numeric
		if (is.character(x))
			{
			cat(x, "\n", sep="")
			}
		if (is.numeric(x))
			{
			cat(x, "\n", sep="")
			} else {
			print(x)
			}
		}
	else
		{
		pass="BLAH"
		}
	}




#######################################################
# np
#######################################################
#' normalizePath shortcut
#' 
#' Utility function that runs \code{\link[base]{normalizePath}}. Useful for
#' running on Mac vs. Windows.
#' 
#' @param path The path to run \code{\link[base]{normalizePath}} on.
#' @param ... Additional arguments to \code{\link[base]{normalizePath}}.
#' @return \code{path} The path that was normalized.
#' @export
#' @seealso \code{\link[base]{normalizePath}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' # Get a path
#' extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
#' extdata_dir
#'
#' path = paste(extdata_dir, "//", "Psychotria_5.2.newick", sep="")
#' path
#'
#' path = np(path)
#' path
#' 
np <- function(path=path, ...)
	{
	path = normalizePath(path, ...)
	return(path)
	}





#######################################################
# strsplit_whitespace
#######################################################
#' Split strings on whitespace
#' 
#' This function splits strings on whitespace (spaces and tabs), so you don't have
#' to remember the \code{regexp}/\code{grep} format codes.
#' 
#' @param tmpline A string containing text.
#' @return \code{list_of_strs} 
#' @export
#' @seealso \code{\link[base]{strsplit}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' tmpline = "Hello world see	my	tabs."
#' strsplit_whitespace(tmpline)
#' 
strsplit_whitespace <- function(tmpline)
	{
	# split on 1 or more whitespaces
	temp = strsplit(tmpline, "[ \t]+")
	
	# get the list
	list_of_strs = temp[[1]]
	
	# remove any leading/trailing ""
	list_of_strs = list_of_strs[list_of_strs != ""]
	
	return(list_of_strs)
	}


#######################################################
# moref
#######################################################
#' print to screen the header of a file
#' 
#' This does the rough equivalent of the \code{UNIX} function \code{more}, but within R.
#' 
#' @param fn A filename.
#' @param printnotcat If \code{TRUE}, use \code{\link[base]{print}} instead of \code{\link[base]{cat}}. Default \code{FALSE}.
#' @return Nothing returned.
#' @export
#' @seealso \code{\link[base]{scan}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
#' 
moref <- function(fn, printnotcat = FALSE)
	{
	lines = scan(file=fn, what="character", sep="\n")
	
	if (printnotcat == TRUE)
		{
		for (i in 1:length(lines))
			{
			print(lines[i])
			}
		}
	else
		{
		for (i in 1:length(lines))
			{
			cat(paste(lines[i], "\n", sep=""))
			}
		}
	}





# return matching TRUE/FALSE values
# list1 (.e.g. a big list) TRUE if it is found in list2 (e.g. a smaller list)

#######################################################
# match_list1_in_list2
#######################################################
#' Return TRUE for list1 items when they occur in list2
#' 
#' Return matching TRUE/FALSE values.  E.g. list1 (e.g. a big list) TRUE if it is found
#' in list2 (e.g. a smaller list)
#'
#' Utility function for %in%, when one's brain gets confused.
#' 
#' @param list1 The list of things you want to check
#' @param list2 The list of things you want to check against
#' @return \code{matchlist} The TRUE/FALSE list for list1
#' @export
#' @seealso \code{\link[base]{match}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
match_list1_in_list2 <- function(list1, list2)
	{
	matchlist = list1 %in% list2
	return(matchlist)
	}



#######################################################
# unlist_dtf_cols
#######################################################
#' Unlist the columns in a data.frame
#' 
#' Utility function. What it says.
#' 
#' @param dtf Input \code{\link[base]{data.frame}}
#' @param printflag Print the results if TRUE.
#' @return \code{dtf} The data.frame, hopefully without lists for columns
#' @export
#' @seealso \code{\link[base]{unlist}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
unlist_dtf_cols <- function(dtf, printflag=FALSE)
	{
	# Sometimes cbind makes each column a list, this can screw up use/searching of
	#  the column later on.  
	# Unlist each column...
	for (i in 1:ncol(dtf))
		{
		tmpstr = paste("unlisting col: ", names(dtf)[i], "...", sep="")
		prflag(tmpstr, printflag=printflag)		
		
		# catch a possible error from unlisting
		# the number of rows needs to stay the same!!
		tmpcol = unlist(dtf[, i])
		if (length(tmpcol) != length(dtf[, i]))
			{
			tmpstr = paste("...failed! unlist(col) length=", length(tmpcol), "; nrow(dtf) = ", nrow(dtf), sep="")
			prflag(tmpstr, printflag=printflag)
			} 
		else
			{
			dtf[, i] = tmpcol
			tmpstr = paste(" ", " ", sep="")
			prflag(tmpstr, printflag=printflag)
			}
		}
	
	#dtf2 = adf(dtf)
	
	return(dtf)
	}


# NOTE!!! THESE MATCH FUNCTIONS JUST RETURN THE *FIRST* MATCH, *NOT* ALL MATCHES
# (argh)
# return indices in 2nd list matching the first list
# It WILL return one match for each item in the list, though...

#######################################################
# get_indices_where_list1_occurs_in_list2
#######################################################
#' Return (first!) indices in second list matching the first list
#' 
#' This function will return one match (the first) for each item in the list; i.e. the second-list
#' index for each item in the first list.  Only the first hit in the second list is returned.
#' 
#' This is used by \code{\link{prt}}.
#'
#' @param list1 The first list. 
#' @param list2 The second list list.
#' @return \code{match_indices} The match indices.
#' @export
#' @seealso \code{\link{prt}}, \code{\link[base]{LETTERS}}, \code{\link{get_indices_where_list1_occurs_in_list2_noNA}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' list1 = c("N", "I", "C", "K")
#' list2 = LETTERS
#' get_indices_where_list1_occurs_in_list2(list1, list2)
get_indices_where_list1_occurs_in_list2 <- function(list1, list2)
	{
	match_indices = match(list1, list2)
	return(match_indices)
	}


# return indices in 2nd list matching the first list
#######################################################
# get_indices_where_list1_occurs_in_list2_noNA
#######################################################
#' Return (first!) indices in second list matching the first list, excluding NAs
#' 
#' This function will return one match (the first) for each item in the list; i.e. the second-list
#' index for each item in the first list.  Only the first hit in the second list is returned.  Unlike 
#' \code{\link{get_indices_where_list1_occurs_in_list2}}, non-hits (NAs) are excluded.
#' 
#' This is used by get_indices_of_branches_under_tips, which is used by \code{\link{extend_tips_to_ultrametricize}}, which can be used by section_the_tree.
#'
#' @param list1 The first list. 
#' @param list2 The second list list.
#' @return \code{match_indices} The match indices.
#' @export
#' @seealso \code{\link{prt}}, \code{\link[base]{LETTERS}}, \code{\link{get_indices_where_list1_occurs_in_list2}}, 
#' \code{\link{extend_tips_to_ultrametricize}}, \code{\link{section_the_tree}}, \code{\link{return_items_not_NA}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' list1 = c("N", "I", "C", "K")
#' list2 = LETTERS
#' get_indices_where_list1_occurs_in_list2_noNA(list1, list2)
#' 
get_indices_where_list1_occurs_in_list2_noNA <- function(list1, list2)
	{
	match_indices = match(list1, list2)
	match_indices = return_items_not_NA(match_indices)
	return(match_indices)
	}


#######################################################
# return_items_not_NA
#######################################################
#' Remove NAs from a vector/list
#' 
#' Utility function. This function returns the non-NA values from a vector.
#' 
#' This is used by \code{\link{get_indices_where_list1_occurs_in_list2_noNA}}, which is used 
#' by \code{\link{get_indices_of_branches_under_tips}}, which is used by 
#' \code{\link{extend_tips_to_ultrametricize}}, which can be used by \code{\link{section_the_tree}}.
#'
#' @param x The vector of items to check for being not NA.
#' @return \code{y} The surviving, non-NA cells of a vector.
#' @export
#' @seealso \code{\link{prt}}, \code{\link[base]{LETTERS}}, \code{\link{get_indices_where_list1_occurs_in_list2_noNA}},  
#' \code{\link{get_indices_where_list1_occurs_in_list2}}, \code{\link{extend_tips_to_ultrametricize}}, \code{\link{section_the_tree}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' list1 = c("N", "I", NA, "C", "K")
#' return_items_not_NA(list1)
#' 
return_items_not_NA <- function(x)
	{
	y = x[!is.na(x)]
	return(y)
	}


#######################################################
# order_tipranges_by_tr
#######################################################
#' Order the tipranges in a tipranges object so they match the order of tips in a tree
#' 
#' Utility function. What it says.  Life can get very confusing if you don't do this before plotting.
#' 
#' @param tipranges A tipranges object.
#' @param tr An ape tree object.
#' @return \code{tipranges} The reordered data.frame
#' @export
#' @seealso \code{\link[base]{unlist}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
order_tipranges_by_tr <- function(tipranges, tr)
	{
	tipranges_names = rownames(tipranges@df)
	tr_names = tr$tip.label
	
	match_indices = get_indices_where_list1_occurs_in_list2(list1=tr_names, list2=tipranges_names)
	
	tmpdf = tipranges@df[match_indices, ]
	tipranges@df = tmpdf
	
	return(tipranges)
	}



# Extract just the numbers from a string
# just the numbers INCLUDING THOSE CONNECTED BY DECIMAL POINTS!!!!
#######################################################
# extract_numbers
#######################################################
#' Extract just the numbers from a string, including decimal points
#' 
#' This function extracts numbers from a string.  Contiguous digits, including
#' decimal points, are made into a single number. A list of numbers is returned.
#' 
#' This saves you having to remember the \code{regexp}/\code{\link[base]{gregexpr}} code for this sort of thing, and
#' makes it much easier to parse numbers out of the text output of various programs.
#' 
#' @param tmpstr An input string.
#' @return \code{x2} The list of numbers
#' @export
#' @seealso \code{\link[base]{gregexpr}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' tmpstr = "190Ma - 65Ma"
#' extract_numbers(tmpstr)
#' 
#' tmpstr = "190.1Ma - 65.5Ma"
#' extract_numbers(tmpstr)
#' 
extract_numbers <- function(tmpstr)
	{
	defaults ='
	tmpstr = "190Ma - 65Ma"
	'
	# pull out the numbers / extract the numbers / extract numbers
	# just the numbers INCLUDING THOSE CONNECTED BY DECIMAL POINTS!!!!
	# (but not negative symbols)
	matches = gregexpr("(?:([0-9\\.]+))+", tmpstr)[[1]]
	
	# Get the ending points of the matches
	matches_end = matches-1+attr(matches,"match.length")
	
	# Extract the numbers from the string
	x = mapply(substr, tmpstr, matches, matches_end)

	# Convert to numeric
	x2 = as.numeric(x)
	return(x2)
	}


#######################################################
# list2str
#######################################################
#' Convert a list of items to a string
#' 
#' This is a shortcut to save time when converting a list of items to a string.
#' 
#' @param list1 The list to convert.
#' @param spacer The space between each item. Default " ".
#' @return \code{tmpstr} The output string.
#' @export
#' @seealso \code{\link[base]{paste}}, \code{\link[base]{as.character}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#'	 @cite FosterIdiots
#' @examples
#' test=1
#' 
list2str <- function(list1, spacer=" ")
	{
	
	for (i in 1:length(list1))
		{
		if (i == 1)
			{
			tmpstr = as.character(list1[1])
			if (length(list1) == 1)
				{
				return(tmpstr)
				}
			next
			}
		addstr = as.character(list1[i])
		tmpstr = paste(tmpstr, addstr, sep=spacer)
		}
	return(tmpstr)
	}





# Get the classes of the columns in a data frame
#######################################################
# cls.df
#######################################################
#' Get the class for each column in a list
#' 
#' This function returns the \code{\link[base]{class}} of each column in a \code{\link[base]{data.frame}}.
#' 
#' R does lots of weird and unpredictable things when you build up tables/matrices/data.frames
#' by e.g. \code{\link[base]{cbind}} and \code{\link[base]{rbind}} on vectors of results.  The major problems 
#' are (1) columns get made into class \code{\link[base]{list}}; (2) \code{\link[base]{numeric}}
#' columns are converted to class \code{\link[base]{factor}}; (3) \code{\link[base]{numeric}} columns
#' are converted to class \code{\link[base]{character}}; (4) you have a \code{\link[base]{matrix}} when
#' you think you have a \code{\link[base]{data.frame}}.
#' 
#' All of this could be taken care of by detailed understanding and tracking of when R recasts values in 
#' vectors, matrices, and data frames...but this is a huge pain, it is easier to just have a function
#' that jams everything back to a \code{\link[base]{data.frame}} with no lists, no factors, and with columns being numeric
#' where possible.  See \code{\link{dfnums_to_numeric}} and \code{\link{unlist_df4}} for these options.
#' 
#' @param dtf Input \code{\link[base]{data.frame}}.
#' @param printout Print the results to screen, if desired.
#' @return \code{dtf_classes} A \code{\link[base]{data.frame}} showing the column, column name, and column class.
#' @export
#' @seealso \code{\link{dfnums_to_numeric}}, \code{\link{unlist_df4}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
#' cls.df(x)
#' dfnums_to_numeric(adf(x))
#' unlist_df4(x)
#'
#' x = matrix(c(1,2,3,4,5,"A"), nrow=3, ncol=2)
#' cls.df(x)
#' dfnums_to_numeric(adf(x))
#' unlist_df4(x)
#'
#' x = adf(matrix(c(1,2,3,4,5,"A"), nrow=3, ncol=2))
#' names(x) = c("A","B")
#' cls.df(x)
#' dfnums_to_numeric(adf(x))
#' unlist_df4(x)
#' 
cls.df <- function(dtf, printout=FALSE)
	{
	# Convert to data.frame if needed
	if (class(dtf) == "matrix")
		{
		dtf = as.data.frame(dtf, stringsAsFactors=FALSE)
		}
	
	dtf_names = names(dtf)
	numcols = ncol(dtf)
	
	cls_col_list = c()
	for (i in 1:numcols)
		{
		# Initialize cls_col
		cls_col = NA
		
		# Get one column:
		cmdstr = paste("cls_col = class(dtf$'", dtf_names[i], "')", sep="")
		eval(parse(text = cmdstr))
		
		#cat(i, ": ", dtf_names[i], "	=	", cls_col, "\n", sep="")
		cls_col_list[i] = cls_col
		}
	
	# row names
	colnum = 1:numcols
	
	dtf_classes = cbind(colnum, dtf_names, cls_col_list)
	dtf_classes = data.frame(dtf_classes, row.names=colnum)
	
	# Print the output if true
	if (printout)
		{
		cat("\n")
		cat("cls.df(dtf) reports: dataframe 'dtf' has ", nrow(dtf), " rows, ", numcols, " columns.\n", sep="")
		cat("...names() and classes() of each column below...\n", sep="")
		cat("\n")
		print(dtf_classes)
		cat("\n")
		}	
	return(dtf_classes)
	}



# Get the classes of the columns in a data frame
#######################################################
# dfnums_to_numeric
#######################################################
#' Get the class for each column in a list
#' 
#' This function converts each column to class \code{\link[base]{numeric}} where possible, and
#' class \code{\link[base]{character}} otherwise.
#' 
#' R does lots of weird and unpredictable things when you build up tables/matrices/data.frames
#' by e.g. \code{\link[base]{cbind}} and \code{\link[base]{rbind}} on vectors of results.  The major problems 
#' are (1) columns get made into class \code{\link[base]{list}}; (2) \code{\link[base]{numeric}}
#' columns are converted to class \code{\link[base]{factor}}; (3) \code{\link[base]{numeric}} columns
#' are converted to class \code{\link[base]{character}}; (4) you have a \code{\link[base]{matrix}} when
#' you think you have a \code{\link[base]{data.frame}}.
#' 
#' All of this could be taken care of by detailed understanding and tracking of when R recasts values in 
#' vectors, matrices, and data frames...but this is a huge pain, it is easier to just have a function
#' that jams everything back to a \code{\link[base]{data.frame}} with no lists, no factors, and with columns being numeric
#' where possible.  See \code{\link{unlist_df4}} for more, and \code{\link{cls.df}} to see the class of each column.
#' 
#' \bold{WARNING: IF A COLUMN IS A MIX OF NUMBERS AND NON-NUMBERS, THE NON-NUMBERS WILL BE CONVERTED TO NA IF 
#' THE COLUMN IS MAJORITY NUMBERS (on default; see \code{max_NAs}).}
#' 
#' @param dtf Input \code{\link[base]{data.frame}}.
#' @param max_NAs Non-numeric cells will get converted to NA, up to the fraction of cells specified by \code{max_NAs}.  Above this
#' fraction, the column is converted to class \code{character}.
#' @param printout Print the results to screen, if desired.
#' @param roundval If not NULL, \code{\link[base]{round}} will be run using this for the number of digits.
#' @return \code{dtf} The output \code{\link[base]{data.frame}}.
#' @export
#' @seealso \code{\link{cls.df}}, \code{\link{unlist_df4}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
#' cls.df(x)
#' dfnums_to_numeric(adf(x))
#' unlist_df4(x)
#'
#' x = matrix(c(1,2,3,4,5,"A"), nrow=3, ncol=2)
#' cls.df(x)
#' dfnums_to_numeric(adf(x))
#' unlist_df4(x)
#'
#' x = adf(matrix(c(1,2,3,4,5,"A"), nrow=3, ncol=2))
#' names(x) = c("A","B")
#' cls.df(x)
#' dfnums_to_numeric(adf(x))
#' unlist_df4(x)
#' 
dfnums_to_numeric <- function(dtf, max_NAs=0.5, printout=FALSE, roundval=NULL)
	{
	dtf_classes = cls.df(dtf, printout=FALSE)
	
	dtf_names = names(dtf)
	numcols = ncol(dtf)
	
	cls_col_list = c()
	for (i in 1:numcols)
		{
		# Initialize cls_col
		cls_col = NA
		
		# Get one column:
		cmdstr = paste("cls_col = class(dtf$'", dtf_names[i], "')", sep="")
		eval(parse(text = cmdstr))
		
		#cat(i, ": ", dtf_names[i], "	=	", cls_col, "\n", sep="")
		cls_col_list[i] = cls_col
		}
	
	for (i in 1:numcols)
		{
		if (cls_col_list[i] == "list")
			{
			next()	# skip the lists
			}
		if (cls_col_list[i] != "numeric")
			{
			# Initialize cls_col
			newcol = NA
	
			# Get one column, convert to numeric:
			cmdstr = paste("newcol = as.numeric(as.character(dtf$'", dtf_names[i], "'))", sep="")
			suppressWarnings(eval(parse(text = cmdstr)))
			
			# If it's less than 50% NAs (or max_NA NAs), then convert to numeric
			#print(newcol)
			#print(max_NAs * length(newcol))
			#print(sum(is.na(newcol)))
			if (sum(is.na(newcol)) < (max_NAs * length(newcol)))
				{
				# Get the column, convert to numeric:
				# (if it's a factor, you have to convert to character, then to numeric)
				# (if it's a character, you can convert to character anyway...)
				#cmdstr = paste("dtf$", dtf_names[i], " = as.numeric(as.character(dtf$", dtf_names[i], "))", sep="")
				cmdstr = paste("dtf$'", dtf_names[i], "' = newcol", sep="")
				suppressWarnings(eval(parse(text = cmdstr)))
				
				
				# If a rounding val is specified, do the rounding
				if (!is.null(roundval))
					{
					cmdstr = paste("dtf$'", dtf_names[i], "' = round(dtf$'", dtf_names[i], "', digits=roundval)", sep="")
					suppressWarnings(eval(parse(text = cmdstr)))
					}
				
				}
			}
		}
	tmp_classes = cls.df(dtf)
	dtf_classes$newclasses = tmp_classes[,ncol(tmp_classes)]
	
	if(printout)
		{
		cat("\n")
		cat("dfnums_to_numeric(dtf, max_NAs=", max_NAs, ") reports: dataframe 'dtf_classes' has ", nrow(dtf_classes), " rows, ", ncol(dtf_classes), " columns.\n", sep="")
		cat("...names() and classes() of each column below...\n", sep="")
		cat("\n")
		print(dtf_classes)
		}
	return(dtf)
	}



#######################################################
# A few generic functions from genericR_v1.R, used in BioGeoBEARS
#######################################################
#' Convert to data.frame, without factors
#' 
#' Shortcut for: \code{as.data.frame(x, row.names=NULL, stringsAsFactors=FALSE)}
#' 
#' This function, and \code{\link{adf2}}, are useful for dealing with errors due to 
#' automatic conversion of some columns to factors.  Another solution may be to prepend
#' \code{options(stringsAsFactors = FALSE)} at the start of one's script, to turn off all default stringsAsFactors silliness.
#' 
#' @param x matrix or other object transformable to data.frame
#' @return data.frame
#' @export
#' @seealso \code{\link{adf2}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
#' adf(x)
#' 
adf <- function(x)
	{
	return(as.data.frame(x, row.names=NULL, stringsAsFactors=FALSE))
	}



#' Convert to data.frame, without factors
#' 
#' Shortcut for: \code{tmp_rownames = 1:nrow(x); as.data.frame(x, row.names=tmp_rownames, stringsAsFactors=FALSE)}
#' 
#' This function, and \code{\link{adf2}}, are useful for dealing with errors due to 
#' automatic conversion of some columns to factors.  Another solution may be to prepend
#' \code{options(stringsAsFactors = FALSE)} at the start of one's script, to turn off all default stringsAsFactors silliness.
#'
#' In adf2, rownames are forced to be numbers; this can prevent errors due to e.g. repeated rownames
#' after an \code{rbind} operation.
#'
#' @param x matrix or other object transformable to data.frame
#' @return data.frame
#' @export
#' @seealso \code{\link{adf}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
#' adf2(x)
#' 
adf2 <- function(x)
	{
	# Deals with the problem of repeated row names
	rownames = 1:nrow(x)
	return(as.data.frame(x, row.names=rownames, stringsAsFactors=FALSE))
	}


#' String splitting shortcut
#' 
#' \code{\link[base]{strsplit}} returns the results inside a list, which is annoying. \code{strsplit2} shortens the process.
#'
#' @param x A string to split
#' @param ... Other arguments to \code{\link[base]{strsplit}}.  The argument \code{split} is \emph{required}.
#' @return \code{out} The output from inside the list.
#' @export
#' @seealso \code{\link[base]{strsplit}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
#' 
#' # strsplit returns the results inside a list element
#' out = strsplit("ABC", split="")
#' out
#' # I.e....
#' out[[1]]
#' 
#' # If this is annoying/ugly in the code, use strsplit2:
#' out = strsplit2("ABC", split="")
#' out
#' 
strsplit2 <- function(x, ...)
	{
	out = strsplit(x, ...)[[1]]
	return(out)
	}


#######################################################
# unlist_df:
#######################################################
#' Unlist the columns in a data.frame
#' 
#' Sometimes, matrices or data.frames will malfunction due to their having lists as columns
#' and other weirdness.  This is a shortcut for \code{data.frame(lapply(df, function(x) unlist(x)))}.
#' 
#' @param df matrix or other object transformable to data.frame
#' @return data.frame
#' @export
#' @seealso \code{\link{unlist_df2}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
#' unlist_df2(x)
#'
unlist_df <- function(df)
	{
	outdf <- data.frame(lapply(df, function(x) unlist(x)))
	}



#######################################################
# unlist_df2:
#######################################################
#' Unlist the columns in a data.frame, with more checks
#' 
#' Sometimes, matrices or data.frames will malfunction due to their having lists as columns
#' and other weirdness. This runs \code{\link{unlist}} and additional checks.
#' 
#' @param df matrix or other object transformable to data.frame
#' @return \code{outdf} A \code{\link[base]{matrix}}.
#' @export
#' @seealso \code{\link{unlist_df}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
#' unlist_df2(x)
#'
unlist_df2 <- function(df)
	{
	store_colnames = names(df)
	
	outdf = NULL
	
	numrows = dim(df)[1]
	
	for (i in 1:ncol(df))
		{
		#print(names(df)[i])
		tmpcol = unlist(df[, i])
		#print(length(tmpcol))
		
		# Error check; e.g. blank cells might screw it up
		if (length(tmpcol) < numrows)
			{
			tmpcol2 = df[,i]
			tmpcol = as.character(tmpcol2)
			}
		
		outdf = cbind(outdf, tmpcol)
		}
	
	#outdf = adf2(outdf)
	
	names(outdf) = store_colnames
	return(outdf)
	}



#######################################################
# unlist_df3:
#######################################################
#' Unlist the columns in a data.frame, with more checks and adf
#' 
#' Sometimes, matrices or data.frames will malfunction due to their having lists as columns
#' and other weirdness. This runs \code{\link[base]{unlist}} and additional checks, and
#' forces conversion to a \code{\link[base]{data.frame}} at the end.
#' 
#' @param df matrix or other object transformable to data.frame
#' @return \code{outdf} data.frame
#' @export
#' @seealso \code{\link{unlist_df}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
#' unlist_df3(x)
#'
unlist_df3 <- function(df)
	{
	store_colnames = names(df)
	store_rownames = rownames(df)
	
	outdf = NULL
	
	numrows = dim(df)[1]
	
	for (i in 1:ncol(df))
		{
		#print(names(df)[i])
		tmpcol = unlist(df[, i])
		#print(length(tmpcol))
		
		# Error check; e.g. blank cells might screw it up
		if (length(tmpcol) < numrows)
			{
			tmpcol2 = df[,i]
			tmpcol = as.character(tmpcol2)
			}
		
		outdf = cbind(outdf, tmpcol)
		}
	
	outdf = adf2(outdf)
	
	names(outdf) = store_colnames
	rownames(outdf) = store_rownames
	return(outdf)
	}


#######################################################
# unlist_df4:
#######################################################
#' Unlist the columns in a data.frame, with more checks, adf, and dfnums_to_numeric
#' 
#' Sometimes, matrices or data.frames will malfunction due to their having lists as columns
#' and other weirdness. This runs \code{\link[base]{unlist}} and additional checks, and
#' forces conversion to a \code{\link[base]{data.frame}} at the end.  It also adds
#' \code{\link{dfnums_to_numeric}} which should remove the problem of numbers columns being of
#' class \code{\link[base]{character}}.
#' 
#' See especially  \code{\link[base]{data.matrix}} for a possibly simpler alternative.
#' 
#' @param df matrix or other object transformable to data.frame
#' @param ... Additional options passed to \code{\link{dfnums_to_numeric}}.
#' @return \code{outdf} data.frame
#' @export
#' @seealso \code{\link{unlist_df}}, \code{\link{dfnums_to_numeric}}, \code{\link{cls.df}}, \code{\link[base]{data.matrix}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
#' cls.df(x)
#' unlist_df4(x)
#'
#' x = matrix(c(1,2,3,4,5,"A"), nrow=3, ncol=2)
#' cls.df(x)
#' unlist_df4(x)
#'
#' x = adf(matrix(c(1,2,3,4,5,"A"), nrow=3, ncol=2))
#' names(x) = c("A","B")
#' cls.df(x)
#' unlist_df4(x)
#' 
unlist_df4 <- function(df, ...)
	{
	store_colnames = names(df)
	store_rownames = rownames(df)
	
	outdf = NULL
	
	numrows = dim(df)[1]
	numcols = dim(df)[2]
	
	for (i in 1:ncol(df))
		{
		#print(names(df)[i])
		tmpcol = unlist(df[, i])
		#print(length(tmpcol))
		
		# Error check; e.g. blank cells might screw it up
		if (length(tmpcol) < numrows)
			{
			tmpcol2 = df[,i]
			tmpcol = as.character(unlist(tmpcol2))
			}
		
		outdf = cbind(outdf, tmpcol)
		}

	# Unlist each row
# 	outdf2 = NULL
# 	for (i in 1:nrow(df))
# 		{
# 		#print(names(df)[i])
# 		tmprow = unlist(df[i, ])
# 		#print(length(tmpcol))
# 		
# 		# Error check; e.g. blank cells might screw it up
# 		if (length(tmprow) < numcols)
# 			{
# 			tmprow2 = df[i,]
# 			tmprow = as.character(unlist(tmprow2))
# 			}
# 		
# 		outdf2 = rbind(outdf2, tmprow)
# 		}
# 	outdf = outdf2
	
	outdf_tmp = adf2(outdf)
	
	# Remove factors and character silliness from numbers
	outdf = dfnums_to_numeric(outdf_tmp, ...)
	
	names(outdf) = store_colnames
	rownames(outdf) = store_rownames
	return(outdf)
	}




#######################################################
# slashslash:
#######################################################
#' Remove double slash (slash a slash)
#' 
#' Shortcut for: \code{gsub(pattern="//", replacement="/", x=tmpstr)}
#' 
#' This function is useful for removing double slashes that can
#' appear in full pathnames due to inconsistencies in trailing
#' slashes in working directories etc.
#'
#' @param tmpstr a path that you want to remove double slashes from
#' @return outstr a string of the fixed path
#' @export
#' @seealso \code{\link[base]{getwd}}, \code{\link[base]{setwd}}, \code{\link[base]{gsub}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' tmpstr = "/Library/Frameworks//R.framework/Versions/"
#'
#' outstr = slashslash(tmpstr)
#' outstr
#'
slashslash <- function(tmpstr)
	{
	outstr = np(gsub(pattern="//", replacement="/", x=tmpstr))
	return(outstr)
	}


# Add a slash to a directory name if needed
#######################################################
# addslash:
#######################################################
#' Add a slash to a directory name if needed
#' 
#' This function adds a slash to the end of the string, if one is not present. Handy for standardizing paths.
#'
#' @param tmpstr a path that you want to possibly add a slash to 
#' @return outstr a string of the fixed path
#' @export
#' @seealso \code{\link[base]{getwd}}, \code{\link[base]{setwd}}, \code{\link[base]{gsub}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' tmpstr = "/Dropbox/_njm/__packages"
#' tmpstr
#' outstr = addslash(tmpstr)
#' outstr
#'
#' # Annoying, getwd() often doesn't return the ending slash, which 
#' # can make life hard for paste() later on
#' tmpstr = getwd()
#' tmpstr
#' outstr = addslash(tmpstr)
#' outstr
#'
addslash <- function(tmpstr)
	{
	tmpchars = strsplit(tmpstr, split="")[[1]]
	last_char = tmpchars[length(tmpchars)]
	
	if (last_char != "/")
		{
		outstr = paste(tmpstr, "/", sep="")
		} else {
		outstr = tmpstr
		}
		
	return(outstr)
	}





#######################################################
# Functions for calculating model summary statistics,
# and displaying nice tables
#######################################################

# Nicely print tables
# Manually format conditional formatting table...

#######################################################
# conditional_format_cell
#######################################################
#' Conditionally format a number (mostly)
#' 
#' When a table has numbers that range over many orders of magnitude, it can be
#' very distracting if the display program forces each column to the same format.  This
#' function formats a cell much like Excel would.
#' 
#' The defaults seem to work well, but could be modified. The current function also extracts just the filename, if a
#' full path is given.
#' 
#' @param cellval The cell value to format.
#' @param numbers_below_this_get_scientific When the absolute value of a number is below this value, scientific notation is used.
#' @param numdigits_for_superlow_scientific Number of digits after the '.' for scientific notation of small numbers.
#' @param numbers_above_this_get_scientific When the absolute value of a number is above this value, scientific notation is used.
#' @param numdigits_for_superhigh_scientific Number of digits after the '.' for scientific notation of large numbers.
#' @param numdigits_inbetween_have_fixed_digits Numbers of medium size have this many fixed digits.  Note that other cutoffs are specified
#' in the code, and \code{\link[base]{signif}} is used to make e.g. integers appear as 0, 1, 2..
#' @return \code{cellval} The value, reformatted and of class \code{\link[base]{character}}.
#' @export
#' @seealso \code{\link[base]{signif}}, \code{\link[base]{sprintf}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
#' 
#' cellval = 143514514514532
#' conditional_format_cell(cellval)
#' 
#' cellval = -42.235235
#' conditional_format_cell(cellval)
#' 
#' cellval = -42.0000000
#' conditional_format_cell(cellval)
#' 
#' cellval = 0.0000
#' conditional_format_cell(cellval)
#' 
#' cellval = 0.0001
#' conditional_format_cell(cellval)
#' 
#' cellval = 0.00001
#' conditional_format_cell(cellval)
#' 
#' cellval = 0.0000111
#' conditional_format_cell(cellval)
#' 
conditional_format_cell <- function(cellval, numbers_below_this_get_scientific=0.0001, numdigits_for_superlow_scientific=1, numbers_above_this_get_scientific=10000000, numdigits_for_superhigh_scientific=2, numdigits_inbetween_have_fixed_digits=4)
 	{
	defaults='
	# Basically, these are our categories
	numbers_below_this_get_scientific = 0.0001
	numdigits_for_superlow_scientific = 2
	
	numbers_above_this_get_scientific = 10000000
	numdigits_for_superhigh_scientific = 2
	
	# numdigits inbetween get fixed to a certain total number of digits
	numdigits_inbetween_have_fixed_digits = 6
	
	
	results_summary_table3
	
	cellval = results_summary_table3$LnLalt[3]
	cellval = results_summary_table3$AIC_wt_vBest[1]
	cellval = 143514514514532
	'
	
	
	# Is it a number, or convertable to a number?
	tmp_cellval = suppressWarnings(as.numeric(cellval))
	if ( !is.na(tmp_cellval)==TRUE )
		{
		# Also, don't do 0.0000, 0.00e+00, etc.
		if (tmp_cellval==0)
			{
			cellval = signif(as.numeric(cellval), 8)
			cellval = as.character(cellval)
			return(cellval)
			}

		
		# Number just right
		if ( (abs(tmp_cellval) >= numbers_below_this_get_scientific) &&  (abs(tmp_cellval) <= numbers_above_this_get_scientific) )
			{
			fmt_str = paste('"%.', numdigits_inbetween_have_fixed_digits, 'f"', sep="")
			sprintf_cmd = paste("cellval = sprintf(fmt=", fmt_str, ", ", tmp_cellval, ")", sep="")
			eval(parse(text=sprintf_cmd))

			if (abs(tmp_cellval) >= 1000)
				{
				fmt_str = paste('"%.', 0, 'f"', sep="")
				sprintf_cmd = paste("cellval = sprintf(fmt=", fmt_str, ", ", tmp_cellval, ")", sep="")
				eval(parse(text=sprintf_cmd))
				}

			if ( (abs(tmp_cellval) >= 100) && (abs(tmp_cellval) < 1000) )
				{
				fmt_str = paste('"%.', 1, 'f"', sep="")
				sprintf_cmd = paste("cellval = sprintf(fmt=", fmt_str, ", ", tmp_cellval, ")", sep="")
				eval(parse(text=sprintf_cmd))
				}

			if ( (abs(tmp_cellval) >= 0.1) && (abs(tmp_cellval) < 100) )
				{
				fmt_str = paste('"%.', 2, 'f"', sep="")
				sprintf_cmd = paste("cellval = sprintf(fmt=", fmt_str, ", ", tmp_cellval, ")", sep="")
				eval(parse(text=sprintf_cmd))
				}

			if ( (abs(tmp_cellval) >= 0.01) && (abs(tmp_cellval) < 0.1) )
				{
				fmt_str = paste('"%.', 3, 'f"', sep="")
				sprintf_cmd = paste("cellval = sprintf(fmt=", fmt_str, ", ", tmp_cellval, ")", sep="")
				eval(parse(text=sprintf_cmd))
				}

			
			# If above 0, don't do 1.00000, 2.0000
			if (abs(tmp_cellval) >= 1)
				{
				cellval = signif(as.numeric(cellval), 8)
				}
			
			cellval = as.character(cellval)
			return(cellval)
			}	

		# Number too small
		if (abs(tmp_cellval) < numbers_below_this_get_scientific)
			{
			fmt_str = paste('"%0.', numdigits_for_superlow_scientific, 'e"', sep="")
			sprintf_cmd = paste("cellval = sprintf(fmt=", fmt_str, ", ", tmp_cellval, ")", sep="")
			eval(parse(text=sprintf_cmd))
			cellval = as.character(cellval)
			return(cellval)
			}
		
		# Number too big
		if (abs(tmp_cellval) > numbers_above_this_get_scientific)
			{
			fmt_str = paste('"%0.', numdigits_for_superhigh_scientific, 'e"', sep="")
			sprintf_cmd = paste("cellval = sprintf(fmt=", fmt_str, ", ", tmp_cellval, ")", sep="")
			eval(parse(text=sprintf_cmd))
			cellval = as.character(cellval)
			return(cellval)
			}
		
		}

	# Check for NA, convert to ""
	if (is.na(cellval) == TRUE)
		{
		cellval = ""
		return(cellval)
		}

	# If character
	if (is.character(cellval) == TRUE)
		{
		# Just get the filename, if it's a path with slashes
		cellval = get_path_last(path=cellval)
		return(cellval)
		}
	return(cellval)
	}



#######################################################
# conditional_format_table
#######################################################
#' Conditionally format the numbers (mostly) in a table
#' 
#' When a table has numbers that range over many orders of magnitude, it can be
#' very distracting if the display program forces each column to the same format.  This
#' function uses \code{\link{conditional_format_cell}} via \code{\link[base]{sapply}} to format a cell much like Excel would.
#' 
#' The defaults seem to work well, but could be modified. The current function also extracts just the filename, if a
#' full path is given.
#' 
#' @param input_table The table to format.
#' @param numbers_below_this_get_scientific When the absolute value of a number is below this value, scientific notation is used.
#' @param numdigits_for_superlow_scientific Number of digits after the '.' for scientific notation of small numbers.
#' @param numbers_above_this_get_scientific When the absolute value of a number is above this value, scientific notation is used.
#' @param numdigits_for_superhigh_scientific Number of digits after the '.' for scientific notation of large numbers.
#' @param numdigits_inbetween_have_fixed_digits Numbers of medium size have this many fixed digits.  Note that other cutoffs are specified
#' in the code, and \code{\link[base]{signif}} is used to make e.g. integers appear as 0, 1, 2..
#' @return \code{output_table} The table, reformatted with cells of class \code{\link[base]{character}}.
#' @export
#' @seealso \code{\link[base]{signif}}, \code{\link[base]{sprintf}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
#' 
#' input_table = adf(c(143514514514532, -42.235235, -42.0000000, 
#' 0.0000, 0.0001, 0.00001, 0.0000111))
#' conditional_format_table(input_table=input_table)
#' 
conditional_format_table <- function(input_table, numbers_below_this_get_scientific=0.0001, numdigits_for_superlow_scientific=1, numbers_above_this_get_scientific=10000000, numdigits_for_superhigh_scientific=2, numdigits_inbetween_have_fixed_digits=4)
	{
	defaults='
	input_table = adf(c(143514514514532, -42.235235, -42.0000000, 0.0000, 0.0001, 0.00001, 0.0000111))
	conditional_format_table(input_table=input_table)
	'
	
	# Unlist table and run through sapply
	chardata = sapply(unlist(input_table), FUN=conditional_format_cell, numbers_below_this_get_scientific, numdigits_for_superlow_scientific, numbers_above_this_get_scientific, numdigits_for_superhigh_scientific, numdigits_inbetween_have_fixed_digits)
	names(chardata) = NULL
	
	output_table = adf(matrix(data=chardata, nrow=nrow(input_table), ncol=ncol(input_table), byrow=FALSE))
	names(output_table) = names(input_table)
	rownames(output_table) = rownames(input_table)
	
	return(output_table)
	}
	









#######################################################
# get_path_last
#######################################################
#' Get the text that comes after the last slash
#' 
#' Extracts the filename from a full path.
#' 
#' @param path A string of class \code{\link[base]{character}}.
#' @return \code{lastword} A string with the filename, without the path.
#' @export
#' @seealso \code{\link{get_path_first}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#'	 @cite FosterIdiots
#' @examples
#' get_path_last("/Users/nickm/Psychotria_geog.data")
#' 
get_path_last <- function(path)
	{
	words = strsplit(path, split="/")[[1]]
	lastword = words[length(words)]
	
	return(lastword)
	}


#######################################################
# get_path_first
#######################################################
#' Get the text that comes before the last slash
#' 
#' Extracts the path from a full path, removing the filename.
#' 
#' @param inpath A string of class \code{\link[base]{character}}.
#' @param addslash If \code{TRUE}, add a slash at the end of the path.
#' @return \code{outpath} A string with the full path, without the file.
#' @export
#' @seealso \code{\link{get_path_last}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' get_path_first("/Users/nickm/Library/Psychotria_geog.data")
#' 
get_path_first <- function(inpath, addslash="FALSE")
	{
	words = strsplit(inpath, split="/")[[1]]
	pathwords = words[1:(length(words)-1)]
	outpath = paste(pathwords, sep="", collapse="/")
	
	if (addslash == TRUE)
		{
		
		}

	# outpath = paste("/", outpath, sep="")	# not needed
	outpath
	
	return(outpath)
	}




# Get everything BEFORE the last suffix (.nex or whatever)
#######################################################
# get_fn_prefix
#######################################################
#' Get everything BEFORE the last suffix (.nex or whatever)
#' 
#' Extracts the string from before the last suffix.  I.e., "filename.nex" becomes "filename".
#' 
#' @param fn The input filename.
#' @return \code{prefix} The output string.
#' @export
#' @seealso \code{\link{get_path_last}}, \code{\link{get_path_first}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' get_fn_prefix("/Users/nickm/Library/R/Psychotria_geog.data")
#' get_fn_prefix("Psychotria_geog.data")
#' 
get_fn_prefix <- function(fn)
	{
	defaults='
	fn = "example-data/aP6.fas"
	'# end defaults
	
	words = strsplit(fn, split="\\.")[[1]]
	words_to_keep = words[1:(length(words)-1)]
	prefix = paste(words_to_keep, collapse=".", sep="")
	
	return(prefix)
	}





#######################################################
# Statistical and likelihood-related functions
#######################################################


#######################################################
# calc_AIC_vals
#######################################################
#' Calculate AICs for a list of models
#' 
#' A list of AICs (Akaike Information Criterion) is calculated from two input lists.  Lower values
#' of AIC indicate some combination of better fit to the data and more parsimony in the model 
#' (fewer free parameters).
#' 
#' The two input lists are:
#' 
#' \bold{1.} A list of data likelihoods under a variety of models.\cr
#' \bold{2.} A list of the number of free parameters under each model.\cr
#' 
#' See \cite{Burnham_Anderson_2002} and \url{http://www.brianomeara.info/tutorials/aic} for 
#' discussion of AIC and its uses.
#' 
#' @param LnL_vals A vector of log-likelihoods (typically negative, but may not be for continuous data).
#' @param nparam_vals A vector of the number of parameters for each model.
#' @return \code{AIC_vals} A vector of AIC results.
#' @export
#' @seealso \code{\link{calc_AIC_column}}, \code{\link{calc_AICc_column}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Burnham_Anderson_2002
#' @examples
#' LnL_vals = c(-34.5, -20.9)
#' nparam_vals = c(2, 3)
#' calc_AIC_vals(LnL_vals, nparam_vals)
#' 
#' LnL_vals = c(-20.9, -20.9, -20.9, -20.9)
#' nparam_vals = c(3, 4, 5, 6)
#' calc_AIC_vals(LnL_vals, nparam_vals)
#' 
calc_AIC_vals <- function(LnL_vals, nparam_vals)
	{
	AIC_vals = mapply(FUN=getAIC, LnL=LnL_vals, numparams=nparam_vals)
	return(AIC_vals)
	}

#######################################################
# calc_AIC_column
#######################################################
#' Calculate AICs to make a column in a table
#' 
#' A list of AICs (Akaike Information Criterion) is calculated from two input lists.  Lower values
#' of AIC indicate some combination of better fit to the data and more parsimony in the model 
#' (fewer free parameters).
#' 
#' The two input lists are:
#' 
#' \bold{1.} A list of data likelihoods under a variety of models.\cr
#' \bold{2.} A list of the number of free parameters under each model.\cr
#' 
#' See \cite{Burnham_Anderson_2002} and \url{http://www.brianomeara.info/tutorials/aic} for 
#' discussion of AIC and its uses.
#' 
#' @param LnL_vals A vector of log-likelihoods (typically negative, but may not be for continuous data).
#' @param nparam_vals A vector of the number of parameters for each model.
#' @return \code{AIC_col} A \code{\link[base]{data.frame}} column of AIC results.
#' @export
#' @seealso \code{\link{calc_AIC_vals}}, \code{\link{calc_AICc_vals}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Burnham_Anderson_2002
#' @examples
#' test=1
#' 
#' LnL_vals = c(-20.9, -20.9, -20.9, -20.9)
#' nparam_vals = c(3, 4, 5, 6)
#' calc_AIC_column(LnL_vals, nparam_vals)
#' 
calc_AIC_column <- function(LnL_vals, nparam_vals)
	{
	AIC_vals = calc_AIC_vals(LnL_vals, nparam_vals)
	
	# Make a column
	AIC_col = matrix(AIC_vals, nrow=length(AIC_vals), ncol=1)
	AIC_col = adf2(AIC_col)
	names(AIC_col) = "AIC"
	
	return(AIC_col)
	}



#######################################################
# calc_AICc_vals
#######################################################
#' Calculate AICc values for a list of models
#' 
#' A list of AICc values (second order Akaike Information Criterion) is calculated from two input lists.  Lower values
#' of AICc indicate some combination of better fit to the data and more parsimony in the model 
#' (fewer free parameters). AICc contains a correction for sample size.
#' 
#' The two input lists are:
#' 
#' \bold{1.} A list of data likelihoods under a variety of models.\cr
#' \bold{2.} A list of the number of free parameters under each model.\cr
#' 
#' \code{samplesize} can be a scalar or vector; but see below.
#' 
#' See \cite{Burnham_Anderson_2002} and \url{http://www.brianomeara.info/tutorials/aic} for 
#' discussion of AIC, AICc and their uses.
#' 
#' @param LnL_vals A vector of log-likelihoods (typically negative, but may not be for continuous data).
#' @param nparam_vals A vector of the number of parameters for each model.
#' @param samplesize A single samplesize, or a vector of the samplesizes each model.  However, samplesize should always be
#' the same for all comparisons, since maximum likelihood and AIC/AICc model-selection methods are always comparing different models
#' on the \emph{same} data, not different data on the same mode.
#' @return \code{AICc_vals} A vector of AICc results.
#' @export
#' @seealso \code{\link{calc_AIC_vals}}, \code{\link{calc_AICc_column}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Burnham_Anderson_2002
#' @examples
#' LnL_vals = c(-34.5, -20.9)
#' nparam_vals = c(2, 3)
#' calc_AICc_vals(LnL_vals, nparam_vals, samplesize=20)
#' 
#' LnL_vals = c(-20.9, -20.9, -20.9, -20.9)
#' nparam_vals = c(3, 4, 5, 6)
#' calc_AICc_vals(LnL_vals, nparam_vals, samplesize=20)
#' 
calc_AICc_vals <- function(LnL_vals, nparam_vals, samplesize)
	{
	AICc_vals = mapply(FUN=getAICc, LnL=LnL_vals, numparams=nparam_vals, samplesize)
	return(AICc_vals)
	}


#######################################################
# calc_AICc_column
#######################################################
#' Calculate AICc values for a list of models
#' 
#' A list of AICc values (second order Akaike Information Criterion) is calculated from two input lists.  Lower values
#' of AICc indicate some combination of better fit to the data and more parsimony in the model 
#' (fewer free parameters). AICc contains a correction for sample size.
#' 
#' The two input lists are:
#' 
#' \bold{1.} A list of data likelihoods under a variety of models.\cr
#' \bold{2.} A list of the number of free parameters under each model.\cr
#' 
#' \code{samplesize} can be a scalar or vector; but see below.
#' 
#' See \cite{Burnham_Anderson_2002} and \url{http://www.brianomeara.info/tutorials/aic} for 
#' discussion of AIC, AICc and their uses.
#' 
#' @param LnL_vals A vector of log-likelihoods (typically negative, but may not be for continuous data).
#' @param nparam_vals A vector of the number of parameters for each model.
#' @param samplesize A single samplesize, or a vector of the samplesizes each model.  However, samplesize should always be
#' the same for all comparisons, since maximum likelihood and AIC/AICc model-selection methods are always comparing different models
#' on the \emph{same} data, not different data on the same mode.
#' @return \code{AICc_col} A \code{\link[base]{data.frame}} column of AICc results.
#' @export
#' @seealso \code{\link{calc_AICc_vals}}, \code{\link{calc_AIC_column}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Burnham_Anderson_2002
#' @examples
#' LnL_vals = c(-34.5, -20.9)
#' nparam_vals = c(2, 3)
#' calc_AICc_column(LnL_vals, nparam_vals, samplesize=20)
#' 
#' LnL_vals = c(-20.9, -20.9, -20.9, -20.9)
#' nparam_vals = c(3, 4, 5, 6)
#' calc_AICc_column(LnL_vals, nparam_vals, samplesize=20)
#' 
calc_AICc_column <- function(LnL_vals, nparam_vals, samplesize)
	{
	AICc_vals = calc_AICc_vals(LnL_vals, nparam_vals, samplesize)
	
	# Make a column
	AICc_col = matrix(AICc_vals, nrow=length(AICc_vals), ncol=1)
	AICc_col = adf2(AICc_col)
	names(AICc_col) = "AICc"
	
	return(AICc_col)
	}


#######################################################
# lrttest_on_summary_table
#######################################################
#' Calculate Likelihood Ratio Test (LRT) results, and add to table
#' 
#' The Likelihood Ratio Test (LRT) is a standard method for testing whether or not the data likelihood
#' conferred by a more complex is significantly better than the data likelihood conferred by
#' the simpler model, given a certain number of extra free parameters for the complex model.  The null hypothesis
#' is that there is no difference; rejection means that there is a statistically significant improvement in 
#' the more complex model.
#' 
#' The LRT only works for situations in which the simpler model is nested within the more complex model (i.e., by 
#' taking some parameters of the more complex model and forcing them to be fixed to a specific value).  In addition,
#' the LRT may be unreliable in data-poor situations, and inherits whatever difficulties there may be
#' in ML searches.  See \cite{Burnham_Anderson_2002} for discussion.
#'
#' This function assumes that the log-likelihoods are in the column "LnL", and the number of parameters is specified in "nparams"
#' 
#' @param restable A \code{\link[base]{data.frame}} with at least columns named "LnL" and "nparams".
#' @param row_to_use_as_null This is the row specifying the model to which the others will be compared in pairwise fashion.
#' @param rows_to_exclude Some rows may have models that the simpler model cannot nest within.  These should be excluded.
#' @param returnwhat If "pval", just return the p-value.  If "all", return all of the intermediate outputs.
#' @param add_to_table If TRUE, add to the main table and return the main table. If FALSE, return just the Akaike Weights results.
#' @return \code{pval} or \code{LRTrow}, both \code{\link[base]{data.frame}}.  Depends on \code{returnwhat}.
#' @export
#' @seealso \code{\link{lrttest}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Burnham_Anderson_2002
#' @examples
#' test=1
#' 
lrttest_on_summary_table <- function(restable, row_to_use_as_null, rows_to_exclude, returnwhat="pval", add_to_table=TRUE)
	{
	defaults='
	restable = results_summary_table
	row_to_use_as_null=2
	rows_to_exclude=c(1)
	returnwhat="all"
	add_to_table=FALSE
	'
	
	# This function assumes that the log-likelihoods are in the column "LnL", and the number of parameters is specified in "nparams"
	
	LRT_table = NULL
	for (j in 1:nrow(restable))
		{
		LnL_1 = restable$LnL[j]
		LnL_2 = restable$LnL[row_to_use_as_null]
		numparams1 = restable$nparams[j]
		numparams2 = restable$nparams[row_to_use_as_null]
		returnwhat=returnwhat		
		LRTrow = lrttest(LnL_1, LnL_2, numparams1, numparams2, returnwhat=returnwhat)
		length_LRTrow = length(LRTrow)
		
		# Check for rows_to_exclude
		if (  ((j %in% rows_to_exclude) == TRUE) || (j==row_to_use_as_null))
			{
			if (length(LRTrow) == 1)
				{
				LRTrow = adf(NA)
				names(LRTrow) = c("pval")
				} else {
				LRTrow = adf(t(rep(NA, times=length_LRTrow)))
				names(LRTrow) = c("alt", "null", "LnLalt", "LnLnull", "DFalt", "DFnull", "DF", "Dstatistic", "pval", "test", "tail")
				}
			}

		if (length(LRTrow) != 1)
			{
			LRTrow$null = row_to_use_as_null
			LRTrow$alt = j
			}
		
		# Add to LRT table
		LRT_table = rbind(LRT_table, LRTrow)		
		}
	
	if (add_to_table == TRUE)
		{
		restable = cbind(restable, LRT_table)
		return(restable)
		} else {
		return(LRT_table)
		}
	return("ERROR!")
	}



#######################################################
# AkaikeWeights_on_summary_table
#######################################################
#' Calculate Akaike Weights, and add to table
#' 
#' This calculates Akaike Weights (relative probabilities on models explaining the same data) for the models
#' in a column in a table.
#' 
#' @param restable A \code{\link[base]{data.frame}} with at least a column named as in \code{add_to_table}.
#' @param colname_to_use The name of the column containing AIC values.
#' @param add_to_table If TRUE, add to the main table and return the main table. If FALSE, return just the Akaike Weights results.
#' @return \code{restable}, the modified table, or \code{wt_vBest}, the Akaike Weights results.
#' @export
#' @seealso \code{\link{calc_AIC_column}}, \code{\link{calc_AICc_column}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Burnham_Anderson_2002
#' @examples
#' test=1
#' 
AkaikeWeights_on_summary_table <- function(restable, colname_to_use="AIC", add_to_table=TRUE)
	{
	# Initialize AICvals, AICval_1
	AICvals = NA
	AICval_1 = NA
	
	# Assumes you will use ALL rows (unlike LRT, where some rows might not be nested)

	# Get the list of AIC values
	cmdstr = paste("AICvals = restable$", colname_to_use, sep="")
	eval(parse(text=cmdstr))
	
	AkaikeWeights = rep(NA, times=nrow(restable))
	
	for (j in 1:nrow(restable))
		{
		# Get the AIC value for this row
		cmdstr = paste("AICval_1 = restable$", colname_to_use, "[j]", sep="")
		eval(parse(text=cmdstr))
		
		AkaikeWeights[j] = getAIC_weight_for_model1(AICval_1, AICvals)
		#print(cmdstr)
		#print(AkaikeWeights)
		}
	
	wt_vBest = adf(AkaikeWeights)
	names(wt_vBest) = paste(colname_to_use, "_wt_vBest", sep="")
	wt_vBest
	
	if (add_to_table == TRUE)
		{
		restable = cbind(restable, wt_vBest)
		return(restable)
		} else {
		return(wt_vBest)
		}
	return("ERROR!")	
	}






#######################################################
# AICstats_2models
#######################################################
#' Calculate all the AIC and LRT stats between two models
#' 
#' The Likelihood Ratio Test (LRT) is a standard method for testing whether or not the data likelihood
#' conferred by a more complex is significantly better than the data likelihood conferred by
#' the simpler model.  See \code{\link{lrttest}} and \code{\link{lrttest_on_summary_table}} for more discussion.
#' 
#' See \cite{Burnham_Anderson_2002} and \url{http://www.brianomeara.info/tutorials/aic} for 
#' discussion of AIC and its uses.
#'
#' This function assumes that \code{LnL_1} and \code{numparams1} refer to the more complex model, and that \code{LnL_2}
#' and \code{numparams2} refer to the simpler model nested within the more complex one.
#' 
#' @param LnL_1 Log-likelihood of more complex model.
#' @param LnL_2 Log-likelihood of simpler complex model.
#' @param numparams1 Number of free parameters of the more complex model.
#' @param numparams2 Number of free parameters of the less complex model.
#' @return \code{LRT_AIC_results} A table of LRT and AIC results.
#' @export
#' @seealso \code{\link{lrttest}}, \code{\link{lrttest_on_summary_table}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Burnham_Anderson_2002
#' @examples
#' test=1
#' 
AICstats_2models <- function(LnL_1, LnL_2, numparams1, numparams2)
	{
	defaults='
	LnL_1=alt_LnL; LnL_2=null_LnL; numparams1=num_free_branches_nonclock; numparams2=num_free_branches_clock
	'
	
	# Run Likelihood Ratio Test (LRT)	
	LRT_result = lrttest(LnL_1, LnL_2, numparams1, numparams2, returnwhat="all")
	
	# AIC results
	AIC_results = NULL
	AIC_results_names = NULL
	
	# Combine into list
	LnLs = c(LnL_1, LnL_2)
	numparams_list = c(numparams1, numparams2)
	
	num_models = length(LnLs)
	
	# Make AIC list
	AIC_list = rep(NA, num_models)
	
	for (i in 1:num_models)
		{
		AIC_list[i] = getAIC(LnL=LnLs[i], numparams_list[i])
		}
	
	# Get relative weight of each model
	AICweight_list = rep(NA, num_models)
	for (i in 1:num_models)
		{
		AICweight_list[i] = getAIC_weight_for_model1(AICval_1=AIC_list[i], AICvals=AIC_list)
		}
	
	AICweight_ratio_model1 = get_AICweight_ratio_model1_over_model2(AICval_1=AIC_list[1], AICval_2=AIC_list[2])
	AICweight_ratio_model2 = get_AICweight_ratio_model1_over_model2(AICval_1=AIC_list[2], AICval_2=AIC_list[1])
	
	nums = seq(1, num_models, 1)
	
	
	# Put together the results
	AIC_results = c(AIC_results, AIC_list)
	AIC_names = paste("AIC", nums, sep="")
	
	AIC_results = c(AIC_results, AICweight_list)
	AICweight_list_names = paste("AICwt", nums, sep="")
	
	AIC_results = c(AIC_results, AICweight_ratio_model1, AICweight_ratio_model2)
	relprob_names = c("AICweight_ratio_model1", "AICweight_ratio_model2")
	
	# Names
	AIC_results_names = c(AIC_names, AICweight_list_names, relprob_names)
	
	AIC_results2 = matrix(AIC_results, nrow=1)
	AIC_results3 = adf2(AIC_results2)
	names(AIC_results3) = AIC_results_names
	
	AIC_results = AIC_results3
	
	
	# Merge in LRT
	LRT_AIC_results = cbind(LRT_result, AIC_results)
	
	return(LRT_AIC_results)
	}




# Get the ratio between the pairwise weights
#######################################################
# AkaikeWeights_and_Ratios_pairwise_on_summary_table_compared_to_ref
#######################################################
#' Get the ratio between the pairwise Akaike Weights
#' 
#' Given the relative likelihoods of the models, calculate the Akaike weight of the models. Akaike
#' weights sum to 1.
#' 
#' See \cite{Burnham_Anderson_2002} and \url{http://www.brianomeara.info/tutorials/aic} for 
#' discussion of AIC and its uses.
#' 
#' @param restable A \code{\link[base]{data.frame}} with at least columns named "LnL" and "nparams".
#' @param colname_to_use The name of the column containing AIC values.
#' @param ref_model What is the row of the reference model?  "best", "worst", or a row number.
#' @param add_to_table If TRUE, add to the main table and return the main table. If FALSE, return just the Akaike Weights results.
#' @return \code{restable}, the modified table, or \code{AICstats_pairwise}, the pairwise Akaike statistics.
#' @export
#' @seealso \code{\link{get_Akaike_weights_from_rel_likes_pairwise}}, \code{\link{get_Akaike_weights_from_rel_likes}}, 
#' \code{\link{rel_likes_from_deltaAICs}}, \code{\link{getAIC}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Burnham_Anderson_2002
#' @examples
#' test=1
#' 
#' tmptable = adf(c(40, 50, 60))
#' names(tmptable) = "AIC"
#' AkaikeWeights_and_Ratios_pairwise_on_summary_table_compared_to_ref(
#' restable=tmptable, colname_to_use="AIC", ref_model="best", add_to_table=TRUE)
#' 
AkaikeWeights_and_Ratios_pairwise_on_summary_table_compared_to_ref <- function(restable, colname_to_use="AIC", ref_model="best", add_to_table=TRUE)
	{
	# Initialize AICvals
	AICvals = NA
	
	# Extract the AIC vals
	cmdstr = paste("AICvals = restable$", colname_to_use, sep="")
	eval(parse(text=cmdstr))
	
	# Get the pairwise delta AICs
	get_rownum_ref_model(AICvals, ref_model="best")
	deltaAICs_pairwise = get_deltaAIC_pairwise_w_ref_model(AICvals, ref_model)

	# Get the pairwise relative likelihoods (likelihoods of the models given the list of models)
	# n x 2 matrix
	rel_likes_AIC_pairwise = rel_likes_from_deltaAICs_pairwise(deltaAICs_pairwise)
	
	# Calculate the pairwise Akaike weights
	Akaike_weights_pairwise = get_Akaike_weights_from_rel_likes_pairwise(rel_likes_AIC_pairwise)
	
	# Calculate the pairwise ratio (for only the row; the named row is obvious)
	Akaike_weight_ratios_pairwise = get_Akaike_weight_ratio_from_Akaike_pairwise_weights(Akaike_weights_pairwise)
	
	# Make the AICstats_pairwise table	
	AICstats_pairwise = cbind(deltaAICs_pairwise, rel_likes_AIC_pairwise, Akaike_weights_pairwise, Akaike_weight_ratios_pairwise)
	
	# Change the names to specify AIC, AICc, etc.
	names(AICstats_pairwise) = paste(colname_to_use, "_", names(AICstats_pairwise), sep="")
	
	
	if (add_to_table == TRUE)
		{
		restable = cbind(restable, AICstats_pairwise)
		return(restable)
		} else {
		return(AICstats_pairwise)
		}
	return("ERROR!")	
	}




#######################################################
# Some good functions relating to AIC
#######################################################

# Do LRT function
# Do the more-parameter model first as #1
# model #1 = alternative model (more params)
# model #2 = null model (fewer params, some of #1 fixed)


#######################################################
# lrttest
#######################################################
#' Calculate Likelihood Ratio Test (LRT)
#' 
#' The Likelihood Ratio Test (LRT) is a standard method for testing whether or not the data likelihood
#' conferred by a more complex is significantly better than the data likelihood conferred by
#' the simpler model, given a certain number of extra free parameters for the complex model.  The null hypothesis
#' is that there is no difference; rejection means that there is a statistically significant improvement in 
#' the more complex model.
#' 
#' The LRT only works for situations in which the simpler model is nested within the more complex model (i.e., by 
#' taking some parameters of the more complex model and forcing them to be fixed to a specific value).  In addition,
#' the LRT may be unreliable in data-poor situations, and inherits whatever difficulties there may be
#' in ML searches.  See \cite{Burnham_Anderson_2002} for discussion.
#'
#' This function assumes that \code{LnL_1} and \code{numparams1} refer to the more complex model, and that \code{LnL_2}
#' and \code{numparams2} refer to the simpler model nested within the more complex one.
#' 
#' @param LnL_1 Log-likelihood of more complex model.
#' @param LnL_2 Log-likelihood of simpler complex model.
#' @param numparams1 Number of free parameters of the more complex model.
#' @param numparams2 Number of free parameters of the less complex model.
#' @param returnwhat If "pval", just return the p-value.  If "all", return all of the intermediate outputs.
#' @return \code{pval} or \code{LRT_result2}. Depends on \code{returnwhat}.
#' @export
#' @seealso \code{\link{lrttest_on_summary_table}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Burnham_Anderson_2002
#' @examples
#' test=1
#' 
lrttest <- function(LnL_1, LnL_2, numparams1, numparams2, returnwhat="pval")
	{
	defaults='
	LnL_1 = restable$LnL[j]
	LnL_2 = restable$LnL[row_to_use_as_null]
	numparams1 = restable$nparams[j]
	numparams2 = restable$nparams[row_to_use_as_null]
	returnwhat=returnwhat
	'
	
	# Calculate delta LnL (difference in log-likelihoods
	delta_LnL = (LnL_1 - LnL_2)
	
	# D-statistic -- difference statistic
	# D-statistic defined at:
	# http://en.wikipedia.org/wiki/Likelihood-ratio_test
	Dstatistic = -2*LnL_2 + 2*LnL_1
	
	# Difference in number of parameters
	# (this is the number of degrees of freedom)
	delta_numparameters = numparams1 - numparams2
	
	# Original (before Feb. 15)
	#pval = 1-pchisq(q=2*delta_LnL, df=delta_numparameters)
	
	# Calculate probability of observed difference being bigger than expected
	# H0 = null hypothesis = D no bigger than expected by chance
	# This is the prob D > nullD
	# pval = pchisq(q=Dstatistic, df=delta_numparameters, lower.tail=TRUE)
	
	# This is the prob D <= nullD, i.e. the prob of the null hypothesis
	pval = pchisq(q=Dstatistic, df=delta_numparameters, lower.tail=FALSE)
	pval
	
	if (returnwhat == "pval")
		{
		return(pval)
		}

	if (returnwhat == "all")
		{
		LRT_result_names = c("alt", "null", "LnLalt", "LnLnull", "DFalt", "DFnull", "DF", "Dstatistic", "pval", "test", "tail")
		LRT_result = c("", "", LnL_1, LnL_2, numparams1, numparams2, delta_numparameters, Dstatistic, pval, "chi-squared", "one-tailed")
		
		LRT_result2 = matrix(LRT_result, nrow=1)
		LRT_result2 = adf2(LRT_result2)
		names(LRT_result2) = LRT_result_names
		
		LRT_result = LRT_result2
		return(LRT_result)
		}

	return(pval)
	}




#######################################################
# getAIC
#######################################################
#' Calculate AIC
#' 
#' Calculate AIC (Akaike Information Criterion).  Lower values
#' of AIC indicate some combination of better fit to the data and more parsimony in the model 
#' (fewer free parameters).
#' 
#' See \cite{Burnham_Anderson_2002} and \url{http://www.brianomeara.info/tutorials/aic} for 
#' discussion of AIC and its uses.
#' 
#' @param LnL The log-likelihood (typically negative, but may not be for continuous data).
#' @param numparams The number of parameters for each model.
#' @return \code{AICval} A vector of AIC results.
#' @export
#' @seealso \code{\link{calc_AIC_column}}, \code{\link{calc_AIC_column}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Burnham_Anderson_2002
#' @examples
#' LnL = -34.5
#' numparams = 2
#' getAIC(LnL, numparams)
#' 
#' LnL = -20.9
#' numparams = 3
#' getAIC(LnL, numparams)
#' 
#' # It turns out to work on lists, also
#' LnL = c(-34.5, -20.9)
#' numparams = c(2, 3)
#' getAIC(LnL, numparams)
#' 
getAIC <- function(LnL, numparams)
	{
	# Calculate AIC
	AICval = 2*numparams - 2*LnL
	return(AICval)
	}



#######################################################
# getAICc
#######################################################
#' Calculate AICc
#' 
#' Calculate AICc (Akaike Information Criterion).  Lower values
#' of AICc indicate some combination of better fit to the data and more parsimony in the model 
#' (fewer free parameters).
#' 
#' See \cite{Burnham_Anderson_2002} and \url{http://www.brianomeara.info/tutorials/AICc} for 
#' discussion of AICc and its uses.
#' 
#' @param LnL The log-likelihood (typically negative, but may not be for continuous data).
#' @param numparams The number of parameters for each model.
#' @param samplesize The number of data on which the model conferred likelihood.
#' @return \code{AICcval} A vector of AICc results.
#' @export
#' @seealso \code{\link{calc_AICc_column}}, \code{\link{calc_AICc_column}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/AICc}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Burnham_Anderson_2002
#' @examples
#' LnL = -34.5
#' numparams = 2
#' samplesize = 20
#' getAICc(LnL, numparams, samplesize)
#' 
#' LnL = -20.9
#' numparams = 3
#' samplesize = 20
#' getAICc(LnL, numparams, samplesize)
#' 
#' LnL = -34.5
#' numparams = 2
#' samplesize = 5
#' getAICc(LnL, numparams, samplesize)
#' 
#' LnL = -20.9
#' numparams = 3
#' samplesize = 5
#' getAICc(LnL, numparams, samplesize)
#' 
getAICc <- function(LnL, numparams, samplesize)
	{
	# http://en.wikipedia.org/wiki/Akaike_information_criterion#AICc
	# Sample size is typically number of tips...
	
	if (numparams >= samplesize)
		{
		stop("ERROR!  You cannot have more parameters than samples in AICc, you get bizarre results.")
		}
	
	# Calculate AIC
	AICval = 2*numparams - 2*LnL
	
	# Correction for finite sample size
	correction_val = (2*numparams*(numparams+1)) / (samplesize - numparams - 1)
	
	# 
	AICc_val = AICval + correction_val
	
	return(AICc_val)
	}




#######################################################
# AIC stats for all models in a column
#######################################################

#######################################################
# get_deltaAIC
#######################################################
#' Calculate deltaAIC
#' 
#' Calculate deltaAIC (Akaike Information Criterion), the absolute difference between the 
#' best model (lowest AIC) and other models.
#' 
#' See \cite{Burnham_Anderson_2002} and \url{http://www.brianomeara.info/tutorials/aic} for 
#' discussion of AIC and its uses.
#' 
#' @param AICvals A vector of AIC values.
#' @return \code{deltaAICs} A vector of deltaAICs.
#' @export
#' @seealso \code{\link{rel_likes_from_deltaAICs}}, \code{\link{getAIC}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Burnham_Anderson_2002
#' @examples
#' test=1
#' 
#' AICvals = c(40, 50, 60)
#' get_deltaAIC(AICvals)
#' 
get_deltaAIC <- function(AICvals)
	{
	# Take the first, if there is a tie
	bestAIC = min(AICvals, na.rm=TRUE)[1]
	
	deltaAICs = AICvals - bestAIC
	
	return(deltaAICs)
	}



# Get the relative likelihoods
# l. To calculate them, for each model first calculate the relative likelihood of the model,
# which is just exp( -0.5 * deltaAIC score for that model).
# http://www.brianomeara.info/tutorials/aic
# rel_likes_from_deltaAICs(deltaAICs)

#######################################################
# rel_likes_from_deltaAICs
#######################################################
#' Calculate the relative likelihoods of the models, from the deltaAIC
#' 
#' Given deltaAIC (Akaike Information Criterion), the absolute difference between the 
#' best model (lowest AIC) and other models, calculate the relative likelihoods of the 
#' models.
#' 
#' See \cite{Burnham_Anderson_2002} and \url{http://www.brianomeara.info/tutorials/aic} for 
#' discussion of AIC and its uses.
#' 
#' @param deltaAICs A vector of deltaAIC values.
#' @return \code{rel_likes_AIC} A vector of relative likelihoods.
#' @export
#' @seealso \code{\link{get_Akaike_weights_from_rel_likes}}, \code{\link{rel_likes_from_deltaAICs}}, \code{\link{getAIC}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Burnham_Anderson_2002
#' @examples
#' test=1
#' 
#' AICvals = c(40, 50, 60)
#' deltaAICs = get_deltaAIC(AICvals)
#' deltaAICs
#' 
#' rel_likes_AIC = rel_likes_from_deltaAICs(deltaAICs)
#' rel_likes_AIC
#' 
rel_likes_from_deltaAICs <- function(deltaAICs)
	{
	rel_likes_AIC = exp(-0.5 * deltaAICs)
	return(rel_likes_AIC)
	}



# http://www.brianomeara.info/tutorials/aic
# Akaike weights are can be used in model averaging. They represent the relative likelihood of a model. 
# To calculate them, for each model first calculate the relative likelihood of the model, which is just 
# exp( -0.5 * ∆AIC score for that model). The Akaike weight for a model is this value divided by the sum 
# of these values across all models.

#######################################################
# get_Akaike_weights_from_rel_likes
#######################################################
#' Calculate the Akaike Weights, from the relative likelihoods of the models
#' 
#' Given the relative likelihoods of the models, calculate the Akaike weight of the models. Akaike
#' weights sum to 1.
#' 
#' See \cite{Burnham_Anderson_2002} and \url{http://www.brianomeara.info/tutorials/aic} for 
#' discussion of AIC and its uses.
#' 
#' @param rel_likes_AIC A vector of relative likelihoods.
#' @return \code{Akaike_weights} A vector of Akaike Weights.
#' @export
#' @seealso \code{\link{get_Akaike_weights_from_rel_likes}}, \code{\link{rel_likes_from_deltaAICs}}, \code{\link{getAIC}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Burnham_Anderson_2002
#' @examples
#' test=1
#' 
#' AICvals = c(40, 50, 60)
#' deltaAICs = get_deltaAIC(AICvals)
#' deltaAICs
#' 
#' Akaike_weights = rel_likes_from_deltaAICs(deltaAICs)
#' Akaike_weights
#' 
get_Akaike_weights_from_rel_likes <- function(rel_likes_AIC)
	{
	# Take the sum of the relative likelihoods
	sum_rel_likes = sum(rel_likes_AIC)
	
	# Initialize Akaike weights
	Akaike_weights = rep(NA, length(rel_likes_AIC))
	
	# Divide each rel_like by the sum
	Akaike_weights = rel_likes_AIC / sum_rel_likes
	
	return(Akaike_weights)
	}


# Get rownum of named model
#######################################################
# get_rownum_ref_model
#######################################################
#' Get rownum of named model
#' 
#' Find the row number of the best model according to AIC, the worst model according to AIC,
#' or just takes the row number if that is what was input.
#' 
#' @param AICvals A vector of AIC values.
#' @param ref_model What is the row of the reference model?  "best", "worst", or a row number.
#' @return \code{ref_model_num} The 
#' @export
#' @seealso \code{\link[stats]{convolve}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#'	 @cite FosterIdiots
#' @examples
#' test=1
#' 
get_rownum_ref_model <- function(AICvals, ref_model="best")
	{
	# Get AIC of named model
	nums = 1:length(AICvals)
	if ( (ref_model == "best") || (ref_model == "worst"))
		{
		if (ref_model == "best")
			{
			# Best AIC is the lowest
			best_AICval = min(AICvals, na.rm=TRUE)
			best_AICnum = nums[AICvals == best_AICval][1]
			ref_model_num = best_AICnum
			}
		if (ref_model == "worst")
			{
			# Worst AIC is the highest
			worst_AICval = max(AICvals, na.rm=TRUE)
			worst_AICnum = nums[AICvals == worst_AICval][1]
			ref_model_num = worst_AICnum
			}
		} else {
		# Just use the number given
		ref_model_num = ref_model
		}
	return(ref_model_num)
	}




#######################################################
# AIC-based stats for pairwise comparisons to a reference model ("ref")
#######################################################

#######################################################
# get_deltaAIC_pairwise_w_ref_model
#######################################################
#' Calculate deltaAIC
#' 
#' Calculate deltaAIC (Akaike Information Criterion), the absolute difference between the 
#' best model (lowest AIC) and other models.  This function does it pairwise only, with a reference model.
#' 
#' See \cite{Burnham_Anderson_2002} and \url{http://www.brianomeara.info/tutorials/aic} for 
#' discussion of AIC and its uses.
#' 
#' @param AICvals A vector of AIC values.
#' @param ref_model What is the row of the reference model?  "best", "worst", or a row number.
#' @return \code{deltaAICs_pairwise} A 2-column \code{\link[base]{data.frame}} of pairwise deltaAICs for each row (column 1) and the reference model (column 2).
#' @export
#' @seealso \code{\link{get_deltaAIC}}, \code{\link{rel_likes_from_deltaAICs}}, \code{\link{getAIC}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Burnham_Anderson_2002
#' @examples
#' test=1
#' 
#' AICvals = c(40, 50, 60)
#' get_deltaAIC(AICvals)
#' get_deltaAIC_pairwise_w_ref_model(AICvals, ref_model="best")
#' 
#' 
get_deltaAIC_pairwise_w_ref_model <- function(AICvals, ref_model="best")
	{
	# Here, just assume there are TWO rows in each case
	deltaAICs_pairwise = matrix(data=NA, nrow=length(AICvals), ncol=2)
	
	ref_model_num = get_rownum_ref_model(AICvals, ref_model)
	
	# Do deltaAIC for each pair
	for (j in 1:length(AICvals))
		{
		tmp_AICvals = c(AICvals[j], AICvals[ref_model_num])
		tmp_deltaAICs = get_deltaAIC(tmp_AICvals)
		deltaAIC_for_this_row = tmp_deltaAICs[j]
		deltaAICs_pairwise[j,] = tmp_deltaAICs
		}
	
	deltaAICs_pairwise = adf(deltaAICs_pairwise)
	names(deltaAICs_pairwise) = c("deltaAIC_row", "deltaAIC_ref")
	
	# n x 2 matrix
	return(deltaAICs_pairwise)
	}



# n x 2 matrix
#######################################################
# rel_likes_from_deltaAICs_pairwise
#######################################################
#' Calculate the relative likelihoods of the models, from the deltaAICs, pairwise
#' 
#' Given deltaAIC (Akaike Information Criterion), the absolute difference between the 
#' best model (lowest AIC) and other models, calculate the relative likelihoods of the 
#' models.
#' 
#' See \cite{Burnham_Anderson_2002} and \url{http://www.brianomeara.info/tutorials/aic} for 
#' discussion of AIC and its uses.
#' 
#' @param deltaAICs_pairwise A vector of AIC values.
#' @return \code{rel_likes_AIC_pairwise} A \code{\link[base]{data.frame}} of relative likelihoods for each row (column 1) and the reference model (column 2).
#' @export
#' @seealso \code{\link{get_Akaike_weights_from_rel_likes}}, \code{\link{rel_likes_from_deltaAICs}}, \code{\link{getAIC}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Burnham_Anderson_2002
#' @examples
#' test=1
#' 
#' AICvals = c(40, 50, 60)
#' deltaAICs = get_deltaAIC_pairwise_w_ref_model(AICvals, ref_model="best")
#' deltaAICs
#' 
#' rel_likes_AIC_pairwise = rel_likes_from_deltaAICs_pairwise(deltaAICs)
#' rel_likes_AIC_pairwise
#' 
rel_likes_from_deltaAICs_pairwise <- function(deltaAICs_pairwise)
	{
	rel_likes_AIC_pairwise = exp(-0.5 * deltaAICs_pairwise)
	rel_likes_AIC_pairwise = adf(rel_likes_AIC_pairwise)
	names(rel_likes_AIC_pairwise) = c("rel_like_pairwise_row", "rel_like_pairwise_ref")
	return(rel_likes_AIC_pairwise)
	}

	

#######################################################
# get_Akaike_weights_from_rel_likes_pairwise
#######################################################
#' Calculate the Akaike Weights, from the relative likelihoods of the models
#' 
#' Given the relative likelihoods of the models, calculate the Akaike weight of the models. Akaike
#' weights sum to 1.
#' 
#' See \cite{Burnham_Anderson_2002} and \url{http://www.brianomeara.info/tutorials/aic} for 
#' discussion of AIC and its uses.
#' 
#' @param rel_likes_AIC_pairwise A 2-column \code{data.frame} of relative likelihoods of each pair of models.
#' @return \code{Akaike_weights_pairwise} A \code{\link[base]{data.frame}} of Akaike Weights for each row (column 1) and the reference
#' model (column 2). Note that only 2 models are being compared in each row, not all of them, as in \code{\link{get_Akaike_weights_from_rel_likes}}.
#' @export
#' @seealso \code{\link{get_Akaike_weights_from_rel_likes}}, \code{\link{rel_likes_from_deltaAICs}}, \code{\link{getAIC}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Burnham_Anderson_2002
#' @examples
#' test=1
#' 
#' AICvals = c(40, 50, 60)
#' deltaAICs = get_deltaAIC_pairwise_w_ref_model(AICvals, ref_model="best")
#' deltaAICs
#' 
#' rel_likes_AIC_pairwise = rel_likes_from_deltaAICs_pairwise(deltaAICs)
#' rel_likes_AIC_pairwise
#' 
#' Akaike_weights_pairwise = get_Akaike_weights_from_rel_likes_pairwise(rel_likes_AIC_pairwise)
#' Akaike_weights_pairwise
#' 
get_Akaike_weights_from_rel_likes_pairwise <- function(rel_likes_AIC_pairwise)
	{
	# Take the sum of the relative likelihoods
	sums_rel_likes = rowSums(rel_likes_AIC_pairwise)

	# Divide each row by the rowsums
	Akaike_weights_pairwise = rel_likes_AIC_pairwise / sums_rel_likes
	Akaike_weights_pairwise = adf(Akaike_weights_pairwise)
	
	names(Akaike_weights_pairwise) = c("Akaike_pairwt_row", "Akaike_pairwt_ref")
	
	return(Akaike_weights_pairwise)
	}



# Get the ratio between the pairwise weights
#######################################################
# get_Akaike_weight_ratio_from_Akaike_pairwise_weights
#######################################################
#' Get the ratio between the pairwise Akaike Weights
#' 
#' Given the relative likelihoods of the models, calculate the Akaike weight of the models. Akaike
#' weights sum to 1.
#' 
#' See \cite{Burnham_Anderson_2002} and \url{http://www.brianomeara.info/tutorials/aic} for 
#' discussion of AIC and its uses.
#' 
#' @param Akaike_weights_pairwise A 2-column \code{data.frame} of Akaike Weights for each pair of models.
#' @return \code{Akaike_weight_ratios_pairwise} A \code{\link[base]{data.frame}} of Akaike Weight Ratios for each row (column 1) and the reference
#' model (column 2). Note that only 2 models are being compared in each row, not all of them.
#' @export
#' @seealso \code{\link{get_Akaike_weights_from_rel_likes_pairwise}}, \code{\link{get_Akaike_weights_from_rel_likes}}, 
#' \code{\link{rel_likes_from_deltaAICs}}, \code{\link{getAIC}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Burnham_Anderson_2002
#' @examples
#' test=1
#' 
#' AICvals = c(40, 50, 60)
#' deltaAICs = get_deltaAIC_pairwise_w_ref_model(AICvals, ref_model="best")
#' deltaAICs
#' 
#' rel_likes_AIC_pairwise = rel_likes_from_deltaAICs_pairwise(deltaAICs)
#' rel_likes_AIC_pairwise
#' 
#' Akaike_weights_pairwise = get_Akaike_weights_from_rel_likes_pairwise(
#' rel_likes_AIC_pairwise)
#' Akaike_weights_pairwise
#' 
#' Akaike_weight_ratios_pairwise = get_Akaike_weight_ratio_from_Akaike_pairwise_weights(
#' Akaike_weights_pairwise)
#' Akaike_weight_ratios_pairwise
#' 
get_Akaike_weight_ratio_from_Akaike_pairwise_weights <- function(Akaike_weights_pairwise)
	{
	# 2nd column is the null/named row; divid by this number
	
	Akaike_weight_ratios_pairwise = Akaike_weights_pairwise[,1] / Akaike_weights_pairwise[,2]
	#Akaike_weight_ratios_pairwise = Akaike_weights_pairwise$Akaike_pairwt_row / Akaike_weights_pairwise$Akaike_pairwt_ref
	
	Akaike_weight_ratios_pairwise = adf(Akaike_weight_ratios_pairwise)
	names(Akaike_weight_ratios_pairwise) = "Akaike_weight_ratio_row_to_ref"
	
	return(Akaike_weight_ratios_pairwise)
	}






# Older methods
# NOTE: AIC weight = the model probability, between 0 and 1.  The weights ALWAYS sum to 1

# AIC weight for model1
# based on: http://www.brianomeara.info/tutorials/aic

#######################################################
# getAIC_weight_for_model1
#######################################################
#' Calculate Akaike Weight
#' 
#' See \cite{Burnham_Anderson_2002} and \url{http://www.brianomeara.info/tutorials/aic} for 
#' discussion of AIC and its uses.
#'
#' @param AICval_1 The AIC of the model of interest.
#' @param AICvals The AICs of all the models being compared.
#' @return \code{AICweight} AICweight for the models.
#' @export
#' @seealso \code{\link{AkaikeWeights_on_summary_table}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Burnham_Anderson_2002
#' @examples
#' test=1
#' 
#' AICval_1 = 20
#' AICvals = c(20,30,40)
#' getAIC_weight_for_model1(AICval_1, AICvals)
#' 
getAIC_weight_for_model1 <- function(AICval_1, AICvals)
	{
	if (AICval_1 < 0) {stop("ERROR: AIC cannot be negative, you probably put in log-likelihood (LnL)")}
	
	# find better AIC model (i.e., the lowest AIC
	# use first, if two match
	better_model_index = (1:length(AICvals))[AICvals == min(AICvals)][1]
	
	deltaAIC = AICval_1 - AICvals[better_model_index]
	deltaAICs = AICvals - AICvals[better_model_index]

	# relative weight of model1 (lower # of paramters)
	# (i.e. likelihood of the model, among others)
	relative_likelihood_model1 = exp( -0.5 * (deltaAIC) )
	relative_likelihood_all_models = exp( -0.5 * (deltaAICs) )

	sum_likes = sum(relative_likelihood_all_models)
	AICweight = relative_likelihood_model1 / sum_likes

	return(AICweight)
	}




#######################################################
# get_AICweight_ratio_model1_over_model2
#######################################################
#' Calculate ratio of Akaike Weights
#' 
#' See \cite{Burnham_Anderson_2002} and \url{http://www.brianomeara.info/tutorials/aic} for 
#' discussion of AIC and its uses.
#'
#' @param AICval_1 The AIC of the model of interest.
#' @param AICval_2 The AIC of another model of interest, for a pairwise comparison.
#' @return \code{AICweight_ratio_model1} Ratio of Akaike Weights.
#' @export
#' @seealso \code{\link{AkaikeWeights_on_summary_table}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Burnham_Anderson_2002
#' @examples
#' test=1
#' 
#' AICval_1 = 20
#' AICval_2 = 30
#' get_AICweight_ratio_model1_over_model2(AICval_1, AICval_2)
#' 
get_AICweight_ratio_model1_over_model2 <- function(AICval_1, AICval_2)
	{
	if (AICval_1 < 0) {stop("ERROR: AIC cannot be negative, you probably put in log-likelihood (LnL)")}
	if (AICval_2 < 0) {stop("ERROR: AIC cannot be negative, you probably put in log-likelihood (LnL)")}

	AICvals = c(AICval_1, AICval_2)
	
	
	# relative weight of model2 (lower # of paramters
	AIC_weight_model2 = getAIC_weight_for_model1(AICval_2, AICvals)

	# relative weight of model1 (higher # of paramters
	AIC_weight_model1 = getAIC_weight_for_model1(AICval_1, AICvals)
	
	AICweight_ratio_model1 = AIC_weight_model1 / AIC_weight_model2

	return(AICweight_ratio_model1)
	}


#######################################################
# get_AICweight_ratio_model1_over_model2
#######################################################
#' Calculate relative probability of model 1 (=Akaike Weight)
#' 
#' See \cite{Burnham_Anderson_2002} and \url{http://www.brianomeara.info/tutorials/aic} for 
#' discussion of AIC and its uses.
#'
#' @param AICval_1 The AIC of the model of interest.
#' @param AICval_2 The AIC of another model of interest, for a pairwise comparison.
#' @return \code{relative_prob_model1} Akaike Weight of model 1.
#' @export
#' @seealso \code{\link{AkaikeWeights_on_summary_table}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Burnham_Anderson_2002
#' @examples
#' test=1
#' 
#' AICval_1 = 20
#' AICval_2 = 30
#' get_relative_prob_model1old(AICval_1, AICval_2)
#' 
get_relative_prob_model1old <- function(AICval_1, AICval_2)
	{
	if (AICval_1 < 0) {stop("ERROR: AIC cannot be negative, you probably put in log-likelihood (LnL)")}
	if (AICval_2 < 0) {stop("ERROR: AIC cannot be negative, you probably put in log-likelihood (LnL)")}

	AICvals = c(AICval_1, AICval_2)
	
	# find better AIC model (take first, if the same)
	better_model_index = (1:2)[AICvals == min(AICvals)][1]
	
	# relative weight of model2 (lower # of paramters
	relative_weight_model2 = exp(-0.5 * (AICval_2-AICvals[better_model_index]))


	# relative weight of model1 (higher # of paramters
	relative_weight_model1 = exp(-0.5 * (AICval_1-AICvals[better_model_index]))
	
	relative_prob_model1 = relative_weight_model1 / (relative_weight_model1 + relative_weight_model2)

	return(relative_prob_model1)
	}


#######################################################
# get_relative_prob_model2old
#######################################################
#' Calculate relative probability of model 1 (Akaike Weight)
#' 
#' See \cite{Burnham_Anderson_2002} and \url{http://www.brianomeara.info/tutorials/aic} for 
#' discussion of AIC and its uses.
#'
#' This is an older version of \code{\link{get_relative_prob_model1old}}, kept for back-compatibility.
#'
#' @param AICval_1 The AIC of the model of interest.
#' @param AICval_2 The AIC of another model of interest, for a pairwise comparison.
#' @return \code{relative_prob_model1} Akaike Weight of model 1.
#' @export
#' @seealso \code{\link{AkaikeWeights_on_summary_table}}, \code{\link{get_relative_prob_model1old}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Burnham_Anderson_2002
#' @examples
#' test=1
#' 
#' AICval_1 = 20
#' AICval_2 = 30
#' get_relative_prob_model1old(AICval_1, AICval_2)
#' 
get_relative_prob_model2old <- function(AICval_1, AICval_2)
	{
	if (AICval_1 < 0) {stop("ERROR: AIC cannot be negative, you probably put in log-likelihood (LnL)")}
	if (AICval_2 < 0) {stop("ERROR: AIC cannot be negative, you probably put in log-likelihood (LnL)")}

	AICvals = c(AICval_1, AICval_2)
	
	# find better AIC model (take first, if the same)
	better_model_index = (1:2)[AICvals == min(AICvals)][1]
	
	# relative weight of model2 (lower # of paramters
	relative_weight_model2 = exp(-0.5 * (AICval_2-AICvals[better_model_index]))


	# relative weight of model1 (higher # of paramters
	relative_weight_model1 = exp(-0.5 * (AICval_1-AICvals[better_model_index]))
	
	relative_prob_model2 = relative_weight_model2 / (relative_weight_model1 + relative_weight_model2)
	
	return(relative_prob_model2)
	}




# Do the more-parameter model first as #1
# model #1 = alternative model (more params)
# model #2 = null model (fewer params, some of #1 fixed)

#######################################################
# AICstats_2models
#######################################################









# Print tables to PDF

# Print to PDF via xtable
# source: http://tex.stackexchange.com/questions/15013/generate-a-pdf-containing-r-output-inside-latex-table
# 
#######################################################
# pdfit
#######################################################
#' Print a table to LaTeX format
#' 
#' This function prints a table to PDF via \code{\link[xtable]{xtable}} and the LaTeX \code{pdflatex} function.  It will only 
#' work if you have command-line LaTeX installed.
#' 
#' This function was inspired by \url{http://tex.stackexchange.com/questions/15013/generate-a-pdf-containing-r-output-inside-latex-table}.
#' 
#' @param table_vals A table, hopefully produced by \code{\link{conditional_format_table}}.
#' @param file_prefix The prefix for the output PDF and the intermediate files.
#' @param size Font size, overriding \code{getOption("xtable.size")}. Default is "tiny" (with backslashes).  
#' You can also try "small".  Input \code{NULL} (without quotes or backslashes) for medium.  (\code{NULL} is the options default.) 
#' @param gettex If TRUE, the \code{tex} code for the table is returned.
#' @param caption A caption, if desired.
#' @return \code{texfile} The filename of the \code{tex} file.
#' @export
#' @seealso \code{\link{pdftable}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
#' 
#' # Setup data
#' \dontrun{
#' data = c(2.768443, 1.869964, 5.303702, 4.733483,  2.123816,  
#' 18.551051, 5.483625,  3.590745,  18.772389)
#' result = matrix(data, nrow=3, byrow=TRUE)
#' result = as.data.frame(result)
#' names(result) = c("CV", "LCB", "UCB")
#' rownames(result) = c("within", "between", "total")
#' result
#' pdfit(table_vals=result)#' }
#' 
pdfit <- function(table_vals, file_prefix="tmptable", size="\\tiny", gettex=FALSE, caption=NULL)
	{
	#require(xtable)
	
	cat("\nNOTE: the pdfit() function will only work if you have 'latex' installed on your machine.\n", sep="")
	
	#tab = table_vals
	
	# Make a temporary tex file
	texfile <- paste(file_prefix, ".tex", sep="")
	cat("\\documentclass{article}\n\\usepackage{changepage}\n\\begin{document}\n\\begin{figure}\n\\begin{adjustwidth}{-2cm}{}\n", file=texfile)
	
	# Stick in the table
	getOption("xtable.size")
	options(xtable.size=size)
	getOption("xtable.size")
	
	print(xtable(table_vals, size, caption), include.rownames=FALSE, floating=FALSE, 
		 file=texfile, append=TRUE)
	cat("\\end{adjustwidth}\n\\end{figure}\n\\end{document}\n", file=texfile, append=TRUE)

	# Print to PDF
	system(paste("pdflatex", texfile))

	options(xtable.size=NULL)
	
	if (gettex == TRUE)
		{
		return(texfile)
		}
	return(texfile)
	}




# Default directory for temp files is ~
#######################################################
# pdfit
#######################################################
#' Print a table to LaTeX format
#' 
#' This function prints a table to PDF via \code{\link{pdfit}}, which calls \code{\link[xtable]{xtable}}
#' and the LaTeX \code{pdflatex} function.  It will only work if you have command-line LaTeX installed.
#' 
#' This function was inspired by \url{http://tex.stackexchange.com/questions/15013/generate-a-pdf-containing-r-output-inside-latex-table}.
#' 
#' @param table_vals A table, hopefully produced by \code{\link{conditional_format_table}}.
#' @param pdffn The filename for the output PDF (and the prefix for the intermediate files).
#' @param size Font size, overriding \code{getOption("xtable.size")}. Default is "tiny" (with backslashes).  
#' You can also try "small".  Input \code{NULL} (without quotes or backslashes) for medium.  (\code{NULL} is the options default.) 
#' @param tmpdir The location for the temporary files.
#' @param openPDF If \code{TRUE}, open the PDF via a \code{\link[base]{system}} command.
#' @param caption A caption, if desired.
#' @return \code{pdffn} The filename of the PDF file.
#' @export
#' @seealso \code{\link{pdfit}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
#' 
#' # Setup data
#' \dontrun{
#' data = c(2.768443, 1.869964, 5.303702, 4.733483,  2.123816,  
#' 18.551051, 5.483625,  3.590745,  18.772389)
#' result = matrix(data, nrow=3, byrow=TRUE)
#' result = as.data.frame(result)
#' names(result) = c("CV", "LCB", "UCB")
#' rownames(result) = c("within", "between", "total")
#' result
#' pdftable(table_vals=result)#' }
#' 
pdftable <- function(table_vals, pdffn="tmptable.pdf", size="\\tiny", tmpdir="~", openPDF=TRUE, caption=NULL)
	{
	#require(xtable)
	
	# pdfif no like it if columns are lists
	#table_vals = unlist_df3(table_vals)
	
	# Set up tmp filenames
	pdf_prefix = get_path_last(pdffn)
	pdf_prefix = get_fn_prefix(pdf_prefix)
	#pdf_prefix = paste(tmpdir, "/", pdf_prefix, sep="")

	# Generate table in .tex form in the tmpdir
	pdfit(table_vals, file_prefix=pdf_prefix, size=size, gettex=FALSE, caption)
	
	# Copy the PDF to pdffn
	cmdstr = paste("cp ", pdf_prefix, ".pdf ", pdffn, sep="")
	system(cmdstr)
	
	if (openPDF)
		{
		cmdstr = paste("open ", pdffn, sep="")
		system(cmdstr)
		}
	
	return(pdffn)
	}


















#######################################################
# TREE FUNCTIONS
#######################################################





#######################################################
# get_daughters
#######################################################
#' Get all the direct daughters nodes of a node
#' 
#' @param nodenum The node number to get the daughters of
#' @param t An ape phylo object
#' @return \code{daughter_nodenums} List of the daughter node numbers
#' @export
#' @seealso \code{\link{findall}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_daughters <- function(nodenum, t)
	{
	daughter_edgenums = findall(nodenum, t$edge[,1])
	daughter_nodenums = t$edge[,2][daughter_edgenums]
	return(daughter_nodenums)
	}




# Get indices of all matches to a list
#######################################################
# findall
#######################################################
#' Get indices of all matches to a list
#'
#' Just a handy shortcut function
#' 
#' @param what The item to find
#' @param inlist The list to search in 
#' @return \code{matching_indices} List of the matching indices
#' @export
#' @seealso \code{\link{get_daughters}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
findall <- function(what, inlist)
	{
	TFmatches = inlist == what
	indices = 1:length(inlist)
	matching_indices = indices[TFmatches]
	return(matching_indices)
	}





#######################################################
# get_parent
#######################################################
#' Get the direct parent node of a node
#' 
#' @param nodenum The node number to get the parent of
#' @param t An ape phylo object
#' @return \code{parent_nodenum}The parent node number
#' @export
#' @seealso \code{\link{findall}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_parent <- function(nodenum, t)
	{
	matching_edges = findall(nodenum, t$edge[,2])
	parent_nodenum = t$edge[,1][matching_edges][1]
	return(parent_nodenum)
	}


#######################################################
# get_level
#######################################################
#' Get a node's level in the tree
#'
#' Finds how many nodes deep a node is.
#' 
#' @param nodenum The node number to get the parent of
#' @param t An ape phylo object
#' @param tmplevel A starting level (the function is recursive)
#' @return \code{tmplevel} The level of the node.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_level <- function(nodenum, t, tmplevel=0)
	{
	parent_nodenum = get_parent(nodenum, t)
	if (is.na(parent_nodenum))
		{
		#tmplevel = 0
		return(tmplevel)
		}
	else
		{
		#print(paste("parent_nodenum: ", parent_nodenum, " level: ", tmplevel, sep=""))
		tmplevel = tmplevel + 1
		tmplevel = get_level(parent_nodenum, t, tmplevel)
		return(tmplevel)
		}
	# If an error occurs
	return(NA)
	}


#######################################################
# get_TF_tips
#######################################################
#' Get TRUE/FALSE for nodes being tips
#'
#' A utility function that returns \code{TRUE}/\code{FALSE} for whether or not each node is a tip.
#' 
#' @param obj An ape phylo object
#' @return \code{TF_tips} The \code{TRUE}/\code{FALSE} list for each tip.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}, \code{\link{match_list1_in_list2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_TF_tips <- function(obj)
	{
	# Get TF for nodes being tips
	
	# BIG CHANGE?
	#TF_tips = match_list1_in_list2(1:length(dists_from_root), obj$tip.label)
	TF_tips = match_list1_in_list2(1:length(obj$edge), 1:length(obj$tip.label))
	#TF_tips = obj$tip.label[TF_tips_indices]
	return(TF_tips)
	}



#######################################################
# get_TF_tips
#######################################################
#' Get TRUE/FALSE for nodes being tips
#'
#' A utility function that returns indices (node numbers) of the tips. This mostly saves typing.
#' 
#' @param obj An ape phylo object
#' @return \code{tip_indices} The node numbers of the tips.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}, \code{\link[ape]{phylo}}, \code{\link{get_indices_of_branches_under_tips}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_indices_of_tip_nodes <- function(obj)
	{
	tip_indices = 1:length(obj$tip.label)
	return(tip_indices)
	}

#######################################################
# get_indices_of_branches_under_tips
#######################################################
#' Get the indices of the branches (row number in edge matrix) below each tip
#'
#' A utility function. Gets the indices of the branches (row number in edge matrix) below each tip.
#' 
#' @param obj An \code{\link[ape]{ape}} \code{\link[ape]{phylo}} object
#' @return \code{branchnums_under_tips} The indices of the branches (row number in edge matrix) below each tip.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}, \code{\link{get_indices_of_tip_nodes}}, \code{\link{get_indices_where_list1_occurs_in_list2_noNA}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_indices_of_branches_under_tips <- function(obj)
	{
	tip_indices = get_indices_of_tip_nodes(obj)
	branchnums_under_tips = get_indices_where_list1_occurs_in_list2_noNA(tip_indices, obj$edge[, 2])
	return(branchnums_under_tips)
	}




#######################################################
# get_node_ages_of_tips
#######################################################
#' Get the ages of each tip above the root
#'
#' A utility function.
#' 
#' @param obj An ape phylo object
#' @return \code{TF_tips} The age (from the root) of each tip.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_node_ages_of_tips <- function(obj)
	{
	TF_tips = get_TF_tips(obj)
	root_node_num = get_nodenum_structural_root(obj)
	dists_from_root = dist.nodes(obj)[root_node_num, ]
	node_ages_of_tips = dists_from_root[TF_tips]
	return(node_ages_of_tips)
	}


#######################################################
# get_all_node_ages
#######################################################
#' Get the ages of all the nodes in the tree (above the root)
#'
#' A utility function. Use of \code{\link[ape]{dist.nodes}} may be slow.
#' 
#' @param obj An ape phylo object
#' @return \code{TF_tips} The age (from the root) of each node.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_all_node_ages <- function(obj)
	{
	node_ages = dist.nodes(obj)[get_nodenum_structural_root(obj), ]
	return(node_ages)
	}


#######################################################
# get_max_height_tree
#######################################################
#' Get the maximum age of all the nodes (above the root)
#'
#' I.e., the distance of the highest node above the root.  A utility function. 
#' Use of \code{\link[ape]{dist.nodes}} may be slow.
#' 
#' @param obj An ape phylo object
#' @return \code{max_height} The age (from the root) of the highest node.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_max_height_tree <- function(obj)
	{
	max_height = max(get_node_ages_of_tips(obj))
	return(max_height)
	}



#######################################################
# get_edge_times_before_present
#######################################################
#' Get the times of the top and bottom of each edge
#'
#' A utility function. 
#' 
#' @param t An ape phylo object
#' @return \code{edge_times_bp} A 2-column matrix with the age (from the present) of the top and bottom of each edge.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
get_edge_times_before_present <- function(t)
	{
	#height above root
	hts_at_end_of_branches_aka_at_nodes = t$edge.length
	hts_at_end_of_branches_aka_at_nodes = get_all_node_ages(t)
	h = hts_at_end_of_branches_aka_at_nodes

	# times before present, below (ultrametric!) tips
	# numbers are positive, i.e. in millions of years before present
	#                       i.e. mybp, Ma
	times_before_present = get_max_height_tree(t) - h

	
	# fill in the ages of each node for the edges
	edge_ages = t$edge
	edge_ages[,1] = h[t$edge[,1]]	# bottom of branch
	edge_ages[,2] = h[t$edge[,2]]	# top of branch

	# fill in the times before present of each node for the edges
	edge_times_bp = t$edge
	edge_times_bp[,1] = times_before_present[t$edge[,1]]	# bottom of branch
	edge_times_bp[,2] = times_before_present[t$edge[,2]]	# top of branch
	
	return(edge_times_bp)
	}








#######################################################
# extend_tips_to_ultrametricize
#######################################################
#' Take a tree, extend all tips (including fossils) up to 0.0 my before present
#' 
#' Makes tree precisely ultrametric by extending the terminal branches up to the highest tip (which is treated as 0 my before present).
#'
#' This function ADDS the time_before_present to everything, including fossils.  You have been warned.
#' 
#' @param obj An \code{\link[ape]{ape}} \code{\link[ape]{phylo}} object.
#' @param age_of_root The length of the branch below the root. Default 0.
#' @param tips_end_at_this_date The tips can be set to something other than 0, if desired.  (This could produce negative branclengths, however.)
#' @return \code{obj} The corrected phylogeny
#' @export
#' @seealso \code{\link[ape]{read.tree}}, \code{\link{prt}}, \code{\link{average_tr_tips}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
extend_tips_to_ultrametricize <- function(obj, age_of_root=0, tips_end_at_this_date=NA)
	{
	#print("node ages of tips:")
	tip_ages = age_of_root + get_node_ages_of_tips(obj)
	#print(tip_ages)
	
	
	if (is.na(tips_end_at_this_date))
		{
		tips_end_at_this_date = max(tip_ages)
		}
	
	nums_to_add_to_tip_to_ultrametricize = tips_end_at_this_date - tip_ages
	
	indices_of_branches_under_tips = get_indices_of_branches_under_tips(obj)

	obj$edge.length[indices_of_branches_under_tips] = obj$edge.length[indices_of_branches_under_tips] + nums_to_add_to_tip_to_ultrametricize
	
	return(obj)
	}




#######################################################
# average_tr_tips
#######################################################
#' Average the heights of (non-fossil) tips to make ultrametric-ish.
#'
#' When you have a digitized tree, or other slightly uneven source tree, average the tips 
#' to get them all to line up at 0 my before present.  This makes an ultrametric tree
#' if and only if there are no fossil tips in the tree.
#'
#' If the user includes fossils accidentally, this function can easily lead to pathological
#' results (negative branch lengths etc.), so use with care!!
#' 
#' @param tr An ape phylo object
#' @param fossils_older_than Tips that are older than \code{fossils_older_than} will be excluded from the tips that 
#' are going to be averaged. This is not currently set to 0, because Newick files can have slight precision issues etc.
#' that mean not all tips quite come to zero (which is why you need \code{\link{average_tr_tips}} in the first place!).  
#' Obviously you should be cautious about the value of , depending on the absolute timescale of your tree. Make sure you do
#' not inappropriately average in fossils!!
#' @return \code{edge_times_bp} A 2-column matrix with the age (from the present) of the top and bottom of each edge.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}, \code{\link{extend_tips_to_ultrametricize}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
average_tr_tips <- function(tr, fossils_older_than=0.6)
	{
	#require(BioGeoBEARS)	# for prt()

	# Check for negative branchlengths
	blren_equal_below_0_TF = tr$edge.length <= 0
	if (sum(blren_equal_below_0_TF) > 0)
		{
		tmptxt = paste(tr$edge.length[blren_equal_below_0_TF], collapse=", ", sep="")
		stoptxt = paste("\nFATAL ERROR in average_tr_tips(): the INPUT tree has branchlengths <= 0:\n", tmptxt, 
		"\nThis can sometimes happen if you (A) are accidentally including fossil tips (change 'fossils_older_than'), or\n",
		"if your input tree was e.g. an MCC (majority clade consensus) trees output by BEAST's TreeAnnotator.\n", 
		"In that case, you must fix the input Newick file. See ?check_BioGeoBEARS_run for comments.\n", sep="")
		cat(stoptxt)
		stop(stoptxt)
		}



	tr_table = prt(tr, printflag=FALSE)

	tipnums_tmp = 1:length(tr$tip.label)
	
	# Cut out fossils!!
	tips_are_fossils_TF = tr_table$time_bp[tipnums_tmp] > fossils_older_than
	tipnums_tmp = tipnums_tmp[tips_are_fossils_TF == FALSE]
	
	edges_w_tips_TF = tr$edge[,2] %in% tipnums_tmp
	edges_rownums = (1:nrow(tr$edge))
	edges_w_tips_rownums = edges_rownums[edges_w_tips_TF]
	
	tipnums = tr$edge[edges_w_tips_rownums, 2]
	

	#tipnums = tr$edge[edges_w_tips_rownums, 2]
	#tipnums
	meanval = mean(tr_table$node_ht[tipnums])
	meanval
	diffs_from_mean = tr_table$node_ht[tipnums] - meanval
	diffs_from_mean

	tr5 = tr
	tr5$edge.length[edges_w_tips_rownums] = tr5$edge.length[edges_w_tips_rownums] - diffs_from_mean
	
	min(tr5$edge.length)
	min(tr$edge.length)


	# Check the output; if IT has negative branchlengths, return NA!!
	# Check for negative branchlengths
	blren_equal_below_0_TF = tr5$edge.length <= 0
	if (sum(blren_equal_below_0_TF) > 0)
		{
		tmptxt = paste(tr$edge.length[blren_equal_below_0_TF], collapse=", ", sep="")
		stoptxt = paste("\nFATAL ERROR in average_tr_tips(): the OUTPUT tree has branchlengths <= 0:\n", tmptxt, 
		"\nThis can sometimes happen if you (A) are accidentally including fossil tips (change 'fossils_older_than'), or\n",
		"(B) if average_tr_tips() introduced more negative branches (especially can happen with shallow branches).\n", 
		"Returning nada\n", sep="")
		cat(stoptxt)
		stop(stoptxt)
		return(NA)
		}
	
	#tr5_table = prt(tr5)
	#tr5_table
		
	return(tr5)
	}






#######################################################
# is.not.na
#######################################################
#' Check for not NA
#'
#' A utility function. 
#' 
#' @param x Thing to check for NA
#' @return \code{TRUE} or \code{FALSE}
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
is.not.na <- function(x)
	{
	return(is.na(x) == FALSE)
	}






