#require("ape")
#require("rexpokit")
#require("cladoRcpp")


# Source ALL of the R files in a particular directory
# (so we don't have to have a fully-functioning package)

# Author: Nick Matzke, nmatzke
# License: GPL-3
#' Source all of the .R files in a directory of
#' the master of a GitHub package, except "compile" and "package" files.
#'
#' \code{sourceall_git} sources all of the .R files in the 
#' R/ directory of
#' a GitHub package. This can be useful for code that 
#' isn't developed into a package yet, or
#' isn't compiling yet, it can be 
#' useful to source all of the R files in a particular 
#' online directory.  
#' 
#' See \code{sourceall} for sourcing .R files in a
#' locally-saved directory.
#'
#' @param repo Should be username/repository_name, e.g. 
#'             "nmatzke/BioGeoBEARS" 
#'
#' @param subdir A desired subdirectory of the GitHub
#'               repository.  E.g. "R" or "R/" will work
#'               on an R package repository. Default "".
#' @param continue_recursion Should .R files in subdirectories of
#'               the specified directory be included? Default
#'               is FALSE.
#' @return Returns a vector of the .R files sourced.
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' \dontrun{
#' 
#'   repo = "nmatzke/BioGeoBEARS"
#'   subdir = "R"
#'   continue_recursion = FALSE
#'   sourceall_git(repo, subdir)
#' 
#' }
sourceall_git <- function(repo, subdir="", continue_recursion=FALSE)
  {
  junk='
  repo = "nmatzke/BioGeoBEARS"
  subdir = "R"
  continue_recursion = FALSE
  sourceall_git(repo, subdir)
  '
  
  
  #library(httr)     # for GET
  #library(devtools) # for source_url

  
  # Based on:
  # http://stackoverflow.com/questions/25485216/how-to-get-list-files-from-a-github-repository-folder-using-r
  # 
  git_api = paste0("https://api.github.com/repos/", repo, "/git/trees/master?recursive=2")

  req <- httr::GET(git_api)
  httr::stop_for_status(req)
  
  # Find the .R files
  filelist <- unlist(lapply(httr::content(req)$tree, "[", "path"), use.names = F)

  # Filter by subdirectory
  if (subdir != "")
  	{
	# Add a slash at the end of subdir, if needed
	string_to_startwith = addslash(subdir)

	# Remove first slash, if needed
	if (startsWith(string_to_startwith, prefix="/") == TRUE)
		{
		string_to_startwith = substr(x=string_to_startwith, start=2, stop=nchar(string_to_startwith))
		}
  	
  	# Go through filelist and subset
  	TF = startsWith(x=filelist, prefix=string_to_startwith)
  	filelist = filelist[TF]
  	}
  
  # Remove deeper subdirectories, if desired
  if (continue_recursion == FALSE)
  	{
  	number_of_slashes = stringr::str_count(string=filelist, pattern="/")
  	if (subdir == "")
  		{
	  	TF = number_of_slashes == 0
	  	} else {
	  	TF = number_of_slashes == 1
	  	}
		filelist = filelist[TF]
  	}
  
	# Keep only files that end with .R
	#Rfiles = grep("\\.R", filelist, value = TRUE, fixed = TRUE)

	# Files to remove
	Rfiles_remove_TF1 = grepl(pattern="compile", x=filelist)
	Rfiles_remove_TF2 = grepl(pattern="package", x=filelist)
	Rfiles_remove_TF3 = grepl(pattern="\\.Rproj", x=filelist)
	Rfiles_remove_TF = (Rfiles_remove_TF1 + Rfiles_remove_TF2 + Rfiles_remove_TF3) >= 1
	filelist = filelist[Rfiles_remove_TF==FALSE]

  
  # Remove .Rproj, add 
  #keepTF = grepl(pattern=".Rproj", x=Rfiles) == FALSE
  keepTF = endsWith(x=filelist, suffix=".R")
  filelist = filelist[keepTF]
  rawurl = paste0("https://raw.githubusercontent.com/", repo, "/master/")
  Rfiles = paste0(rawurl, filelist)
  cat("\nAttempting to source all *.R files in ", rawurl, ":\n", sep="")
  for (i in 1:length(Rfiles))
    {
    cat("Sourcing ", Rfiles[i], " ...\n", sep="")
    devtools::source_url(Rfiles[i])
    }
  cat("\n...done.")
  return(Rfiles)
  } # END sourceall_git


#######################################################
# sourceall
####################################ex###################
#' Source all .R files in a directory, except "compile" and "package" files
#' 
#' Utility function.
#' 
#' @param path The path to source
#' @param pattern Default is .R
#' @param ... Additional arguments to \code{\link{source}}.
#' @return \code{path} The path that was sourced.
#' @export
#' @seealso \code{\link[base]{source}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' test=1
#' 
#'
#' \dontrun{
#' 
#' sourceall("/GitHub/BioGeoBEARS/R/")
#' 
#' }
#' 
sourceall <- function(path=path, pattern="\\.R", ...)
	{
	tmppath = np(addslash(path))
	Rfiles = list.files(path=tmppath, pattern="\\.R", ...)
	
	# Keep only if it ends in .R
	TF = endsWith(x=Rfiles, suffix=".R")
	Rfiles = Rfiles[TF]
	
	# Files to remove
	Rfiles_remove_TF1 = grepl("compile", Rfiles)
	Rfiles_remove_TF2 = grepl("package", Rfiles)
	Rfiles_remove_TF3 = grepl("\\.Rproj", Rfiles)
	Rfiles_remove_TF = (Rfiles_remove_TF1 + Rfiles_remove_TF2 + Rfiles_remove_TF3) >= 1
	
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
	} # END sourceall



#' Shortcut for names()
n <- function(x)
	{
	return(names(x))
	}

#' Shortcut for names()
rn <- function(x)
	{
	return(names(x))
	}



#######################################################
# catdf
#######################################################
# Cat a data frame to screen, tab-delimited
# (handy for pasting into Excel)

#' Print a data.frame to screen with tab-delimination
#'
#' Uses \code{cat} to print a data.frame to screen with tab-
#' delimited columns. Useful for pasting into
#' e.g. Excel.
#'
#' @param dtf An input data.fram
#' @return NULL
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' 
#' # Set up a data.frame
#' student_names = c("Nick", "Tom", "Bob")
#' grade1 = c(37, 100, 60)
#' grade2 = c(43, 80, 70)
#' grade3 = c(100, 90, 100)
#' 
#' # convert to data frame
#' temp_table = cbind(student_names, grade1, grade2, grade3)
#' grade_dtf = as.data.frame(temp_table, stringsAsFactors=FALSE)
#' col_headers = c("names", "test1", "test2", "test3")
#' names(grade_dtf) = col_headers
#' 
#' # This might not paste well into Excel
#' grade_dtf
#' 
#' # But this will
#' catdf(grade_dtf)
#' 
catdf <- function(dtf)
	{
	cat(names(dtf), sep="\t")
	cat("\n")
	for (i in 1:nrow(dtf))
		{
		cat(unlist(dtf[i,]), sep="\t")
		cat("\n")
		} # END for (i in 1:nrow(dtf))
		
	} # END catdf


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
#' @param chunksize_toprint Number of lines to print (1 line is added, to 
#'        repeat between printing pages). Default 40.
#' @param printflag For optional printing. Passed to \code{\link{prflag}}.
#' @return NULL
#' @export
#' @seealso \code{\link[base]{print}}, \code{\link{prflag}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' # Set up a data.frame
#' student_names = c("Nick", "Tom", "Bob")
#' grade1 = c(37, 100, 60)
#' grade2 = c(43, 80, 70)
#' grade3 = c(100, 90, 100)
#' 
#' # convert to data frame
#' temp_table = cbind(student_names, grade1, grade2, grade3)
#' grade_dtf = as.data.frame(temp_table, stringsAsFactors=FALSE)
#' col_headers = c("names", "test1", "test2", "test3")
#' names(grade_dtf) = col_headers
#' 
#' # Print
#' x=printall(dtf=grade_dtf, chunksize_toprint=40, printflag=TRUE)
#' x=printall(dtf=grade_dtf, chunksize_toprint=12, printflag=TRUE)
#' 
#' # Don't print
#' x=printall(dtf=grade_dtf, chunksize_toprint=40, printflag=FALSE)
#' 
printall <- function(dtf, chunksize_toprint=40, printflag=TRUE)
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
#' @seealso \code{get_daughters}, \code{chainsaw2}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' test=1
#' prflag(x=test, printflag=TRUE)
#' 
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
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
np <- function(path=path, mustWork=FALSE, ...)
	{
	path = normalizePath(path, mustWork=mustWork, ...)
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' test=1
#' 
#' \dontrun{
#' fn = "filename"  # put the name of a text file here
#' moref(fn, printnotcat=FALSE)
#' moref(fn, printnotcat=TRUE)
#' 
#' }

#' 
moref <- function(fn, printnotcat=FALSE)
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
#' Utility function for \code{\link{%in% }}and \code{\link[base]{match}}, when one's brain gets confused.
#' 
#' @param list1 The list of things you want to check
#' @param list2 The list of things you want to check against
#' @return \code{matchlist} The TRUE/FALSE list for list1
#' @export
#' @seealso \code{\link[base]{match}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
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
#' This function takes \code{\link[base]{data.frame}} and applies \code{\link{unlist}}
#' to each column, returning a \code{\link[base]{data.frame}} that hopefully
#' can be used without list-related errors.
#' 
#' Because R is not a strongly typed language, odd things 
#' can happen when e.g. sticking together columns of data
#' to make a \code{data.frame} object (a data table). 
#' 
#' After many years of encountering this error at 
#' inconvenient times, I found it was easier to just
#' write a function to nuke the issue when it appeared. 
#' 
#' @param dtf Input \code{\link[base]{data.frame}}
#' @param printflag Print the results if TRUE (default: FALSE).
#' @return \code{dtf} The data.frame, hopefully without lists for columns
#' @export
#' @seealso \code{\link[base]{unlist}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' test=1
#' 
#' # Set up a data.frame
#' student_names = c("Nick", "Tom", "Bob")
#' grade1 = c(37, 100, 60)
#' grade2 = c(43, 80, 70)
#' grade3 = c(100, 90, 100)
#' 
#' # convert to data frame
#' temp_table = cbind(student_names, grade1, grade2, grade3)
#' grade_dtf = as.data.frame(temp_table, stringsAsFactors=FALSE)
#' col_headers = c("names", "test1", "test2", "test3")
#' names(grade_dtf) = col_headers
#' 
#' unlist_dtf_cols(grade_dtf, printflag=FALSE)
#' unlist_dtf_cols(grade_dtf, printflag=TRUE)
#' 
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
#' This function will return one 2nd-list index (the first match) for each item in the 1st list.
#' (In other words, the second-list index for each item in the first list.  Only the 
#' first hit in the second list is returned.
#' 
#' This is used by \code{\link{prt}}.
#'
#' @param list1 The first list. 
#' @param list2 The second list list.
#' @return \code{match_indices} The match indices.
#' @export
#' @seealso \code{\link{prt}}, \code{\link[base]{LETTERS}}, \code{\link{get_indices_where_list1_occurs_in_list2_noNA}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
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
#' This function will return one 2nd-list index (the first match) for each item in the 1st list.
#' (In other words, the second-list index for each item in the first list.  Only the 
#' first hit in the second list is returned.  Unlike 
#' \code{\link{get_indices_where_list1_occurs_in_list2}}, non-hits (NAs) are excluded.
#' 
#' This is used by \code{\link{get_indices_of_branches_under_tips}}, which is used by 
#' \code{\link{extend_tips_to_ultrametricize}}, which can be used by \code{\link{section_the_tree}}.
#'
#' @param list1 The first list. 
#' @param list2 The second list list.
#' @return \code{match_indices} The match indices.
#' @export
#' @seealso \code{\link{prt}}, \code{\link[base]{LETTERS}}, \code{\link{get_indices_where_list1_occurs_in_list2}}, 
#' \code{\link{extend_tips_to_ultrametricize}}, \code{\link{section_the_tree}}, \code{\link{return_items_not_NA}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' list1 = c("N", "I", NA, "C", "K")
#' return_items_not_NA(list1)
#' 
return_items_not_NA <- function(x)
	{
	y = x[!is.na(x)]
	return(y)
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' test=1
#' list1 = c("one", "two", "three")
#' list2str(list1)
#' list2 = c(1,2,3)
#' list2str(list2)

list2str <- function(list1, spacer=" ")
	{
	junk='
	list1 = c("one", "two", "three")
	list2str(list1)
	list2 = c(1,2,3)
	list2str(list2)
	'
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
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




# dfnums_to_numeric: convert each column to a non-list of type numeric, where possible
#######################################################
# dfn2n -- shortcut for dfnums_to_numeric
#######################################################
#' dfn2n/dfnums_to_numeric: convert each column to a non-list of type numeric, where possible
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
#' @param disallow_duplicate_colnames If TRUE (default), duplicate column names are disallowed. Duplicate column names cause \code{dfnums_to_numeric()} to fail to change the class of the second column with the same name.
#' @return \code{dtf} The output \code{\link[base]{data.frame}}.
#' @export
#' @seealso \code{\link{cls.df}}, \code{\link{unlist_df4}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' x = matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
#' cls.df(x)
#' dfn2n(adf(x))
#' unlist_df4(x)
#'
#' x = matrix(c(1,2,3,4,5,"A"), nrow=3, ncol=2)
#' cls.df(x)
#' dfn2n(adf(x))
#' unlist_df4(x)
#'
#' x = adf(matrix(c(1,2,3,4,5,"A"), nrow=3, ncol=2))
#' names(x) = c("A","B")
#' cls.df(x)
#' dfn2n(adf(x))
#' unlist_df4(x)
#' 
dfn2n <- function(dtf, max_NAs=0.5, printout=FALSE, roundval=NULL, disallow_duplicate_colnames=FALSE)
	{
	dfnums_to_numeric(dtf=dtf, max_NAs=max_NAs, printout=printout, roundval=roundval, disallow_duplicate_colnames=disallow_duplicate_colnames)
	}

# Force everything to numeric
df_to_numeric <- function(dtf, max_NAs=1.0, printout=FALSE, roundval=NULL, disallow_duplicate_colnames=FALSE)
	{
	dfnums_to_numeric(dtf=dtf, max_NAs=max_NAs, printout=printout, roundval=roundval, disallow_duplicate_colnames=disallow_duplicate_colnames)
	}


# dfnums_to_numeric: convert each column to a non-list of type numeric, where possible
#######################################################
# dfn2n -- shortcut for dfnums_to_numeric
#######################################################
#' dfn2n/dfnums_to_numeric: convert each column to a non-list of type numeric, where possible
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
#' @param disallow_duplicate_colnames If TRUE (default), duplicate column names are disallowed. Duplicate column names cause \code{dfnums_to_numeric()} to fail to change the class of the second column with the same name.
#' @return \code{dtf} The output \code{\link[base]{data.frame}}.
#' @export
#' @seealso \code{\link{cls.df}}, \code{\link{unlist_df4}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
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
dfnums_to_numeric <- function(dtf, max_NAs=0.5, printout=FALSE, roundval=NULL, disallow_duplicate_colnames=FALSE)
	{
	dtf_classes = cls.df(dtf, printout=FALSE)
	
	dtf_names = names(dtf)
	numcols = ncol(dtf)

	# 2014-11-18_NJM Error check: do NOT allow duplicated area names
	if (disallow_duplicate_colnames == TRUE)
		{
		uniq_dtf_names = unique(dtf_names)
		if (length(uniq_dtf_names) != length(dtf_names))
			{
			# Count area names
			countnames = rep(0, times=length(uniq_dtf_names))
			for (i in 1:length(uniq_dtf_names))
				{
				TF = uniq_dtf_names[i] == dtf_names
				countnames[i] = sum(TF)
				}
			TF_gt_1 = countnames > 1
		
			error_txt = paste("\n\nSTOP ERROR in dfnums_to_numeric(..., disallow_duplicate_colnames=TRUE): Duplicate column names in your data.frame.\n\nEach column name should be used only once for dfnums_to_numeric() to function properly. The counts of each name are printed below:\n\n", sep="")
			cat(error_txt)
		
			dup_names = adf2(cbind(uniq_dtf_names, countnames))
			names(dup_names) = c("column_name", "count")
			print(dup_names)
		
			cat("\n\n")
			error_txt2 = paste("STOP ERROR in dfnums_to_numeric(..., disallow_duplicate_colnames=TRUE): Duplicate column names in your data.frame. Each column name should be used only once for dfnums_to_numeric() to function properly.", sep="")
			stop(error_txt2)
			} # END if (length(uniq_dtf_names) != length(dtf_names))
		} # END if (disallow_duplicate_colnames == TRUE)
	

	
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
	
	#print(tmpstr)
	
	if (length(last_char) > 0)
		{
		if (last_char != "/")
			{
			outstr = paste(tmpstr, "/", sep="")
			} else {
			outstr = tmpstr
			} # END if (last_char != "/")
		} else {
		outstr = tmpstr
		} # END if (length(last_char) > 0)
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
#' NOTE: "numdigits_inbetween_have_fixed_digits"
#' is the closest to "round()" -- i.e., to round to 3 decimals, set "numdigits_inbetween_have_fixed_digits" 
#' to 3 + 1 (the dot) + typical number of integer digits.
#' NOTE: "numdigits_inbetween_have_fixed_digits"
#' is the closest to "round()" -- i.e., to round to 3 decimals, set "numdigits_inbetween_have_fixed_digits" 
#' to 3 + 1 (the dot) + typical number of integer digits.
#' @return \code{cellval} The value, reformatted and of class \code{\link[base]{character}}.
#' @export
#' @seealso \code{\link[base]{signif}}, \code{\link[base]{sprintf}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
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
	cellval = results_summary_table3$AIC_wt[1]
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
				tmp_num_digits = max(c(numdigits_inbetween_have_fixed_digits - 4, 0))
				fmt_str = paste('"%.', tmp_num_digits, 'f"', sep="")
				sprintf_cmd = paste("cellval = sprintf(fmt=", fmt_str, ", ", tmp_cellval, ")", sep="")
				eval(parse(text=sprintf_cmd))
				}

			if ( (abs(tmp_cellval) >= 100) && (abs(tmp_cellval) < 1000) )
				{
				tmp_num_digits = max(c(numdigits_inbetween_have_fixed_digits - 3, 1))
				fmt_str = paste('"%.', tmp_num_digits, 'f"', sep="")
				sprintf_cmd = paste("cellval = sprintf(fmt=", fmt_str, ", ", tmp_cellval, ")", sep="")
				eval(parse(text=sprintf_cmd))
				}

			if ( (abs(tmp_cellval) >= 0.1) && (abs(tmp_cellval) < 100) )
				{
				tmp_num_digits = max(c(numdigits_inbetween_have_fixed_digits - 2, 2))
				fmt_str = paste('"%.', tmp_num_digits, 'f"', sep="")
				sprintf_cmd = paste("cellval = sprintf(fmt=", fmt_str, ", ", tmp_cellval, ")", sep="")
				eval(parse(text=sprintf_cmd))
				}

			if ( (abs(tmp_cellval) >= 0.01) && (abs(tmp_cellval) < 0.1) )
				{
				tmp_num_digits = max(c(numdigits_inbetween_have_fixed_digits - 1, 3))
				fmt_str = paste('"%.', tmp_num_digits, 'f"', sep="")
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
#' NOTE: "numdigits_inbetween_have_fixed_digits"
#' is the closest to "round()" -- i.e., to round to 3 decimals, set "numdigits_inbetween_have_fixed_digits" 
#' to 3 + 1 (the dot) + typical number of integer digits.
#' @digits The "digits" argument overrides the "numdigits_inbetween_have_fixed_digits" argument, 
#' replacing it with digits+2 (e.g. digits=4 means numdigits_inbetween_have_fixed_digits=6
#' which means 0.12 will print 0.01234 instead. Default is NULL.
#' @return \code{output_table} The table, reformatted with cells of class \code{\link[base]{character}}.
#' @export
#' @seealso \code{\link[base]{signif}}, \code{\link[base]{sprintf}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' test=1
#' 
#' input_table = adf(c(143514514514532, -42.235235, -42.0000000, 
#' 0.0000, 0.0001, 0.00001, 0.0000111))
#' conditional_format_table(input_table=input_table)
#' 
conditional_format_table <- function(input_table, numbers_below_this_get_scientific=0.0001, numdigits_for_superlow_scientific=1, numbers_above_this_get_scientific=10000000, numdigits_for_superhigh_scientific=2, numdigits_inbetween_have_fixed_digits=4, digits=NULL)
	{
	defaults='
	input_table = adf(c(143514514514532, -42.235235, -42.0000000, 0.0000, 0.0001, 0.00001, 0.0000111))
	conditional_format_table(input_table=input_table)
	'
	
	# The "digits" argument overrides the "numdigits_inbetween_have_fixed_digits" argument, 
	# replacing it with digits+2 (e.g. digits=4 means numdigits_inbetween_have_fixed_digits=6
	# which means 0.12 will print 0.01234 instead.
	if (is.null(digits) == FALSE)
		{
		numdigits_inbetween_have_fixed_digits = digits + 2
		}
	
	# Error check for row names
	uniq_rownames = unique(rownames(input_table))
	if (length(uniq_rownames) < nrow(input_table))
		{
		numrows = nrow(input_table)
		nums = 1:numrows
		rownames(input_table) = paste0("row", nums)
		}
	
	# Unlist table and run through sapply
	chardata = sapply(unlist(input_table), FUN=conditional_format_cell,
		numbers_below_this_get_scientific=numbers_below_this_get_scientific, 
		numdigits_for_superlow_scientific=numdigits_for_superlow_scientific, 
		numbers_above_this_get_scientific=numbers_above_this_get_scientific, 
		numdigits_for_superhigh_scientific=numdigits_for_superhigh_scientific, 
		numdigits_inbetween_have_fixed_digits=numdigits_inbetween_have_fixed_digits)
	names(chardata) = NULL
	
	
	
	output_table = adf(matrix(data=chardata, nrow=nrow(input_table), ncol=ncol(input_table), byrow=FALSE))
	names(output_table) = names(input_table)
	rownames(output_table) = rownames(input_table)
	
	if (is.null(names(output_table)))
		{
		numcols = ncol(input_table)
		nums = 1:numcols
		names(output_table) = paste0("col", nums)
		}
	
	return(output_table)
	} # END conditional_format_table
	

#######################################################
# cft
#######################################################
#' Conditionally format the numbers (mostly) in a table
#' 
#' Shortcut for \code{\link{conditional_format_table}}.
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
#' NOTE: "numdigits_inbetween_have_fixed_digits"
#' is the closest to "round()" -- i.e., to round to 3 decimals, set "numdigits_inbetween_have_fixed_digits" 
#' to 3 + 1 (the dot) + typical number of integer digits.
#' @digits The "digits" argument overrides the "numdigits_inbetween_have_fixed_digits" argument, 
#' replacing it with digits+2 (e.g. digits=4 means numdigits_inbetween_have_fixed_digits=6
#' which means 0.12 will print 0.01234 instead. Default is digits=NULL.#' @return \code{output_table} The table, reformatted with cells of class \code{\link[base]{character}}.
#' @export
#' @seealso \code{\link[base]{signif}}, \code{\link[base]{sprintf}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' test=1
#' 
#' input_table = adf(c(143514514514532, -42.235235, -42.0000000, 
#' 0.0000, 0.0001, 0.00001, 0.0000111))
#' conditional_format_table(input_table=input_table)
#' 
cft <- function(input_table, numbers_below_this_get_scientific=0.0001, numdigits_for_superlow_scientific=1, numbers_above_this_get_scientific=10000000, numdigits_for_superhigh_scientific=2, numdigits_inbetween_have_fixed_digits=4, digits=NULL)
	{	
	return(conditional_format_table(input_table, 
		numbers_below_this_get_scientific=numbers_below_this_get_scientific,
		numdigits_for_superlow_scientific=numdigits_for_superlow_scientific,
		numbers_above_this_get_scientific=numbers_above_this_get_scientific,
		numdigits_for_superhigh_scientific=numdigits_for_superhigh_scientific,
		numdigits_inbetween_have_fixed_digits=numdigits_inbetween_have_fixed_digits, digits=digits))
	} # END cft










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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
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
#' This function assumes that the log-likelihoods are in the column "LnL", and the number of parameters is specified in "nparams".
#' 
#' @param restable A \code{\link[base]{data.frame}} with at least columns named "LnL" and "nparams".
#' @param row_to_use_as_null This is the row specifying the model to which the others will be compared in pairwise fashion.
#' @param rows_to_exclude Some rows may have models that the simpler model cannot nest within.  These should be excluded.
#' @param returnwhat If "pval", just return the p-value.  If "all", return all of the intermediate outputs.
#' @param add_to_table If TRUE, add to the main table and return the main table. If FALSE, return just the Akaike Weights results.
#' @return \code{pval} or \code{LRTrow}, both \code{\link[base]{data.frame}}.  Depends on \code{returnwhat}.
#' @export
#' @seealso \code{\link{lrttest}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
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
	
	names(wt_vBest) = paste(colname_to_use, "_wt", sep="")
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/AICc}
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
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
# exp( -0.5 * deltaAIC score for that model). The Akaike weight for a model is this value divided by the sum 
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' test=1
#' 
#' AICvals = c(40, 50, 60)
#' get_rownum_ref_model(AICvals, ref_model="best")
#' get_rownum_ref_model(AICvals, ref_model="worst")
#' 
get_rownum_ref_model <- function(AICvals, ref_model="best")
	{
	ex='
	AICvals = c(40, 50, 60)
	get_rownum_ref_model(AICvals, ref_model="best")
	get_rownum_ref_model(AICvals, ref_model="worst")
	'
	
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
#' @examples
#' test=1
#' 
#' AICval_1 = 20
#' AICval_2 = 30
#' get_AICweight_ratio_model1_over_model2(AICval_1, AICval_2)
#' 
get_AICweight_ratio_model1_over_model2 <- function(AICval_1, AICval_2)
	{
	if (AICval_1 < 0) {stop("ERROR: AIC cannot be negative (when the data are discrete), you probably put in log-likelihood (LnL)")}
	if (AICval_2 < 0) {stop("ERROR: AIC cannot be negative (when the data are discrete), you probably put in log-likelihood (LnL)")}

	AICvals = c(AICval_1, AICval_2)
	
	
	# relative weight of model2 (lower # of paramters
	AIC_weight_model2 = getAIC_weight_for_model1(AICval_2, AICvals)

	# relative weight of model1 (higher # of paramters
	AIC_weight_model1 = getAIC_weight_for_model1(AICval_1, AICvals)
	
	AICweight_ratio_model1 = AIC_weight_model1 / AIC_weight_model2

	return(AICweight_ratio_model1)
	}



# Open a directory
#######################################################
# opd
#######################################################
#' Open a directory
#' 
#' This function opens an e.g. Finder window. This is merely a shortcut for 
#' e.g. \code{system(paste0("open ", getwd()))}.
#'
#' This is meant to be a within-R equivalent to the Terminal command "\code{open .}".
#' 
#' @param wd The directory to open. Defaults to \code{getwd()}.
#' @return NULL
#' @export
#' @seealso \code{\link{system}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' # Use system() to open the current working directory.
#' system(paste0("open ", getwd()))
#' 
#' # Use opd() to open the current working directory.
#' opd()
#' 
opd <- function(wd=getwd())
	{
	cat("\nopd() or openwd() is using the system('open ", wd, "') command to open a Finder window in Mac. May or may not not work on Windows.\n", sep="")
	wd = getwd()
	system(paste0("open ", wd))
	return(NULL)
	} # END opd <- function(wd=getwd()



#######################################################
# openwd
#######################################################
#' Open a directory
#' 
#' This function opens an e.g. Finder window. This is merely a shortcut for 
#' e.g. \code{system(paste0("open ", getwd()))}.
#'
#' This is meant to be a within-R equivalent to the Terminal command "\code{open .}".
#' 
#' @param wd The directory to open. Defaults to \code{getwd()}.
#' @return NULL
#' @export
#' @seealso \code{\link{system}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' # Use system() to open the current working directory.
#' system(paste0("open ", getwd()))
#' 
#' # Use opd() to open the current working directory.
#' opd()
#' 
openwd <- function(wd=getwd())
	{
	cat("\nopd() or openwd() is using the system('open ", wd, "') command to open a Finder window in Mac. May or may not not work on Windows.\n", sep="")
	wd = getwd()
	system(paste0("open ", wd))
	return(NULL)
	} # END opd <- function(wd=getwd()























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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' test=1
#' 
#' # Load hard-coded tree
#' tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
#' get_daughters(nodenum=4, t=tr)
#' get_daughters(nodenum=5, t=tr)
#' 
get_daughters <- function(nodenum, t)
	{
	ex='
	tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
	get_daughters(nodenum=4, t=tr)
	get_daughters(nodenum=5, t=tr)
	'
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
#' Get the indices of \code{what}, in list \code{inlist}.
#' 
#' @param what The item to find
#' @param inlist The list to search in 
#' @return \code{matching_indices} List of the matching indices
#' @export
#' @seealso \code{\link{get_daughters}}, \code{\link{chainsaw2}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' test=1
#' 
#' list_items = c("a", "b", "c", "b")
#' findall(what="b", inlist=list_items)
#' 
findall <- function(what, inlist)
	{
	ex='
	list_items = c("a", "b", "c", "b")
	findall(what="b", inlist=list_items)	
	'
	
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' test=1
#' tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
#' get_parent(nodenum=1, t=tr)
#' get_parent(nodenum=2, t=tr)
#' get_parent(nodenum=3, t=tr)
#' get_parent(nodenum=4, t=tr)
#' get_parent(nodenum=5, t=tr)
#' 
get_parent <- function(nodenum, t)
	{
	ex='
	tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
	get_parent(nodenum=1, t=tr)
	get_parent(nodenum=2, t=tr)
	get_parent(nodenum=3, t=tr)
	get_parent(nodenum=4, t=tr)
	get_parent(nodenum=5, t=tr)
	'
	matching_edges = findall(nodenum, t$edge[,2])
	parent_nodenum = t$edge[,1][matching_edges][1]
	return(parent_nodenum)
	}


#######################################################
# get_level
#######################################################
#' Get a node's level in the tree
#'
#' Finds how many nodes above the root a node is.
#' 
#' @param nodenum The node number to get the parent of
#' @param t An ape phylo object
#' @param tmplevel A starting level (the function is recursive)
#' @return \code{tmplevel} The level of the node.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' test=1
#' tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
#' get_level(nodenum=1, t=tr)
#' get_level(nodenum=2, t=tr)
#' get_level(nodenum=3, t=tr)
#' get_level(nodenum=4, t=tr)
#' get_level(nodenum=5, t=tr)
#' 
get_level <- function(nodenum, t, tmplevel=0)
	{
	ex='
	tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
	get_level(nodenum=1, t=tr)
	get_level(nodenum=2, t=tr)
	get_level(nodenum=3, t=tr)
	get_level(nodenum=4, t=tr)
	get_level(nodenum=5, t=tr)
	'
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' test=1
#' tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
#' get_TF_tips(obj=tr)
#' 
get_TF_tips <- function(obj)
	{
	ex='
	tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
	get_TF_tips(obj=tr)
	'
	
	# Get TF for nodes being tips
	
	# BIG CHANGE?
	#TF_tips = match_list1_in_list2(1:length(dists_from_root), obj$tip.label)
	ntips = length(obj$tip.label)
	numnodes = obj$Nnode
	total_numnodes = ntips+numnodes
	
	# 2019-01-24_change
	#TF_tips = match_list1_in_list2(list1=1:length(obj$edge), list2=1:length(obj$tip.label))
	TF_tips = match_list1_in_list2(list1=1:total_numnodes, list2=1:length(obj$tip.label))
	total_numnodes
	#TF_tips = obj$tip.label[TF_tips_indices]
	return(TF_tips)
	}



#######################################################
# get_indices_of_tip_nodes
#######################################################
#' Get TRUE/FALSE for nodes being tips
#'
#' A utility function that returns indices (node numbers) of the tips. This mostly saves typing.
#' 
#' @param obj An ape phylo object
#' @return \code{tip_indices} The node numbers of the tips.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}, \code{\link[ape]{phylo}}, \code{\link{get_indices_of_branches_under_tips}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' test=1
#' tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
#' get_indices_of_tip_nodes(obj=tr)
#' 
get_indices_of_tip_nodes <- function(obj)
	{
	ex='
	tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
	get_indices_of_tip_nodes(obj=tr)
	'

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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' test=1
#' tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
#' get_indices_of_branches_under_tips(obj=tr)
#' 
get_indices_of_branches_under_tips <- function(obj)
	{
	ex='
	tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
	get_indices_of_branches_under_tips(obj=tr)
	'

	tip_indices = get_indices_of_tip_nodes(obj)
	branchnums_under_tips = get_indices_where_list1_occurs_in_list2_noNA(list1=tip_indices, list2=obj$edge[, 2])
	return(branchnums_under_tips)
	}




#######################################################
# get_node_ages_of_tips
#######################################################
#' Get the ages of each tip above the root
#'
#' A utility function to get the ages of each tip above the root.
#' Uses \code{dist.nodes}, which may get slow with large trees.
#' 
#' @param obj An ape phylo object
#' @return \code{TF_tips} The age (from the root) of each tip.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' test=1
#' tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
#' get_node_ages_of_tips(obj=tr)
#' 
get_node_ages_of_tips <- function(obj)
	{
	ex='
	tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
	get_node_ages_of_tips(obj=tr)
	'
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
#' A utility function. \code{get_all_node_ages} uses of \code{\link[ape]{dist.nodes}}
#' internally, which may be slow for large trees.
#' 
#' @param obj An ape phylo object
#' @return \code{TF_tips} The age (from the root) of each node.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' test=1
#' tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
#' get_all_node_ages(obj=tr)
#' 
get_all_node_ages <- function(obj)
	{
	ex='
	tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
	get_all_node_ages(obj=tr)
	'
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
#' Same function as \code{\link{get_root_age}}.
#' 
#' @param obj An ape phylo object
#' @return \code{max_height} The age (from the root) of the highest node.
#' @export
#' @seealso \code{\link{get_root_age}}, \code{\link{prt}}, \code{\link{chainsaw2}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' test=1
#' tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
#' get_all_node_ages(obj=tr)
#' max(get_all_node_ages(obj=tr))
#' get_max_height_tree(obj=tr)
#' 
get_max_height_tree <- function(obj)
	{
	ex='
	tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
	get_max_height_tree(obj=tr)
	'
	max_height = max(get_node_ages_of_tips(obj))
	return(max_height)
	}

get_tree_height <- function(obj)
	{
	get_max_height_tree(obj)
	}

get_tr_height <- function(obj)
	{
	get_max_height_tree(obj)
	}

get_treeheight <- function(obj)
	{
	get_max_height_tree(obj)
	}

get_trht <- function(obj)
	{
	get_max_height_tree(obj)
	}


#######################################################
# get_root_age
#######################################################
#' Get the age of the root node.
#'
#' I.e., the distance of the highest node above the root.  A utility function. 
#' Use of \code{\link[ape]{dist.nodes}} may be slow.
#'
#' Same function as \code{\link{get_max_height_tree}}.
#' 
#' @param obj An ape phylo object
#' @return \code{root_age} The age of the root
#' @export
#' @seealso \code{\link{get_max_height_tree}}, \code{\link{prt}}, \code{\link{chainsaw2}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' test=1
#' tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
#' get_all_node_ages(obj=tr)
#' max(get_all_node_ages(obj=tr))
#' get_root_age(obj=tr)
#'  
get_root_age <- function(tr)
	{
	tip_ages_above_root = get_node_ages_of_tips(tr)
	root_age = max(tip_ages_above_root)
	return(root_age)
	}
	



#######################################################
# get_edge_times_before_present
#######################################################
#' Get the times of the top and bottom of each edge
#'
#' A utility function that gets the times of the top and 
#' bottom of each edge of the input tree.
#' 
#' @param t An ape phylo object
#' @return \code{edge_times_bp} A 2-column matrix with the age (from the present) of the top and bottom of each edge.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' test=1
#' tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
#' get_edge_times_before_present(t=tr)
#' 
get_edge_times_before_present <- function(t)
	{
	ex='
	tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
	get_edge_times_before_present(t=tr)
	'

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
# get_nodenums
#######################################################
#' Get the unique node numbers in a tree
#' 
#' This is a utility function for \code{\link{get_nodenum_structural_root}}.
#' It returns the the NUMBERS identifying each node.
#'
#' @param t A tree object in \code{\link[ape]{phylo}} format.
#' @return \code{ordered_nodenames} The node numbers, in order.
#' @export
#' @seealso \code{\link[ape]{phylo}}, \code{\link{get_nodenum_structural_root}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' test = 1
#' tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
#' get_nodenums(t=tr)
#'  
get_nodenums <- function(t)
	{
	ex='
	tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
	get_nodenums(t=tr)
	'
	# get just the unique node numbers from the edge list (left column: start node; right column: end node):
	nodenames = unique(c(t$edge))
	ordered_nodenames = nodenames[order(nodenames)]
	return(ordered_nodenames)
	}

#######################################################
# get_nodenum_structural_root
#######################################################
#' Gets the root node 
#' 
#' This function gets the root node by finding the node not 
#' in the descendants list (\code{edge[,2]}). This
#' may be more reliable than e.g. assuming \code{length(tr$tip.label)+1}.
#'
#' @param t A tree object in \code{\link[ape]{phylo}} format.
#' @param print_nodenum Print the node numbers as you go through the list? Default FALSE.
#' @return \code{root_nodenums_list} 
#' @export
#' @seealso \code{\link[ape]{phylo}}, \code{\link{get_nodenums}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' blah=1
#' 
#' tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
#' x = get_nodenum_structural_root(t=tr, print_nodenum=FALSE)
#' x
#' x = get_nodenum_structural_root(t=tr, print_nodenum=TRUE)
#' x
get_nodenum_structural_root <- function(t, print_nodenum=FALSE)
	{
	ex='
	tr = read.tree(file="", text="((human:1,chimp:1):1,gorilla:2);")
	x = get_nodenum_structural_root(t=tr, print_nodenum=FALSE)
	x
	x = get_nodenum_structural_root(t=tr, print_nodenum=TRUE)
	x
	'

	#numnodes = length(t$tip.label) + length(t$node.label)
	#ordered_nodes = 1:length(numnodes)
	
	ordered_nodes = get_nodenums(t)

	root_nodenums_list = c()
	for (n in 1:length(ordered_nodes))
		{
		tmpnode = ordered_nodes[n]
		if (tmpnode %in% t$edge[,2])
			{
			blah = TRUE
			}
		else
			{
			if (print_nodenum == TRUE)
				{
				cat("get_nodenum_structural_root(): Root nodenum = ", tmpnode, sep="")
				}
			root_nodenums_list = c(root_nodenums_list, tmpnode)
			}
		}
	return(root_nodenums_list)
	}






# Fix a tree where tips don't come to zero
#######################################################
# level_tree_tips
#######################################################
#' Fix a tree where tips don't all have age of precisely 0.0
#'
#' Due to rounding errors or other issues, sometimes newick 
#' strings encode phylogenies that are supposed to be ultrametric,
#' but where all of the tips do not have precisely the same 
#' height above the root. This function attempts to correct this by
#' adjusting all of the tip branches so that all tips have the same 
#' height. The methods are described in the method argument.
#'
#' NOTE: Use with caution!  If the tree is a non-dated tree, or
#' is a dated tree with fossils, pathological behavior could
#' easily result from applying this function (e.g., negative 
#' branchlengths, or a radical alteration of the tree). As with 
#' everything, you need to think before applying an particular
#' function to a tree. This function is intended for fixing 
#' very small discrepancies in tip height.
#'
#' That said, sometimes you want to equalize the tips that are 
#' supposed to have age 0.0 Ma (mega-annum, million years before 
#' present), but 
#' leave the fossils unchanged. Any tips older than 
#' \code{fossils_older_than} will be treated as fossils, and 
#' left unchanged. (By default, tips older than 0.6 Ma will
#' be considered fossils.)
#'
#' @param tr An ape \code{phylo} object.
#' @param method The option "mean" (default) takes the average 
#'               of the tip heights and adjusts all tips to 
#'               match. The option "highest" brings all 
#'               tips up to the heighest tip. The option 
#'               "lowest" decreases all tips to the lowest tipl
#' @param printflag \code{TRUE} or \code{FALSE}, passed to \code{\link{prt}}.
#'                  Default is \code{FALSE}.
#' @param fossils_older_than Tips older than this are left out of the 
#'                adjustment of tip ages. Default 0.6.
#' @return tr, the adjusted tree with leveled tips.
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
#' tr = read.tree(file="", text="((human:0.9,chimp:1):1,gorilla:2);")
#' method="mean"
#' printflag=TRUE
#' fossils_older_than=0.6
#' tr2 = level_tree_tips(tr, method="mean", printflag=FALSE, fossils_older_than=0.6)
#' write.tree(tr2, file="")
#' tr2 = level_tree_tips(tr, method="highest", printflag=FALSE, fossils_older_than=0.6)
#' write.tree(tr2, file="")
#' tr2 = level_tree_tips(tr, method="lowest", printflag=FALSE, fossils_older_than=0.6)
#' write.tree(tr2, file="")
#' 
level_tree_tips <- function(tr, method="mean", printflag=FALSE, fossils_older_than=0.6)
	{
	ex='
	tr = read.tree(file="", text="((human:0.9,chimp:1):1,gorilla:2);")
	method="mean"
	printflag=TRUE
	fossils_older_than=0.6
	tr2 = level_tree_tips(tr, method="mean", printflag=FALSE, fossils_older_than=0.6)
	write.tree(tr2, file="")
	tr2 = level_tree_tips(tr, method="highest", printflag=FALSE, fossils_older_than=0.6)
	write.tree(tr2, file="")
	tr2 = level_tree_tips(tr, method="lowest", printflag=FALSE, fossils_older_than=0.6)
	write.tree(tr2, file="")
	'


	# Look at the tree table:
	trtable = prt(tr, printflag=printflag, fossils_older_than=fossils_older_than)
	trtable$time_bp

	ntips = length(tr$tip.label)
	fossils_TF = trtable$fossils[1:ntips]

	tipnums_to_change = (1:ntips)[fossils_TF == FALSE]
	tip_ages = trtable$time_bp[tipnums_to_change]
	
	# Get the mean of the tipages, add difference
	if (method == "mean")
		{
		mean_tipage = mean(tip_ages)
		}
	# Adjust everything to match the highest tip
	if (method == "highest")
		{
		mean_tipage = 0
		}
	# Adjust everything to match the lowest tip
	# (could introduce errors)
	if (method == "lowest")
		{
		mean_tipage = max(tip_ages)
		}


	change_branchlength_by = tip_ages - mean_tipage

	# Edit the tree:
	edgenums = trtable$parent_br[tipnums_to_change]
	tr$edge.length[edgenums] = tr$edge.length[edgenums] + change_branchlength_by
	
	# Check for negative branchlengths
	brlen_negative_count = sum( tr$edge.length[edgenums] < 0)
	
	if (brlen_negative_count > 0)
		{
		error_txt = paste0("STOP ERROR in level_tree_tips(): correcting uneven tip ages with method='", method, "', and fossils_older_than=", fossils_older_than, " has resulted in a tree with ", brlen_negative_count, " branches with negative branchlengths. This is Very Bad. You should correct your tree with another method or by hand.  See e.g. impose_min_brlen().")
		cat("\n\n")
		cat(error_txt)
		cat("\n\n")
		stop(error_txt)
		} # END if (brlen_negative_count > 0)
	
	return(tr)
	}


# print tree in hierarchical format
#######################################################
# prt
#######################################################
#' Print tree in table format
#' 
#' Learning and using APE's tree structure can be difficult and confusing because much of the information is
#' implicit.  This function prints the entire
#' tree to a table, and makes much of the implicit information explicit.  It is not particularly fast, but
#' it is useful, especially for learning about ape's tree structure, and/or for easily accessing various
#' information from the tree.
#'
#' See \url{http://ape.mpl.ird.fr/ape_development.html} for the official documentation of R tree objects.
#' 
#' @param t A \code{\link[ape]{phylo}} tree object.
#' @param printflag Should the table be printed to screen?  Default TRUE.
#' @param relabel_nodes Manually re-number the internal nodes, if desired. Default FALSE.
#' @param time_bp_digits The number of digits to print in the time_bp (time before present) column. Default=7.
#' @param add_root_edge Should a root edge be added?  Default \code{TRUE}.
#' @param get_tipnames Should the list of tipnames descending from each node be printed as a string in another column?  
#' This is slow-ish, but useful for matching up nodes between differing trees. Default \code{FALSE}.
#' @param fossils_older_than Tips that are older than \code{fossils_older_than} will be marked as \code{TRUE} in a column called \code{fossil}. Default 0.001.
#' @param silence_warnings Suppress warnings about missing branchlengths (prt makes each branchlength equal 1)
#' This is not currently set to 0, because Newick files can have slight precision issues etc. that mean not all tips quite come to zero.  You 
#' can attempt to fix this with \code{\link{average_tr_tips}} (but make sure you do not inappropriately average in fossils!!).
#' @return \code{dtf} A \code{\link[base]{data.frame}} holding the table. (Similar to the printout of a \code{\link[phylobase]{phylo4}} object.)
#' @export
#' @seealso \code{\link[ape]{phylo}}, \code{\link{average_tr_tips}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://ape.mpl.ird.fr/ape_development.html}
#' @examples
#' test=1
#' 	tr = read.tree(file="", text="((human:0.9,chimp:1):1,gorilla:2);")
#' 	trtable1 = prt(tr)
#' 	trtable2 = prt(t=tr, printflag=FALSE, relabel_nodes=TRUE, get_tipnames=TRUE, fossils_older_than=0.000001)
#' 
prt <- function(t, printflag=FALSE, relabel_nodes=FALSE, time_bp_digits=7, add_root_edge=TRUE, get_tipnames=TRUE, fossils_older_than=0.001, silence_warnings=FALSE)
	{
	defaults='
	#wd = "/drives/GDrive/__classes/BIOSCI210/lab3_genome_size/"
	#setwd(wd)

	#library(ape)

	# Read a Newick-formatted phylogeny file (which has been subset to birds and mammals 
	# found in the table) to an APE tree object 
	#trfn = "birds_mammals_subset_tree.newick"
	#tr = read.tree(trfn)
	
	tr = read.tree(file="", text="((human:0.9,chimp:1):1,gorilla:2);")
	trtable1 = prt(tr)
	trtable2 = prt(t=tr, printflag=FALSE, relabel_nodes=TRUE, get_tipnames=TRUE, fossils_older_than=0.000001)

	t = tr
	printflag=TRUE;
	relabel_nodes = FALSE;
	time_bp_digits=7;
	add_root_edge=TRUE;
	get_tipnames=FALSE;
	fossils_older_than=0.6
	silence_warnings=FALSE
	'
	
	if (class(t) != "phylo")
		{
		txt = paste0("STOP ERROR IN prt(): the input 't' must be of class 'phylo'. You had class(t)='", class(t), "'.")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		} # END if (class(t) != "phylo")
	

	# assemble beginning table
	
	# check if internal node labels exist
	if (("node.label" %in% attributes(t)$names == FALSE) || (is.null(t$node.label) == TRUE))
		{
		rootnum = get_nodenum_structural_root(t)
		
		new_node_labels = paste("inNode", rootnum:(rootnum+t$Nnode-1), sep="")
		t$node.label = new_node_labels
		}
	
	# or manually relabel the internal nodes, if desired
	if (relabel_nodes == TRUE)
		{
		rootnum = get_nodenum_structural_root(t)
		
		new_node_labels = paste("inNode", rootnum:(rootnum+t$Nnode-1), sep="")
		t$node.label = new_node_labels
		}
	
	labels = c(t$tip.label, t$node.label)
	ordered_nodenames = get_nodenums(t)
	#nodenums = 1:length(labels)
	node.types1 = rep("tip", length(t$tip.label))
	node.types2 = rep("internal", length(t$node.label))
	node.types2[1] = "root"
	node.types = c(node.types1, node.types2)
	

	# These are the index numbers of the edges below each node
	parent_branches = get_indices_where_list1_occurs_in_list2(ordered_nodenames, t$edge[,2])


	# Error check for branchlengths
	edgelength_NULL_warning = FALSE
	if (is.null(t$edge.length))
		{
		edgelength_NULL_warning = TRUE
		txt = paste0("\n\nWARNING: 'brlen_to_parent = t$edge.length[parent_branches]'...produced NULL as a result. This is probably a parsimony tree/cladogram with no branchlengths.\nprt() is inserting '1' for each branchlength. You may or may not want this. THESE ARE NOT REAL BRANCHLENGTHS AND NEITHER THE BRANCHLENGTHS NOR THE DERIVED 'TIMES' ETC SHOULD BE USED!!! You have been warned!\n\n")
		if (silence_warnings == FALSE)
			{
			cat(txt)
			}
		
		t$edge.length =rep(1, times=length(parent_branches))
		} # END if (is.null(brlen_to_parent))


	#parent_edges = parent_branches
	brlen_to_parent = t$edge.length[parent_branches]
	
	
	
	parent_nodes = t$edge[,1][parent_branches]
	daughter_nodes = lapply(ordered_nodenames, get_daughters, t)
	
	# print out the structural root, if desired
	root_nodenum = get_nodenum_structural_root(t)
	tmpstr = paste("prt(t): root_nodenum=", root_nodenum, "\n", sep="")
	prflag(tmpstr, printflag=printflag)
	
	levels_for_nodes = unlist(lapply(ordered_nodenames, get_level, t))
	#tmplevel = get_level(23, t)
	#print(tmplevel)
	
	
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
	
	
	# If desired, get the list of all tipnames descended from a node, in alphabetical order
	if (get_tipnames == TRUE)
		{
		# Make the empty list
		list_of_clade_members_lists = rep(list(NA), length(ordered_nodenames))
		
		# Tips have only one descendant
		list_of_clade_members_lists[1:length(t$tip.label)] = t$tip.label
		list_of_clade_members_lists
		
		
		nontip_nodenums = (length(t$tip.label)+1) : length(ordered_nodenames)
		if (length(nontip_nodenums) > 1)
			{
			# More than 1 node
			nontip_nodenames = ordered_nodenames[nontip_nodenums]
			nontip_cladelists = sapply(X=nontip_nodenames, FUN=get_all_daughter_tips_of_a_node, t=t)
			nontip_cladelists
			
			nontip_cladelists_alphabetical = sapply(X=nontip_cladelists, FUN=sort)
			nontip_cladelists_alphabetical
			
			nontip_cladelists_alphabetical_str = sapply(X=nontip_cladelists_alphabetical, FUN=paste, collapse=",")
			nontip_cladelists_alphabetical_str
			
			# Store the results
			list_of_clade_members_lists[nontip_nodenums] = nontip_cladelists_alphabetical_str
			list_of_clade_members_lists
			} else {
			# Just one node
			nontip_nodenames = ordered_nodenames[nontip_nodenums]
			nontip_cladelists = sapply(X=nontip_nodenames, FUN=get_all_daughter_tips_of_a_node, t=t)
			nontip_cladewords = unlist(sapply(X=nontip_cladelists, FUN=strsplit, split=","))
			
			nontip_cladelists_alphabetical = sort(nontip_cladewords)
			nontip_cladelists_alphabetical
			
			nontip_cladelists_alphabetical_str = paste(nontip_cladelists_alphabetical, collapse=",", sep="")
			nontip_cladelists_alphabetical_str
			
			# Store the results
			list_of_clade_members_lists[nontip_nodenums] = nontip_cladelists_alphabetical_str
			list_of_clade_members_lists			
			}
			
		}

	
	# Add fossils TRUE/FALSE column.  You can turn this off with fossils_older_than=NULL.
	fossils = times_before_present > fossils_older_than

	# Obviously, internal nodes are irrelevant and should be NA
	tmpnodenums = (length(t$tip.label)+1) : ( length(t$tip.label) + t$Nnode )
	fossils[tmpnodenums] = NA
	
	if (get_tipnames == FALSE)
		{
		# Don't put in the list of clade names
		tmpdtf = cbind(1:length(ordered_nodenames), ordered_nodenames, levels_for_nodes, node.types, parent_branches, brlen_to_parent, parent_nodes, daughter_nodes, h, times_before_present, fossils, labels)
		
		dtf = as.data.frame(tmpdtf, row.names=NULL)
		# nd = node
		
		# edge.length is the same as brlen_2_parent
		names(dtf) = c("node", "ord_ndname", "node_lvl", "node.type", "parent_br", "edge.length", "ancestor", "daughter_nds", "node_ht", "time_bp", "fossils", "label")
		
		# convert the cols from class "list" to some natural class
		dtf = unlist_dtf_cols(dtf, printflag=FALSE)
		} else {
		# Put in the list of clade names
		tmpdtf = cbind(1:length(ordered_nodenames), ordered_nodenames, levels_for_nodes, node.types, parent_branches, brlen_to_parent, parent_nodes, daughter_nodes, h, round(times_before_present, digits=time_bp_digits), fossils, labels, list_of_clade_members_lists)
		
		dtf = as.data.frame(tmpdtf, row.names=NULL)
		# nd = node
		
		# edge.length is the same as brlen_2_parent
		#print(names(names(dtf)))
		#print(c("node", "ord_ndname", "node_lvl", "node.type", "parent_br", "edge.length", "ancestor", "daughter_nds", "node_ht", "time_bp", "fossils", "label", "tipnames"))
		names(dtf) = c("node", "ord_ndname", "node_lvl", "node.type", "parent_br", "edge.length", "ancestor", "daughter_nds", "node_ht", "time_bp", "fossils", "label", "tipnames")
		
		# convert the cols from class "list" to some natural class
		dtf = unlist_dtf_cols(dtf, printflag=FALSE)		
		}
	

	
	
	
	
	# Add the root edge, if desired
	# (AND, only if t$root.edge exists)
	if ( (add_root_edge == TRUE) && (!is.null(t$root.edge)) )
		{
		root_row_TF = dtf$node.type == "root"
		root_edge_length = t$root.edge
		
		# Stick in this edge length
		dtf$edge.length[root_row_TF] = root_edge_length
		
		# Add the root edge length to all node heights
		dtf$node_ht = dtf$node_ht + root_edge_length
		}
	
	# print if desired
	prflag(dtf, printflag=printflag)

	if ((edgelength_NULL_warning == TRUE) && (printflag == TRUE))
		{
		txt = paste0("\n\nWARNING: 'brlen_to_parent = t$edge.length[parent_branches]'...produced NULL as a result. This is probably a parsimony tree/cladogram with no branchlengths.\nprt() is inserting '1' for each branchlength. You may or may not want this. THESE ARE NOT REAL BRANCHLENGTHS AND NEITHER THE BRANCHLENGTHS NOR THE DERIVED 'TIMES' ETC SHOULD BE USED!!! You have been warned!\n\n")
		if (silence_warnings == FALSE)
			{
			cat(txt)
			}
		} # END if (edgelength_NULL_warning == TRUE)
	
	#tree_strings = c()
	#root_str = get_node_info(root_nodenum, t)
	return(dtf)
	} # END prt()











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
#' @param tips_end_at_this_date The tips can be set to something other than 0, if desired.
#' (This could produce negative branchlengths, however.)
#' @return \code{obj} The corrected phylogeny
#' @export
#' @seealso \code{\link[ape]{read.tree}}, \code{\link{prt}}, \code{\link{average_tr_tips}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' test=1
#' tr = read.tree(file="", text="((human:0.9,chimp:1):1,gorilla:2);")
#' tr2 = extend_tips_to_ultrametricize(obj=tr, age_of_root=0, tips_end_at_this_date=NA, fossils_older_than=NULL)
#' write.tree(tr2, file="")
#' 
#' tr2 = extend_tips_to_ultrametricize(obj=tr, age_of_root=0, tips_end_at_this_date=NA, fossils_older_than=0.05)
#' write.tree(tr2, file="")
#' 
extend_tips_to_ultrametricize <- function(obj, age_of_root=0, tips_end_at_this_date=NA, fossils_older_than=NULL)
	{
	ex='
	tr = read.tree(file="", text="((human:0.9,chimp:1):1,gorilla:2);")
	tr2 = extend_tips_to_ultrametricize(obj=tr, age_of_root=0, tips_end_at_this_date=NA, fossils_older_than=NULL)
	write.tree(tr2, file="")

	tr2 = extend_tips_to_ultrametricize(obj=tr, age_of_root=0, tips_end_at_this_date=NA, fossils_older_than=0.05)
	write.tree(tr2, file="")
	'

	#print("node ages of tips:")
	tip_ages = age_of_root + get_node_ages_of_tips(obj)
	#print(tip_ages)
	
	# Change ONLY things YOUNGER than fossils
	if (is.null(fossils_older_than) == FALSE)
		{
		tips_must_be_higher_than_this = age_of_root + max(tip_ages) - fossils_older_than
		change_TF = tip_ages > tips_must_be_higher_than_this
		} else {
		change_TF = rep(TRUE, length(tip_ages))
		}
	
	
	if (is.na(tips_end_at_this_date))
		{
		tips_end_at_this_date = max(tip_ages)
		}
	
	nums_to_add_to_tip_to_ultrametricize = tips_end_at_this_date - tip_ages[change_TF]
	
	indices_of_branches_under_tips = get_indices_of_branches_under_tips(obj)

	obj$edge.length[indices_of_branches_under_tips][change_TF] = obj$edge.length[indices_of_branches_under_tips][change_TF] + nums_to_add_to_tip_to_ultrametricize
	
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
#' results (negative branch lengths etc.), so use with care!! (See e.g.
#' impose_min_brlen() for a method to remove negative branchlengths)
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' test=1
#' tr = read.tree(file="", text="((human:0.9,chimp:1):1,gorilla:2);")
#' tr2 = average_tr_tips(tr, fossils_older_than=0.6)
#' write.tree(tr2, file="")
#' tr2 = average_tr_tips(tr, fossils_older_than=0.05)
#' write.tree(tr2, file="")
#' 
average_tr_tips <- function(tr, fossils_older_than=0.6)
	{
	ex='
	tr = read.tree(file="", text="((human:0.9,chimp:1):1,gorilla:2);")
	tr2 = average_tr_tips(tr, fossils_older_than=0.6)
	write.tree(tr2, file="")
	tr2 = average_tr_tips(tr, fossils_older_than=0.05)
	write.tree(tr2, file="")
	
	'
	#require(BioGeoBEARS)	# for prt()

	# Check for negative branchlengths
	brlen_equal_below_0_TF = tr$edge.length <= 0
	if (sum(brlen_equal_below_0_TF) > 0)
		{
		tmptxt = paste(tr$edge.length[brlen_equal_below_0_TF], collapse=", ", sep="")
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

	tr6 = extend_tips_to_ultrametricize(obj=tr5, age_of_root=0, tips_end_at_this_date=NA, fossils_older_than=fossils_older_than)
	min(tr6$edge.length)

	# Check the output; if IT has negative branchlengths, return NA!!
	# Check for negative branchlengths
	brlen_equal_below_0_TF = tr5$edge.length <= 0
	if (sum(brlen_equal_below_0_TF) > 0)
		{
		tmptxt = paste(tr5$edge.length[brlen_equal_below_0_TF], collapse=", ", sep="")
		stoptxt = paste("\nFATAL ERROR in average_tr_tips(): the OUTPUT tree has branchlengths <= 0:\n", tmptxt, 
		"\nThis can sometimes happen if you (A) are accidentally including fossil tips (change 'fossils_older_than'), or\n",
		"(B) if average_tr_tips() introduced more negative branches (especially can happen with shallow branches). (Note: see e.g. impose_min_brlen() for a method to remove negative branchlengths.)\n", 
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
# Detect blank items
# (from e.g. reading in an Excel spreadsheet with
#  XLConnect::readWorksheetFromFile
#######################################################
#' Return TRUE for any "blank" items in list
#'
#' This function is especially useful for processing e.g. 
#' data read in or copied from Excel or other derived 
#' text files, where apparently 
#' blank cells may convert to "", \code{\link[base]{NA}}, 
#' \code{\link[base]{NaN}}, " ", "\t" (tab),
#' "NA", "NaN" (character versions of NA and NaN) etc. 
#' 
#' isblank_TF runs a series of tests for these various forms
#' of blank, and returns TRUE for each cell that matches 
#' any of the tests.
#'
#' @param items A vector (e.g., a column of a data.frame).
#' @return \code{blank_TF} A vector of \code{TRUE}/\code{FALSE}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
#' items = c(NA, NaN, "NA", "NaN", "na", "n/a", "nan", 0, " ", "", "\t", "\t\t")
#' isblank_TF(items)
#' 
isblank_TF <- function(items)
	{
	ex='
	# Showing which things currently match
	items = c(NA, NaN, "NA", "NaN", "na", "n/a", "nan", 0, " ", "", "\t", "\t\t")
	isblank_TF(items)
	'
	TF0 = is.null(items)
	if (TF0 == TRUE)
		{
		blank_TF = TRUE
		return(blank_TF)
		} # END if (TF0 == TRUE)
	TF1 = items == ""
	TF2 = items == " "
	TF3 = items == "\t"
	TF4 = is.na(items)
	TF5 = is.nan(items)
	TF6 = items == "NaN"
	TF7 = items == "NA"
	
	# Keep only the items where none of the above occur
	TFall = TF1 + TF2 + TF3 + TF4 + TF5 + TF6 + TF7
	
	# Correct for NA, NaNs
	TFall[is.na(TFall)] = 1
	TFall[is.nan(TFall)] = 1

	blank_TF = TFall > 0
	return(blank_TF)
	}




#######################################################
# is.not.na
#######################################################
#' Check for not NA
#'
#' A utility function, equivalent to \code{!is.na(x)}. 
#' 
#' @param x Item to check for NA
#' @return \code{TRUE} or \code{FALSE}
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}, \code{\link[base]{is.na}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' test=1
is.not.na <- function(x)
	{
	return(is.na(x) == FALSE)
	}


# Avoid negative branchlengths
#' Fix negative branchlengths, e.g. in a BEAST MCC tree
#'
#' Sometimes the MCC operation from TreeAnnotator can take
#' a posterior sample of trees from e.g. \code{BEAST} or \code{Beast2},
#' and produce a summary MCC tree where some of the branchlengths are
#' negative. 
#'
#' This function iterates through the tree from tips to root, and
#' if a negative branchlength (technically, a branchlength less 
#' than \code{min_brlen}) is encountered, it is changed to 
#' have the branchlength of \code{min_brlen}.  By default \code{min_brlen} is a
#' low value (0.01), indicating a branchlength close to 0. Branches
#' below are changed slightly to compensate for this change, keeping
#' the tips of the output tree at the same height as the tips in the
#' input tree.
#'
#' Advantage: this should produce a tree that can be run in 
#' BioGeoBEARS or other comparative methods without producing
#' the math errors that result from negative branchlengths.
#' (BEAST starting trees also have to have positive branchlengths.)
#'
#' Disadvantage: This is effectively a manual modification of your
#' tree, although a small one.  Branches with very short branchlengths
#' are likely statistically unresolved, so the actual topology and 
#' branchlength you are using are being arbitrarily determined. For a 
#' few branches in a large tree this likely doesn't effect downstream
#' analyses very much, but it might. An alternative is to run your 
#' analysis over one or a series of trees sampled from the posterior
#' sample of trees. 
#' 
#' MCC = Maximum Clade Credibility tree. This is the tree topology
#'       that contains the set of clades that have the highest 
#'       clade credibility. Once the topology is decided, the MCC
#'       algorithm then calculates node heights (e.g. mean or median)
#'       for the clades contained in the MCC tree, and also calculates
#'       the posterior probabilities of each clade (the branch supports),
#'       i.e. the frequency of each clade in the posterior sample.
#' 
#'       This procedure will not necessarily pick the absolute optimal
#'       tree (e.g. the Maximum Likelihood tree), but it usually is a 
#'       decent representation of the central tendency of the posterior
#'       sample of trees.  I suspect that the reason for negative 
#'       branchlengths is clades with very low support, which can 
#'       happen especially in large phylogenies or with gaps in the
#'       data matrix for some taxa.
#' 
#' @param phy An ape \code{\link[ape]{phylo}} object. 
#' @param min_brlen Any negative branches (or technically, branches below 
#'                  \code{min_brlen}) are converted to have length \code{min_brlen}.
#'                  Default is 0.01 (reasonable if the time tree's scale
#'                  is millions of years. Obviously, \code{min_brlen} should be adjusted
#'                  for different timescales or time units.
#' @param leave_BL0_terminals In "sampled ancestor" trees, some tips in 
#'                  sampled trees can have 0 branchlengths, indicating
#'                  that, according to the assumptions of the model, the
#'                  sampling rate, etc., it is somewhat probable that the
#'                  fossil is a member of a lineage that is directly ancestral
#'                  to later taxa. If \code{leave_BL0_terminals=TRUE} (default), then
#'                  these branches are left to have length zero (0.0). 
#' @param direct_ancestor_brlen Branches that are between 0 and \code{direct_ancestor_brlen} 
#'                  (inclusive) are treated as direct ancestors, and not modified.
#' @param printlevel Default is 2 (prints all messages). 1 prints messages if something is
#'                  edited, 0 prints nothing.
#' @return phy3 The edited phylogeny.
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
#' tr = read.tree(file="", text="((human:0.9,chimp:1):-0.1,gorilla:2);")
#' tr2 = impose_min_brlen(phy=tr, min_brlen=0.01, leave_BL0_terminals=TRUE, direct_ancestor_brlen=1e-07, printlevel=2)
#' write.tree(tr2, file="")
#' 
impose_min_brlen <- function(phy, min_brlen=0.01, leave_BL0_terminals=TRUE, direct_ancestor_brlen=1e-07, printlevel=2)
	{
	defaults='
	tr = read.tree(file="", text="((human:0.9,chimp:1):-0.1,gorilla:2);")
	tr2 = impose_min_brlen(phy=tr, min_brlen=0.01, leave_BL0_terminals=TRUE, direct_ancestor_brlen=1e-07, printlevel=2)
	write.tree(tr2, file="")
	
	min_brlen = 1e-6
	leave_BL0_terminals=TRUE
	direct_ancestor_brlen=1e-07
	'
	
	num_internal_nodes = phy$Nnode
	numtips = length(phy$tip.label)
	rootnodenum = numtips+1
	# likelihoods are computed at all nodes
	# make a list to store the 
	numnodes = numtips + num_internal_nodes
	
	# Reorder the edge matrix into pruningwise order
	# This is CRUCIAL!!
	phy2 <- reorder(phy, "pruningwise")
	phy2_orig = phy2
	
	observed_min_orig = min(phy2_orig$edge.length)
	TF = phy2_orig$edge.length < min_brlen
	num_branches_below_min = sum(TF)
	if (observed_min_orig >= min_brlen)
		{
		if (printlevel >= 2)
			{
			txt = paste0("No branches found =< min_brlen (", min_brlen, "). Returning original tree, pruningwise reordered.")
			cat("\n")
			cat(txt)
			cat("\n")
			}
		return(phy2)
		} else {
		if (printlevel >= 1)
			{
			txt = paste0(num_branches_below_min, " branches found < min_brlen (", min_brlen, "). Running downpass to edit branchlengths by pushing nodes down so that minimum branchlength=", min_brlen, ". Adjusting i:Leftnode,j:Rightnode,anc:Ancnode;...")
			cat("\n")
			cat(txt)
			cat("\n\n")
			}
		} 
		
		# END if (observed_min_orig > min_brlen)



	# DEFINE DOWNPASS THROUGH THE BRANCHES	
	tipnums <- 1:numtips
	i = 1
	edges_to_visit = seq(from=1, by=2, length.out=num_internal_nodes)

	#######################################################
	#######################################################
	# THIS IS A DOWNPASS FROM THE TIPS TO THE ROOT
	#######################################################
	#######################################################
	direct_ancestor_count = 0
	for (i in edges_to_visit)
		{
		# First edge visited is i
		#print(i)
		
		# Its sister is j 
		j <- i + 1
		#print(j)

		# Get the node numbers at the tips of these two edges		
		left_desc_nodenum <- phy2$edge[i, 2]
		right_desc_nodenum <- phy2$edge[j, 2]
		
		# And for the ancestor edge (i or j shouldn't matter, should produce the same result!!!)
		ancnode <- phy2$edge[i, 1]
		anc_edge_TF = phy2$edge[,2] == ancnode
		anc_edgenum = (1:nrow(phy2$edge))[anc_edge_TF]
		if (length(anc_edgenum) == 0)
			{
			anc_edgenum = NA
			}

		# Get the edge length (left, right, anc)
		Llength = phy2$edge.length[i]
		Rlength = phy2$edge.length[j]
		Alength = phy2$edge.length[anc_edgenum]
		
		txt = paste("ancnode:", ancnode, " left:", left_desc_nodenum, " right:", right_desc_nodenum, sep="")
		#print(txt)
		
		# If 0-length terminal branches are allowed (e.g. a 
		# sampled-ancestor tree), and
		# If a node is a terminal node,
		# and if it is of length 0 and the other branch > 0,
		# then don't change this branch length
		edit_brlength = TRUE

		if (leave_BL0_terminals == TRUE)
			{
			if (i <= numtips)
				{
				TF1  = (Llength >= 0) && (Llength <= direct_ancestor_brlen)
				TF2  = (Rlength >= 0) && (Rlength <= direct_ancestor_brlen)
				#cat(i, Llength, Rlength, sep="\t")
				#cat("\n")
				if ( TF1 && (Rlength > 0) )
					{
					edit_brlength = FALSE
					direct_ancestor_count = direct_ancestor_count + 1
					}
				} # END if (i <= numtips)
			if (j <= numtips)
				{
				if ( TF2 && (Llength > 0) )
					{
					edit_brlength = FALSE
					direct_ancestor_count = direct_ancestor_count + 1
					}
				} # END if (i <= numtips)
			} # END if (leave_BL0_terminals == TRUE)
		
		
		# If either the left or right branch is < min_brlen, change both
		if (edit_brlength == TRUE) {
		if ( (Llength < min_brlen) || (Rlength < min_brlen) )
			{
			cat(i, ":", left_desc_nodenum, ",", j, ":", right_desc_nodenum, ",anc:", ancnode, "; ", sep="")
			
			amount_to_add_L = min_brlen - Llength
			amount_to_add_R = min_brlen - Rlength
			amount_to_add = max(c(amount_to_add_L, amount_to_add_R))
			
			# Change both branchlengths by the same amount
			phy2$edge.length[i] = phy2$edge.length[i] + amount_to_add
			phy2$edge.length[j] = phy2$edge.length[j] + amount_to_add
			
			# Subtract that amount from the branch below
			# (may result in negative, which will be fixed further in the downpass)
			if (is.na(anc_edgenum) == FALSE)
				{
				phy2$edge.length[anc_edgenum] = phy2$edge.length[anc_edgenum] - amount_to_add
				}
			}
			} # END if (edit_brlength == TRUE)
		
		} # end downpass
	#######################################################
	# END DOWNPASS FROM THE TIPS TO THE ROOT
	#######################################################
	rootnode = ancnode
	phy3 = phy2

	num_branches_below_min2 = num_branches_below_min - direct_ancestor_count
	observed_min_brlen_new = min(phy3$edge.length)

	if (printlevel >= 1)
		{
		cat("\n")
		cat("\n")
		txt = paste0("Originally, ", num_branches_below_min2, " (non-ancestor) branches were found < min_brlen (", min_brlen, "). Running downpass to edit branchlengths by pushing nodes down so that minimum branchlength=", min_brlen, ". Adjusting i:Leftnode,j:Rightnode,anc:Ancnode;...")
		cat("\n")
		cat(txt)
		cat("\n\n")
		cat("New observed_min_brlen_new=", observed_min_brlen_new, sep="")
		cat("\n")
		}
	return(phy3)
	} # END impose_min_brlen <- function(phy, min_brlen=0.1)



#' Replacement for trim
#'
#' A quick version of \code{\link[gdata]{trim}}, so that we don't have
# to depend on \code{gdata}'s \code{\link[gdata]{trim}}.
#'
#' @param x An string to apply quicktrim to. 
#' @return Returns a string without leading or trailing whitespace.
#' @seealso \code{\link[gdata]{trim}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @references
#' \url{http://stackoverflow.com/questions/2261079/how-to-trim-leading-and-trailing-whitespace-in-r}
#' @examples
#' test = "   hello!\t\t"
#' test
#' quicktrim(test)
#' 
quicktrim <- function (x)
	{
	gsub("^\\s+|\\s+$", "", x)
	} # END quicktrim <- function (x)




#######################################################
# Compare two lists, produce an output table of 
# matches/mismatches
#######################################################
#' Compare two lists, produce an output table of matches/mismatches
#'
#' This function is for purposes like matching a list of tipnames in a phylogeny
#' against a list of taxa names in a traits file. The output table is 
#' more readable -- help for identifying small differences in taxon names
#' (e.g. Homo_sapiens vs. Homo sapiens vs. Homo Sapiens).
#'
#' Note: If you want to make your life easier, as well as that of 
#' your friendly neighborhood computational biologist, just always
#' use taxon names with "_" instead of spaces, and never use periods,
#' commas, quote marks, parentheses, etc. in your taxon names.  Many of
#' these special characters will mess up the reading of the data
#' into other programs (NEXUS and Newick files, etc.)
#'
#' @param names1 The first list of names to match.
#' @param names2 The second list of names to match.
#' @param listdesc1 Description of the first list (e.g., "tree file").
#' @param listdesc2 Description of the first list (e.g., "traits file").
#' @param list_or_file_txt Description for comparison output text (default="lists").  
#' @return something
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
#' names2 = c("sp1", "sp2", "sp3", "sp4")
#' names1 = c("sp1", "sp2", "sp3", "Sp4")
#' listdesc1="file1"
#' listdesc2="file2"
#' list_or_file_txt="lists"
#' compare_two_name_lists(names1, names2, listdesc1="file1", listdesc2="file2", list_or_file_txt="lists")
#' 
compare_two_name_lists <- function(names1, names2, listdesc1="file1", listdesc2="file2", list_or_file_txt="lists")
	{
	defaults='
	names2 = c("sp1", "sp2", "sp3", "sp4")
	names1 = c("sp1", "sp2", "sp3", "Sp4")
	listdesc1="file1"
	listdesc2="file2"
	list_or_file_txt="lists"
	compare_two_name_lists(names1, names2, listdesc1="file1", listdesc2="file2", list_or_file_txt="lists")
	'
	
	# Let's make the geography names the primary ones
	names2_found_TF = names2 %in% names1
	geognames_found_TF = names1 %in% names2
	
	names2_matched = names2[names2_found_TF == TRUE]
	geognames_matched = names1[geognames_found_TF == TRUE]
	names2_not_matched = names2[names2_found_TF == FALSE]
	geognames_not_matched = names1[geognames_found_TF == FALSE]

	# Total list
	all_names_from_both_inputs = sort(unique(c(names1, names2)))
	matched_or_missing = rep("", times=length(all_names_from_both_inputs))
	
	TF = all_names_from_both_inputs %in% names2_matched
	matched_or_missing[TF] = paste0("(found in both ", list_or_file_txt, ")")
	TF = all_names_from_both_inputs %in% geognames_matched
	matched_or_missing[TF] = paste0("(found in both ", list_or_file_txt, ")")
	TF = all_names_from_both_inputs %in% names2_not_matched
	matched_or_missing[TF] = paste0("found only in: ", listdesc2)
	TF = all_names_from_both_inputs %in% geognames_not_matched
	matched_or_missing[TF] = paste0("found only in: ", listdesc1)
	
	match_mismatch_table = cbind(all_names_from_both_inputs, matched_or_missing)
	match_mismatch_table_df = as.data.frame(match_mismatch_table, stringsAsFactors=FALSE)
	names(match_mismatch_table_df) = c("all_names_from_both_inputs", "matched_or_missing")
	row.names(match_mismatch_table_df) = 1:nrow(match_mismatch_table)
	
	return(match_mismatch_table_df)
	} # END compare_two_name_lists <- function(names1, names2, listdesc1="file1", listdesc2="file2")




#######################################################
# allQs
#######################################################
#' Checks a row/vector to see if it is all question marks ("?")
#' 
#' Rows with geography data consisting of all "?" (i.e., no information) 
#' are NOT allowed in standard BioGeoBEARS analysis.  Presumably any 
#' specimen/fossil was observed somewhere, whether or not the complete 
#' range is known.  
#'
#' By changing symbol_to_check, you can search for other values, e.g., all 0s.
#' 
#' This function is used by check_tipranges_for_allQs(). 
#' 
#' @rowvals The values of a row in tipranges@df
#' @symbol_to_check The symbol to check for uniformity in the row. Default "?".
#' @return \code{TRUE} or \code{FALSE}
#' @export
#' @seealso \code{\link{check_tipranges_for_allQs}}, \code{\link{getranges_from_LagrangePHYLIP}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#'	 @cite SmithRee2010_CPPversion
#' @examples
#' 
#' x = c("?", "?", "?", "?")
#' allQs(rowvals=x, symbol_to_check="?")
#' 
#' x = c("0", "0", "0", "0")
#' allQs(rowvals=x, symbol_to_check="0")
#' 
allQs <- function(rowvals, symbol_to_check="?")
	{
	ex='
	x = c("?", "?", "?", "?")
	allQs(rowvals=x, symbol_to_check="?")

	x = c("0", "0", "0", "0")
	allQs(rowvals=x, symbol_to_check="0")
	'
	TF = rowvals == symbol_to_check
	if (sum(TF) == length(rowvals))
		{
		return(TRUE)
		} else {
		return(FALSE)
		}
	}


#' Remove strings of commas from lists of text strings
#'
#' When pruning lists of OTUs, you can end up with lots of ",,,," 
#' strings of various lengths. This function attempts to remove
#' these with multiple \code{\link[base]{gsub}} calls.
#'
#' There is probably a better way, but regex strings with 
#' special characters are annoyingly hard to figure out.
#'
#' @param list_of_strings A list of strings
#' @return OTUs_list_of_lists_reduced A list of strings, with all multiple,
#' starting, and ending commas removed.
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' test=1
#' 
#' list_of_strings = NULL
#' list_of_strings[[1]] = c("A,,B", "B,C,,D,,,E,,,,F")
#' list_of_strings[[2]] = c(",A,,B", "B,C,,D,,,E,,,,F")
#' list_of_strings[[3]] = c(",,,A,,B", "B,C,,D,,,E,,,,F")
#' list_of_strings[[4]] = c(",,,A,,B", "B,C,,D,,,E,,,,F,")
#' list_of_strings[[5]] = c(",,,A,,B", "B,C,,D,,,E,,,,F,,,")
#' list_of_strings
#' 
#' OTUs_list_of_lists_reduced = remove_multiple_commas_from_strings_in_list(list_of_strings)
#' OTUs_list_of_lists_reduced
#' 
remove_multiple_commas_from_strings_in_list <- function(list_of_strings)
	{
	setup='
	list_of_strings = NULL
	list_of_strings[[1]] = c("A,,B", "B,C,,D,,,E,,,,F")
	list_of_strings[[2]] = c(",A,,B", "B,C,,D,,,E,,,,F")
	list_of_strings[[3]] = c(",,,A,,B", "B,C,,D,,,E,,,,F")
	list_of_strings[[4]] = c(",,,A,,B", "B,C,,D,,,E,,,,F,")
	list_of_strings[[5]] = c(",,,A,,B", "B,C,,D,,,E,,,,F,,,")
	list_of_strings
	
	OTUs_list_of_lists_reduced = remove_multiple_commas_from_strings_in_list(list_of_strings)
	OTUs_list_of_lists_reduced
	'
	
	# Convert variable name
	OTUs_list_of_lists_reduced = list_of_strings
	
	# Remove all resulting extraneous commas
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,", replacement=",", x=OTUs_list_of_lists_reduced)

	# Remove commas at the end of strings
	lengths = sapply(X=OTUs_list_of_lists_reduced, FUN=nchar)
	names(lengths) = NULL
	last_chars = mapply(FUN=substr, OTUs_list_of_lists_reduced, start=lengths, stop=lengths)
	names(last_chars) = NULL
	ending_comma_TF = last_chars == ","

	stops = -1 + lengths[ending_comma_TF]
	starts = rep(1, times=length(stops))

	no_ending_commas = mapply(FUN=substr, OTUs_list_of_lists_reduced[ending_comma_TF], start=starts, stop=stops)
	names(no_ending_commas) = NULL
	no_ending_commas
	OTUs_list_of_lists_reduced[ending_comma_TF] = no_ending_commas


	# Remove commas at the beginning of strings
	first_chars = sapply(X=OTUs_list_of_lists_reduced, FUN=substr, start=1, stop=1)
	names(first_chars) = NULL
	starting_comma_TF = first_chars == ","
	sum(starting_comma_TF)
	#cbind(first_chars, starting_comma_TF)[500:750,]

	starts = rep(2, times=sum(starting_comma_TF))
	stops = sapply(X=OTUs_list_of_lists_reduced[starting_comma_TF], FUN=nchar)

	no_starting_commas = mapply(FUN=substr, OTUs_list_of_lists_reduced[starting_comma_TF], start=starts, stop=stops)
	names(no_starting_commas) = NULL
	OTUs_list_of_lists_reduced[starting_comma_TF] = no_starting_commas

	#OTUs_list_of_lists_reduced[500:750]

	#sum(string::str_count(string=OTUs_list_of_lists_reduced, pattern="fossil"))
	
	return(OTUs_list_of_lists_reduced)
	}










# Check if object is a list BUT NOT a data.frame
# At some point, R started saying that is.list() on a data.frame would 
# evaluate to TRUE. This is ANNOYING.
# The function is_list_not_dataframe() function gives the desired behavior.
is_list_not_dataframe <- function(obj)
	{
	example_code='
	# Set up a BioGeoBEARS_run_object
	BioGeoBEARS_run_object = define_BioGeoBEARS_run()

	# Get the dmat (trivial case, non-stratified Psychotria example)
	tmpres = get_dmat_times_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object, numstates=NULL, max_range_size=4)
	tmpres
	class(tmpres)
	length(tmpres)
	
	# Is a matrix a list?
	is.list(tmpres$dmat)
	is_list_not_dataframe(obj=tmpres$dmat)
	
	# Convert the matrix to a data.frame
	tmpdf = as.data.frame(tmpres)
	tmpdf
	class(tmpdf)
	length(tmpdf)
	
	# Is a data.frame a list?
	is.list(tmpdf)
	is_list_not_dataframe(obj=tmpdf)

	ana_events_table = example_ana_events_table()
	is.list(ana_events_table)
	is_list_not_dataframe(obj=ana_events_table)

	is.list(c(NA, NA))
	is_list_not_dataframe(obj=c(NA, NA))
	
	# How about a dmat
	dmat = structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), dim = c(4L, 
4L))
	obj = dmat
	is_list_not_dataframe(obj)
	
	' # END example_code


	TF1 = is.list(obj) # This can evaluate to TRUE for a data.frame
	
	# But this can't
	if ( ("data.frame" %in% class(obj)) == TRUE)
		{
		not_a_df_TF = FALSE
		} else {
		not_a_df_TF = TRUE
		}
	
	TF = (TF1 + not_a_df_TF) == 2
	return(TF)
	} # END is_list_not_dataframe <- function(obj)



