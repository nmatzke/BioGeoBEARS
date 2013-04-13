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

