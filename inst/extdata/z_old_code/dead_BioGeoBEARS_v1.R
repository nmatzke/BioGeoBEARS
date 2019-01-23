#######################################################
# These are "dead" functions -- never used, and 
# so removed from version 1.1.2+
# Nicholas J. Matzke
# 2019-01-22
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











#' Assign rate categories by CRP
#'
#' Assigns rate categories to an input matrix bye CRP (Chinese Restaurant
#' Process, a form of clustering where the number of clusters is estimated
#' and can range from 1 to infinity (or really, the maximum number of clusters
#' is the number of things being clustered).
#'
#' @param nrows number of rows in the matrix (order of the square matrix)
#' @param alpha the clustering (or concentration) parameter of the Chinese Restaurant Process (CRP)
#' @param priorlambda a hyperparameter specifying the mean of the exponential for
#' drawing a rate for each category.  This could itself be drawn from a hyperprior
#' in another function.
#' @return Qmat, the instantaneous transition matrix (Q matrix)
#' @references
#' \url{http://en.wikipedia.org/wiki/Chinese_restaurant_process}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' 
assign_rate_categories_rchinese <- function(nrows, alpha, priorlambda = 0.1)
	{
	# empty list of categories
	Qcategories = ones(nrows)
	Qmat = Qcategories
	diag(Qcategories) = 0
	numrates = nrows * nrows - nrows

	categories = rchinese(numrates, alpha)
	
	# randomly assign the categories to cells
	Qcategories[Qcategories == 1] = categories
	
	# randomly draw rates for each category
	uniq_categories = sort(unique(categories))
	
	# using a uniform as the prior
	numcats = length(uniq_categories)
	cat("# rate categories = ", numcats, "\n", sep="")
	#uniq_rates = runif(numcats, priormin, priormax)
	uniq_rates = rexp(numcats, priorlambda)
	
	for (i in 1:numcats)
		{
		tmpcat = uniq_categories[i]
		Qmat[Qcategories == tmpcat] = uniq_rates[i]
		}
	
	# normalize the diagonals so that the rows sum to 1:
	Qmat = mat2q(Qmat)
	return(Qmat)
	}


run_BioGeoBEARS_reorder <- function()
		{
	#######################################################
	# Reorder Robin Becks tree so that the BioGeoBEARS plots look better
	#######################################################
	# Convert BioGeoBEARS object from one tree to another
	# tr1 = new tree
	# tr2 = original tree, used the BioGeoBEARS analysis
	wd  = "/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/z_help/Robin_Beck/"
	setwd(wd)

	# Original tree
	trfn = "PBG_TE_no_Alo_IGR_topcons_MCC_no_zero.tre"
	tr2 = read.tree(trfn)
	plot(tr2)

	# Ladderize the tree
	# GOOD LADDERIZE
	tr1 = ladderize_and_reorder(phy=tr2, right=FALSE)
	plot(tr1)

	tr1 = rotate_tips(tr=tr1, tipname1="Kulbeckia", tipname2="Paranyctoides")
	plot(tr1)

	tr1 = rotate_tips(tr=tr1, tipname1="Paranyctoides", tipname2="Sheikhdzheilia")
	plot(tr1)




	new_trfn = "tree_ladderized.newick"
	write.tree(phy=tr1, file=new_trfn)



	resfns = c("Mammals_BAYAREALIKE+J_M0_unconstrained_v1.Rdata",
	"Mammals_BAYAREALIKE_M0_unconstrained_v1.Rdata",
	"Mammals_DIVALIKE+J_M0_unconstrained_v1.Rdata",
	"Mammals_DIVALIKE_M0_unconstrained_v1.Rdata",
	"Mammals_DEC+J_M0_unconstrained_v1.Rdata",
	"Mammals_DEC_M0_unconstrained_v1.Rdata")

	for (i in 1:length(resfns))
		{
		# Convert each file
	
		# Loads to res
		load(file=resfns[i])
		res2 = BioGeoBEARS_reorder(res, tr1, tr2, trfn_for_BGB_inputs=new_trfn)
	
		new_resfn = gsub(pattern="\\.Rdata", replacement="_ladderized.Rdata", x=resfns[i])
		res = res2
		# Loads to res
		save(res, file=new_resfn)
		} # END for (i in 1:length(resfns))

	} # END run_BioGeoBEARS_reorder









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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
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
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://www.brianomeara.info/tutorials/aic}
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




