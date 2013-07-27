# source("/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_basics_v1.R")

require("ape")
require("rexpokit")
require("cladoRcpp")




# axisPhylo() with more flexibility in labeling

#######################################################
# axisPhylo2
#######################################################
#' axisPhylo with more flexibility in labeling
#' 
#' Hacking axisPhylo to make it more flexible
#' 
#' @param side The side to plot on (default 1, bottom)
#' @param roundlabels Number of digits to round to, if desired
#' @param minage Starting age, if desired
#' @param ... Additional arguments to standard functions
#' @return nothing
#' @export
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'	 @cite FosterIdiots
#' @examples
#' testval=1
#' 
axisPhylo2 <- function (side = 1, roundlabels=FALSE, minage=0, ...) 
	{
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    if (lastPP$type %in% c("phylogram", "cladogram")) {
        if (lastPP$direction %in% c("rightwards", "leftwards")) {
            x <- pretty(lastPP$xx)
            if (lastPP$direction == "rightwards") 
                maxi <- max(lastPP$xx)
            else {
                maxi <- min(lastPP$xx)
                x <- -x
            }
        }
        else {
            x <- pretty(lastPP$yy)
            if (lastPP$direction == "upwards") 
                maxi <- max(lastPP$yy)
            else {
                maxi <- min(lastPP$yy)
                x <- -x
            }
        }
    }
    if (roundlabels)
		{
		axis(side = side, at = c(maxi - x), labels = abs(minage+x), ...)
		}
	else
		{
		axis(side = side, at = c(maxi - x), labels = round(abs(minage+x), digits=roundlabels), ...)
		}
	}



#######################################################
# Matrix functions
#' Utility functions to help deal with matrices
#######################################################

#######################################################
# normat
#######################################################
#' Normalize a transition matrix
#' 
#' \code{normat} normalizes a square transition matrix, such that each row sums to 0, and 
#' the diagonal equals the negative of the sum of the rest of the cells in the row.  This 
#' matrix can then be exponentiated by values of \code{t} (time or another measure of branch
#' length) to produce transition probabilities for any given value of \code{t}.
#' 
#' See \cite{FosterIdiots} for a succinct summary of transition matrices and their exponentiation.
#' 
#' @param relative_matrix A square matrix giving the relative probabilities/weights of transitions.
#' @return \code{m} A Q matrix, i.e. normalized transition matrix (Qmat)
#' @export
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'	 @cite FosterIdiots

#' @examples
#' testval=1
#' 
normat <- function(relative_matrix)
	{
	# make the rows sum to zero (not 1), as in here:
	# http://en.wikipedia.org/wiki/Substitution_model
	m = as.matrix(relative_matrix)
	
	diag(m) = 0
	rowsums = rowSums(m)
	diag(m) = -rowsums
	
	return(m)
	}



#######################################################
# prune_specimens_to_species
#######################################################
#' Take a tree and species names/geography table and produce a pruned tree and tipranges object
#' 
#' This function takes a tree and species names/geography table and produces a pruned tree 
#' and (optionally) a tipranges object.
#' 
#' Often, users will have an phylogeny where the tips/OTUs (operational taxonomic units) are 
#' specimens rather than species.  The analyses done by models like DEC, DEC+J, etc., in programs 
#' like LAGRANGE and BioGeoBEARS, assume as a core part of the model that species might occupy
#' more than one areas.  A phylogeny of specimens, then, would not be an appropriate input 
#' to these programs, as each single specimen can only be found in one region. The exception 
#' would occur when the researcher is confident that each species lives in only one region; in
#' that case, the specimen geography is representative of the species geography.
#' 
#' This function requires a table containing 
#'
#' (1) Column "OTUs": all tipnames in the input tree (often, 
#' original specimen/original OTU names) ); 
#' 
#' (2) Column "species": the corresponding species names;
#' 
#' (3) optionally, the geographic range inhabited by each specimen (column "region").  If 
#' an OTU has more than one geographic range in the original table, these should be 
#' split by "|".
#' 
#' When the pruning occurs, all tips belonging to the same species are cut, except the first.
#'
#' NOTE: Tips that should be cut because they are outgroups, or because they are geographically 
#' outside of your domain of analysis, should be represented in xls$region by "out_group" or "Out". These
#' will be cut from the final tree/geography table.
#' 
#' @param original_tr The input tree (an \code{\link[ape]{ape}} \code{\link[ape]{phylo}} object).
#' @param xls The input table (a \code{\link[base]{data.frame}})
#' @param group_name The name of the clade in the tree.  For use in plots and output files. Default="default".
#' @param titletxt Additional text for the plots. Default "".
#' @param areas_abbr An optional table, containing the abbreviations (e.g. letters) corresponding to each
#' region in xls$region. Default is NULL, in which case the program imposes A, B, C, D, etc. \code{areas_abbr}
#' must have column headings \code{abbr} and \code{letter}.
#' @param plot_intermediate If TRUE, the starting, ending, and intermediate stages of tree pruning are plotted.
#' @return The outputs are a \code{\link[base]{list}} with a pruned tree and, optionally, a tipranges object.
#' @export
#' @seealso \code{\link[ape]{drop.tip}}, \code{\link{define_tipranges_object}}, 
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' testval=1
#' tipranges_object = define_tipranges_object()
#' tipranges_object
#'
#' areanames = getareas_from_tipranges_object(tipranges_object)
#' areanames
#'
prune_specimens_to_species <- function(original_tr, xls, group_name="default", titletxt="", areas_abbr=NULL,  plot_intermediate=TRUE)#, inputs="Robjects", outputs="Robjects")
	{
	runjunk='
	group_name="default"; titletxt=""; areas_abbr=NULL;  plot_intermediate=TRUE
	group_name="default"; titletxt=""; areas_abbr=NULL;  plot_intermediate=FALSE
	'
	
	# Check that input table has the correct column headings
	if ("OTUs" %in% names(xls) == FALSE)
		{
		stoptxt = "\nFATAL ERROR: input table 'xls' must have a column named 'OTUs'\n"
		cat(stoptxt)
		stop(stoptxt)
		}
	if ("species" %in% names(xls) == FALSE)
		{
		stoptxt = "\nFATAL ERROR: input table 'xls' must have a column named 'species'\n"
		cat(stoptxt)
		stop(stoptxt)
		}	
	
	#################################################################
	# CLEAN UP NAMES
	# Biologists love Excel but are often very sloppy about making sure their names in the 
	# table match the names in their tree.  This section fixes the more common mistakes.
	#################################################################
	
	# Convert all "_" to " "
	xls$OTUs = gsub(pattern="_", replacement=" ", x=xls$OTUs)
	xls$species = gsub(pattern="_", replacement=" ", x=xls$species)

	# Also remove all "'"
	xls$OTUs = gsub(pattern="'", replacement="", x=xls$OTUs)
	xls$species = gsub(pattern="'", replacement="", x=xls$species)


	# Before converting spaces to underscores, trim all beginning and ending spaces
	#require(gdata)	# for trim
	xls$OTUs = trim(xls$OTUs)
	xls$species = trim(xls$species)
	
	# Also in regions. Jesus, MVD!  What kind of crazy person has some cells 
	# say "Out" and others say "Out " ??????!???!?!?!?!?!
	if ("region" %in% names(xls) == TRUE)
		{
		xls$region = trim(xls$region)
		}

	# Convert all spaces to "_"
	xls$OTUs = gsub(pattern=" ", replacement="_", x=xls$OTUs)
	xls$species = gsub(pattern=" ", replacement="_", x=xls$species)

	# Load the tree file
	#nexfn = nexus_trfns[i]
	#original_tr = read.nexus(nexfn)

	# Cut trailing spaces and underscores from tip.labels
	# Also cut "'", which e.g. Mesquite inserts whenever an OTU name has spaces, but which
	# e.g. MrBayes won't read. (NEVER USE SPACES IN COMPUTER STUFF, PEOPLE!!)
	tmp_tipnames = original_tr$tip.label
	tmp_tipnames = gsub(pattern="'", replacement="", x=tmp_tipnames)
	tmp_tipnames = gsub(pattern="_", replacement=" ", x=tmp_tipnames)
	tmp_tipnames = trim(tmp_tipnames)
	tmp_tipnames = gsub(pattern=" ", replacement="_", x=tmp_tipnames)
	original_tr$tip.label = tmp_tipnames
	
	
	#######################################################
	# Now, go through the unique species in the XLS file,
	# for each, keep just one OTU, remove the rest
	#######################################################
	tr = original_tr	# copy the original tree

	if (plot_intermediate == TRUE)
		{
		plot(tr, cex=0.5)
		axisPhylo()
		#titletxt = paste(group_name, ": BEAST MCC tree", "\nNEXUS file='", nexfn, "'", sep="")
		titletxt = paste(group_name, titletxt, sep="")
		title(titletxt)
		}

	#######################################################
	# Error check
	#######################################################

	# Check for names matching tips in tree
	tipnames = sort(tr$tip.label)
	tablenames = sort(xls$OTUs)

	length(tipnames)
	length(tablenames)
	tipnames == tablenames

	# Find phylogeny tips not in tablenames
	tipnames_in_tablenames_TF = tipnames %in% tablenames
	tipnames_not_in_table = tipnames[tipnames_in_tablenames_TF == FALSE]

	# Find phylogeny tips not in tablenames
	tablenames_not_in_tipnames_TF = tablenames %in% tipnames
	tablenames_not_in_tips = tablenames[tablenames_not_in_tipnames_TF == FALSE]

	num_tips_nomatch = sum(tipnames_in_tablenames_TF == FALSE)
	num_table_nomatch = sum(tablenames_not_in_tipnames_TF == FALSE)

	if (num_tips_nomatch > 0)
		{
		cat("\nWARNING: Tip names not in table:\n", sep="")

		cat("length(tipnames)=", length(tipnames), "\n", sep="")
		cat("length(tablenames)=", length(tablenames), "\n", sep="")
		
		TF = tipnames == tablenames
		tmptxt = paste(as.character(TF), collapse=" ", sep="")
		cat("tipnames == tablenames:\n", tmptxt, "\n", sep="")

		cat("\ntipnames_not_in_table:\n", sep="")
		cat(tipnames_not_in_table, sep="\n")

		#cat(sort(tipnames), sep="\n")
		#cat(sort(tablenames), sep="\n")
		}

	if (num_table_nomatch > 0)
		{
		cat("\nWARNING: Table names not in tips:\n")
		cat("length(tipnames)=", length(tipnames), "\n", sep="")
		cat("length(tablenames)=", length(tablenames), "\n", sep="")
		
		TF = tipnames == tablenames
		tmptxt = paste(as.character(TF), collapse=" ", sep="")
		cat("tipnames == tablenames:\n", tmptxt, "\n", sep="")

		cat("\ntablenames_not_in_tips:\n", sep="")
		cat(tablenames_not_in_tips, sep="\n")
	
		#cat(sort(tipnames), sep="\n")
		#cat(sort(tablenames), sep="\n")
		}
	
	
	#######################################################
	# Cut the outgroups and taxa outside of your analysis 
	# (represented by "out_group" or "Out")
	#######################################################
	if ("region" %in% names(xls) == TRUE)
		{
		# Get outgroups and cut; also cut blank ranges
		outgroups_TF = ((xls$region == "Out") + (xls$region == "out_group") + (xls$region == "") + (xls$region == " ")) == 1
		outgroup_tips = xls$OTUs[outgroups_TF]
		outgroup_species = unique(xls$species[outgroups_TF])
		tmptxt = paste(outgroup_species, collapse=", ", sep="")
		cat("\nRemoving outgroups: ", tmptxt, sep="")

		# Drop the tips
		tr = drop.tip(phy=tr, tip=outgroup_tips)

		if (plot_intermediate == TRUE)
			{
			plot(tr, cex=0.5)
			title("Cutting outgroups")
			axisPhylo()
			}

		cat("...done.", sep="")
		
		# Also, remove from xls
		xls = xls[outgroups_TF==FALSE, ]
		}	
	
	# Get unique species
	species_uniq = unique(xls$species)
	species_uniq


	# Assemble geography data
	#dispersal_xlsfn = "__dispersal_matrix_v01.xlsx"
	#areas_abbr = read.xls(dispersal_xlsfn, sheet="area_codes")
	if ("region" %in% names(xls) == TRUE)
		{
		if (is.null(areas_abbr))
			{
			areas_split = sapply(X=xls$region, FUN=strsplit2, split="\\|")
			areas_unlisted = unlist(areas_split)
			abbr = sort(unique(areas_unlisted))
			letters = LETTERS[1:length(abbr)]
			areas_abbr = adf2(data.matrix(cbind(abbr, letters)))
			}
		}
	areas_abbr

	# Go through and drop all non-uniq tips
	tmp_species_ranges_table = NULL
	for (j in 1:length(species_uniq))
		{
		#######################################################
		# Prune down the tree
		#######################################################
		# Collapsing species
		species_name = species_uniq[j]
		cat("\nCollapsing OTUs to species: ", species_name, sep="")

		OTU_names_for_this_species_TF = xls$species == species_name
		OTU_names_for_this_species = xls$OTUs[OTU_names_for_this_species_TF]
	
		# The OTU to keep/rename
		OTU_to_keep = OTU_names_for_this_species[1]
	
		# Leave the first OTU, drop the rest
		OTU_names_for_this_species = OTU_names_for_this_species[-1]
	
		# Find the tips that match this species
		tr = drop.tip(phy=tr, tip=OTU_names_for_this_species)
	
		# Rename the kept tip
		tr$tip.label[tr$tip.label == OTU_to_keep] = species_name
	
		if (plot_intermediate == TRUE)
			{
			plot(tr, cex=0.5)
			title("Cutting redundant OTUs")
			axisPhylo()
			}
		cat("		...done.", sep="")



		if ("region" %in% names(xls) == TRUE)
			{
			#######################################################
			# Assemble geography data
			#######################################################
			region_codes = xls$region[OTU_names_for_this_species_TF]
			# Some tips have multiple, find with "|"
			region_codes_multiple_TF = grepl(pattern="\\|", x=region_codes)
	
	
			if (sum(region_codes_multiple_TF) > 0)
				{
				list_o_multiple_regions = NULL
				for (k in 1:sum(region_codes_multiple_TF))
					{
					tmp_regions_before_split = region_codes[region_codes_multiple_TF][k]
					tmp_regions_after_split = strsplit(tmp_regions_before_split, split="\\|")[[1]]
					list_o_multiple_regions = c(list_o_multiple_regions, c(tmp_regions_after_split))
					}
				region_codes = region_codes[region_codes_multiple_TF == FALSE]	# cut the ones that were split
				all_regions_after_splitting = c(region_codes, list_o_multiple_regions)		# add them, post-parsing
				} else {
				all_regions_after_splitting = region_codes
				}
	
			areas_for_this_sp = sort(unique(all_regions_after_splitting))
			tmp_letter_codes_for_this_species = NULL
			for (k in 1:length(areas_for_this_sp))
				{
				# Abbreviation for this unique area
				abbr_TF = areas_abbr$abbr == areas_for_this_sp[k]
		
				# Error check
				if (sum(abbr_TF) == 0)
					{
					stoptxt = cat("\nMAJOR WARNING!: the areaname '", areas_for_this_sp[k], "' of species '", species_uniq, "' is not found in your abbreviation table!\n", sep="")
					cat(stoptxt)
					#stop(stoptxt)
					}
		
				tmp_letter_code = areas_abbr$letter[abbr_TF]
				tmp_letter_codes_for_this_species = c(tmp_letter_codes_for_this_species, tmp_letter_code)
				}
	
			# Add to table
			range_letter_strs = paste(tmp_letter_codes_for_this_species, collapse="", sep="")
			tmprow = c(species_name, range_letter_strs)
			tmp_species_ranges_table = rbind(tmp_species_ranges_table, tmprow)
			} # End if regions
		} # end forloop

	if ("region" %in% names(xls) == TRUE)
		{
		tmp_species_ranges_table = dfnums_to_numeric(adf2(tmp_species_ranges_table))
		names(tmp_species_ranges_table) = c("species", "letterrange")

		tipranges = letter_strings_to_tipranges_df(letter_strings=tmp_species_ranges_table$letterrange, letter_codes_in_desired_order=areas_abbr$letter, tipnames_in_order=tmp_species_ranges_table$species)
		tipranges

		# Then, check the output
		tipnames1 = sort(tr$tip.label)
		tipnames2 = sort(rownames(tipranges@df))
		
		if (length(tipnames1) != length(tipnames2))
			{
			stoptxt = paste("\nFATAL ERROR: Your output tree and your output tipranges have different lengths:\nlength(tr$tip.label)=", length(tr$tip.label), ", length(rownames(tipranges@df))=", length(rownames(tipranges@df)), "\n", sep="")
			cat(stoptxt)
			stop(stoptxt)
			}
		
		TF = tipnames1 == tipnames2
		if (sum(TF) != length(TF))
			{
			stoptxt = paste("\nFATAL ERROR: Your output tree tips, and your output tipranges names, don't match:\n", sep="")
			cat(stoptxt)
			
			cat(tipnames1, ", ")
			cat(tipnames2, ", ")
			cat(TF, ", ")
			
			stop(stoptxt)
			}
		
		
		TF = rowSums(tipranges@df) == 0
		if (sum(TF) > 0)
			{
			stoptxt = paste("\nFATAL ERROR: tipranges contains ", sum(TF), " rows with 0 ranges:\n", sep="")
			cat(stoptxt)
			
			print(tipranges@df[TF, ])
			stop(stoptxt)
			}



		}

	is.ultrametric(tr)













	if (plot_intermediate == TRUE)
		{
		# Plot final tree
		titletxt = paste(group_name, ": Final pruned tree", sep="")
		plot(tr, cex=0.5)
		axisPhylo()
		title(titletxt)
		}


	pruning_results = list()
	if ("region" %in% names(xls) == TRUE)
		{
		pruning_results$tr = tr
		pruning_results$tipranges = tipranges
		} else {
		pruning_results$tr = tr
		}

	return(pruning_results)
	}









#######################################################
# Lagrange-related functions
#######################################################





#######################################################
# getareas_from_tipranges_object
#######################################################
#' Get the names of the areas in a tipranges object
#' 
#' This function extracts the names of the areas in a \code{tipranges} object.  Just a shortcut for 
#' \code{names(tipranges@@df)}.
#' 
#' @param tipranges An object of class \code{tipranges}.
#' @return \code{areanames}, a list of the names of the areas
#' @export
#' @seealso \code{\link{define_tipranges_object}}, \code{\link[cladoRcpp]{areas_list_to_states_list_old}},
#' \code{\link{areas_list_to_states_list_new}}, \code{\link{tipranges_to_tip_condlikes_of_data_on_each_state}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' testval=1
#' tipranges_object = define_tipranges_object()
#' tipranges_object
#'
#' areanames = getareas_from_tipranges_object(tipranges_object)
#' areanames
#'
getareas_from_tipranges_object <- function(tipranges)
	{
	areanames = names(tipranges@df)
	return(areanames)
	}


#######################################################
# get_tiplabel_ranges
#######################################################
#' For each tip, get a text string of the areas in a tipranges object.
#' 
#' This function extracts the names of the areas in a \code{tipranges} object.  Just a shortcut for 
#' \code{names(tipranges@@df)}.
#' 
#' @param tipranges An object of class \code{tipranges}.
#' @param tr An ape phylo object.
#' @param sep Input to \code{\link[base]{paste}}.
#' @return \code{areanames}, a list of the names of the areas
#' @export
#' @seealso \code{\link{define_tipranges_object}}, \code{\link[cladoRcpp]{areas_list_to_states_list_old}},
#' \code{\link{areas_list_to_states_list_new}}, \code{\link{tipranges_to_tip_condlikes_of_data_on_each_state}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @examples
#' testval=1
#' tipranges_object = define_tipranges_object()
#' tipranges_object
#'
#' areanames = getareas_from_tipranges_object(tipranges_object)
#' areanames
#'
get_tiplabel_ranges <- function(tipranges, tr, sep="")
	{
	# Make sure the order matches the tr order
	tipranges = order_tipranges_by_tr(tipranges, tr)

	tmpdf = tipranges@df
	tmpnames = names(tmpdf)
	tmpnames_list = list()
	
	# TRUE/FALSE occupied
	tmpdf_TF = (tmpdf != 0)
	
	for (i in 1:nrow(tmpdf))
		{
		tmpnames_list[[i]] = c(tmpnames[tmpdf_TF[i,]])
		}
	
	tipranges_txt = sapply(FUN=paste, tmpnames_list, collapse="", sep=sep)
	return(tipranges_txt)
	}



#######################################################
# tipranges_to_tip_condlikes_of_data_on_each_state
#######################################################
#' Convert a tipranges object to the tip likelihoods
#' 
#' This function takes a tipranges object, and converts it to tip likelihoods for 
#' input into the likelihood calculations of \code{\link{calc_loglike_sp}}.
#' 
#' This (like LAGRANGE (\cite{ReeSmith2008}) and every other available program) assumes that the geographic
#' ranges at the tips are known with certainty.  Reality may be different, particularly
#' for sparsely-studied, scarce, or fossil taxa.  In such a case, a detection model 
#' is needed to specify the likelihood of the observation data under each possible
#' geographic range at the tips.
#'
#' Note that data likelihoods under this or that hypothesis are not the same thing as
#' probabilities.  E.g., with DNA, if sequencing machine says that the base could be either
#' A or C, but not G or T, then the likelihood of the data for that nucleotide position for
#' that species would be 1 1 0 0, not 0.5 0.5 0 0.  See Felsenstein (2004), p. 255, for more.
#' 
#' @param tipranges An object of class \code{tipranges}.
#' @param phy A phylogenetic tree (\code{\link[ape]{ape}} object of class \code{\link[ape]{phylo}})
#' @param states_list A complete list of the different states, of class \code{\link{list}} form
#' @param maxareas The maximum number of areas in a geographic range, if the user does 
#' @return \code{tip_condlikes_of_data_on_each_state} For each tip/row, likelihood of that tip's data under each possible true geographic range (columns)
#' @export
#' @seealso \code{\link{define_tipranges_object}}, \code{\link{getareas_from_tipranges_object}}, 
#' \code{\link{areas_list_to_states_list_new}}, 
#' \code{\link[cladoRcpp]{areas_list_to_states_list_old}}, \code{\link{binary_ranges_to_letter_codes}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Felsenstein2004
#' @examples
#' testval=1
#' # Define a tipranges object
#' tipranges_object = define_tipranges_object()
#' tipranges_object
#'
#' areanames = getareas_from_tipranges_object(tipranges_object)
#' areanames
#'
#' # Specify phylogeny to go with default tipranges object
#' newick_str = "((tip1:1,tip2:1):1,tip3:2):1;"
#' phy = read.tree(file="", text=newick_str)
#' 
#' # Here, we will assume the maximum range size is all areas, but it could be smaller
#' maxareas = length(areanames)
#' \dontrun{
#' states_list = areas_list_to_states_list_old(areas=areanames, include_null_range=TRUE, 
#' maxareas=maxareas)
#' states_list
#' }
#' 
#' states_list = areas_list_to_states_list_new(areas=areanames, include_null_range=TRUE, 
#' maxareas=maxareas)
#' states_list
#' 
#' tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(
#' tipranges=tipranges_object, phy=phy, states_list=states_list, maxareas=maxareas )
#' tip_condlikes_of_data_on_each_state
#' 
tipranges_to_tip_condlikes_of_data_on_each_state <- function(tipranges, phy, states_list=NULL, maxareas=length(getareas_from_tipranges_object(tipranges)) )
	{
	default='
	tipranges
	phy
	maxareas=length(getareas_from_tipranges_object(tipranges))
	'
	require(cladoRcpp)
	
	# Reorder the edge matrix into pruningwise order
	# This (may be) CRUCIAL!!
	#phy2 <- reorder(phy, "pruningwise")
	# NOTE: phy and phy2 have the tips in the same order, it's just the
	# edge matrix that gets reordered...
	
	areanames = getareas_from_tipranges_object(tipranges)
	
	if (is.null(states_list))
		{
		cat("Note: tipranges_to_tip_condlikes_of_data_on_each_state() is\n")
		cat("      creating 'states_list' automatically.\n")
		states_list = areas_list_to_states_list_new(areas=areanames, include_null_range=TRUE, maxareas=maxareas)
		}
	
	# Check for ranges greater than the maximum number of areas in states_list
	ranges_greater_than_maxareas_TF = sapply(states_list, length) > maxareas
	# Reduce states_list accordingly, if any are found
	if (sum(ranges_greater_than_maxareas_TF) > 0)
		{
		states_list = states_list[ranges_greater_than_maxareas_TF]
		}
	
	# Get the letter code ranges for each binary state
	letter_code_ranges = binary_ranges_to_letter_codes(tipranges, areanames)
	names(letter_code_ranges) = NULL
	letter_code_ranges
	
	# Empty matrix of tip relative likelihoods
	tip_condlikes_of_data_on_each_state = matrix(data=0, nrow=length(phy$tip.label), ncol=length(states_list))

	# Make sure the letter code ranges are in the same order as the 
	# phylogeny tips
	tipranges_df_order = match(phy$tip.label, rownames(tipranges@df))
	letter_code_ranges = letter_code_ranges[tipranges_df_order]
	
	# Collapse the states to text fields and compare
	states_list_txt = sapply(X=states_list, FUN=paste, collapse="", sep="")
	for (rownum in 1:nrow(tip_condlikes_of_data_on_each_state))
		{
		temp_state = paste(letter_code_ranges[rownum][[1]], collapse="", sep="")
		state_match_TF = temp_state == states_list_txt
		tip_condlikes_of_data_on_each_state[rownum, state_match_TF] = 1
		}
	
	return(tip_condlikes_of_data_on_each_state)
	}





#######################################################
# areas_list_to_states_list_new
#######################################################
#' Convert a list of areas to a list of geographic ranges (states); R version
#' 
#' R version of areas_list_to_states_list_old, which makes use of \code{\link[cladoRcpp]{cladoRcpp}}'s
#' \code{\link[cladoRcpp]{rcpp_areas_list_to_states_list}}. 
#' 
#' This is the original R version of the function which converts a list of possible areas to
#' a list of all possible states (geographic ranges).  This gets slow for large numbers of areas.
#' 
#' The function is mostly replaced by \code{\link[cladoRcpp]{rcpp_areas_list_to_states_list}} in optimized code, but is still used in some places
#' for display purposes.
#' 
#' @param areas a list of areas (character or number; the function converts these to numbers, starting with 0)
#' @param maxareas maximum number of areas in this analyses
#' @param include_null_range \code{TRUE} or \code{FALSE}, should the \code{NULL} range be included in the possible states? (e.g., LAGRANGE default is yes)
#' @param split_ABC \code{TRUE} or \code{FALSE} If \code{TRUE} the output will consist of a list of lists (c("A","B","C"), c("A","B"), c("A","D"), etc.); 
#' if \code{FALSE}, the list of areas will be collapsed ("ABC", "AB", "AD", etc.).
#' @return \code{states_list} A list of the states.
#' @export
#' @seealso \code{\link{numstates_from_numareas}}, \code{\link{rcpp_areas_list_to_states_list}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @examples
#' areas = c("A","B","C")
#' areas_list_to_states_list_new(areas=areas, maxareas=length(areas), 
#' include_null_range=TRUE, split_ABC=TRUE)
#' areas_list_to_states_list_new(areas=areas, maxareas=length(areas), 
#' include_null_range=TRUE, split_ABC=FALSE)
#' areas_list_to_states_list_new(areas=areas, maxareas=length(areas), 
#' include_null_range=FALSE, split_ABC=TRUE)
#' areas_list_to_states_list_new(areas=areas, maxareas=length(areas), 
#' include_null_range=FALSE, split_ABC=FALSE)
#' areas_list_to_states_list_new(areas=areas, maxareas=2, 
#' include_null_range=TRUE, split_ABC=TRUE)
#' areas_list_to_states_list_new(areas=areas, maxareas=2, 
#' include_null_range=TRUE, split_ABC=FALSE)
#' areas_list_to_states_list_new(areas=areas, maxareas=2, 
#' include_null_range=FALSE, split_ABC=TRUE)
#' areas_list_to_states_list_new(areas=areas, maxareas=2, 
#' include_null_range=FALSE, split_ABC=FALSE)
#' areas_list_to_states_list_new(areas=areas, maxareas=1, 
#' include_null_range=TRUE, split_ABC=TRUE)
#' areas_list_to_states_list_new(areas=areas, maxareas=1, 
#' include_null_range=TRUE, split_ABC=FALSE)
#' areas_list_to_states_list_new(areas=areas, maxareas=1, 
#' include_null_range=FALSE, split_ABC=TRUE)
#' areas_list_to_states_list_new(areas=areas, maxareas=1, 
#' include_null_range=FALSE, split_ABC=FALSE)
#' 
areas_list_to_states_list_new <- function(areas=c("A","B","C"), maxareas=length(areas), include_null_range=TRUE, split_ABC=TRUE)
	{
	runjunk='
	areas = c("K", "O", "M", "H")
	maxareas=length(areas)
	include_null_range=FALSE
	split_ABC=FALSE
	'
	
	# Error trap
	if (maxareas > length(areas))
		{
		maxareas = length(areas)
		}
	
	
	# Initialize the states_list to the correct size
	nstates = numstates_from_numareas(numareas=length(areas), maxareas=maxareas, include_null_range=include_null_range)
	states_list = rep(NA, times=nstates)

	if (split_ABC == TRUE)
		{
		# Option #1: Don't split states
		# Add range combinations to the list
		states_list_area_indexes = rcpp_areas_list_to_states_list(areas=areas, include_null_range=include_null_range, maxareas=maxareas)
		
		# +1 to indexes
		states_list_area_indexes = lapply(X=states_list_area_indexes, FUN="+", 1)
		
		# convert to letters
		tmpfun <- function(x, abbr) { abbr[x] }
		states_list_areas = lapply(X=states_list_area_indexes, FUN=tmpfun, abbr=areas)
		states_list_areas
		
		# convert NA to "_"
		states_list_areas[is.na(states_list_areas)] = "_"
		return(states_list_areas)
		}


	if (split_ABC == FALSE)
		{
		# Option #1: Don't split states
		# Add range combinations to the list
		states_list_area_indexes = rcpp_areas_list_to_states_list(areas=areas, include_null_range=include_null_range, maxareas=maxareas)
		
		# +1 to indexes
		states_list_area_indexes = lapply(X=states_list_area_indexes, FUN="+", 1)
		
		# convert to letters
		tmpfun <- function(x, abbr) { abbr[x] }
		states_list_areas = lapply(X=states_list_area_indexes, FUN=tmpfun, abbr=areas)
		states_list_areas

	
		# convert NA to "_"
		states_list_areas[is.na(states_list_areas)] = "_"
		
		# Collapse the ranges
		states_list_areas = lapply(X=states_list_areas, FUN=paste, collapse="")
		states_list_areas
		
		return(states_list_areas)
		}

	return(stop("areas_list_to_states_list_new(): ERROR, you shouldn't reach this."))
	}
	




#######################################################
# Binary ranges to letter codes
#######################################################

#######################################################
# binary_range_to_letter_code_txt
#######################################################
#' Convert binary presence/absence codes (1/0) to text area names
#' 
#' Given a row of a \code{tipranges} object, converts to the corresponding name(s), collapsed
#' into a string.  E.g., if the areas were \code{(A,B,C,D)}, and the tipranges row had \code{(1 0 1 0)}, the 
#' output statename would be "AC".
#' 
#' @param tipranges_row row of a \code{tipranges} object.
#' @param areanames a list of the names of the areas
#' @return \code{statename} The corresponding name(s), collapsed into a string
#' @export
#' @seealso \code{\link{binary_range_to_letter_code_list}}, \code{\link{tipranges_to_tip_condlikes_of_data_on_each_state}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' tipranges_row = c(1, 0, 1, 0)
#' areanames = c("A", "B", "C", "D")
#' statename = binary_range_to_letter_code_txt(tipranges_row, areanames)
#' statename
#' 
binary_range_to_letter_code_txt <- function(tipranges_row, areanames)
	{
	present_TF = tipranges_row == 1
	statename = paste(areanames[present_TF], collapse="", sep="")
	
	return(statename)
	}

#######################################################
# binary_range_to_letter_code_list
#######################################################
#' Convert binary presence/absence codes (1/0) to a list of text area names
#' 
#' Given a row of a \code{tipranges} object, converts to a list of the corresponding name(s).
#' E.g., if the areas were \code{(A,B,C,D)}, and the tipranges row had \code{(1 0 1 0)}, the 
#' output statename would be ("A","C").
#' 
#' @param tipranges_row row of a \code{tipranges} object.
#' @param areanames a list of the names of the areas
#' @return \code{list_of_areas_in_the_state} A list of the name(s) of the areas corresponding to the presence/absence coding in the row
#' @export
#' @seealso \code{\link{binary_ranges_to_letter_codes}}, \code{\link{letter_string_to_binary}}, \code{\link{letter_strings_to_tipranges_df}}, \code{\link{binary_range_to_letter_code_txt}}, \code{\link{tipranges_to_tip_condlikes_of_data_on_each_state}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' tipranges_row = c(1, 0, 1, 0)
#' areanames = c("A", "B", "C", "D")
#' list_of_areas_in_the_state = binary_range_to_letter_code_list(tipranges_row, 
#' areanames)
#' list_of_areas_in_the_state
#' 
binary_range_to_letter_code_list <- function(tipranges_row, areanames)
	{
	present_TF = tipranges_row == 1
	list_of_areas_in_the_state = c(areanames[present_TF])
	
	return(list_of_areas_in_the_state)
	}


#######################################################
# binary_ranges_to_letter_codes
#######################################################
#' Convert binary presence/absence codes (1/0) to a list of text area names
#' 
#' Given a row of a \code{tipranges} object, converts to a list of the corresponding statenames
#' for each row.
#' 
#' @param tipranges a \code{tipranges} object.
#' @param areanames a list of the names of the areas
#' @return \code{letter_code_ranges} A list of the states -- there will be as many states as there are rows/tips in \code{tipranges}.
#' Each state will be a list of area names.
#' @export
#' @seealso \code{\link{binary_range_to_letter_code_list}}, \code{\link{letter_string_to_binary}}, \code{\link{letter_strings_to_tipranges_df}}, \code{\link{tipranges_to_tip_condlikes_of_data_on_each_state}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' # Define a tipranges object
#' tipranges_object = define_tipranges_object()
#' tipranges_object
#'
#' areanames = getareas_from_tipranges_object(tipranges_object)
#' areanames
#'
#' letter_code_ranges = binary_ranges_to_letter_codes(tipranges=tipranges_object, 
#' areanames)
#' letter_code_ranges
#' 
binary_ranges_to_letter_codes <- function(tipranges, areanames)
	{
	letter_code_ranges = apply(X=tipranges@df, MARGIN=1, FUN=binary_range_to_letter_code_list, areanames=areanames)
	
	return(letter_code_ranges)
	}


#######################################################
# letter_string_to_binary
#######################################################
#' Convert ranges in the form of letters (A, AB, BFG, etc.) to binary state number codes
#' 
#' This function takes a letter string (e.g. ABD) and converts to binary encoding (e.g. 1101).
#' 
#' @param letter_string A string of letters (e.g. "ABD")
#' @param letter_codes_in_desired_order The letter codes in the desired order. The default keyword, "alphabet", uses the standard 26 capital letters;
#' the output binary codes will thus have 26 positions.  If the user inputs fewer letters here, or puts them in another order, those will be used.
#' @return \code{numcodes} A list with the binary codes.
#' @export
#' @seealso \code{\link{binary_ranges_to_letter_codes}}, \code{\link{binary_range_to_letter_code_list}}, \code{\link{letter_strings_to_tipranges_df}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' letter_string = "ABD"
#' letter_string_to_binary(letter_string, letter_codes_in_desired_order="alphabet")
#' 
#' letter_string = "ABD"
#' letter_string_to_binary(letter_string, 
#' letter_codes_in_desired_order=c("A","B","C","D","E","F"))
#' 
#' letter_string = "ABD"
#' letter_string_to_binary(letter_string, 
#' letter_codes_in_desired_order=strsplit("ABCDEF", split="")[[1]])
#' 
letter_string_to_binary <- function(letter_string, letter_codes_in_desired_order="alphabet")
	{
	tests='
	letter_string = "BGD"
	alphabet = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z")
	letter_codes_in_desired_order = "alphabet"
	'
	
	if ( (length(letter_codes_in_desired_order)==1 ) && (letter_codes_in_desired_order == "alphabet"))
		{
		# Alphabet, just so we have it
		alphabet = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z")
		letter_codes_in_desired_order = alphabet
		} 
	
	nums = seq(1, length(letter_codes_in_desired_order), by=1)
	
	letter_codes = strsplit(letter_string, split="")[[1]]
	
	# Start the numcodes as all zeros
	numcodes = rep(0, times=length(letter_codes_in_desired_order))
	
	# Check for which are in the full list
	codes_found_TF = letter_codes_in_desired_order %in% letter_codes
	numcodes[codes_found_TF] = 1
	
	return(numcodes)
	}
	

#######################################################
# Convert ranges in the form of letters (A, AB, BFG, etc.) to binary state number codes
# Apply to a vector of such strings, output a tipranges object
#######################################################
#######################################################
# letter_strings_to_tipranges_df
#######################################################
#' Convert ranges in the form of letters (A, AB, BFG, etc.) to a \code{tipranges} object
#' 
#' This function converts ranges in the form of concatenated letters (A, AB, BFG, etc.) to binary state number codes. 
#' Via \code{\link{apply}}, this is done to each member of the entire input vector of strings.  It outputs \code{\link{tipranges}} object.
#' 
#' @param letter_strings A list of ranges in concatenated letter form ("A", "AB", "BFG", etc.)
#' @param letter_codes_in_desired_order The letter codes in the desired order. The default keyword, "alphabet", uses the standard 26 capital letters;
#' the output binary codes will thus have 26 positions.  If the user inputs fewer letters here, or puts them in another order, those will be used.
#' @param tipnames_in_order If given, the input tipnames will be applied as rownames in the tipranges object.  Default is \code{NULL}, which results in numbering the rows.
#' @return \code{tipranges} An object of class \code{tipranges}.
#' @export
#' @seealso \code{\link{letter_string_to_binary}}, \code{\link{binary_range_to_letter_code_list}}, \code{\link{binary_ranges_to_letter_codes}},  \code{\link{getranges_from_LagrangePHYLIP}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' letter_strings = c("A", "B", "C", "AB", "AC", "BC", "ABC")
#' letter_strings_to_tipranges_df(letter_strings)
#' 
#' letter_strings = c("A", "B", "C", "AB", "AC", "BC", "ABC")
#' letter_strings_to_tipranges_df(letter_strings, 
#' tipnames_in_order=paste("tip", seq(1,7), sep=""))
#' 
letter_strings_to_tipranges_df <- function(letter_strings, letter_codes_in_desired_order="alphabet", tipnames_in_order=NULL)
	{
	all_letters = paste(letter_strings, collapse="")
	all_letters_list = strsplit(all_letters, split="")[[1]]
	unique_letters = unique(all_letters_list)
	
	
	# Check if letter_codes should be default
	if ( (length(letter_codes_in_desired_order)==1 ) && (letter_codes_in_desired_order == "alphabet"))
		{
		# Alphabet, just so we have it
		alphabet = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z")
		letter_codes_in_desired_order = alphabet[1: length(unique_letters)]
		} 
	
	# Make the binary ranges (columns = each species, rows = presence/absence in that region)
	#binary_ranges = sapply(X=xls$tipranges, FUN=letter_string_to_binary, letter_codes_in_desired_order=Webb_Ree_2011_letter_codes)
	binary_ranges = sapply(X=letter_strings, FUN=letter_string_to_binary, letter_codes_in_desired_order=letter_codes_in_desired_order)
	
	binary_ranges = base::t(binary_ranges)
	binary_ranges_df = adf2(binary_ranges)
	
	# Apply tipnames, if available
	if (!is.null(tipnames_in_order))
		{
		rownames(binary_ranges_df) = tipnames_in_order
		}
	
	# Column names are the letter codes for each region
	colnames(binary_ranges_df) = letter_codes_in_desired_order
	
	# Use function to put in well-defined object
	tipranges = define_tipranges_object(tmpdf=binary_ranges_df)
	
	return(tipranges)
	}


#######################################################
# getranges_from_LagrangePHYLIP
#######################################################
#' Read a LAGRANGE PHYLIP-style file containing geographic ranges into a \code{tipranges} object
#' 
#' Given some geographic range data for tips in the Lagrange C++/PHYLIP format (\cite{SmithRee2010_CPPversion}), this function imports
#' the range data into a \code{tipranges}-class \code{data.frame} structure.
#' 
#' LAGRANGE C++ geographic range files are ASCII text files with the format:
#' 
#' \code{19	4 (A B C D)}\cr
#' \code{P_mariniana_Kokee2	1000}\cr
#' \code{P_mariniana_Oahu	0100}\cr
#' \code{P_mariniana_MauiNui	0010}\cr
#' \code{P_hawaiiensis_Makaopuhi	0001}\cr
#' \code{P_wawraeDL7428	1000}\cr
#' [...]\cr
#' \cr
#'
#' The first row specifies the number of taxa (here, 19), the number of areas (here, 4), and finally, the names/abbreviations of the areas.
#' The rest of the rows give the taxon names, followed by a tab and then the presence/absence in each range with 1s/0s.
#'
#' The file above is part of the geographic range data for the Hawaiian \emph{Psychotria} dataset used by \cite{ReeSmith2008}.
#'
#' @param lgdata_fn The LAGRANGE geographic data file to be read.
#' @return \code{tipranges_object} An object of class \code{tipranges}
#' @export
#' @seealso \code{\link{define_tipranges_object}}, \code{\link{save_tipranges_to_LagrangePHYLIP}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#'	 @cite SmithRee2010_CPPversion
#' @examples
#' testval=1
#' 
#' # Get the example files directory
#' extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
#' # tmp hard code: 
#' # extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"
#' # Set the filename (Hawaiian Psychotria from Ree & Smith 2008)
#' fn = np(paste(extdata_dir, "/Psychotria_geog.data", sep=""))
#' getranges_from_LagrangePHYLIP(lgdata_fn=fn)
#' 
getranges_from_LagrangePHYLIP <- function(lgdata_fn="lagrange_area_data_file.data")
	{
	# Make a tipranges instance
	#tipranges_object = define_tipranges_object()
	#tipranges_object@df = dfnums_to_numeric(tipranges_object@df)
	
	# Read the Lagrange geographic ranges data file
	tmp_blah = read_PHYLIP_data(lgdata_fn)
	tmp_input = adf2(data.matrix(tmp_blah))
	#nums_as_char = as.numeric(unlist(tmp_input))
	#tmpdf = adf2(matrix(data=nums_as_char, nrow=nrow(tmp_input), ncol=ncol(tmp_input), byrow=TRUE))
	#names(tmpdf) = names(tmp_input)
	#rownames(tmpdf) = rownames(tmp_input)
	#sum(tmpdf)
	#tmpdf = dfnums_to_numeric(tmp_input, printout=TRUE)
	#cls.df(tmpdf)
	
	# Put into a tipranges object
	tipranges_object = define_tipranges_object(tmpdf=tmp_input)
	# you can get the dataframe with
	# tipranges_object@df
	
	tipranges_object@df = adf2(data.matrix(tipranges_object@df))
	
	rownames(tipranges_object@df) = rownames(tmp_blah)
	
	return(tipranges_object)
	}



#######################################################
# order_tipranges_by_tree_tips
#######################################################
#' Reorder the rows in a \code{tipranges} object, to correspond to tree tips
#' 
#' The tipranges object, as read from a LAGRANGE/PHYLIP-style geography file, may not
#' have the species names as the same order as they are in the tips of the tree. This 
#' function allows the user to reorder them to match the tree
#'
#' @param tipranges An object of class \code{\link{tipranges}}.
#' @param tr A \code{\link[ape]{phylo}} tree object.
#' @return \code{tipranges} An object of class \code{tipranges}
#' @export
#' @seealso \code{\link{tipranges_to_area_strings}}, \code{\link{define_tipranges_object}}, \code{\link{save_tipranges_to_LagrangePHYLIP}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#'	 @cite SmithRee2010_CPPversion
#' @examples
#' testval=1
#' 
#' # Get the example files directory
#' extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
#' # tmp hard code: 
#' # extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"
#' # Set the filename (Hawaiian Psychotria from Ree & Smith 2008)
#' 
#' trfn = np(paste(extdata_dir, "/Psychotria_5.2.newick", sep=""))
#' tr = read.tree(trfn)
#' 
#' fn = np(paste(extdata_dir, "/Psychotria_geog.data", sep=""))
#' tipranges1 = getranges_from_LagrangePHYLIP(lgdata_fn=fn)
#' tipranges1
#' 
#' # Reorder the tipranges object
#' tipranges2 = order_tipranges_by_tree_tips(tipranges1, tr)
#' tipranges2
#' 
order_tipranges_by_tree_tips <- function(tipranges, tr)
	{
	# Get the order of the tips
	matchnums = get_indices_where_list1_occurs_in_list2(list1=tr$tip.label, list2=rownames(tipranges@df))
	
	# Reorder the tipranges object
	tipranges@df = tipranges@df[matchnums,]
	
	return(tipranges)
	}




#######################################################
# tipranges_to_area_strings
#######################################################
#' Convert tipranges binary coding to range strings
#' 
#' This function converts the 0110-type format of the tipranges object into a list of strings
#' describing the geographic ranges.  E.g., 1100 becomes AB, 0111 become BCD (assuming the regions
#' are abbreviated A, B, C...).  Users can input their preferred abbreviations with \code{areaabbr}.
#' 
#' Note that you will HAVE to use \code{\link{order_tipranges_by_tree_tips}} on the tipranges object first,
#' to make sure the tipranges are in the correct order on the tree tips.
#'
#' @param tipranges An object of class \code{\link{tipranges}}.
#' @param areaabbr A vector of the abbreviations (preferably 1 character each).
#' @return \code{tiprange_names} A vector of strings.
#' @export
#' @seealso \code{\link{order_tipranges_by_tree_tips}}, \code{\link{define_tipranges_object}}, \code{\link{save_tipranges_to_LagrangePHYLIP}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#'	 @cite SmithRee2010_CPPversion
#' @examples
#' 
#' # Get the example files directory
#' extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
#' # tmp hard code: 
#' # extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"
#' # Set the filename (Hawaiian Psychotria from Ree & Smith 2008)
#' 
#' trfn = np(paste(extdata_dir, "/Psychotria_5.2.newick", sep=""))
#' tr = read.tree(trfn)
#' 
#' fn = np(paste(extdata_dir, "/Psychotria_geog.data", sep=""))
#' tipranges1 = getranges_from_LagrangePHYLIP(lgdata_fn=fn)
#' tipranges1
#' tipranges_to_area_strings(tipranges=tipranges1, areaabbr=NULL)
#' tipranges_to_area_strings(tipranges=tipranges1, areaabbr=c("K", "O", "M", "H"))
#' 
#' # Reorder the tipranges object
#' tipranges2 = order_tipranges_by_tree_tips(tipranges1, tr)
#' tipranges2
#' tipranges_to_area_strings(tipranges=tipranges2, areaabbr=NULL)
#' tipranges_to_area_strings(tipranges=tipranges2, areaabbr=c("K", "O", "M", "H"))
#' 
tipranges_to_area_strings <- function(tipranges, areaabbr=NULL)
	{
	if (is.null(areaabbr))
		{
		areaabbr = names(tipranges@df)
		} else {
		names(tipranges@df) = areaabbr
		}
	
	tiprangesTF = tipranges@df == 1
	tiprange_names = apply(X=tiprangesTF, 1, getname, tiparea_names=areaabbr)
	
	
	return(tiprange_names)
	}



#######################################################
# getname
#######################################################
#' Collapse range abbreviations to strings
#' 
#' This is a utility function used by apply in \code{\link{tipranges_to_area_strings}}. It extracts
#' the present areas and concatenates the abbreviations for one row.
#'
#' @param TFrow A list of TRUE and FALSE
#' @param tiparea_names The names of each area
#' @param concat If TRUE (default), merge the areas in a state into a single string.
#' @param sep The sep argument for \code{\link[base]{paste}}.
#' @return \code{tiparea} A string.
#' @export
#' @seealso \code{\link{states_list_indexes_to_areastxt}}, \code{\link{order_tipranges_by_tree_tips}}, \code{\link{define_tipranges_object}}, \code{\link{save_tipranges_to_LagrangePHYLIP}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' getname(TFrow=c(FALSE, TRUE, TRUE, FALSE), 
#' tiparea_names=c("K", "O", "M", "H"), sep="")
#' getname(TFrow=c(FALSE, TRUE, TRUE, FALSE), 
#' tiparea_names=c("K", "O", "M", "H"), sep="_")
#' 
getname <- function(TFrow, tiparea_names, concat=TRUE, sep="")
	{
	words = tiparea_names[TFrow]
	
	if (concat == TRUE)
		{
		tiparea = paste(words, sep=sep, collapse="")
		} else {
		tiparea = paste(words, sep=sep)
		}
	return(tiparea)
	}


#######################################################
# states_list_indexes_to_areastxt
#######################################################
#' States (ranges) lists to txt string of the areas
#' 
#' This is a utility function.
#'
#' @param states_list A list of states, where each state consists of a list of areas.
#' @param areanames A list of areanamess.
#' @param counting_base Does states_list start indexing areas from 0 (default) or 1?
#' @param concat If TRUE (default), merge the areas in a state into a single string.
#' @param sep Character to merge on, as in \code{\link[base]{paste}}. Default "".
#' @return \code{tiparea} A string.
#' @export
#' @seealso \code{\link{getname}}, \code{\link{order_tipranges_by_tree_tips}}, \code{\link{define_tipranges_object}}, \code{\link{save_tipranges_to_LagrangePHYLIP}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
#' 
states_list_indexes_to_areastxt <- function(states_list, areanames, counting_base=0, concat=TRUE, sep="")
	{
	defaults='
	areas = c("K", "O", "M", "H")
	states2 = rcpp_areas_list_to_states_list(areas=areas)
	areanames = c("K", "O", "M", "H")
	counting_base=0
	states_list_indexes_to_areanames(states_list=states2, areanames, counting_base=0)
	'
	# Increment by 1 if your state indexes are 0-based
	if (counting_base == 0)
		{
		states_list = sapply(X=states_list, FUN="+", 1)
		}
	
	# sapply getname through the list of lists of indexes
	areastxt = sapply(X=states_list, FUN=getname, tiparea_names=areanames, concat=concat, sep=sep)
	
	return(areastxt)
	}


#######################################################
# save_tipranges_to_LagrangePHYLIP
#######################################################
#' Save a tipranges object to a LAGRANGE PHYLIP-style file containing binary-encoded geographic ranges
#' 
#' Given some geographic range data for tips in the \code{tipranges} object, this function
#' exports them to an ASCII text file in the Lagrange C++/PHYLIP format (\cite{SmithRee2010_CPPversion}).  This file can then be
#' read by \code{\link{getranges_from_LagrangePHYLIP}}.
#' 
#' LAGRANGE C++ geographic range files are ASCII text files with the format:
#' 
#' \code{19	4 (A B C D)}\cr
#' \code{P_mariniana_Kokee2	1000}\cr
#' \code{P_mariniana_Oahu	0100}\cr
#' \code{P_mariniana_MauiNui	0010}\cr
#' \code{P_hawaiiensis_Makaopuhi	0001}\cr
#' \code{P_wawraeDL7428	1000}\cr
#' [...]\cr
#' \cr
#'
#' The first row specifies the number of taxa (here, 19), the number of areas (here, 4), and finally, the names/abbreviations of the areas.
#' The rest of the rows give the taxon names, followed by a tab and then the presence/absence in each range with 1s/0s.
#'
#' The file above is part of the geographic range data for the Hawaiian \emph{Psychotria} dataset used by \cite{ReeSmith2008}.
#' @param tipranges_object An object of class \code{tipranges}.
#' @param lgdata_fn The LAGRANGE geographic data file to be output.
#' @param areanames A list of the names of the areas.
#' @return \code{tipranges_object} An object of class \code{tipranges}
#' @export
#' @seealso \code{\link{define_tipranges_object}}, \code{\link{getranges_from_LagrangePHYLIP}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#'	 @cite SmithRee2010_CPPversion
#' @examples
#' testval=1
#' 
#' # Create an example tipranges object
#' tipranges = define_tipranges_object()
#' 
#' # See current directory
#' getwd()
#' 
#' \dontrun{
#' # Save the file
#' # Set the filename
#' fn = "example_tipranges.data"
#' save_tipranges_to_LagrangePHYLIP(tipranges_object=tipranges, lgdata_fn=fn)
#'
#' # Show the file
#' tmplines = scan(file=fn, what="character", sep="\n")
#' cat(tmplines, sep="\n")
#' 
#' # Again, with areanames
#' save_tipranges_to_LagrangePHYLIP(tipranges_object=tipranges, 
#' lgdata_fn=fn, areanames=c("area1","area2","area3"))
#'
#' # Show the file
#' tmplines = scan(file=fn, what="character", sep="\n")
#' cat(tmplines, sep="\n")
#' } # End dontrun
#' 
save_tipranges_to_LagrangePHYLIP <- function(tipranges_object, lgdata_fn="lagrange_area_data_file.data", areanames=colnames(tipranges_object@df))
	{
	# Extract the tipranges data
	tipranges_df = tipranges_object@df

	# If there are no area names, or if they are the wrong length, assign A B C D E etc...
	if ( (areanames == "replace") || (length(areanames) != ncol(tipranges_df)))
		{
		new_areanames = LETTERS[1:ncol(tipranges_df)]
		areanames_txt = paste(new_areanames, sep=" ")

		cat("\nNote: assigning '", areanames_txt, "' as area names.\n", sep="")

		} else {
		areanames_txt = paste(areanames, collapse=" ")
		}


	
	# Collapse ranges to strings
	ranges_strings = apply(X=tipranges_df, MARGIN=1, FUN=paste, collapse="")
	
	# Get taxon names
	taxon_names = row.names(tipranges_df)
	
	ranges_table = cbind(taxon_names, ranges_strings)
	ranges_table_txt = apply(X=ranges_table, MARGIN=1, FUN=paste, collapse="\t")
	
	header_string = paste(nrow(tipranges_df), "\t", ncol(tipranges_df), "\t(", areanames_txt, ")", sep="")
	
	# Write the header string
	#write_lines_good(dtf=header_string, outfn=lgdata_fn, sepval="\n", tmpappend=FALSE)
	write.table(x=header_string, file=lgdata_fn, append=FALSE, quote=FALSE, sep="\n", row.names=FALSE, col.names=FALSE)

	# Write the rest
	#write_lines_good(dtf=ranges_table_txt, outfn=lgdata_fn, sepval="\n", tmpappend=TRUE)
	write.table(x=ranges_table_txt, file=lgdata_fn, append=TRUE, quote=FALSE, sep="\n", row.names=FALSE, col.names=FALSE)
	
	return(lgdata_fn)
	}



#######################################
# Read PHYLIP data
# This assumes data are interleaved, and that names are separated from data
# by a tab character; none of the 10-character limit stuff.
#
# Also assumes no ambiguity in characters
#######################################
#######################################################
# read_PHYLIP_data
#######################################################
#' Read a PHYLIP-format file
#' 
#' This assumes data are interleaved, and that names are separated from data
#' by a tab character; there is no 10-character limit on names.
#'
#' This function is a precursor to \code{\link{getranges_from_LagrangePHYLIP}}.
#' 
#' @param lgdata_fn The filename to read.
#' @param regionnames A list of the names of the areas. Only used if the names are NOT specified in the file.
#' @return \code{tmpdf} A \code{\link[base]{data.frame}} containing the data.
#' @export
#' @seealso \code{\link{getranges_from_LagrangePHYLIP}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' 
#' # Get the example files directory
#' extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
#' # tmp hard code: extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"
#' # Set the filename (Hawaiian Psychotria from Ree & Smith 2008)
#' fn = np(paste(extdata_dir, "/Psychotria_geog.data", sep=""))
#'
#' # Read in the file
#' tmpdf = read_PHYLIP_data(lgdata_fn=fn, regionnames=NULL)
#' tmpdf
#' 
#' # Read in the file
#' tmpdf = read_PHYLIP_data(lgdata_fn=fn, 
#' regionnames=c("Kauai", "Oahu", "Maui-Nui","Big Island"))
#' tmpdf	# Note that regionnames are only 
#' # used if they are NOT specified in the file.
#' # But, you could put them on manually
#' names(tmpdf) = c("Kauai", "Oahu", "Maui-Nui","Big Island")
#' tmpdf
#' 
#' # This one has no area names
#' fn = np(paste(extdata_dir, "/Psychotria_geog_noAreaNames.data", sep=""))
#' tmpdf = read_PHYLIP_data(lgdata_fn=fn, 
#' regionnames=c("Kauai", "Oahu", "Maui-Nui","Big Island"))
#' tmpdf	# Note that regionnames are only 
#' # used if they are NOT specified in the file.
#' 
read_PHYLIP_data <- function(lgdata_fn="lagrange_area_data_file.data", regionnames=NULL)
	{
	setup='
	lgdata_fn = "/Users/nickm/Desktop/__projects/_2011-07-15_Hannah_spider_fossils/_data/lagrange_for_nick2/palp_no_Lacun.data"
	'

	# Read the 1st line, split on whitespaces
	firstline = scan(lgdata_fn, what="character", nlines=1)
	
	# Parse the firstline
	ntips = as.numeric(firstline[1])
	nareas = as.numeric(firstline[2])
	
	#######################################################
	# If the length of the first line is > 2, parse that information 
	#######################################################
	if (length(firstline) > 2)
		{
		# Get region names	
		regionnames = firstline[3:length(firstline)]
		
		# Remove "(" and ")"
		regionnames = mapply(sub, pattern="\\(", replacement="", x=regionnames)
		regionnames = mapply(sub, pattern="\\)", replacement="", x=regionnames)
		names(regionnames) = NULL
		}
	
	# Make the data.frame
	tmpdf = matrix(data=NA, nrow=ntips, ncol=nareas)
	
	# Parse the remaining rows
	tmplines = scan(lgdata_fn, what="character", sep="\t", skip=1)
	tmplines = matrix(data=tmplines, ncol=2, byrow=TRUE)
	tmplines

	tipnames = tmplines[,1]
	areas_char = tmplines[,2]
	
	areas_char2 = unlist(mapply(strsplit, x=areas_char, split=""))
	areas_char3 = matrix(data=areas_char2, nrow=ntips, byrow=TRUE)


	# Store in a data.frame and return
	tmpdf = areas_char3
	tmpdf = adf(tmpdf)
	
	if ( is.null(regionnames) )
		{
		alphabet = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z")
		regionnames	= alphabet[1:ncol(tmpdf)]
		colnames(tmpdf) = regionnames
		} else {
		names(tmpdf) = regionnames
		}
	row.names(tmpdf) = tipnames
	
	return(tmpdf)	
	}







#######################################################
# Biogeographic transition matrices
#######################################################


# Get the maximum rangesize for a given ancestral rangesize
# Now, go through, and make a list of the max minsize for each decsize
# e.g.
# 
# 	# Set the parameter controlling the size distribution of 
# 	# the smaller descendant species
# 	maxent01s_param = maxent_constraint_01
# 	maxent01v_param = maxent_constraint_01v
# 	maxent01j_param = maxent_constraint_01
# 	maxent01y_param = maxent_constraint_01
# 
# 	# Cladogenesis model inputs
# 	spPmat_inputs = NULL
# 	states_indices = states_list
# 	states_indices[1] = NULL	# shorten the states_indices by 1 (cutting the null range state from the speciation matrix)
# 	l = states_indices
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
# 	maxent01s = relative_probabilities_of_subsets(max_numareas=numareas, maxent_constraint_01=maxent01s_param, NA_val=0)
# 	maxent01v = relative_probabilities_of_vicariants(max_numareas=numareas, maxent_constraint_01=maxent01v_param, NA_val=0)
# 	maxent01j = relative_probabilities_of_subsets(max_numareas=numareas, maxent_constraint_01=maxent01j_param, NA_val=0)
# 	maxent01y = relative_probabilities_of_subsets(max_numareas=numareas, maxent_constraint_01=maxent01y_param, NA_val=0)
# 
# 	
# 	# Matrix of probs for each ancsize
# 	maxprob_as_function_of_ancsize_and_decsize = mapply(FUN=max, maxent01s, maxent01v, maxent01j, maxent01y, MoreArgs=list(na.rm=TRUE))
# 	maxprob_as_function_of_ancsize_and_decsize = matrix(data=maxprob_as_function_of_ancsize_and_decsize, nrow=nrow(maxent01s), ncol=ncol(maxent01s))
# 	maxprob_as_function_of_ancsize_and_decsize[maxprob_as_function_of_ancsize_and_decsize > 0] = 1
# 	maxprob_as_function_of_ancsize_and_decsize[maxprob_as_function_of_ancsize_and_decsize <= 0] = 0
# 	
# 	# Now, go through, and make a list of the max minsize for each decsize
# 	maxsize <- function (areasizes_possible_01)
# 		{
# 		max(1:sum(areasizes_possible_01 > 0, na.rm=TRUE))
# 		}
# 	# Now, go through, and make a list of the max minsize for each decsize
# 	max_minsize_as_function_of_ancsize = apply(X=maxprob_as_function_of_ancsize_and_decsize, MARGIN=1, FUN=maxsize)
# 	max_minsize_as_function_of_ancsize

#######################################################
# maxsize
#######################################################
#' Get the maximum rangesize for a given ancestral rangesize
#' 
#' This function returns the maximum descendant rangesize for a given ancestral rangesize, given a list of
#' 0/1 values specifying the possibility of each descendant rangesizes.
#' 
#' This is mostly a utility function used within \code{\link[base]{apply}} within other functions.
#' 
#' @param areasizes_possible_01 A list of 0/1 values, indicating whether an range of that size (rangesize = 1-based index = 1, 2, 3...) is possible (1) or not (0).
#' @return \code{max_number_of_areas} The maximum number of areas 
#' @export
#' @seealso \code{\link[base]{apply}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' areasizes_possible_01 = c(1,1,1,0,0)
#' maxsize(areasizes_possible_01)
#' 
maxsize <- function (areasizes_possible_01)
	{
	indexes_possible = 1:length(areasizes_possible_01)
	
	areasizes_possible_01_v2 = areasizes_possible_01 * indexes_possible
	max_number_of_areas = max(areasizes_possible_01_v2, na.rm=TRUE)
	
	return(max_number_of_areas)
	}




#######################################################
# default_states_list
#######################################################
#' Default input for a states_list
#' 
#' R CMD check limits the length of inputs to variables for functions; this is a 
#' workaround.
#' 
#' @return \code{states_list} The list of states
#' @export
#' @seealso \code{\link{make_dispersal_multiplier_matrix}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' states_list = default_states_list()
#' 
default_states_list <- function ()
	{
	states_list=list("_", c("A"), c("B"), c("C"), c("A","B"), c("B","C"), c("A","C"), c("A","B","C"))
	
	return(states_list)
	}




#######################################################
# make_dispersal_multiplier_matrix
#######################################################
#' Make a default matrix of relative dispersal probabilities between areas
#' 
#' Given either a list of areas, or a list of states, this function provides a square
#' dispersal matrix giving the relative probability of dispersal between areas.  The function
#' fills in these dispersals probabilities with the value 1.  The user can then modify this
#' as desired.
#' dispersal_multipliers_matrix Default NULL
#' distances_mat Default NULL
#' x_exponent Default 0
#' 
#' If only a states list is given, the list of areas is calculated by getting \code{\link[base]{unique}} values from
#' the concatenated states list.
#' 
#' @param areas A list of areas; if \code{NULL}, the states list will be used.
#' @param states_list A list of states, where each state consists of a list of areas. A default example list is provided.
#' @param dispersal_multipliers_matrix Default NULL.
#' @param distances_mat Default NULL.
#' @param x_exponent Default 0.
#' @return \code{dispersal_multiplier_matrix} A square matrix, with 1s for all cells.
#' @export
#' @seealso \code{\link{make_relprob_matrix_de}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite FosterIdiots
#' @examples
#' testval=1
#' make_dispersal_multiplier_matrix(areas=NULL, 
#' states_list=list("_", c("A"), c("B"), c("C"), 
#' c("A","B"), c("B","C"), c("A","C"), c("A","B","C")))
#' make_dispersal_multiplier_matrix(areas=c("A","B","C","D"))
#' 
make_dispersal_multiplier_matrix <- function(areas=NULL, states_list=default_states_list(), dispersal_multipliers_matrix=NULL, distances_mat=NULL, x_exponent=0)
	{
	defaults='
	areas=NULL
	states_list=list("_", c("A"), c("B"), c("C"), c("A","B"), c("B","C"), c("A","C"), c("A","B","C"))
	make_dispersal_multiplier_matrix()
	'
	
	if (is.null(areas))
		{
		if (is.null(states_list))
			{
			stop("ERROR: for make_dispersal_multiplier_matrix() you must specify either an areas list or a states_list")
			}
		
		tmpareas = unlist(states_list)
		areas = unique(tmpareas)
		# Remove null range
		areas = areas[areas != "_"]
		}
	
	# Make a dispersal matrix giving the relative
	# probabilities of dispersal between
	# each area (areas, not ranges)
	if (is.null(dispersal_multipliers_matrix))
		{
		numareas = length(areas)
		dispersal_multipliers_matrix = as.data.frame(matrix(data=rep(1, numareas*numareas), nrow=numareas, ncol=numareas, byrow=TRUE))
	
		names(dispersal_multipliers_matrix) = areas
		rownames(dispersal_multipliers_matrix) = areas
		}
	
	# If the distances matrix is NOT null, multiply it by the distance^-x_exponent matrix
	if (!is.null(distances_mat))
		{
		# By default, x_exponent is 0, i.e., no effect
		prob_dispersal_by_distance = distances_mat ^ (1 * x_exponent)
		dispersal_multipliers_matrix = dispersal_multipliers_matrix * prob_dispersal_by_distance
		}
	
	return(dispersal_multipliers_matrix)
	}


#######################################################
# make_relprob_matrix_de
#######################################################
#' Make a relative dispersal probability matrix (in text form)
#' 
#' This function takes a list of states/geographic ranges, and makes a relative probability matrix
#' describing the probability of transition between each state.  These probabilities are described 
#' in terms of d, "dispersal" (actually range expansion)
#' and "extinction" (actually local extirpation, or range contraction), as done in the program 
#' \code{LAGRANGE} (\cite{ReeSmith2008}, \cite{SmithRee2010_CPPversion}).
#'
#' The output \code{\link[base]{data.frame}}, termed \code{dedf} (dedf=dispersal-extinction data.frame), contains the 
#' actual text of the formulas by which the transition probability matrix would be calculated.  E.g., the example calculates
#' the matrix corresponding to Equation 1 on p. 6 of Ree & Smith (2008).
#'
#' Note that the geographic range-change process described here is a continuous-time process, where the probability of change
#' is a function of branch length, and all transitions occur because of dispersal and extinction.  LAGRANGE also implements a 
#' cladogenesis model (thus DEC -- dispersal-extinction-cladogenesis) which describes an "instantaneous" process of geographic
#' range change at speciation/lineage-splitting events. \code{BioGeoBEARS} allows users to turn on, turn off, or otherwise customize
#' both the continuous-time model and the cladogenesis model.
#'
#' @param states_list A list of states, where each state consists of a list of areas. A default example list is provided.
#' @param split_ABC \code{TRUE} or \code{FALSE} If \code{TRUE} then each state/range in the input geographic ranges (\code{states_list}) 
#' will be split on the argument contained in \code{split}.
#' @param split The character to split on.
#' @param remove_simultaneous_events If \code{TRUE} (default, as in \code{LAGRANGE} and almost all phylogenetic Markov models), then it is assumed that 
#' all changes in geographic range along branches must happen one event at a time.  If \code{FALSE}, simultaneous events are not excluded;
#' this is not recommended. However, notably, a commonly-used biogeographic model (treating biogeography as a multistate discrete character in an ML
#' framework, where every species/lineage inhabits one and only one area at any point in time) effectively is invoking a simultaneous event: e.g., A->B is a 
#' simultaneous range gain and range loss, from the perspective of the dispersal-extinction framework.
#' @param add_multiple_Ds If \code{TRUE} (default, as in \code{LAGRANGE}), the probabilities of dispersal from each possible source area are added together.
#' @param dispersal_multiplier_matrix A user-provided dispersal multiplier matrix; the default is a matrix of 1s from \code{\link{make_dispersal_multiplier_matrix}}(states_list=states_list).
#' @return \code{dedf} The output \code{\link[base]{data.frame}}, termed \code{dedf} (dedf=dispersal-extinction data.frame), contains the 
#' actual text of the formulas by which the transition probability matrix would be calculated.
#' @export
#' @seealso \code{\link{make_dispersal_multiplier_matrix}}
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
#' testval=1
#' 
#' states_list = list("_", c("A"), c("B"), c("C"), c
#' ("A","B"), c("B","C"), c("A","C"), c("A","B","C"))
#' 
#' states_list = areas_list_to_states_list_new(
#' areas=c("A","B","C"), include_null_range=TRUE, split_ABC=TRUE)
#' states_list
#' 
#' dedf = make_relprob_matrix_de(states_list=states_list, 
#' split_ABC=FALSE, split="", remove_simultaneous_events=TRUE, 
#' add_multiple_Ds=TRUE, 
#' dispersal_multiplier_matrix=make_dispersal_multiplier_matrix(states_list=states_list))
#' 
#' dedf
#' 
make_relprob_matrix_de <- function(states_list=default_states_list(), split_ABC=FALSE, split="", remove_simultaneous_events=TRUE, add_multiple_Ds=TRUE, 
dispersal_multiplier_matrix=make_dispersal_multiplier_matrix(states_list))
	{
	defaults='
	states_list=list("_", c("A"), c("B"), c("C"), c("A","B"), c("B","C"), c("A","C"), c("A","B","C"))
	split_ABC=FALSE
	split=""
	add_multiple_Ds=TRUE
	remove_simultaneous_events=TRUE
	# Default is a matrix of 1s
	dispersal_multiplier_matrix=make_dispersal_multiplier_matrix()
	'
	
	
	# If split_ABC=TRUE, assume areas can be split on ""
	if (split_ABC == TRUE)
		{
		states_list = mapply(FUN=strsplit2, split=split, states_list)
		names(states_list) = NULL		
		}	
	
	# Make an empty matrix
	tmpmat = matrix(0, nrow=length(states_list), ncol=length(states_list))
	
	# Go through the rows and columns
	for (i in 1:length(states_list))
		{
		ancstates = states_list[[i]]
		for (j in 1:length(states_list))
			{
			decstates = states_list[[j]]
			
			# If the ancestor is null, the descendant prob is 0
			if ( (length(ancstates)==1) && (ancstates == "_") )
				{
				tmpmat[i,j] = 0
				next()
				}
			
			# If the states are all the same, "-"
			if ( (length(ancstates)==length(decstates)) && (all(decstates == ancstates)==TRUE) )
				{
				tmpmat[i,j] = "-"
				next()
				}
			
			# If the descendant state is null (global extinction)
			if ( (length(decstates)==1) && (decstates == "_") )
				{
				# This many local extinction events (e) are required
				events = rep("e", times=length(ancstates))
				tmpmat[i,j] = paste(events, collapse="")
				next()
				}
			
			
			# If none of the above apply (ending the loop with next()) then
			# put in ancestral & descendant states
			tmplist = NULL
			
			# If the ancestral state is missing from descendants -- local extinction
			
			# Make sure that ALL the decstates are within the ancstates
			# No new states when an extinction event is being considered!
			if ( all(decstates %in% ancstates) == TRUE)
				{
				for (m in ancstates)
					{
					if ( (m %in% decstates) == FALSE)
						{
						tmplist = c(tmplist, "e")
						}
					}
				tmpout = paste(tmplist, collapse="")
				tmpmat[i,j] = tmpout
				next()
				}

			# If the descendant state is missing from ancestors -- local dispersal
			# THE DESCENDANT MUST BE LONGER THAN THE ANCESTOR, 
			# AND (FOR NOW) ONLY ONE AREA CAN BE ADDED
			# SO ALL AREAS IN ANCESTOR MUST BE IN THE DESCENDANT
			max_rangesize_increase = 1
			if ( (all(decstates %in% ancstates) == FALSE) && (all(ancstates %in% decstates)) && (length(decstates) == (length(ancstates)+max_rangesize_increase)) )
				{
				tmplist = NULL
				for (m in decstates)
					{
					if ( (m %in% ancstates) == FALSE)
						{
						if (add_multiple_Ds == TRUE)
							{
							# Go through each ancestral state, multiply by the relevant dispersal probability
							
							these_dispersal_probs = dispersal_multiplier_matrix[ancstates,m]
							ds_with_these_dispersal_probs = sapply(X=these_dispersal_probs, FUN=paste, "*d", sep="")
							#ds_with_these_dispersal_probs
							
							ds_to_add = paste(ds_with_these_dispersal_probs, collapse="+", sep="")			
							tmplist = c(tmplist, ds_to_add)
							} else {
							# (Note: A distance matrix has little meaning in this context)
							tmplist = c(tmplist, "d")
							}
						}
					}
				tmpout = paste(tmplist, collapse="+")
				tmpmat[i,j] = tmpout
				}
			} # end j, columns, descendant states
		} # end i, rows, ancestral states
	
	
	# Remove multiple simultaneous events if desired (default: TRUE)
	if (remove_simultaneous_events == TRUE)
		{	
		for (i in 1:nrow(tmpmat))
			{
			for (j in 1:ncol(tmpmat))
				{
				# Remove e.g. "dd"
				word = tmpmat[i,j]
				if (grepl(pattern="dd", x=word) == TRUE)
					{
					tmpmat[i,j] = 0
					}
				if (grepl(pattern="ee", x=word) == TRUE)
					{
					tmpmat[i,j] = 0
					}
				}
			}

		#character_counts = nchar(tmpmat)
		#character_counts
		#tmpmat[character_counts > 1] = 0
		}
	
	
	# Show null range as "()"
	tmpnames = unlist(lapply(X=states_list, FUN=paste, collapse=""))
	tmpnames[tmpnames == ""] = "()"
	
	
	# Dispersal-extinction matrix as data.frame
	dedf = adf(tmpmat)
	names(dedf) = tmpnames
	row.names(dedf) = tmpnames
	dedf
	
	return(dedf)
	}



#######################################################
# size_species_matrix
#######################################################
#' Calculate the dimensions of the cladogenesis/speciation matrix
#' 
#' This function calculates the dimensions of the cladogenesis/speciation matrix describing the 
#' transition probabilities between ancestral geographic ranges and descendant geographic range pairs
#' on Left (L) and Right (R) branches.
#' 
#' Under a cladogenesis model of geographic range change, the model will
#' give the conditional probability of each possible combination of geographic ranges on the 
#' Left (L) and Right (R) descendant branches, conditional on a particular ancestral state.  A matrix
#' representing these transitions will have \code{numstates} ancestral states, and \code{numstates*numstates} possible 
#' descendant pairs.  Many of these will have 0 conditional probability under the model, but,
#' for visualization or experimental purposes it can be useful to display them all.
#' 
#' However, because \code{numstates = 2^numareas} under default conditions, and the number of cells the processor has to 
#' consider (without optimization tricks) is \code{numstates^3}, this transition matrix can very quickly become cumbersome
#' to explicitly calculate or display.  \code{size_species_matrix} allows the user to check this ahead of time.
#' 
#' See \code{\link[cladoRcpp]{numstates_from_numareas}} for the details of calculating \code{numstates}.
#'
#' At various points in \code{BioGeoBEARS} code, the text and numeric versions of the cladogenesis matrix are named \code{spmat} and \code{spPmat}, respectively.
#'
#' @param states_list A list of states, where each state consists of a list of areas. A default example list is provided.
#' @param printwarn If \code{printwarn>0} (\code{printwarn=1} by default), then print to screen a message describing the size of the cladogenesis matrix.
#' @return \code{spmat_dimensions} The dimensions of the cladogenesis matrix.
#' @export
#' @seealso \code{\link{make_relprob_matrix_de}}, \code{\link{make_spmat_row}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @examples
#' testval=1
#' spmat_dimensions = size_species_matrix(
#' states_list=list("_", c("A"), c("B"), c("C"), c("A","B"), 
#' c("B","C"), c("A","C"), c("A","B","C")), printwarn=1)
#' spmat_dimensions
#'
size_species_matrix <- function(states_list=default_states_list(), printwarn=1)
	{
	states_list = states_list[states_list != "_"]
	areas = unique(unlist(states_list))
	
	num_nonNull_states = length(states_list)
	
	if (printwarn > 0)
		{
		cat("\nNOTE: Your states_list has ", length(areas), " geographic areas (", areas, "), and ", num_nonNull_states, " states.\n")
		cat("        This means that at every speciation event in your tree, you have\n")
		cat("        ", num_nonNull_states, "x", num_nonNull_states, "=", num_nonNull_states*num_nonNull_states, " possible combinations of states in the two branches\n")
		cat("        just above every node (speciation event) in the tree.  This means you have to calculate\n")
		cat("        ", num_nonNull_states^3, " cells in the relative probability matrix describing the various possible scenarios.\n")
		}
	
	spmat_dimensions = c(num_nonNull_states, num_nonNull_states^2)
	
	return(spmat_dimensions)
	}







#######################################################
# expand.grid.alt
#######################################################
#' A faster version of expand.grid
#' 
#' This should be faster than \code{\link[base]{expand.grid}}, which "[c]reate[s] a data frame from all combinations of the supplied vectors or factors" (R documentation).  
#' 
#' The source of this function was this discussion thread: \url{http://stackoverflow.com/questions/10405637/use-outer-instead-of-expand-grid}
#' 
#' @param seq1 A sequence of elements
#' @param seq2 A sequence of elements
#' @return \code{matrix_of_combinations} A matrix of all the possible combinations.
#' @export
#' @seealso \code{\link[stats]{convolve}}, \code{\link{expand.grid}}, \code{\link{expand.grid.jc}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' seq1 = c("A","B","C")
#' seq2 = seq1
#' expand.grid(seq1,seq2)
#' expand.grid.alt(seq1,seq2)
#' expand.grid.jc(seq1,seq2)
#' 
expand.grid.alt <- function(seq1,seq2)
	{
	matrix_of_combinations = cbind(rep.int(seq1, length(seq2)), c(base::t(matrix(rep.int(seq2, length(seq1)), nrow=length(seq2)))))
	return(matrix_of_combinations)
	}

#######################################################
# expand.grid.jc
#######################################################
#' An even faster version of expand.grid
#' 
#' This should be faster than \code{\link[base]{expand.grid}}, which "[c]reate[s] a data frame from all combinations of the supplied vectors or factors" (R documentation).  
#' 
#' The source of this function was this discussion thread: \url{http://stackoverflow.com/questions/10405637/use-outer-instead-of-expand-grid}
#' 
#' @param seq1 A sequence of elements
#' @param seq2 A sequence of elements
#' @return \code{matrix_of_combinations} A matrix of all the possible combinations.
#' @export
#' @seealso \code{\link[stats]{convolve}}, \code{\link{expand.grid}}, \code{\link{expand.grid.jc}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' seq1 = c("A","B","C")
#' seq2 = seq1
#' expand.grid(seq1,seq2)
#' expand.grid.alt(seq1,seq2)
#' expand.grid.jc(seq1,seq2)
#' 
expand.grid.jc <- function(seq1,seq2)
	{
    matrix_of_combinations = cbind(Var1 = rep.int(seq1, length(seq2)), Var2 = rep.int(seq2, rep.int(length(seq1),length(seq2))))
	return(matrix_of_combinations)
	}







#######################################################
# make_spmat_row
#######################################################
#' Construct a (text) cell of the cladogenesis/speciation matrix
#' 
#' Given the identity of the states/geographic ranges on the left branch (\code{Lstates}), right branch (\code{Rstates}), and ancestral areas (\code{ancareas_txt_tmp}),
#' construct the (text version) of the row of transition probabilities.  This means that each nonzero cell gets a \emph{v} for a vicariance event, a \emph{y}
#' for a sympatric speciation/range-copying event, a \emph{j} for a founder-event/jump speciation event, and an \emph{s} for a sympatric-subset event.
#'
#' This function is utilized by \code{\link[base]{apply}} in other functions (e.g. ) in an attempt to speed up calculation over rows.  However, processing of
#' text formulas via \code{\link[base]{apply}} will
#' never be fast enough for large matrices; see \code{\link[cladoRcpp]{cladoRcpp}} for optimized functions.
#' 
#' This text-based matrix later gets evaluated by other functions to calculate the numerical probabilities.  I.e., if j=0 and the other forms of speciation have weights 
#' equal to each other, this is the \code{LAGRANGE} cladogenesis model.
#'
#' @param Lstates A string listing the possible left states, which will be split by \code{splitval}.
#' @param Rstates A string listing the possible right states, which will be split by \code{splitval}.
#' @param ancareas_txt_tmp A string listing the possible ancestral states, which will be split by \code{splitval}.
#' @param splitval The character to split on.
#' @param code_for_overlapping_subsets Hypothetically, there is no reason that a vicariance event could happen, e.g. ABC-->AB, BC.  This is disallowed in \code{LAGRANGE}
#' BioGeoBEARS defaults, and, if one is going to employ the construct of discrete areas in the first place, overlaps should probably be avoided.  But this parameter
#' will allow experimentation.  Here, \code{code_for_overlapping_subsets=NA} equals the default, and any other value means that overlapping vicariance events are
#' included, with a number describing the number of areas in the overlap.  Users could then manually convert this to a probability according to some function.
#' @return \code{returncell} The text specifying the type of transition.
#' @export
#' @seealso \code{\link{size_species_matrix}}, \code{\link{make_relprob_matrix_bi}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @examples
#' testval=1
#' 
make_spmat_row = function(Lstates, Rstates, ancareas_txt_tmp, splitval="", code_for_overlapping_subsets=NA)
	{
	Lareas = strsplit(Lstates, split=splitval)[[1]]
	Rareas = strsplit(Rstates, split=splitval)[[1]]
	
	#Lareas = strsplit(LR2cells[1], split=splitval)[[1]]
	#Rareas = strsplit(LR2cells[2], split=split)[[1]]
	ancareas = strsplit(ancareas_txt_tmp, split=splitval)[[1]]
	
	# OK, now apply the rules
	# First, rule out impossible cases
	allareas = unique(c(Lareas, Rareas))
	
	# Any cells with "_" states in ancestor or descendant are ruled out
	if ( any(ancareas == "_") == TRUE) 
		{
		returncell = 0
		return(returncell)
		}
	
	# These include any where ranges are lost in both lineages
	if (any( (ancareas %in% allareas) == FALSE) )
		{
		returncell = 0
		return(returncell)
		}
	
	# And any cases where ALL of the descendant ranges are NOT in anc
	if (all( ((allareas %in% ancareas) == FALSE) ))
		{
		returncell = 0
		return(returncell)
		}
	
	# And any cases where >1 descstates are not in anc
	# (The only way to add a range is with jump dispersal, in
	#  this formulation. Which I'm not sure if I believe, colonization
	#  range-addition could be very common compared to background
	#  shifts in range over geological time.  But whatevs.)
	new_states_in_daughters_TF = (allareas %in% ancareas) == FALSE
	num_newareas_in_daughters = sum(new_states_in_daughters_TF)
	
	maxnum_new_areas_at_speciation = 1
	if (num_newareas_in_daughters > maxnum_new_areas_at_speciation)
		{
		returncell = 0
		return(returncell)
		}
	
	# If both descendants are identical to the ancestor,
	# we have "perfect" sympatric speciation
	#
	# Lagrange allows sympatry only when range size = 1
	# 
	# Case: sympatry
	if ( (length(Lareas)==length(ancareas)) && (length(Lareas)==length(Rareas)) && all(Lareas==Rareas) && all(Lareas==ancareas) )
		{
		rangesize = length(ancareas)
		returncell = paste("y", rangesize, sep="")
		return(returncell)
		}
	
	
	# If one of the descendants is identical to the ancestor, 
	# we can have jump dispersal or subset speciation
	LRareas = list(Lareas, Rareas)
	
	
	# FIX LRareas!!!
	
	LR_identical_TF = c( ((length(ancareas)==length(LRareas[[1]])) &&  all(ancareas==LRareas[[1]]) ), ((length(ancareas)==length(LRareas[[2]])) &&  all(ancareas==LRareas[[2]])))
	if ( any(LR_identical_TF) )
		{
		# BUT: We have to make sure that the non-identical daughter
		# is either (a) a single-new-area or (b) an actual subset
		identical_daughter = list(Lareas, Rareas)[LR_identical_TF][[1]]
		nonidentical_daughter = list(Lareas, Rareas)[LR_identical_TF==FALSE][[1]]
		
						
		# Case: jump dispersal ("j"), founding a new lineage in a single region
		# Jump dispersal if 1 area identical and the other with a single new area
		rangesize_nonidentical_daughter = length(nonidentical_daughter)
		if ( (rangesize_nonidentical_daughter == 1) && ( all(nonidentical_daughter %in% identical_daughter) == FALSE) )
			{
			#returncell = paste(returncell, "+j", sep="")
			returncell = "j"
			return(returncell)
			}
		
		# If both are identical, that's sympatry ("y"), we've already done that
		
		# If one new species has a subset of the range, that's sympatry2 (subset of range, "s")
		
		# For now, we only allow subset speciation where the new species has a range of 1
		if ( (all(nonidentical_daughter %in% identical_daughter) == TRUE) )
			{
			returncell = paste("s", length(nonidentical_daughter), "_", length(identical_daughter), sep="")
			return(returncell)		
			} else {
			returncell = 0
			return(returncell)					
			}


		} # end where one branch has an identical subset
	
	
	# Vicariance ("v") can only happen when daughter states cover all of the ancestral states
	if (all(allareas %in% ancareas))
		{
		# Classic vicariance occurs where the descendant ranges add up to the
		# ancestral range, perfectly
		desc_states_in_order = sort(c(Lareas, Rareas))
		
		# If and only if this matches the ancareas can we have vicariance
		# (the length test avoids a warning message)
		if ((length(desc_states_in_order)==length(ancareas)) && all(ancareas == desc_states_in_order))
			{
			# Get the smaller of the two range sizes
			range_sizes = c(length(Lareas), length(Rareas))
			smaller_vicariant_rangesize = min(range_sizes)
			ancestral_rangesize = length(ancareas)
			
			# Add the smaller range size
			#returncell = "v"
			returncell = paste("v", smaller_vicariant_rangesize, "_", ancestral_rangesize, sep="")
			return(returncell)
			} else {
			# The remaining cases are basically vicariance but with overlapping subsets
			if (is.na(code_for_overlapping_subsets))
				{
				returncell = 0
				return(returncell)
				} else {
				# Calculate length of overlap
				numareas_overlapping_in_daughters = sum(Lareas %in% Rareas)
				
				# put e.g. "v?2", vicariance with 2 overlapping areas
				returncell = paste(code_for_overlapping_subsets, numareas_overlapping_in_daughters, sep="")
				return(returncell)
				}
			}
		}
	
	# Remaining cases (hopefully none) are left blank
	returncell = 0			
	} # end tmpfunc




#######################################################
# make_relprob_matrix_bi
#######################################################
#' Make a relative probability matrix for a single speciation (bifurcation) event
#' 
#' Given the identity of the states/geographic ranges on the left branch (\code{Lstates}), right branch (\code{Rstates}), and ancestral areas (\code{ancareas_txt_tmp}),
#' construct the (text version) of the row of transition probabilities.  This means that each nonzero cell gets a \emph{v} for a vicariance event, a \emph{y}
#' for a sympatric speciation/range-copying event, a \emph{j} for a founder-event/jump speciation event, and an \emph{s} for a sympatric-subset event.
#'
#' This function is utilized by \code{\link[base]{apply}} in other functions (e.g. ) in an attempt to speed up calculation over rows.  However, processing of
#' text formulas via \code{\link[base]{apply}} will
#' never be fast enough for large matrices; see \code{\link[cladoRcpp]{cladoRcpp}} for optimized functions.
#' 
#' This text-based matrix later gets evaluated by other functions to calculate the numerical probabilities.  I.e., if j=0 and the other forms of speciation have weights 
#' equal to each other, this is the \code{LAGRANGE} cladogenesis model.
#'
#' NOTE: This function is veeeeeeery slow, even for only 3 areas (i.e. \code{2^3=8} geographic ranges).  It is mostly useful for illustration. See \code{\link[cladoRcpp]{cladoRcpp}} for drastic improvements in calculating cladogenesis models.
#'
#' @param states_list A list of states, where each state consists of a list of areas. A default example list is provided.
#' @param split_ABC \code{TRUE} or \code{FALSE} If \code{TRUE} then each state/range in the input geographic ranges (\code{states_list}) 
#' will be split on the argument contained in \code{split}.
#' @param splitval The character to split on.
#' @param code_for_overlapping_subsets Hypothetically, there is no reason that a vicariance event could happen, e.g. ABC-->AB, BC.  This is disallowed in \code{LAGRANGE}
#' BioGeoBEARS defaults, and, if one is going to employ the construct of discrete areas in the first place, overlaps should probably be avoided.  But this parameter
#' will allow experimentation.  Here, \code{code_for_overlapping_subsets=NA} equals the default, and any other value means that overlapping vicariance events are
#' included, with a number describing the number of areas in the overlap.  Users could then manually convert this to a probability according to some function.
#' @param printwarn If \code{printwarn>0} (\code{printwarn=1} by default), then print to screen a message describing the size of the cladogenesis matrix.
#' @return \code{probmat} A matrix of strings, where each cell contains the parameters describing the conditional probability of that ancestor-->(Left descendant,Right descendant) range inheritance scenario.
#' @export
#' @seealso \code{\link{size_species_matrix}}, \code{\link{make_spmat_row}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @examples
#' testval=1
#' probmat = make_relprob_matrix_bi(states_list=list("_", 
#' c("A"), c("B"), c("C"), c("A","B"), c("B","C"), c("A","C"), 
#' c("A","B","C")), split_ABC=FALSE, splitval="", 
#' code_for_overlapping_subsets=NA, printwarn=1)
#' probmat
#' 
make_relprob_matrix_bi <- function(states_list=default_states_list(), split_ABC=FALSE, splitval="", code_for_overlapping_subsets=NA, printwarn=1)
	{
	defaults='
	states_list=list("_", c("A"), c("B"), c("C"), c("A","B"), c("B","C"), c("A","C"), c("A","B","C"))
	split_ABC=FALSE
	splitval=""
	code_for_overlapping_subsets=NA
	printwarn=1
	'
	
	# If split_ABC=TRUE, assume areas can be split on ""
	if (split_ABC == TRUE)
		{
		states_list = mapply(FUN=strsplit2, split=splitval, states_list)
		names(states_list) = NULL		
		}	
	
	# Remove any null state
	#cat("\n\nNote: make_relprob_matrix_bi() is about range change at the point event of
	#  speciation, therefore the geographic range of NULL is irrelevant, it is
	#  being deleted.\n\n")
	#states_list = states_list[unlist(lapply(X=states_list, FUN=is.null) == FALSE)]
	# Remove the "NULL" range of "_"
	states_list = states_list[states_list != "_"]
	# Instead, set the probability of any range combo with "_" in it to 0


	# Check for the size of this matrix
	if (printwarn > 0)
		{
		# Print the predicted size of the spmat for the user.
		size_species_matrix(states_list, printwarn)
		}


	
	# We need descendant states for each of two branches
	# (this does not include diagonal, same-state for some reason)
	# row1 = states on left branch just above node
	# row2 = states on right branch just above node
	
	# 	JUNK
	# 	if (length(states_list) >= 2)
	# 		{
	# 		desc_states_mat = combn(x=states_list, m=2)
	# 		} else {
	# 		desc_states_mat = matrix
	# 		}
	

	# Insert the sympatric events (e.g., A --> A,A)
	# comma: yes
	# sympatric_desc_names = sapply(X=states_list, FUN=paste, collapse=",")
	# comma: no
	sympatric_desc_names = sapply(X=states_list, FUN=paste, collapse="")

	len_states_list = length(states_list)	
	if (len_states_list > 1)
		{
		# Initialize blank matrices
		RL = matrix(NA, nrow=len_states_list^2, ncol=2)
		RL2 = matrix(NA, nrow=len_states_list^2, ncol=2)
		LR = matrix(NA, nrow=len_states_list^2, ncol=2)
		descstates_txt = rep(NA, len_states_list^2)
		
		# Right and Left descendants
		# RL, RL2, and LR have 2 columns; these are then 
		#RL = expand.grid(R=sympatric_desc_names, L=sympatric_desc_names)
		#RL = expand.grid.jc(seq1=sympatric_desc_names, seq2=sympatric_desc_names)
		#RL2 = sapply(FUN=as.character, RL)
		#LR = cbind(RL2[,2], RL2[,1])

		RL = expand.grid.jc(seq1=sympatric_desc_names, seq2=sympatric_desc_names)
		#RL2 = sapply(FUN=as.character, RL)
		LR = cbind(RL[,2], RL[,1])
		#descstates_txt = paste(RL2[,2], RL2[,1], collapse="|")
		descstates_txt = apply(LR, 1, paste, collapse="|")
	
		# comma: yes
		# ancstates_txt = sapply(states_list, FUN=paste, collapse=",")
		# comma: no
		ancstates_txt = sapply(states_list, FUN=paste, collapse="")
		} else {
		# Sympatry only allowed
		descstates_txt = paste(as.character(sympatric_desc_names), "|", as.character(sympatric_desc_names), sep="")
		ancstates_txt = as.character(sympatric_desc_names)
		LR = cbind(sympatric_desc_names, sympatric_desc_names)
		}
	
	# Now make a matrix describing the probabilities
	probmat = matrix("", nrow=length(ancstates_txt), ncol=length(descstates_txt))
	
	
	# APPLY HERE ALSO
	# Rows/ancestral states
	for (i in 1:length(ancstates_txt))
		{
		#ancstates = strsplit(ancstates_txt[i], split=",")[[1]]
		#ancstates = strsplit(ancstates_txt[i], split=splitval)[[1]]
		
		print(i)
		probmat[i,] = mapply(FUN=make_spmat_row, Lstates=LR[,1], Rstates=LR[,2], ancareas_txt_tmp=ancstates_txt[i], splitval=splitval)

		
		} # end ancstates
	
	probmat = adf2(probmat)
	names(probmat) = descstates_txt
	row.names(probmat) = ancstates_txt
	
	return(probmat)
	}




#######################################################
# divide_probs_by_number_of_options_txt
#######################################################
#' Divide each type of event by its frequency
#' 
#' In a speciation/cladogenesis matrix, the conditional probabilities of each row must sum to 
#' 1.  This function sums the number of events of each category and scales them accordingly.
#'
#' This function returns the strings, which can then be processed in other functions by e.g. find/replace or \code{\link[base]{eval}}.
#'
#' @param probmat A character matrix of probabilities in the form of formulas, not normalized by the sum of each row.
#' @return \code{probmat} A matrix of strings, where each cell contains the parameters describing the conditional probability of
#' that ancestor-->(Left descendant,Right descendant) range inheritance scenario.
#' @export
#' @seealso \code{\link{make_relprob_matrix_bi}}, \code{\link{divide_probs_by_number_of_options_nums}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' probmat = make_relprob_matrix_bi()
#' probmat
#' 
#' probmat2 = divide_probs_by_number_of_options_txt(probmat)
#' probmat2
#' 
#' 
divide_probs_by_number_of_options_txt <- function(probmat)
	{
	defaults='
	probmat = spmat
	'
	
	for (i in 1:nrow(probmat))
		{
		tmprow = unlist(probmat[i, ])
		
		# Group everything by literal match (commented out)
		#unique_transitions = unlist(unique(tmprow))

		# Just take the first letter, and uniquify that
		# (so that e.g. s1, s2, or v1, v2, are in the same category)
		# Group everything with the same first letter
		tmprow_letters = sapply(X=tmprow, FUN=substr, start=1, stop=1)
		unique_transitions = unlist(unique(tmprow_letters))
		
		# Remove zeros
		unique_transitions = unique_transitions[unique_transitions != "0"]
		unique_transitions
		
		# Go through the types of transitions, and weight by their frequency
		for (j in 1:length(unique_transitions))
			{
			# Group everything by literal match (commented out)
			#count_of_this_transition = sum(tmprow == unique_transitions[j])
			
			# Group everything with the same first letter
			count_of_this_transition = sum(tmprow_letters == unique_transitions[j])
			
			# If there is only one of this type of transition, do nothing
			if (count_of_this_transition == 1)
				{
				next()
				} else{
				# Otherwise, divide by its frequency in the row
				cells_to_change_TF = tmprow_letters == unique_transitions[j]
				
				# Group everything by literal match (commented out)
				#newtxt = paste(unique_transitions[j], "/", count_of_this_transition, sep="")
				# Group everything with the same first letter
				newtxt = paste(tmprow[cells_to_change_TF], "/", count_of_this_transition, sep="")
				tmprow[cells_to_change_TF] = newtxt
				}
			}		
		# Input the revised cells into the original probmat	
		probmat[i, ] = tmprow
		}
	return(probmat)
	}




#######################################################
# divide_probs_by_number_of_options_nums
#######################################################
#' Divide each type of event by its frequency, return calculated probabilities
#' 
#' In a speciation/cladogenesis matrix, the conditional probabilities of each row must sum to 
#' 1.  This function sums the number of events of each category and scales them accordingly.
#'
#' This function returns the calculated conditional probabilities.
#'
#' @param spPmat A matrix of numbers, where each cell contains the conditional probability of
#' that ancestor-->(Left descendant,Right descendant) range inheritance scenario.
#' @param probmat A matrix of text, describing each of the allowed range-inheritance events.
#' @return \code{spPmat} A matrix of numbers, where each cell contains the conditional probability of
#' that ancestor-->(Left descendant,Right descendant) range inheritance scenario.
#' @export
#' @seealso \code{\link{make_relprob_matrix_bi}}, \code{\link{divide_probs_by_number_of_options_txt}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' spmat = make_relprob_matrix_bi()
#' spmat
#' 
#' spmat1 = divide_probs_by_number_of_options_txt(spmat)
#' spmat1
#' 
#' 
#' probmat = spmat
#' spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", 
#' mergesym="*", ys=1, j=0, v=1, maxent_constraint_01=0.0001, 
#' maxent_constraint_01v=0.0001, max_numareas=3)
#' spPmat
#' probmat2 = divide_probs_by_number_of_options_nums(spPmat, probmat)
#' probmat2
#' 
#' probmat = spmat1
#' spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\+", 
#' mergesym="*", ys=1, j=0, v=1, maxent_constraint_01=0.0001, 
#' maxent_constraint_01v=0.0001, max_numareas=3)
#' spPmat
#' probmat3 = divide_probs_by_number_of_options_nums(spPmat, probmat)
#' probmat3
#' 
divide_probs_by_number_of_options_nums <- function(spPmat, probmat)
	{
	defaults='
	probmat = spmat
	'
	
	for (i in 1:nrow(probmat))
		{
		tmprow = unlist(probmat[i, ])

		# Group everything by literal match (commented out)
		#unique_transitions = unlist(unique(tmprow))

		# Just take the first letter, and uniquify that
		# (so that e.g. s1, s2, or v1, v2, are in the same category)
		# Group everything with the same first letter
		tmprow_letters = sapply(X=tmprow, FUN=substr, start=1, stop=1)
		unique_transitions = unlist(unique(tmprow_letters))
	
		# Remove zeros
		unique_transitions = unique_transitions[unique_transitions != "0"]
		unique_transitions
		
		# Go through the types of transitions, and weight by their frequency
		for (j in 1:length(unique_transitions))
			{
			# Group everything by literal match (commented out)
			#count_of_this_transition = sum(tmprow == unique_transitions[j])
			
			# Group everything with the same first letter
			count_of_this_transition = sum(tmprow_letters == unique_transitions[j])
			
			# If there is only one of this type of transition, do nothing
			if (count_of_this_transition == 1)
				{
				next()
				} else{
				# Otherwise, divide by its frequency in the row
				#cells_to_change_TF = tmprow == unique_transitions[j]

				# Otherwise, divide by its frequency in the row
				cells_to_change_TF = tmprow_letters == unique_transitions[j]

				
				#newtxt = paste(unique_transitions[j], "/", count_of_this_transition, sep="")
				#tmprow[cells_to_change_TF] = newtxt
				spPmat[i, cells_to_change_TF] = spPmat[i, cells_to_change_TF] / count_of_this_transition
				}
			}		
		# Input the revised cells into the original probmat	
		#probmat[i, ] = tmprow
		}
	
	# Re-normalize so rows sum to 1
	spPmat = spPmat / rowSums(spPmat)
	return(spPmat)
	}




# Convert a observed-speciation transition matrix to an
# unobserved-speciation transition matrix
#######################################################
# make_relprob_txtmatrix_sp1
#######################################################
#' Convert a observed-speciation transition matrix to an unobserved-speciation transition matrix (text version)
#'
#' Convert a cladogenesis/speciation transition matrix (specifying the probability of each Left/Right descendant range pair,
#' conditional on each ancestral state) of dimensions \code{numstates} by \code{numstates^2} to a square transition matrix of dimensions
#' \code{numstates} by \code{numstates}, representing the probability of a transition when only one daughter survives in the tree.
#' 
#' This matrix could be used to quantify the probability of range-change along a branch due to unobserved speciation events; all that would 
#' be required would be an estimate of the number of unobserved speciation events on the branch, and treating this as a Poisson process.
#' (Note: this assumes that the probability of either branch surviving is identical, which might not be the case. See the GeoSSE
#' (\cite{Goldberg_etal_2011_GeoSSE}) and ClaSSE (\cite{Goldberg_Igic_2012_ClaSSE}) for the beginnings of work on this, with 2 and 3
#' geographic areas, respectively.
#' 
#' @param probmat A matrix of text, describing each of the allowed range-inheritance events.  Assumes that column names are in the "A|B" format.
#' @param split The value to split Left/Right pairs on (e.g., "A|B" --> "A", "B")
#' @return \code{newmat} A new square matrix.
#' @export
#' @seealso \code{\link{make_relprob_matrix_bi}}, \code{\link{make_relprob_nummatrix_sp1}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' \url{http://tigger.uic.edu/~eeg/code/code.html}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#'	 @cite Goldberg_etal_2011_GeoSSE
#'	 @cite Goldberg_Igic_2012_ClaSSE
#' @examples
#' testval=1
#' probmat = make_relprob_matrix_bi(states_list=list("_", c("A"), 
#' c("B"), c("C"), c("A","B"), c("B","C"), c("A","C"), c("A","B","C")), 
#' split_ABC=FALSE, splitval="", code_for_overlapping_subsets=NA, printwarn=1)
#' probmat
#' 
#' newmat = make_relprob_txtmatrix_sp1(probmat=probmat, split="\\\\|")
#' newmat
#' 
make_relprob_txtmatrix_sp1 <- function(probmat, split="\\|")
	{
	# In this case, at the unobserved speciation event, there is only one 
	# known descendant.  This could (hypothetically) be 
	# either the "left" or "right" branch of a speciation
	# event -- obviously the labels are arbitrary.
	#
	# The probability of any particular ancestor-descendant transition can
	# therefore be found by looking at the transition probabilities for the
	# species bifurcation matrix, and summing the probabilities resulting
	# in the same state for the left branch.
	
	# Get unique states
	tmpnames = names(probmat)
	states = unlist(lapply(tmpnames, strsplit2, split="\\|"))
	uniq_states = unique(states)
	uniq_states
	
	nr = dim(probmat)[1]
	nc = dim(probmat)[2]
	
	newmat = matrix(NA, nrow=nr, ncol=nr)
	
	for (i in 1:length(uniq_states))
		{
		u = uniq_states[i]
		startcol = (i-1)*nr + 1
		endcol = (i)*nr
		
		# Subset to just the columns for a particular state on 
		# the left branch
		tmpmat = probmat[,startcol:endcol]
		tmpmat
		
		# Paste each row together
		tmpcol = paste_rows_without_zeros(tmpmat)
		
		# Insert into matrix
		newmat[,i] = tmpcol
		}
	newmat = adf2(newmat)
	names(newmat) = uniq_states
	row.names(newmat) = uniq_states
	
	return(newmat)
	}


# Convert a observed-speciation transition matrix to an
# unobserved-speciation transition matrix
#######################################################
# make_relprob_nummatrix_sp1
#######################################################
#' Convert a observed-speciation transition matrix to an unobserved-speciation transition matrix (numeric version)
#'
#' Convert a cladogenesis/speciation transition matrix (specifying the probability of each Left/Right descendant range pair,
#' conditional on each ancestral state) of dimensions \code{numstates} by \code{numstates^2} to a square transition matrix of dimensions
#' \code{numstates} by \code{numstates}, representing the probability of a transition when only one daughter survives in the tree.
#' 
#' This matrix could be used to quantify the probability of range-change along a branch due to unobserved speciation events; all that would 
#' be required would be an estimate of the number of unobserved speciation events on the branch, and treating this as a Poisson process.
#' (Note: this assumes that the probability of either branch surviving is identical, which might not be the case. See the GeoSSE
#' (\cite{Goldberg_etal_2011_GeoSSE}) and ClaSSE (\cite{Goldberg_Igic_2012_ClaSSE}) for the beginnings of work on this, with 2 and 3
#' geographic areas, respectively.
#' 
#' @param probmat A matrix of text, describing each of the allowed range-inheritance events.  Assumes that column names are in the "A|B" format.
#' @param spPmat A matrix of numbers, where each cell contains the conditional probability of
#' that ancestor-->(Left descendant,Right descendant) range inheritance scenario.
#' @param split The value to split Left/Right pairs on (e.g., "A|B" --> "A", "B")
#' @return \code{newmat} A new square matrix.
#' @export
#' @seealso \code{\link{make_relprob_matrix_bi}}, \code{\link{make_relprob_txtmatrix_sp1}}, \code{\link{paste_rows_without_zeros}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' \url{http://tigger.uic.edu/~eeg/code/code.html}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#'	 @cite Goldberg_etal_2011_GeoSSE
#'	 @cite Goldberg_Igic_2012_ClaSSE
#' @examples
#' testval=1
#' spmat = make_relprob_matrix_bi(states_list=list("_", c("A"), 
#' c("B"), c("C"), c("A","B"), c("B","C"), c("A","C"), c("A","B","C")), 
#' split_ABC=FALSE, splitval="", code_for_overlapping_subsets=NA, printwarn=1)
#' spmat
#' 
#' spPmat = symbolic_to_relprob_matrix_sp(spmat, cellsplit="\\\\+", 
#' mergesym="*", ys=1, j=0, v=1, maxent_constraint_01=0.0001, 
#' maxent_constraint_01v=0.0001, max_numareas=3)
#' spPmat
#' 
#' newmat = make_relprob_nummatrix_sp1(probmat=spmat, spPmat=spPmat, split="\\\\|")
#' newmat
#'
make_relprob_nummatrix_sp1 <- function(probmat, spPmat, split="\\|")
	{
	# In this case, at the unobserved speciation event, there is only one 
	# known descendant.  This could (hypothetically) be 
	# either the "left" or "right" branch of a speciation
	# event -- obviously the labels are arbitrary.
	#
	# The probability of any particular ancestor-descendant transition can
	# therefore be found by looking at the transition probabilities for the
	# species bifurcation matrix, and summing the probabilities resulting
	# in the same state for the left branch.
	
	# Get unique states
	tmpnames = names(probmat)
	states = unlist(lapply(tmpnames, strsplit2, split="\\|"))
	uniq_states = unique(states)
	uniq_states
	
	nr = dim(probmat)[1]
	nc = dim(probmat)[2]
	
	newmat = matrix(NA, nrow=nr, ncol=nr)
	
	for (i in 1:length(uniq_states))
		{
		u = uniq_states[i]
		startcol = (i-1)*nr + 1
		endcol = (i)*nr
		
		# Subset to just the columns for a particular state on 
		# the left branch
		#tmpmat = probmat[,startcol:endcol]
		tmpmat = spPmat[,startcol:endcol]
		tmpmat
		
		# Paste each row together
		tmpcol = apply(tmpmat, 1, paste, collapse="+")
		
		tmpcol_nums = sapply(X=tmpcol, FUN=meval)
		
		# Insert into matrix
		newmat[,i] = tmpcol_nums
		}
	newmat = adf2(newmat)
	names(newmat) = uniq_states
	row.names(newmat) = uniq_states
	
	return(newmat)
	}



#######################################################
# meval
#######################################################
#' \code{eval()} function for use in \code{sapply}
#' 
#' \code{meval} is a wrapper for \code{\link{eval}}, to allow use in \code{sapply}.
#' 
#' This is an attempt to speed up the use of \code{\link{eval}}; in general use of \code{\link{eval}} to 
#' convert a text version of a transition matrix to a numeric version with probabilities
#' is a poor, slow choice; but it can be useful for examples and display purposes.
#' 
#' See \code{\link[cladoRcpp]{cladoRcpp}} for fast C++ implementations of transition matrix setup.
#' 
#' @param equation_txt The text of the equation to run \code{\link{eval}} on -- e.g., from a cell of a text-based transition matrix.
#' @return \code{outval} The numeric result of \code{\link{eval}}.
#' @export
#' @seealso \code{\link[stats]{convolve}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' d = 0.1
#' equation_txt = "1*d+1*d"
#' meval(equation_txt)
meval <- function(equation_txt)
	{
	outval = eval(parse(text = equation_txt))
	return(outval)
	}

	
# Convert e.g.
#
#       A|A A|B A|C A|A,B A|B,C A|A,C A|A,B,C
# A       s   j   j     0     0     0       0
# B       0   j   0     0     0     0       0
# C       0   0   j     0     0     0       0
# A,B     0   v   0    b1     0     0       0
# B,C     0   0   0     0     j     0       0
# A,C     0   0   v     0     0    b1       0
# A,B,C   0   0   0     0     v     0      b1
#
# ...to...
# 
#    A         B       C     A,B     B,C     A,C   A,B,C 
# "s+j+j"     "j"     "j"  "v+b1"     "j"  "v+b1"  "v+b1" 
# 


#######################################################
# paste_rows_without_zeros
#######################################################
#' Concatenate cells in each row of a text-based transition matrix, excluding zeros
#' 
#' This is a utility function for \code{\link{make_relprob_txtmatrix_sp1}}.
#' 
#' Convert e.g.:\cr
#' \cr
#' \code{       A|A A|B A|C A|A,B A|B,C A|A,C A|A,B,C}\cr
#' \code{ A       s   j   j     0     0     0       0}\cr
#' \code{ B       0   j   0     0     0     0       0}\cr
#' \code{ C       0   0   j     0     0     0       0}\cr
#' \code{ A,B     0   v   0    b1     0     0       0}\cr
#' \code{ B,C     0   0   0     0     j     0       0}\cr
#' \code{ A,C     0   0   v     0     0    b1       0}\cr
#' \code{ A,B,C   0   0   0     0     v     0      b1}\cr
#' \cr
#' ...to...\cr
#' \cr
#' \code{  A         B       C     A,B     B,C     A,C   A,B,C }\cr
#' \code{"s+j+j"     "j"     "j"  "v+b1"     "j"  "v+b1"  "v+b1" }\cr
#' \cr
#' 
#' @param tmpmat A cladogenesis/speciation probability matrix (text-based) to collapse each row of.
#' @return \code{tmpcol} A list containing each row, concatenated
#' @export
#' @seealso \code{\link{make_relprob_txtmatrix_sp1}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' 
#' spmat = make_relprob_matrix_bi(states_list=list("_", c("A"), 
#' c("B"), c("C"), c("A","B"), c("B","C"), c("A","C"), c("A","B","C")), 
#' split_ABC=FALSE, splitval="", code_for_overlapping_subsets=NA, printwarn=1)
#' spmat
#' tmpcol = paste_rows_without_zeros(tmpmat=spmat)
#' tmpcol
#' 
paste_rows_without_zeros <- function(tmpmat)
	{
	tmpcol = apply(tmpmat, 1, paste, collapse="+")
	tmpcol = gsub(pattern="\\+0", replacement="", tmpcol)
	tmpcol = gsub(pattern="0\\+", replacement="", tmpcol)
	return(tmpcol)
	}



#######################################################
# remove_null_rowcols_from_mat
#######################################################
#' Remove rows or columns representing a null geographic range from a matrix
#' 
#' This function removes rows or columns representing a null geographic range from a matrix. 
#' 
#' LAGRANGE (\cite{ReeSmith2008}) and other models often assume that a null geographic range (the lineage inhabits no
#' areas, i.e. is extinct) is a possible state.  However, this is never a possible ancestral state (since an extinct lineage will
#' never have descendants) so sometimes we must remove it.
#'
#' @param tmpmat The matrix to check for null ranges.  Function will only work if rows and columns have names, and one of the names matches \code{null_sym}.
#' @param null_sym The character(s) denoting a null range.
#' @return \code{tmpmat3} The revised matrix.
#' @export
#' @seealso \code{\link{areas_list_to_states_list_new}}, \code{\link[cladoRcpp]{areas_list_to_states_list_old}},
#' \code{\link{make_relprob_matrix_de}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @examples
#' testval=1
#' states_list = list("_", c("A"), c("B"), c("C"), c("A","B"), 
#' c("B","C"), c("A","C"), c("A","B","C"))
#' 
#' states_list = areas_list_to_states_list_new(areas=c("A","B","C"), 
#' include_null_range=TRUE, split_ABC=TRUE)
#' states_list
#' 
#' dedf = make_relprob_matrix_de(states_list=states_list, 
#' split_ABC=FALSE, split="", remove_simultaneous_events=TRUE, 
#' add_multiple_Ds=TRUE, 
#' dispersal_multiplier_matrix=make_dispersal_multiplier_matrix(states_list=states_list))
#' 
#' spmat_noNulls = remove_null_rowcols_from_mat(tmpmat=dedf, null_sym="()")
#' spmat_noNulls
#' 
#' spmat_noNulls = remove_null_rowcols_from_mat(tmpmat=dedf, null_sym="_")
#' spmat_noNulls
#' 
remove_null_rowcols_from_mat <- function(tmpmat, null_sym="()")
	{
	# Remove column with null state
	tmpnames = names(tmpmat)
	TF = tmpnames != null_sym
	tmpmat2 = tmpmat[,TF]
	
	# Remove row with null state
	tmpnames = row.names(tmpmat)
	TF = tmpnames != null_sym
	tmpmat3 = tmpmat2[TF,]
	
	return(tmpmat3)
	}




# Convert symbolic matrix to relprob matrix
#######################################################
# symbolic_to_P_matrix
#######################################################
#' Convert symbolic matrix to relprob matrix
#' 
#' This function takes a transition probability matrix (in text form) and converts to numeric form, given
#' values for \emph{d}, \emph{e}, or other parameters in the text formulas.
#' 
#' This is not particularly fast, but good for illustrative purposes.
#'
#' @param dedf The transition matrix or dispersal-extinction data.frame (dedf), contains the 
#' actual text of the formulas by which the transition probability matrix would be calculated.
#' @param cellsplit The symbol to split the formulas on. Default "\\+" (plus symbol, with escape code).
#' @param mergesym The symbol to merge the formulas with. Default "+".
#' @param diags_sum_to_1 Calculate the diagonals such that, when added to the sum of the off-diagonals in a row, the entire row sums to 1. This creates a transition
#' probability matrix where each row sums to 1, i.e. each cell represents the conditional probability of the column state, given the ancestral
#' row state.  The diagonal values represent the probability of staying the same.
#' @param d The dispersal/range expansion rate. Default \code{d=0.1}.
#' @param e The extinction/range contraction rate. Default \code{e=0.01}.
#' @param ... Additional arguments to pass to \code{\link{symbolic_cell_to_relprob_cell}} via \code{\link[base]{sapply}}, and thence to cell\code{\link[base]{strsplit}}.
#' @return \code{dedf_vals} The output \code{\link[base]{data.frame}}, contains the 
#' numeric results of the formulas calculating the transition probability matrix.
#' @export
#' @seealso \code{\link{areas_list_to_states_list_new}}, \code{\link[cladoRcpp]{areas_list_to_states_list_old}},
#' \code{\link{make_relprob_matrix_de}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' 
#' states_list = list("_", c("A"), c("B"), c("C"), c("A","B"), 
#' c("B","C"), c("A","C"), c("A","B","C"))
#' 
#' states_list = areas_list_to_states_list_new(areas=c("A","B","C"), 
#' include_null_range=TRUE, split_ABC=TRUE)
#' states_list
#' 
#' dedf = make_relprob_matrix_de(states_list=states_list,
#'  split_ABC=FALSE, split="", remove_simultaneous_events=TRUE, 
#' add_multiple_Ds=TRUE, 
#' dispersal_multiplier_matrix=make_dispersal_multiplier_matrix(states_list=states_list))
#' dedf
#' 
#' # Defaults
#' Pmat = symbolic_to_P_matrix(dedf, cellsplit="\\\\+", mergesym="+", 
#' diags_sum_to_1=FALSE, d=0.1, e=0.01)
#' Pmat
#' 
#' # Calculate diagonal
#' Pmat = symbolic_to_P_matrix(dedf, cellsplit="\\\\+", mergesym="+", 
#' diags_sum_to_1=TRUE, d=0.1, e=0.01)
#' Pmat
#' 
#' # You don't have to split, if the formulas are directly parsable
#' Pmat = symbolic_to_P_matrix(dedf, cellsplit="yadda", mergesym="", 
#' diags_sum_to_1=FALSE, d=0.1, e=0.01)
#' Pmat
#' 
symbolic_to_P_matrix <- function(dedf, cellsplit="\\+", mergesym="+", diags_sum_to_1=FALSE, d=0.1, e=0.01, ...)
	{
	defaults='
	cellsplit=""
	mergesym="*"
	diags_sum_to_1=TRUE
	d=0.1
	e=0.01
	'
	
	dedf[dedf == "-"] = NA
	
	dedf_vals = matrix(data=0, nrow=nrow(dedf), ncol=ncol(dedf))
	
	dedf_vals[!is.na(dedf_vals)] = sapply(X=dedf[!is.na(dedf_vals)], FUN=symbolic_cell_to_relprob_cell, cellsplit, mergesym, d, e, ...)
	
	# DON'T Normalize matrix into an appropriate Q matrix
	#diag_TF = diag(x=1, nrow=nrow(dedf), ncol(dedf)) == 1
	#dedf_vals[diag_TF] = -1 * rowSums(dedf_vals, na.rm=TRUE)
	
	# But, if desired, DO sum to 1:
	diag_vals = rep(NA, nrow(dedf_vals))
	if (diags_sum_to_1 == TRUE)
		{
		tmp_rowsums = rowSums(dedf_vals, na.rm=TRUE)
		diag(dedf_vals) = 1-tmp_rowsums
		}
	
	return(dedf_vals)
	}


#######################################################
# symbolic_cell_to_relprob_cell
#######################################################
#' Convert symbolic cell (a text equation) to relprob matrix (a numeric value).
#' 
#' This is a utility function for \code{\link{symbolic_to_P_matrix}} and \code{\link{symbolic_to_Q_matrix}}.
#'
#' This function can be used in \code{\link[base]{sapply}}.  It still will not be very fast compared to the calculations in 
#' \code{\link[cladoRcpp]{cladoRcpp}}, but can be useful for demonstrative purposes.
#' 
#' @param charcell The text formula.
#' @param cellsplit The symbol to split the formulas on. Default "\\+" (plus symbol, with escape code).
#' @param mergesym The symbol to merge the formulas with. Default "+".
#' @param d The dispersal/range expansion rate. Default \code{d=0.1}.
#' @param e The extinction/range contraction rate. Default \code{e=0.01}.
#' @param ... Additional arguments to pass to \code{\link[base]{strsplit}}.
#' @return \code{cellval} The output cell value.
#' @export
#' @seealso \code{\link{symbolic_to_P_matrix}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' 
#' charcell = "1*d+1*d"
#' 
#' # Right
#' cellval = symbolic_cell_to_relprob_cell(charcell, cellsplit="yadda", 
#' mergesym="", d=0.1, e=0.01)
#' cellval
#' 
#' # Wrong
#' cellval = symbolic_cell_to_relprob_cell(charcell, cellsplit="\\+", 
#' mergesym="*", d=0.1, e=0.01)
#' cellval
#' 
#' # Right
#' cellval = symbolic_cell_to_relprob_cell(charcell, cellsplit="\\+", 
#' mergesym="+", d=0.1, e=0.01)
#' cellval
#' 
symbolic_cell_to_relprob_cell <- function(charcell, cellsplit="", mergesym="*", d=0.1, e=0.01, ...)
	{
	# split the cell
	words = strsplit(charcell, split=cellsplit, ...)[[1]]
	
	# initialize the cell value
	cellval = NA
	equation = paste(words, sep="", collapse=mergesym)
	equation2 = paste("cellval = ", equation, sep="")
	
	eval(parse(text=equation2))
	
	return(cellval)
	}



#######################################################
# symbolic_to_Q_matrix
#######################################################
#' Convert symbolic matrix to an instantaneous rate matrix (Q matrix)
#' 
#' This function takes a transition probability matrix (in text form) and converts it to an 
#' instantaneous rate matrix (Q matrix), given values for \emph{d}, \emph{e}, or other parameters
#' in the text formulas.
#' 
#' This is not particularly fast, but good for illustrative purposes.
#'
#' @param dedf The transition matrix or dispersal-extinction data.frame (dedf), contains the 
#' actual text of the formulas by which the transition probability matrix would be calculated.
#' @param cellsplit The symbol to split the formulas on. Default "\\+" (plus symbol, with escape code).
#' @param mergesym The symbol to merge the formulas with. Default "+".
#' @param d The dispersal/range expansion rate. Default \code{d=0.1}.
#' @param e The extinction/range contraction rate. Default \code{e=0.01}.
#' @param ... Additional arguments to pass to \code{\link{symbolic_cell_to_relprob_cell}} via \code{\link[base]{sapply}}, and thence to cell\code{\link[base]{strsplit}}.
#' @return \code{dedf_vals} The output \code{\link[base]{data.frame}}, contains the Q matrix 
#' @export
#' @seealso \code{\link{areas_list_to_states_list_new}}, \code{\link[cladoRcpp]{areas_list_to_states_list_old}},
#' \code{\link{make_relprob_matrix_de}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite FosterIdiots
#' @examples
#' testval=1
#' 
#' states_list = list("_", c("A"), c("B"), c("C"), c("A","B"), 
#' c("B","C"), c("A","C"), c("A","B","C"))
#' 
#' states_list = areas_list_to_states_list_new(areas=c("A","B","C"), 
#' include_null_range=TRUE, split_ABC=TRUE)
#' states_list
#' 
#' dedf = make_relprob_matrix_de(states_list=states_list, split_ABC=FALSE, 
#' split="", remove_simultaneous_events=TRUE, add_multiple_Ds=TRUE, 
#' dispersal_multiplier_matrix=make_dispersal_multiplier_matrix(states_list=states_list))
#' dedf
#' 
#' # Right
#' Qmat = symbolic_to_Q_matrix(dedf, cellsplit="\\\\+", mergesym="+", d=0.1, e=0.01)
#' Qmat
#' 
#' # Wrong
#' Qmat = symbolic_to_Q_matrix(dedf, cellsplit="\\\\+", mergesym="*", d=0.1, e=0.01)
#' Qmat
#' 
#' # You don't have to split, if the formulas are directly parsable
#' Qmat = symbolic_to_Q_matrix(dedf, cellsplit="yadda", mergesym="", d=0.1, e=0.01)
#' Qmat
#' 
symbolic_to_Q_matrix <- function(dedf, cellsplit="\\+", mergesym="*", d=0.1, e=0.01, ...)
	{
	defaults='
	cellsplit="\\+"
	mergesym="*"
	'

	defaults2='
	d=10.9212
	e=1.73083e-07
	dedf = dedf_M2
	cellsplit="\\+"
	mergesym="+"
	'
	
	dedf[dedf == "-"] = NA
	
	dedf_vals = matrix(data=0, nrow=nrow(dedf), ncol=ncol(dedf))
	
	dedf_vals[!is.na(dedf_vals)] = sapply(X=dedf[!is.na(dedf_vals)], FUN=symbolic_cell_to_relprob_cell, cellsplit, mergesym, d, e, ...)
	dedf_vals
	
	# Normalize matrix into an appropriate Q matrix
	diag_TF = diag(x=1, nrow=nrow(dedf), ncol(dedf)) == 1
	dedf_vals[diag_TF] = -1 * rowSums(dedf_vals, na.rm=TRUE)
	
	return(dedf_vals)
	}


#######################################################
# symbolic_to_Q_matrix_exper
#######################################################
#' Experimental version of \code{symbolic_to_Q_matrix_exper}, including base frequencies
#' 
#' Still experimental.
#'
#' This function takes a transition probability matrix (in text form) and converts it to an 
#' instantaneous rate matrix (Q matrix), given values for \emph{d}, \emph{e}, or other parameters
#' in the text formulas.
#' 
#' This is not particularly fast, but good for illustrative purposes.
#'
#' @param dedf The transition matrix or dispersal-extinction data.frame (dedf), contains the 
#' actual text of the formulas by which the transition probability matrix would be calculated.
#' @param cellsplit The symbol to split the formulas on. Default "\\+" (plus symbol, with escape code).
#' @param mergesym The symbol to merge the formulas with. Default "+".
#' @param d The dispersal/range expansion rate. Default \code{d=0.1}.
#' @param e The extinction/range contraction rate. Default \code{e=0.01}.
#' @param basefreqs Base frequencies, i.e. the equilibrium probabilities of the different states; the meaning of such an idea is debatable
#' in the context of a LAGRANGE-like model where the null range (extinct everywhere) is included in the matrix and is 
#' a nonreversible absorbing state.  Default is \code{rep(1,nrow(dedf))/nrow(dedf)}.
#' @param ... Additional arguments to pass to \code{\link{symbolic_cell_to_relprob_cell}} via \code{\link[base]{sapply}}, and thence to cell\code{\link[base]{strsplit}}.
#' @return \code{dedf_vals} The output \code{\link[base]{data.frame}}, contains the Q matrix 
#' @export
#' @seealso \code{\link[cladoRcpp]{areas_list_to_states_list_old}}, \code{\link{make_relprob_matrix_de}}
#' \code{\link{areas_list_to_states_list_new}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite FosterIdiots
#' @examples
#' testval=1
#' 
#' states_list = list("_", c("A"), c("B"), c("C"), c("A","B"), 
#' c("B","C"), c("A","C"), c("A","B","C"))
#' 
#' states_list = areas_list_to_states_list_new(areas=c("A","B","C"), 
#' include_null_range=TRUE, split_ABC=TRUE)
#' states_list
#' 
#' dedf = make_relprob_matrix_de(states_list=states_list, split_ABC=FALSE, 
#' split="", remove_simultaneous_events=TRUE, add_multiple_Ds=TRUE, 
#' dispersal_multiplier_matrix=make_dispersal_multiplier_matrix(states_list=states_list))
#' dedf
#' 
#' # Right
#' Qmat = symbolic_to_Q_matrix_exper(dedf, cellsplit="\\\\+", mergesym="+", d=0.1, e=0.01)
#' Qmat
#' 
#' # Wrong
#' Qmat = symbolic_to_Q_matrix_exper(dedf, cellsplit="\\\\+", mergesym="*", d=0.1, e=0.01)
#' Qmat
#' 
#' # You don't have to split, if the formulas are directly parsable
#' Qmat = symbolic_to_Q_matrix_exper(dedf, cellsplit="yadda", mergesym="", d=0.1, e=0.01)
#' Qmat
#' 
#' # Compare to symbolic_to_Q_matrix
#' Qmat = symbolic_to_Q_matrix(dedf, cellsplit="yadda", mergesym="", d=0.1, e=0.01)
#' Qmat
#' 
symbolic_to_Q_matrix_exper <- function(dedf, cellsplit="\\+", mergesym="*", d=0.1, e=0.01, basefreqs=rep(1,nrow(dedf))/nrow(dedf), ...)
	{
	defaults='
	cellsplit=""
	mergesym="*"
	basefreqs=(rep(1,nrow(dedf))/nrow(dedf))

	cellsplit="*"
	mergesym="+"
	basefreqs=(rep(1,nrow(dedf))/nrow(dedf))
	'
	
	dedf[dedf == "-"] = NA
	
	dedf_vals = matrix(data=0, nrow=nrow(dedf), ncol=ncol(dedf))
	
	dedf_vals[!is.na(dedf_vals)] = sapply(X=dedf[!is.na(dedf_vals)], FUN=symbolic_cell_to_relprob_cell, cellsplit, mergesym, d, e, ...)
	
	
	# Now you have a relative rate matrix; you multiply the columns by the composition, then scale so they sum to 1
	Rmat = dedf_vals
	#Rmat2 = basefreqs * Rmat
	
	
	
	# Normalize matrix into an appropriate Q matrix
	Qmat = Rmat
	diag_TF = diag(x=1, nrow=nrow(dedf), ncol(dedf)) == 1
	Qmat[diag_TF] = -1 * rowSums(Qmat, na.rm=TRUE)
	
	Qmat = diag(basefreqs) %*% Qmat
	rowSums(Qmat)
	
	return(Qmat)
	}



# These are 1-event probability matrices, not instantaneous rate matrices
# Speciation matrix version
#######################################################
# symbolic_cell_to_relprob_cell_sp
#######################################################
#' Convert symbolic cell (a text equation) to relprob cell (a numeric value) -- speciation matrix version
#' 
#' This does the equivalent of \code{\link{symbolic_to_P_matrix}}, but for a speciation/cladogenesis matrix.
#' 
#' These are 1-event probability matrices, not instantaneous rate matrices.
#'
#' This function can be used in \code{\link[base]{sapply}}.  It still will not be very fast compared to the calculations in 
#' \code{\link[cladoRcpp]{cladoRcpp}}, but can be useful for demonstrative purposes.
#' 
#' @param charcell The text formula.
#' @param cellsplit The symbol to split the formulas on. Default "\\\\+" (plus symbol, with escape code).
#' @param mergesym The symbol to merge the formulas with. Default "+".
#' @param ys Relative weight of fully sympatric speciation (range-copying) and sympatric "subset" speciation. Default \code{s=1} mimics LAGRANGE model.
#' @param v Relative weight of vicariant speciation. Default \code{v=1} mimics LAGRANGE model.
#' @param j Relative weight of "founder event speciation"/jump speciation. Default \code{j=0} mimics LAGRANGE model.
#' @param relprob_subsets_matrix A numeric matrix describing the relative probability of each smaller daughter range, conditional on the ancestral rangesize.
#' @param relprob_vicar_matrix A numeric matrix describing the relative probability of each smaller daughter range, conditional on the ancestral rangesize.
#' @param ... Additional arguments to pass to \code{\link{relative_probabilities_of_subsets}} and \code{\link{relative_probabilities_of_vicariants}}, and
#' thence to \code{\link[base]{strsplit}}.
#' @return \code{cellval} The output cell value.
#' @export
#' @seealso \code{\link{symbolic_to_relprob_matrix_sp}}, \code{\link{make_relprob_matrix_de}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' 
#' charcell = "y1"
#' symbolic_cell_to_relprob_cell_sp(charcell, cellsplit="\\\\+", mergesym="*", ys=1, 
#' j=0, v=1, relprob_subsets_matrix=relative_probabilities_of_subsets(max_numareas=3, 
#' maxent_constraint_01=0.0001), 
#' relprob_vicar_matrix=relative_probabilities_of_vicariants(max_numareas=3, 
#' maxent_constraint_01v=0.0001))
#' 
#' charcell = "y1"
#' symbolic_cell_to_relprob_cell_sp(charcell, cellsplit="\\\\+", mergesym="*", ys=1, 
#' j=1, v=1, relprob_subsets_matrix=relative_probabilities_of_subsets(max_numareas=3, 
#' maxent_constraint_01=0.0001), 
#' relprob_vicar_matrix=relative_probabilities_of_vicariants(max_numareas=3, 
#' maxent_constraint_01v=0.0001))
#' 
#' charcell = "j"
#' symbolic_cell_to_relprob_cell_sp(charcell, cellsplit="\\\\+", mergesym="*", ys=1, 
#' j=0, v=1, relprob_subsets_matrix=relative_probabilities_of_subsets(max_numareas=3, 
#' maxent_constraint_01=0.0001), 
#' relprob_vicar_matrix=relative_probabilities_of_vicariants(max_numareas=3, 
#' maxent_constraint_01v=0.0001))
#' 
#' charcell = "j"
#' symbolic_cell_to_relprob_cell_sp(charcell, cellsplit="\\\\+", mergesym="*", ys=1, 
#' j=1, v=1, relprob_subsets_matrix=relative_probabilities_of_subsets(max_numareas=3, 
#' maxent_constraint_01=0.0001), 
#' relprob_vicar_matrix=relative_probabilities_of_vicariants(max_numareas=3, 
#' maxent_constraint_01v=0.0001))
#' 
#' charcell = "v1_2"
#' symbolic_cell_to_relprob_cell_sp(charcell, cellsplit="\\\\+", mergesym="*", ys=1, 
#' j=0, v=1, relprob_subsets_matrix=relative_probabilities_of_subsets(max_numareas=3, 
#' maxent_constraint_01=0.0001), 
#' relprob_vicar_matrix=relative_probabilities_of_vicariants(max_numareas=3, 
#' maxent_constraint_01v=0.0001))
#' 
#' charcell = "v1_2"
#' symbolic_cell_to_relprob_cell_sp(charcell, cellsplit="\\\\+", mergesym="*", ys=1,
#'  j=1, v=1, relprob_subsets_matrix=relative_probabilities_of_subsets(max_numareas=3, 
#' maxent_constraint_01=0.0001), 
#' relprob_vicar_matrix=relative_probabilities_of_vicariants(max_numareas=3, 
#' maxent_constraint_01v=0.0001))
#' 
#' charcell = "s1_2"
#' symbolic_cell_to_relprob_cell_sp(charcell, cellsplit="\\\\+", mergesym="*", ys=1, 
#' j=0, v=1, relprob_subsets_matrix=relative_probabilities_of_subsets(max_numareas=3, 
#' maxent_constraint_01=0.0001), 
#' relprob_vicar_matrix=relative_probabilities_of_vicariants(
#' max_numareas=3, 
#' maxent_constraint_01v=0.0001))
#' 
#' charcell = "s1_2"
#' symbolic_cell_to_relprob_cell_sp(charcell, cellsplit="\\\\+", mergesym="*", 
#' ys=1, j=1, v=1, relprob_subsets_matrix=relative_probabilities_of_subsets(
#' max_numareas=3, maxent_constraint_01=0.0001), 
#' relprob_vicar_matrix=relative_probabilities_of_vicariants(
#' max_numareas=3, 
#' maxent_constraint_01v=0.0001))
#' 
symbolic_cell_to_relprob_cell_sp <- function(charcell, cellsplit="\\+", mergesym="*", ys=1, j=0, v=1, relprob_subsets_matrix=relative_probabilities_of_subsets(6,0.0001), relprob_vicar_matrix=relative_probabilities_of_vicariants(6,0.0001), ...)
	{
	# 
	word = strsplit(charcell, split=cellsplit)[[1]]
	
	# If word DOES NOT have "s" or "y" or "v", just evaluate straight to a number
	if ( all( (grepl(pattern="s", x=word)==FALSE), (grepl(pattern="y", x=word)==FALSE), (grepl(pattern="v", x=word)==FALSE) ) )
		{
		# Evaluate straight to an answer (I guess, j only, here)
		equation = paste(word, sep="", collapse=mergesym)
		equation2 = paste("cellval = ", equation, sep="")
		eval(parse(text=equation2))
		# evaluation done
		} else {
		# If word includes "y" or "s", use sfunc or yfunc
		if (grepl(pattern="s", x=word)==TRUE)
			{
			cellval = ys * sfunc(charcell, relprob_subsets_matrix)
			}
		if (grepl(pattern="y", x=word)==TRUE)
			{
			cellval = ys * yfunc(charcell, relprob_subsets_matrix)
			}
# 		if (grepl(pattern="j", x=word)==TRUE)
# 			{
# 			cellval = j * jfunc(charcell, relprob_subsets_matrix)
# 			}
		if (grepl(pattern="v", x=word)==TRUE)
			{
			cellval = v * vfunc(charcell, relprob_vicar_matrix)
			}
		# evaluation done
		}

	return(cellval)
	}

# Convert symbolic matrix to relprob matrix
# These are 1-event probability matrices, not instantaneous rate matrices
# Speciation matrix version
#######################################################
# symbolic_to_relprob_matrix_sp
#######################################################
#' Convert symbolic matrix (with text equations) to relprob matrix (numeric values) -- speciation matrix version
#' 
#' This does the equivalent of \code{\link{symbolic_to_P_matrix}}, but for a speciation/cladogenesis matrix.
#' 
#' These are 1-event probability matrices, not instantaneous rate matrices.
#'
#' This function uses \code{\link{symbolic_cell_to_relprob_cell_sp}} in an \code{\link[base]{sapply}} call.  It still will not be very fast compared to the calculations in 
#' \code{\link[cladoRcpp]{cladoRcpp}}, but can be useful for demonstrative purposes.
#' 
#' @param spmat The speciation/cladogenesis matrix, with text formula.
#' @param cellsplit The symbol to split the formulas on. Default "\\\\+" (plus symbol, with escape code).
#' @param mergesym The symbol to merge the formulas with. Default "+".
#' @param ys Relative weight of fully sympatric speciation (range-copying) and sympatric "subset" speciation. Default \code{s=1} mimics LAGRANGE model.
#' @param v Relative weight of vicariant speciation. Default \code{v=1} mimics LAGRANGE model.
#' @param j Relative weight of "founder event speciation"/jump speciation. Default \code{j=0} mimics LAGRANGE model.
#' @param maxent_constraint_01 Parameter which assigns relative probabilities to different descendants range sizes, for the smaller descendant.  
#' Values can range from 0.0001 to 1. If \code{maxent_constraint_01=0.0001}, then the smaller descendant has a range size of 1 with probability 1
#' (i.e., the \code{LAGRANGE} default). If \code{maxent_constraint_01=0.5}, then all range sizes are equally weighted.  If \code{maxent_constraint_01=1},
#' then the largest possible smaller descendant gets probability 1.  The reference to "maxent" derives from the fact that the maxent probability
#' distribution on a multistate, ordered, discrete variable -- e.g. a die roll -- can be calculated given just the mean value.  Here, the
#' \code{maxent_constraint_01} parameter is multiplied by the (maximum rangesize + 1).  Thus, when \code{maxent_constraint_01=0.5}, if there are 6
#' possible states, then the parameter becomes 3.5, which sets equal probabilities of all possible descendant ranges sizes, when range size can range
#' from 1 to 6.
#' @param maxent_constraint_01v Works the same as \code{maxent_constraint_01}, but just for descendants of vicariant events.
#' @param max_numareas The maximum number of areas possible allowed for the smaller-ranged-daughter in either vicariant or sympatric types of cladogenesis/speciation.
#' @param ... Additional arguments to pass to \code{\link{relative_probabilities_of_subsets}} and \code{\link{relative_probabilities_of_vicariants}}, and
#' thence to \code{\link[base]{strsplit}}.
#' @return \code{cellval} The output cell value.
#' @export
#' @seealso \code{\link{symbolic_cell_to_relprob_cell_sp}}, \code{\link{make_relprob_matrix_de}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @examples
#' testval=1
#' # Generate the text version of the speciation/cladogenesis probability matrix 
#' # (actually a relative weights matrix
#' # until the rows are normalized so that each sums to 1).
#' spmat = make_relprob_matrix_bi(states_list=list("_", c("A"), c("B"), c("C"), 
#' c("A","B"), c("B","C"), c("A","C"), c("A","B","C")), split_ABC=FALSE, splitval="", 
#' code_for_overlapping_subsets=NA, printwarn=1)
#' spmat
#' 
#' # Look at the conditional probabilities generated by a variety of models
#' spPmat = symbolic_to_relprob_matrix_sp(spmat=spmat, cellsplit="\\\\+", mergesym="*", 
#' ys=1, j=0, v=1, maxent_constraint_01=0.0001, maxent_constraint_01v=0.0001, 
#' max_numareas=3)
#' spPmat = adf(spPmat); names(spPmat) = names(spmat); rownames(spPmat) = rownames(spmat)
#' spPmat
#' 
#' spPmat = symbolic_to_relprob_matrix_sp(spmat=spmat, cellsplit="\\\\+", 
#' mergesym="*", ys=0.5, j=0, v=0.5, maxent_constraint_01=0.0001, 
#' maxent_constraint_01v=0.0001, max_numareas=3)
#' spPmat = adf(spPmat); names(spPmat) = names(spmat); rownames(spPmat) = rownames(spmat)
#' spPmat
#' 
#' spPmat = symbolic_to_relprob_matrix_sp(spmat=spmat, cellsplit="\\\\+", 
#' mergesym="*", ys=1, j=1, v=1, maxent_constraint_01=0.0001, 
#' maxent_constraint_01v=0.0001, max_numareas=3)
#' spPmat = adf(spPmat); names(spPmat) = names(spmat); rownames(spPmat) = rownames(spmat)
#' spPmat
#' 
#' spPmat = symbolic_to_relprob_matrix_sp(spmat=spmat, cellsplit="\\\\+", 
#' mergesym="*", ys=0.25, j=0.25, v=0.25, maxent_constraint_01=0.0001, 
#' maxent_constraint_01v=0.0001, max_numareas=3)
#' spPmat = adf(spPmat); names(spPmat) = names(spmat); rownames(spPmat) = rownames(spmat)
#' spPmat
#' 
#' spPmat = symbolic_to_relprob_matrix_sp(spmat=spmat, cellsplit="\\\\+", 
#' mergesym="*", ys=1, j=1, v=0, maxent_constraint_01=0.0001, 
#' maxent_constraint_01v=0.0001, max_numareas=3)
#' spPmat = adf(spPmat); names(spPmat) = names(spmat); rownames(spPmat) = rownames(spmat)
#' spPmat
#' 
#' spPmat = symbolic_to_relprob_matrix_sp(spmat=spmat, cellsplit="\\\\+", 
#' mergesym="*", ys=1, j=1, v=0, maxent_constraint_01=0.5, 
#' maxent_constraint_01v=0.0001, max_numareas=3)
#' spPmat = adf(spPmat); names(spPmat) = names(spmat); rownames(spPmat) = rownames(spmat)
#' spPmat
#' 
#' spPmat = symbolic_to_relprob_matrix_sp(spmat=spmat, cellsplit="\\\\+", 
#' mergesym="*", ys=1, j=0, v=0, maxent_constraint_01=0.5, 
#' maxent_constraint_01v=0.0001, max_numareas=3)
#' spPmat = adf(spPmat); names(spPmat) = names(spmat); rownames(spPmat) = rownames(spmat)
#' spPmat
#' 
#' spPmat = symbolic_to_relprob_matrix_sp(spmat=spmat, cellsplit="\\\\+", 
#' mergesym="*", ys=1, j=0, v=1, maxent_constraint_01=0.0001, 
#' maxent_constraint_01v=0.5, max_numareas=3)
#' spPmat = adf(spPmat); names(spPmat) = names(spmat); rownames(spPmat) = rownames(spmat)
#' spPmat
#' 
symbolic_to_relprob_matrix_sp <- function(spmat, cellsplit="\\+", mergesym="*", ys=1, j=0, v=1, maxent_constraint_01=0.0001, maxent_constraint_01v=0.0001, max_numareas=6, ...)
	{
	defaults='
	cellsplit="\\+"
	mergesym="*"
	ys=0.2
	j=0.1
	v=0.03
	maxent_constraint_01=0.5
	max_numareas
	'
	
	# Calculate the relative probability of each rangesize
	relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=max_numareas, maxent_constraint_01=maxent_constraint_01, ...)
	
	if (max_numareas > 1)
		{
		relprob_vicar_matrix = relative_probabilities_of_vicariants(max_numareas=max_numareas, maxent_constraint_01v=maxent_constraint_01v, ...)
		} else {
		relprob_vicar_matrix = matrix(NA, nrow=1, ncol=1)
		}
	
	spmat[spmat == "-"] = NA
	
	spmat_vals = matrix(data=0, nrow=nrow(spmat), ncol=ncol(spmat))
	
	spmat_vals[!is.na(spmat_vals)] = sapply(X=spmat[!is.na(spmat_vals)], FUN=symbolic_cell_to_relprob_cell_sp, cellsplit=cellsplit, mergesym=mergesym, ys=ys, j=j, v=v, relprob_subsets_matrix=relprob_subsets_matrix, relprob_vicar_matrix=relprob_vicar_matrix)
	
	# DO NOT normalize, this is a single-event probability matrix,
	# not an instantaneous rate matrix
	# diag_TF = diag(x=1, nrow=nrow(spmat), ncol(spmat)) == 1
	# spmat_vals[diag_TF] = -1 * rowSums(spmat_vals, na.rm=TRUE)
	
	# But, it SHOULD be normalized so each row adds to 1
	# NOTE: THIS WILL TEND TO MEAN THAT EVENTS THAT OCCUPY MORE 
	# CELLS WILL GET MORE WEIGHT FOR THE SAME OVERALL SB RATE
	
	# The rows with "_" will get 0
# 	rowsums = rowSums(spmat_vals)
# 	rowsums[rowsums == 0] = 1
# 	spmat_vals = spmat_vals/rowsums
	
	spmat_vals = spmat_vals/rowSums(spmat_vals)
	rowSums(spmat_vals)
	
	return(spmat_vals)
	}



# Fun with Maximum Entropy probability distribution for discrete variable with given mean
# (and discrete uniform flat prior)
# http://en.wikipedia.org/wiki/Maximum_entropy_probability_distribution
# 
# Equation 6.4 of Harte 2011

#######################################################
# calcZ_part
#######################################################
#' Calculate Z (equation 6.3 of Harte 2011)
#' 
#' This function is a used by \code{\link{calcP_n}} via \code{\link[base]{apply}}, all within \code{\link{get_probvals}}. \code{\link{get_probvals}}
#' calculates the Maximum Entropy (\cite{Harte2011}) discrete probability distribution of a number of ordered states (e.g., faces of a 6-sided 
#' die) given the mean of many rolls.  Here, this is merely used so that a single parameter can control the probability distribution of small
#' versus large descendant areas during cladogenesis.
#' 
#' See also: Maximum Entropy probability distribution for discrete variable with given mean
#' (and discrete uniform flat prior)
#' \url{http://en.wikipedia.org/wiki/Maximum_entropy_probability_distribution}
#' 
#' @param n Value of the state (e.g., which of a number of faces on a die, or number of different size classes of geographic range)
#' @param lambda1 Lambda parameter (\cite{Harte2011}).
#' @return \code{Z}, numeric value
#' @export
#' @seealso \code{\link{calcP_n}}, \code{\link[FD]{maxent}}, \code{\link{symbolic_to_relprob_matrix_sp}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://en.wikipedia.org/wiki/Maximum_entropy_probability_distribution}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Harte2011
#' @examples
#' testval=1
#' n=6
#' lambda1 = 0.5
#' calcZ_part(n, lambda1)
#'
calcZ_part <- function(n, lambda1)
	{
	return(exp(-1*lambda1*n))
	}



# Harte 2011, equation 6.3
# Calc prob(n)
#######################################################
# calcP_n
#######################################################
#' Calculate Z (part of equation 6.4 of Harte 2011)
#' 
#' This function is a used by \code{\link{get_probvals}}, which 
#' calculates the Maximum Entropy (\cite{Harte2011}) discrete probability distribution of a number of ordered states (e.g., faces of a 6-sided 
#' die) given the mean of many rolls.  Here, this is merely used so that a single parameter can control the probability distribution of small
#' versus large descendant areas during cladogenesis.
#' 
#' See also: Maximum Entropy probability distribution for discrete variable with given mean
#' (and discrete uniform flat prior)
#' \url{http://en.wikipedia.org/wiki/Maximum_entropy_probability_distribution}
#' 
#' @param n Value of the state (e.g., which of a number of faces on a die, or number of different size classes of geographic range).
#' @param lambda1 Lambda parameter (\cite{Harte2011}).
#' @param Z numeric values from \code{\link{calcZ_part}}.
#' @return \code{Prob_n}, numeric value of the probability of state \code{n}.
#' @export
#' @seealso \code{\link{calcZ_part}}, \code{\link[FD]{maxent}}, \code{\link{symbolic_to_relprob_matrix_sp}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://en.wikipedia.org/wiki/Maximum_entropy_probability_distribution}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Harte2011
#' @examples
#' testval=1
#' n = 6
#' lambda1 = 0.5
#' Z = 1
#' calcP_n(n, lambda1, Z)
#'
calcP_n = function(n, lambda1, Z)
	{
	Prob_n = (1/Z) * exp(-1*lambda1*n)
	return(Prob_n)
	}

	
	
	

# Set up the function based on
# Harte (2011), Maximum Entropy Theory of Ecology
#

#######################################################
# get_probvals
#######################################################
#' Calculate probability of ordered discrete states using a maxent distribution (equations 6.3-6.4 of Harte 2011)
#' 
#' This function is calculates the Maximum Entropy (\cite{Harte2011}) discrete probability distribution of a number of ordered states (e.g., faces
#' of a 6-sided die) given the mean of many rolls.  Here, this is merely used so that a single parameter can control the probability distribution
#' of small versus large descendant areas during cladogenesis.  This function could then used by \code{\link{relative_probabilities_of_subsets}} in BioGeoBEARS
#' to weight different descendant range sizes (although, currently, the function \code{\link[FD]{maxent}} from the \code{\link[FD]{FD}} package is used).
#'
#' This calculation is based on Equations 6.3-6.4 of \cite{Harte2011}.
#' 
#' See also: Maximum Entropy probability distribution for discrete variable with given mean
#' (and discrete uniform flat prior)
#' \url{http://en.wikipedia.org/wiki/Maximum_entropy_probability_distribution}
#' 
#' @param die_vals Values of the ordered discrete variable state (e.g., \code{seq(1,6)} for a six-sided die) 
#' @param meanval Mean value (the knowledge supplied to the MaxEnt function).
#' @return \code{Prob_nvals}, numeric values of the probability of each state from \code{die_vals}.
#' @export
#' @seealso \code{\link{calcZ_part}}, \code{\link{calcP_n}}, \code{\link[FD]{maxent}}, \code{\link{symbolic_to_relprob_matrix_sp}}, \code{\link{relative_probabilities_of_subsets}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://en.wikipedia.org/wiki/Maximum_entropy_probability_distribution}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Harte2011
#' @examples
#' testval=1
#' # Examples
#' # Set up subplots
#' par(mfrow=c(3,2))
#' 
#' # Flat distribution (equal prob of any descendent size)
#' N = 6
#' # n = die_vals
#' die_vals = seq(1,N)
#' meanval = 3.5
#' probvals = get_probvals(die_vals, meanval)
#' probvals
#' barplot(height=probvals, width=1, names.arg=die_vals, ylim=c(0,1))
#' title(paste("Probabilities of each state, mean val=", meanval, sep=""))
#' 
#' # Descendents tend to have large ranges
#' N = 6
#' # n = die_vals
#' die_vals = seq(1,N)
#' meanval = 5.999
#' probvals = get_probvals(die_vals, meanval)
#' probvals
#' barplot(height=probvals, width=1, names.arg=die_vals, ylim=c(0,1))
#' title(paste("Probabilities of each state, mean val=", meanval, sep=""))
#' 
#' # Flat distribution (equal prob of any descendent size)
#' N = 6
#' # n = die_vals
#' die_vals = seq(1,N)
#' meanval = 5
#' probvals = get_probvals(die_vals, meanval)
#' probvals
#' barplot(height=probvals, width=1, names.arg=die_vals, ylim=c(0,1))
#' title(paste("Probabilities of each state, mean val=", meanval, sep=""))
#' 
#' # Flat distribution (equal prob of any descendent size)
#' N = 6
#' # n = die_vals
#' die_vals = seq(1,N)
#' meanval = 4
#' probvals = get_probvals(die_vals, meanval)
#' probvals
#' barplot(height=probvals, width=1, names.arg=die_vals, ylim=c(0,1))
#' title(paste("Probabilities of each state, mean val=", meanval, sep=""))
#' 
#' # Flat distribution (equal prob of any descendent size)
#' N = 6
#' # n = die_vals
#' die_vals = seq(1,N)
#' meanval = 2
#' probvals = get_probvals(die_vals, meanval)
#' probvals
#' barplot(height=probvals, width=1, names.arg=die_vals, ylim=c(0,1))
#' title(paste("Probabilities of each state, mean val=", meanval, sep=""))
#' 
#' # This produces the LAGRANGE default
#' # (all smaller descendents are of size 1)
#' N = 6
#' # n = die_vals
#' die_vals = seq(1,N)
#' meanval = 1.0001
#' probvals = get_probvals(die_vals, meanval)
#' probvals
#' barplot(height=probvals, width=1, names.arg=die_vals, ylim=c(0,1))
#' title(paste("LAGRANGE 'default', mean val=", meanval, sep=""))
#' 
#' # This is stopped by the error check
#' # (all smaller descendents are of size 1)
#' # N = 6
#' # # n = die_vals
#' # die_vals = seq(1,N)
#' # meanval = 0.5
#' # probvals = get_probvals(die_vals, meanval)
#' # probvals
#' # barplot(height=probvals, width=1, names.arg=die_vals, ylim=c(0,1))
#' # title(paste("Probabilities of each state, mean val=", meanval, sep=""))
#' 
get_probvals <- function(die_vals, meanval)
	{
	# Take Harte (2011), equation 6.4; set exp(-lambda) to == x
	
	defaults='
	# n = number of integers
	N = 6
	# n = die_vals
	die_vals = seq(1,N)
	meanval = 5.999
	'
	
	# Error check
	if (meanval <= min(die_vals) || meanval >= max(die_vals))
		{
		error_msg = paste("ERROR: meanval (", meanval, ") cannot be below min (", min(die_vals), ") or above max (", max(die_vals), ") die_vals!", sep="")
		stop(error_msg)
		}
	
	
	
	Harte_2011_eqn_6_9_left_coefficients = die_vals
	Harte_2011_eqn_6_9_right_coefficients = rep(meanval, length(die_vals))
	polynomial_coefs_that_EQ_0 = Harte_2011_eqn_6_9_left_coefficients - Harte_2011_eqn_6_9_right_coefficients
	
	solutions_for_x = polyroot(z=polynomial_coefs_that_EQ_0)
	solutions_for_x
	
	# Exclude those with nonzero imaginary values
	keepTF = abs(Im(solutions_for_x)) < 1e-14	
	solutions_for_x = solutions_for_x[keepTF]
	solutions_for_x
	
	# Exclude negative and 0 solutions
	keepTF = Re(solutions_for_x) > 1e-14	
	solutions_for_x = Re(solutions_for_x[keepTF])
	solutions_for_x
	
	# Lagrange multiplier Lambda1
	# lambda = -1 * ln(x)
	lambda1 = -1 * log(solutions_for_x)
	lambda1
	
	# Equation 6.4
	# 		calcZ_part <- function(n, lambda1)
	# 			{
	# 			return(exp(-1*lambda1*n))
	# 			}
	
	Z = sum(sapply(X=die_vals, FUN=calcZ_part, lambda1=lambda1))
	Z
	
	
	# Harte 2011, equation 6.3
	# Calc prob(n)
	# 		calcP_n = function(n, lambda1, Z)
	# 			{
	# 			Prob_n = (1/Z) * exp(-1*lambda1*n)
	# 			return(Prob_n)
	# 			}
	
	Prob_nvals = sapply(X=die_vals, FUN=calcP_n, lambda1=lambda1, Z=Z)
	Prob_nvals
	sum(Prob_nvals)
	
	
	
	return(Prob_nvals)
	}



#######################################################
# relative_probabilities_of_subsets
#######################################################
#' Calculate probability of different descendant rangesizes, for the smaller descendant, in subset speciation
#'
#' "Rangesize" here means "number of areas in a geographic range".  The \code{LAGRANGE} cladogenesis model requires that, during cladogenesis events, one daughter
#' lineage will ALWAYS have a geographic range of size 1.  This is argued for in \cite{ReeSmith2008} on the grounds that new species usually get isolated and start
#' in a new area.  This is a reasonable proposition, but still, it would be nice to test the assumption.  In addition, it could be that 
#' some speciation modes, especially vicariance, obey different rules.  E.g., \code{DIVA} (\cite{Ronquist1996_DIVA}, \cite{Ronquist_1997_DIVA})
#' allows vicariant speciation to divide up the ancestral range in every possible way (e.g., ABCD-->AB|CD, or AC|BD, or A|BCD, or D|ABC, etc.),
#' but \code{LAGRANGE} would only allow vicariance to split off areas of size 1: (ABCD-->A|BCD, B|ACD, etc.) (Ronquist_Sanmartin_2011).
#' 
#' To test different models, the user has to have control of the relative probability of different descendant rangesizes.  The probability of each
#' descendant rangesize could be parameterized individually, but we have a limited amount of observational data (essentially one character), so
#' efficient parameterizations should be sought.  
#'
#' One way to do this is with the Maximum Entropy (\cite{Harte2011}) discrete probability distribution
#' of a number of ordered states.  Normally this is applied (in examples) to the problem of estimation of the relative probability of the different
#' faces of a 6-sided die.  The input "knowledge" is the true mean of the dice rolls.  If the mean value is 3.5, then each face of the die will have 
#' probability 1/6.  If the mean value is close to 1, then the die is severely skewed such that the probability of rolling 1 is 99%+ and the probability
#' of other die rolls is very small.  If the mean value is close to 6, then the probability distribution is skewed towards higher numbers.
#' 
#' Here in \code{BioGeoBEARS}, we use the same Maximum Entropy function to specify the relative probability of geographic ranges of a number of different
#' rangesizes.  This is merely used so that a single parameter can control the probability distribution -- there is no MaxEnt estimation going
#' on here.  The user specifies a value for the parameter \code{maxent_constraint_01} between 0.0001 and 0.9999.  This can then be applied to 
#' all of the different ancestor-descendant range combinations in the cladogenesis/speciation matrix.  
#' 
#' Example values of \code{maxent_constraint_01} would give the following results:
#'
#' \code{maxent_constraint_01 = 0.0001} -- The smaller descendant has rangesize 1 with 100% probability (as in \code{LAGRANGE})\cr
#' \code{maxent_constraint_01 = 0.5} -- The smaller descendant can be any rangesize equal probability. This is effectively what happens in \code{DIVA}'s version of vicariance speciation\cr
#' \code{maxent_constraint_01 = 0.9999} -- The smaller descendant will take the largest possible rangesize for a given type of speciation, and a given ancestral
#' rangesize.  E.g., for sympatric/range-copying speciation (the ancestor is simply copied to both descendants, as in a continuous-time model with no
#' cladogenesis effect), an ancestor of size 3 would product two descendant lineages of size 3. Such a model is implemented in the program \code{BayArea} (\cite{Landis_Matzke_etal_2013_BayArea}). 
#' \code{LAGRANGE}, on the other hand, would only allow range-copying for ancestral ranges of size 1.\cr
#'
#' \bold{Note:} In \code{LAGRANGE}-type models, at speciation/cladogenesis events,
#' one descendant daughter branch ALWAYS has size 1, whereas the other descendant daughter branch either (a) is the same (in sympatric/range-copying
#' speciation), (b) inherits the complete ancestral range (in sympatric/subset speciation) or (c) inherits the remainder of the range (in 
#' vicariant/range-division speciation).  LAGRANGE-type behavior (the smaller descendant has rangesize 1 with 100% probability, and 0% probability for
#' any other rangesize) can be achieved by setting the \code{maxent_constraint_01} parameter to 0.0001.
#' 
#' \code{See also:} Maximum Entropy probability distribution for discrete variable with given mean
#' (and discrete uniform flat prior)
#' \url{http://en.wikipedia.org/wiki/Maximum_entropy_probability_distribution}
#' 
#' Currently, the function \code{\link[FD]{maxent}} from the \code{\link[FD]{FD}} package is used to get the discrete probability distribution, given the number of 
#' states and the \code{maxent_constraint_01} parameter.  This could also be done with \code{\link{get_probvals}}, which uses \code{\link{calcZ_part}},
#' \code{\link{calcP_n}}, following equations 6.3-6.4 of \cite{Harte2011}, although this is not yet implemented.
#' 
#' @param max_numareas The maximum number of areas possible allowed for the smaller-ranged-daughter in this type of cladogenesis/speciation.
#' @param maxent_constraint_01 The parameter describing the probability distribution on descendant rangesizes for the smaller descendant. See above.
#' @param NA_val The output matrix consists of ancestral rangesizes and rangesizes of the smaller descendant.  Some values are disallowed -- e.g. descendant ranges
#' larger than the ancestor; or, in subset speciation, descendant ranges the same size as the ancestor are disallowed.  All disallowed descendant rangesizes get
#' \code{NA_val}.
#' @return \code{\link{relative_probabilities_of_vicariants}}, \code{relprob_subsets_matrix}, a numeric matrix giving the relative probability of each rangesize for the smaller descendant of an ancestral range, 
#' conditional on the ancestral rangesize.
#' @export
#' @seealso \code{\link{symbolic_to_relprob_matrix_sp}}, \code{\link{get_probvals}}, \code{\link[FD]{maxent}}, \code{\link{calcZ_part}}, \code{\link{calcP_n}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://en.wikipedia.org/wiki/Maximum_entropy_probability_distribution}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Harte2011
#'	 @cite ReeSmith2008
#'	 @cite Ronquist1996_DIVA
#'	 @cite Ronquist_1997_DIVA
#'	 @cite Ronquist_Sanmartin_2011
#'	 @cite Landis_Matzke_etal_2013_BayArea
#' @examples
#' testval=1
#' # Examples
#' 
#' # Probabilities of different descendant rangesizes, for the smaller 
#' # descendant, under sympatric/subset speciation
#' # (plus sympatric/range-copying, which is folded in):
#' relative_probabilities_of_subsets(max_numareas=6, maxent_constraint_01=0.0001, 
#' NA_val=NA)
#' relative_probabilities_of_subsets(max_numareas=6, maxent_constraint_01=0.5, 
#' NA_val=NA)
#' relative_probabilities_of_subsets(max_numareas=6, maxent_constraint_01=0.9999, 
#' NA_val=NA)
#' 
#' # Probabilities of different descendant rangesizes, for the smaller descendant, 
#' # under vicariant speciation
#' relative_probabilities_of_vicariants(max_numareas=6, maxent_constraint_01v=0.0001, 
#' NA_val=NA)
#' relative_probabilities_of_vicariants(max_numareas=6, maxent_constraint_01v=0.5, 
#' NA_val=NA)
#' relative_probabilities_of_vicariants(max_numareas=6, maxent_constraint_01v=0.9999, 
#' NA_val=NA)
#' 
relative_probabilities_of_subsets <- function(max_numareas=6, maxent_constraint_01=0.5, NA_val=NA)
	{
	defaults='
	
	# DEFAULTS	
	max_numareas=6
	maxent_constraint=0.01
	
	# The user provides a relative weight between 0 and 1.
	# This goes into quantile(), on a series of values:
	# 0, 1...max_numareas, max_numareas+1
	#
	# This means that 0.5 gives the median as a constraint,
	# and maxent() will produce a flat categorical distribution
	
	# This default puts all the weight on the "subset" species being
	# sympatric
	maxent_constraint_01 = 1
	maxent_constraint = quantile(x=seq(0,max_numareas+1,1), probs=maxent_constraint_01)
	maxent_constraint

	# This default puts all the weight on the 1-area subset species
	maxent_constraint_01 = 0.0001
	maxent_constraint= quantile(x=seq(0,max_numareas+1,1), probs=maxent_constraint_01)
	maxent_constraint

	# This default produces even weights (picks the median)
	maxent_constraint_01 = 0.5
	maxent_constraint = quantile(x=seq(0,max_numareas+1,1), probs=maxent_constraint_01)
	maxent_constraint
	'
	
	# Require FD for maxent function
	# FD::maxent
	#require(FD)

	# rows = number of areas in ancestor
	# cols = number of areas in subset daughter
	
	relprob_subsets_matrix = matrix(NA_val, nrow=max_numareas, ncol=max_numareas)
	
	# rows
	rownum=max_numareas
	for (rownum in 1:max_numareas)
		{
		numcols = rownum
		tmpstates = seq(1, numcols, by=1)

		# Convert the user-specified 0-1 value into a MaxEnt constraint
		# between 0 and max_numareas+1
		# 0.5 = median
		maxent_constraint = quantile(x=seq(0,length(tmpstates)+1,1), probs=maxent_constraint_01)

		# Apply Maxent constraint to weight the different numbers of areas
		maxent_result = maxent(constr=maxent_constraint, states=tmpstates)
		probs_of_subset_ranges = maxent_result$prob
		probs_of_subset_ranges
		
		# Add NAs to get the correct rowlength
		numzeros_to_add = max_numareas-length(probs_of_subset_ranges)
		probs_of_subset_ranges = c(probs_of_subset_ranges, rep.int(NA_val, times=numzeros_to_add))

		relprob_subsets_matrix[rownum, ] = probs_of_subset_ranges
		}
		
	# Round to 3 significant digits
	relprob_subsets_matrix = round(relprob_subsets_matrix, 3)
	relprob_subsets_matrix
	
	return(relprob_subsets_matrix)
	}



#######################################################
# relative_probabilities_of_vicariants
#######################################################
#' Calculate probability of different descendant rangesizes, for the smaller descendant, in vicariant speciation
#'
#' "Rangesize" here means "number of areas in a geographic range".  The \code{LAGRANGE} cladogenesis model requires that, during cladogenesis events, one daughter
#' lineage will ALWAYS have a geographic range of size 1.  This is argued for in \cite{ReeSmith2008} on the grounds that new species usually get isolated and start
#' in a new area.  This is a reasonable proposition, but still, it would be nice to test the assumption.  In addition, it could be that 
#' some speciation modes, especially vicariance, obey different rules.  E.g., \code{DIVA} (\cite{Ronquist1996_DIVA}, \cite{Ronquist_1997_DIVA})
#' allows vicariant speciation to divide up the ancestral range in every possible way (e.g., ABCD-->AB|CD, or AC|BD, or A|BCD, or D|ABC, etc.),
#' but \code{LAGRANGE} would only allow vicariance to split off areas of size 1: (ABCD-->A|BCD, B|ACD, etc.) (Ronquist_Sanmartin_2011).
#' 
#' To test different models, the user has to have control of the relative probability of different descendant rangesizes.  The probability of each
#' descendant rangesize could be parameterized individually, but we have a limited amount of observational data (essentially one character), so
#' efficient parameterizations should be sought.  
#'
#' One way to do this is with the Maximum Entropy (\cite{Harte2011}) discrete probability distribution
#' of a number of ordered states.  Normally this is applied (in examples) to the problem of estimation of the relative probability of the different
#' faces of a 6-sided die.  The input "knowledge" is the true mean of the dice rolls.  If the mean value is 3.5, then each face of the die will have 
#' probability 1/6.  If the mean value is close to 1, then the die is severely skewed such that the probability of rolling 1 is 99%+ and the probability
#' of other die rolls is very small.  If the mean value is close to 6, then the probability distribution is skewed towards higher numbers.
#' 
#' Here in \code{BioGeoBEARS}, we use the same Maximum Entropy function to specify the relative probability of geographic ranges of a number of different
#' rangesizes.  This is merely used so that a single parameter can control the probability distribution -- there is no MaxEnt estimation going
#' on here.  The user specifies a value for the parameter \code{maxent_constraint_01} between 0.0001 and 0.9999.  This can then be applied to 
#' all of the different ancestor-descendant range combinations in the cladogenesis/speciation matrix.  
#' 
#' Example values of \code{maxent_constraint_01} would give the following results:
#'
#' \code{maxent_constraint_01 = 0.0001} -- The smaller descendant has rangesize 1 with 100% probability (as in \code{LAGRANGE})\cr
#' \code{maxent_constraint_01 = 0.5} -- The smaller descendant can be any rangesize equal probability. This is effectively what happens in \code{DIVA}'s version of vicariance speciation\cr
#' \code{maxent_constraint_01 = 0.9999} -- The smaller descendant will take the largest possible rangesize for a given type of speciation, and a given ancestral
#' rangesize.  E.g., for sympatric/range-copying speciation (the ancestor is simply copied to both descendants, as in a continuous-time model with no
#' cladogenesis effect), an ancestor of size 3 would product two descendant lineages of size 3. Such a model is implemented in the program \code{BayArea} (\cite{Landis_Matzke_etal_2013_BayArea}). 
#' \code{LAGRANGE}, on the other hand, would only allow range-copying for ancestral ranges of size 1.\cr
#'
#' \bold{Note:} In \code{LAGRANGE}-type models, at speciation/cladogenesis events,
#' one descendant daughter branch ALWAYS has size 1, whereas the other descendant daughter branch either (a) is the same (in sympatric/range-copying
#' speciation), (b) inherits the complete ancestral range (in sympatric/subset speciation) or (c) inherits the remainder of the range (in 
#' vicariant/range-division speciation).  LAGRANGE-type behavior (the smaller descendant has rangesize 1 with 100% probability, and 0% probability for
#' any other rangesize) can be achieved by setting the \code{maxent_constraint_01} parameter to 0.0001.
#' 
#' \code{See also:} Maximum Entropy probability distribution for discrete variable with given mean
#' (and discrete uniform flat prior)
#' \url{http://en.wikipedia.org/wiki/Maximum_entropy_probability_distribution}
#' 
#' Currently, the function \code{\link[FD]{maxent}} from the \code{\link[FD]{FD}} package is used to get the discrete probability distribution, given the number of 
#' states and the \code{maxent_constraint_01} parameter.  This could also be done with \code{\link{get_probvals}}, which uses \code{\link{calcZ_part}},
#' \code{\link{calcP_n}}, following equations 6.3-6.4 of \cite{Harte2011}, although this is not yet implemented.
#' 
#' @param max_numareas The maximum number of areas possible allowed for the smaller-ranged-daughter in this type of cladogenesis/speciation.
#' @param maxent_constraint_01v The parameter describing the probability distribution on descendant rangesizes for the smaller descendant, in a vicariance event
#' (where the maximum size of the smaller range is numareas/2, rounded down). See above.
#' @param NA_val The output matrix consists of ancestral rangesizes and rangesizes of the smaller descendant.  Some values are disallowed -- e.g. descendant ranges
#' larger than the ancestor; or, in subset speciation, descendant ranges the same size as the ancestor are disallowed.  All disallowed descendant rangesizes get
#' \code{NA_val}.
#' @return \code{relprob_subsets_matrix}, a numeric matrix giving the relative probability of each rangesize for the smaller descendant of an ancestral range, 
#' conditional on the ancestral rangesize.
#' @export
#' @seealso \code{\link{relative_probabilities_of_subsets}}, \code{\link{symbolic_to_relprob_matrix_sp}}, \code{\link{get_probvals}}, \code{\link[FD]{maxent}}, \code{\link{calcZ_part}}, \code{\link{calcP_n}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://en.wikipedia.org/wiki/Maximum_entropy_probability_distribution}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Harte2011
#'	 @cite ReeSmith2008
#'	 @cite Ronquist1996_DIVA
#'	 @cite Ronquist_1997_DIVA
#'	 @cite Ronquist_Sanmartin_2011
#'	 @cite Landis_Matzke_etal_2013_BayArea
#' @examples
#' testval=1
#' # Examples
#' 
#' # Probabilities of different descendant rangesizes, for the smaller descendant, 
#' # under sympatric/subset speciation
#' # (plus sympatric/range-copying, which is folded in):
#' relative_probabilities_of_subsets(max_numareas=6, maxent_constraint_01=0.0001, 
#' NA_val=NA)
#' relative_probabilities_of_subsets(max_numareas=6, maxent_constraint_01=0.5, 
#' NA_val=NA)
#' relative_probabilities_of_subsets(max_numareas=6, maxent_constraint_01=0.9999, 
#' NA_val=NA)
#' 
#' # Probabilities of different descendant rangesizes, for the smaller descendant, 
#' # under vicariant speciation
#' relative_probabilities_of_vicariants(max_numareas=6, maxent_constraint_01v=0.0001, 
#' NA_val=NA)
#' relative_probabilities_of_vicariants(max_numareas=6, maxent_constraint_01v=0.5, 
#' NA_val=NA)
#' relative_probabilities_of_vicariants(max_numareas=6, maxent_constraint_01v=0.9999, 
#' NA_val=NA)
#' 
relative_probabilities_of_vicariants <- function(max_numareas=6, maxent_constraint_01v=0.0001, NA_val=NA)
	{
	defaults='
	max_numareas=6
	maxent_constraint=0.01
	
	# The user provides a relative weight between 0 and 1.
	# This goes into quantile(), on a series of values:
	# 0, 1...max_numareas, max_numareas+1
	#
	# This means that 0.5 gives the median as a constraint,
	# and maxent() will produce a flat categorical distribution
	
	# This default means that the only vicariance event allowed is
	# splitting off a single range
	maxent_constraint_01v = 0.0001
	maxent_constraint_v = quantile(x=seq(0,max_numareas+1,1), probs=maxent_constraint_01v)
	maxent_constraint_v
	'
	
	# Require FD for maxent function
	# FD::maxent
	#require(FD)

	# rows = number of areas in ancestor
	# cols = number of areas in subset daughter
	
	relprob_vicar_matrix = matrix(NA_val, nrow=max_numareas, ncol=max_numareas)
	
	# rows
	rownum=max_numareas

	# Row 1 (ancestor size 1) has to be NA for any vicariance event
	# (redundant with NA declaration of relprob_vicar_matrix;
	#  basically an error check)
	relprob_vicar_matrix[1,] = NA_val

	for (rownum in 2:max_numareas)
		{
		numcols = rownum
		tmpstates = seq(1, numcols, by=1)

		# The maximum smallest rangesize is
		# whatever ranges are SMALLER than the median
		# i.e., for 1-6, median is 3.5, ranges 1-3 are allowed for the smaller range
		# i.e., for 1-5, median is 3, ranges 1-2 are allowed for the smaller range
		max_smaller_rangesize = median(tmpstates)

		possible_vicariance_smaller_rangesizes = tmpstates[tmpstates < max_smaller_rangesize]

		# Convert the user-specified 0-1 value into a MaxEnt constraint
		# between 0 and max_numareas+1
		# 0.5 = median
		# Just do +0, since our state space of possible vicariants is only -- NO!
		# Do +1, so that we get even probabilities at maxent_constraint_01v=0.5
		maxent_constraint = quantile(x=seq(0,length(possible_vicariance_smaller_rangesizes)+1,1), probs=maxent_constraint_01v)
		
		# Apply Maxent constraint to weight the different numbers of areas
		maxent_result = maxent(constr=maxent_constraint, states=possible_vicariance_smaller_rangesizes)
		probs_of_subset_ranges = maxent_result$prob
		probs_of_subset_ranges
		
		# Add NAs to get the correct rowlength
		numzeros_to_add = max_numareas-length(probs_of_subset_ranges)
		probs_of_subset_ranges = c(probs_of_subset_ranges, rep.int(NA_val, times=numzeros_to_add))

		relprob_vicar_matrix[rownum, ] = probs_of_subset_ranges
		}
		
	# Round to 3 significant digits
	relprob_vicar_matrix = round(relprob_vicar_matrix, 3)
	relprob_vicar_matrix
	
	
	return(relprob_vicar_matrix)
	}





# Function giving the relative probability of s (subset speciation)
# given different range sizes of the subset descendent
#######################################################
# sfunc
#######################################################
#' Extract the appropriate probability for a subset speciation event, given text code for rangesize of smaller descendant, and ancestor
#'
#' @param charcell The text in the cell, indicating the type of speciation/cladogenesis range inheritance event.
#' @param relprob_subsets_matrix A numeric matrix describing the relative probability of each smaller daughter range, conditional on the ancestral rangesize.
#' @return \code{prob_of_this_b}, a numeric value giving the relative probability of that descendent-ancestor rangesize pair.
#' @export
#' @seealso \code{\link{yfunc}}, \code{\link{vfunc}}, \code{\link{relative_probabilities_of_subsets}}, \code{\link{symbolic_to_relprob_matrix_sp}}, \code{\link{get_probvals}}, \code{\link[FD]{maxent}}, \code{\link{calcZ_part}}, \code{\link{calcP_n}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://en.wikipedia.org/wiki/Maximum_entropy_probability_distribution}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Harte2011
#'	 @cite ReeSmith2008
#'	 @cite Ronquist1996_DIVA
#'	 @cite Ronquist_1997_DIVA
#'	 @cite Ronquist_Sanmartin_2011
#'	 @cite Landis_Matzke_etal_2013_BayArea
#' @examples
#' testval=1
#' # Examples
#' 
#' # Probabilities of different descendant rangesizes, for the smaller descendant, 
#' # under sympatric/subset speciation
#' # (plus sympatric/range-copying, which is folded in):
#' relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, 
#' maxent_constraint_01=0.0001, NA_val=NA)
#' relprob_subsets_matrix
#' sfunc(charcell="s1_1", relprob_subsets_matrix)
#' sfunc(charcell="s1_2", relprob_subsets_matrix)
#' sfunc(charcell="s1_3", relprob_subsets_matrix)
#' sfunc(charcell="s2_3", relprob_subsets_matrix)
#' 
#' relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, 
#' maxent_constraint_01=0.5, NA_val=NA)
#' relprob_subsets_matrix
#' sfunc(charcell="s1_1", relprob_subsets_matrix)
#' sfunc(charcell="s1_2", relprob_subsets_matrix)
#' sfunc(charcell="s1_3", relprob_subsets_matrix)
#' sfunc(charcell="s2_3", relprob_subsets_matrix)
#' 
#' relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, 
#' maxent_constraint_01=0.9999, NA_val=NA)
#' relprob_subsets_matrix
#' sfunc(charcell="s1_1", relprob_subsets_matrix)
#' sfunc(charcell="s1_2", relprob_subsets_matrix)
#' sfunc(charcell="s1_3", relprob_subsets_matrix)
#' sfunc(charcell="s2_3", relprob_subsets_matrix)
#' 
#' relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, 
#' maxent_constraint_01=0.0001, NA_val=NA)
#' relprob_subsets_matrix
#' yfunc(charcell="y1", relprob_subsets_matrix)
#' yfunc(charcell="y2", relprob_subsets_matrix)
#' yfunc(charcell="y3", relprob_subsets_matrix)
#' yfunc(charcell="y4", relprob_subsets_matrix)
#' 
#' relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, 
#' maxent_constraint_01=0.5, NA_val=NA)
#' relprob_subsets_matrix
#' yfunc(charcell="y1", relprob_subsets_matrix)
#' yfunc(charcell="y2", relprob_subsets_matrix)
#' yfunc(charcell="y3", relprob_subsets_matrix)
#' yfunc(charcell="y4", relprob_subsets_matrix)
#' 
#' relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, 
#' maxent_constraint_01=0.9999, NA_val=NA)
#' relprob_subsets_matrix
#' yfunc(charcell="y1", relprob_subsets_matrix)
#' yfunc(charcell="y2", relprob_subsets_matrix)
#' yfunc(charcell="y3", relprob_subsets_matrix)
#' yfunc(charcell="y4", relprob_subsets_matrix)
#' 
#' 
#' # Probabilities of different descendant rangesizes, for the smaller descendant, 
#' # under vicariant speciation
#' relprob_subsets_matrix = relative_probabilities_of_vicariants(max_numareas=6, 
#' maxent_constraint_01v=0.0001, NA_val=NA)
#' relprob_subsets_matrix
#' vfunc(charcell="v1_1", relprob_subsets_matrix)
#' vfunc(charcell="v1_2", relprob_subsets_matrix)
#' vfunc(charcell="v1_3", relprob_subsets_matrix)
#' vfunc(charcell="v1_4", relprob_subsets_matrix)
#' vfunc(charcell="v2_4", relprob_subsets_matrix)
#' vfunc(charcell="v2_2", relprob_subsets_matrix)
#' vfunc(charcell="v1_6", relprob_subsets_matrix)
#' vfunc(charcell="v2_6", relprob_subsets_matrix)
#' vfunc(charcell="v3_6", relprob_subsets_matrix)
#' 
#' relprob_subsets_matrix = relative_probabilities_of_vicariants(max_numareas=6, 
#' maxent_constraint_01v=0.5, NA_val=NA)
#' relprob_subsets_matrix
#' vfunc(charcell="v1_1", relprob_subsets_matrix)
#' vfunc(charcell="v1_2", relprob_subsets_matrix)
#' vfunc(charcell="v1_3", relprob_subsets_matrix)
#' vfunc(charcell="v1_4", relprob_subsets_matrix)
#' vfunc(charcell="v2_4", relprob_subsets_matrix)
#' vfunc(charcell="v2_2", relprob_subsets_matrix)
#' vfunc(charcell="v1_6", relprob_subsets_matrix)
#' vfunc(charcell="v2_6", relprob_subsets_matrix)
#' vfunc(charcell="v3_6", relprob_subsets_matrix)
#' 
#' relprob_subsets_matrix = relative_probabilities_of_vicariants(max_numareas=6, 
#' maxent_constraint_01v=0.9999, NA_val=NA)
#' relprob_subsets_matrix
#' vfunc(charcell="v1_1", relprob_subsets_matrix)
#' vfunc(charcell="v1_2", relprob_subsets_matrix)
#' vfunc(charcell="v1_3", relprob_subsets_matrix)
#' vfunc(charcell="v1_4", relprob_subsets_matrix)
#' vfunc(charcell="v2_4", relprob_subsets_matrix)
#' vfunc(charcell="v2_2", relprob_subsets_matrix)
#' vfunc(charcell="v1_6", relprob_subsets_matrix)
#' vfunc(charcell="v2_6", relprob_subsets_matrix)
#' vfunc(charcell="v3_6", relprob_subsets_matrix)
#' 
sfunc <- function(charcell, relprob_subsets_matrix)
	{
	defaults='
	charcell="s1_6"
	maxent_constraint_01=1
	relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, maxent_constraint_01)
	'
	# Remove the "s"
	#"s1_3"
	charcell2 = gsub(pattern="s", replacement="", x=charcell)
	
	words = strsplit(charcell2, split="_")[[1]]
	#"1" "3"
	
	subset_daughter_range = as.numeric(words[1])
	widespread_daughter_range = as.numeric(words[2])
	
	
	# The relative probability of the different options
	# subset daughter areas 1, 2, 3, 4...
	#rangesizes = seq(1, widespread_daughter_range, by=1)
	
	prob_of_this_b = relprob_subsets_matrix[widespread_daughter_range, subset_daughter_range]
	prob_of_this_b
	
	return(prob_of_this_b)
	}



# Function giving the relative probability of y (sympatric speciation)
# given different range sizes of the subset descendent
#######################################################
# yfunc
#######################################################
#' Extract the appropriate probability for a sympatric/range-copying speciation event, given text code for rangesize of smaller descendant, and ancestor
#'
#' @param charcell The text in the cell, indicating the type of speciation/cladogenesis range inheritance event.
#' @param relprob_subsets_matrix A numeric matrix describing the relative probability of each smaller daughter range, conditional on the ancestral rangesize.
#' @return \code{prob_of_this_s}, a numeric value giving the relative probability of that descendent-ancestor rangesize pair.
#' @export
#' @seealso \code{\link{sfunc}}, \code{\link{vfunc}}, \code{\link{relative_probabilities_of_subsets}}, \code{\link{symbolic_to_relprob_matrix_sp}}, \code{\link{get_probvals}}, \code{\link[FD]{maxent}}, \code{\link{calcZ_part}}, \code{\link{calcP_n}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://en.wikipedia.org/wiki/Maximum_entropy_probability_distribution}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Harte2011
#'	 @cite ReeSmith2008
#'	 @cite Ronquist1996_DIVA
#'	 @cite Ronquist_1997_DIVA
#'	 @cite Ronquist_Sanmartin_2011
#'	 @cite Landis_Matzke_etal_2013_BayArea
#' @examples
#' testval=1
#' # Examples
#' 
#' # Probabilities of different descendant rangesizes, for the smaller 
#' # descendant, under sympatric/subset speciation
#' # (plus sympatric/range-copying, which is folded in):
#' relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, 
#' maxent_constraint_01=0.0001, NA_val=NA)
#' relprob_subsets_matrix
#' sfunc(charcell="s1_1", relprob_subsets_matrix)
#' sfunc(charcell="s1_2", relprob_subsets_matrix)
#' sfunc(charcell="s1_3", relprob_subsets_matrix)
#' sfunc(charcell="s2_3", relprob_subsets_matrix)
#' 
#' relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, 
#' maxent_constraint_01=0.5, NA_val=NA)
#' relprob_subsets_matrix
#' sfunc(charcell="s1_1", relprob_subsets_matrix)
#' sfunc(charcell="s1_2", relprob_subsets_matrix)
#' sfunc(charcell="s1_3", relprob_subsets_matrix)
#' sfunc(charcell="s2_3", relprob_subsets_matrix)
#' 
#' relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, 
#' maxent_constraint_01=0.9999, NA_val=NA)
#' relprob_subsets_matrix
#' sfunc(charcell="s1_1", relprob_subsets_matrix)
#' sfunc(charcell="s1_2", relprob_subsets_matrix)
#' sfunc(charcell="s1_3", relprob_subsets_matrix)
#' sfunc(charcell="s2_3", relprob_subsets_matrix)
#' 
#' relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, 
#' maxent_constraint_01=0.0001, NA_val=NA)
#' relprob_subsets_matrix
#' yfunc(charcell="y1", relprob_subsets_matrix)
#' yfunc(charcell="y2", relprob_subsets_matrix)
#' yfunc(charcell="y3", relprob_subsets_matrix)
#' yfunc(charcell="y4", relprob_subsets_matrix)
#' 
#' relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, 
#' maxent_constraint_01=0.5, NA_val=NA)
#' relprob_subsets_matrix
#' yfunc(charcell="y1", relprob_subsets_matrix)
#' yfunc(charcell="y2", relprob_subsets_matrix)
#' yfunc(charcell="y3", relprob_subsets_matrix)
#' yfunc(charcell="y4", relprob_subsets_matrix)
#' 
#' relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, 
#' maxent_constraint_01=0.9999, NA_val=NA)
#' relprob_subsets_matrix
#' yfunc(charcell="y1", relprob_subsets_matrix)
#' yfunc(charcell="y2", relprob_subsets_matrix)
#' yfunc(charcell="y3", relprob_subsets_matrix)
#' yfunc(charcell="y4", relprob_subsets_matrix)
#' 
#' 
#' # Probabilities of different descendant rangesizes, for the smaller descendant, 
#' # under vicariant speciation
#' relprob_subsets_matrix = relative_probabilities_of_vicariants(max_numareas=6, 
#' maxent_constraint_01v=0.0001, NA_val=NA)
#' relprob_subsets_matrix
#' vfunc(charcell="v1_1", relprob_subsets_matrix)
#' vfunc(charcell="v1_2", relprob_subsets_matrix)
#' vfunc(charcell="v1_3", relprob_subsets_matrix)
#' vfunc(charcell="v1_4", relprob_subsets_matrix)
#' vfunc(charcell="v2_4", relprob_subsets_matrix)
#' vfunc(charcell="v2_2", relprob_subsets_matrix)
#' vfunc(charcell="v1_6", relprob_subsets_matrix)
#' vfunc(charcell="v2_6", relprob_subsets_matrix)
#' vfunc(charcell="v3_6", relprob_subsets_matrix)
#' 
#' relprob_subsets_matrix = relative_probabilities_of_vicariants(max_numareas=6, 
#' maxent_constraint_01v=0.5, NA_val=NA)
#' relprob_subsets_matrix
#' vfunc(charcell="v1_1", relprob_subsets_matrix)
#' vfunc(charcell="v1_2", relprob_subsets_matrix)
#' vfunc(charcell="v1_3", relprob_subsets_matrix)
#' vfunc(charcell="v1_4", relprob_subsets_matrix)
#' vfunc(charcell="v2_4", relprob_subsets_matrix)
#' vfunc(charcell="v2_2", relprob_subsets_matrix)
#' vfunc(charcell="v1_6", relprob_subsets_matrix)
#' vfunc(charcell="v2_6", relprob_subsets_matrix)
#' vfunc(charcell="v3_6", relprob_subsets_matrix)
#' 
#' relprob_subsets_matrix = relative_probabilities_of_vicariants(max_numareas=6, 
#' maxent_constraint_01v=0.9999, NA_val=NA)
#' relprob_subsets_matrix
#' vfunc(charcell="v1_1", relprob_subsets_matrix)
#' vfunc(charcell="v1_2", relprob_subsets_matrix)
#' vfunc(charcell="v1_3", relprob_subsets_matrix)
#' vfunc(charcell="v1_4", relprob_subsets_matrix)
#' vfunc(charcell="v2_4", relprob_subsets_matrix)
#' vfunc(charcell="v2_2", relprob_subsets_matrix)
#' vfunc(charcell="v1_6", relprob_subsets_matrix)
#' vfunc(charcell="v2_6", relprob_subsets_matrix)
#' vfunc(charcell="v3_6", relprob_subsets_matrix)
#' 
yfunc <- function(charcell, relprob_subsets_matrix)
	{
	defaults='
	charcell="y3"
	maxent_constraint_01=1
	relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, maxent_constraint_01)
	
	'
	# Remove the "y"
	#"y3"
	charcell2 = as.numeric(gsub(pattern="y", replacement="", x=charcell))
	
	
	subset_daughter_range = charcell2
	widespread_daughter_range = charcell2
	
	
	# The relative probability of the different options
	# subset daughter areas 1, 2, 3, 4...
	#rangesizes = seq(1, widespread_daughter_range, by=1)
	
	prob_of_this_s = relprob_subsets_matrix[widespread_daughter_range, subset_daughter_range]
	prob_of_this_s
	
	return(prob_of_this_s)
	}


# Function giving the relative probability of j (jump dispersal/founder event speciation)
# given different range sizes of the subset descendent
# NOT DONE
# jfunc <- function(charcell, relprob_subsets_matrix)
# 	{
# 	defaults='
# 	charcell="j3"
# 	maxent_constraint_01=1
# 	relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, maxent_constraint_01)
# 	
# 	'
# 	# Remove the "j"
# 	#"j3"
# 	charcell2 = as.numeric(gsub(pattern="j", replacement="", x=charcell))
# 	
# 	
# 	subset_daughter_range = charcell2
# 	widespread_daughter_range = charcell2
# 	
# 	
# 	# The relative probability of the different options
# 	# subset daughter areas 1, 2, 3, 4...
# 	#rangesizes = seq(1, widespread_daughter_range, by=1)
# 	
# 	prob_of_this_s = relprob_subsets_matrix[widespread_daughter_range, subset_daughter_range]
# 	prob_of_this_s
# 	
# 	return(prob_of_this_s)
# 	}

# Function giving the relative probability of v (vicariant speciation)
# given different smaller range size and ancestral range size
#######################################################
# vfunc
#######################################################
#' Extract the appropriate probability for a vicariant speciation event, given text code for rangesize of smaller descendant, and ancestor
#'
#' @param charcell The text in the cell, indicating the type of speciation/cladogenesis range inheritance event.
#' @param relprob_vicar_matrix A numeric matrix describing the relative probability of each smaller daughter range, conditional on the ancestral rangesize.
#' @return \code{prob_of_this_v}, a numeric value giving the relative probability of that descendent-ancestor rangesize pair.
#' @export
#' @seealso \code{\link{sfunc}}, \code{\link{vfunc}}, \code{\link{relative_probabilities_of_subsets}}, \code{\link{symbolic_to_relprob_matrix_sp}}, \code{\link{get_probvals}}, \code{\link[FD]{maxent}}, \code{\link{calcZ_part}}, \code{\link{calcP_n}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://en.wikipedia.org/wiki/Maximum_entropy_probability_distribution}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Harte2011
#'	 @cite ReeSmith2008
#'	 @cite Ronquist1996_DIVA
#'	 @cite Ronquist_1997_DIVA
#'	 @cite Ronquist_Sanmartin_2011
#'	 @cite Landis_Matzke_etal_2013_BayArea
#' @examples
#' testval=1
#' # Examples
#' 
#' # Probabilities of different descendant rangesizes, for the smaller descendant, 
#' # under sympatric/subset speciation
#' # (plus sympatric/range-copying, which is folded in):
#' relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, 
#' maxent_constraint_01=0.0001, NA_val=NA)
#' relprob_subsets_matrix
#' sfunc(charcell="s1_1", relprob_subsets_matrix)
#' sfunc(charcell="s1_2", relprob_subsets_matrix)
#' sfunc(charcell="s1_3", relprob_subsets_matrix)
#' sfunc(charcell="s2_3", relprob_subsets_matrix)
#' 
#' relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, 
#' maxent_constraint_01=0.5, NA_val=NA)
#' relprob_subsets_matrix
#' sfunc(charcell="s1_1", relprob_subsets_matrix)
#' sfunc(charcell="s1_2", relprob_subsets_matrix)
#' sfunc(charcell="s1_3", relprob_subsets_matrix)
#' sfunc(charcell="s2_3", relprob_subsets_matrix)
#' 
#' relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, 
#' maxent_constraint_01=0.9999, NA_val=NA)
#' relprob_subsets_matrix
#' sfunc(charcell="s1_1", relprob_subsets_matrix)
#' sfunc(charcell="s1_2", relprob_subsets_matrix)
#' sfunc(charcell="s1_3", relprob_subsets_matrix)
#' sfunc(charcell="s2_3", relprob_subsets_matrix)
#' 
#' relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, 
#' maxent_constraint_01=0.0001, NA_val=NA)
#' relprob_subsets_matrix
#' yfunc(charcell="y1", relprob_subsets_matrix)
#' yfunc(charcell="y2", relprob_subsets_matrix)
#' yfunc(charcell="y3", relprob_subsets_matrix)
#' yfunc(charcell="y4", relprob_subsets_matrix)
#' 
#' relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, 
#' maxent_constraint_01=0.5, NA_val=NA)
#' relprob_subsets_matrix
#' yfunc(charcell="y1", relprob_subsets_matrix)
#' yfunc(charcell="y2", relprob_subsets_matrix)
#' yfunc(charcell="y3", relprob_subsets_matrix)
#' yfunc(charcell="y4", relprob_subsets_matrix)
#' 
#' relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, 
#' maxent_constraint_01=0.9999, NA_val=NA)
#' relprob_subsets_matrix
#' yfunc(charcell="y1", relprob_subsets_matrix)
#' yfunc(charcell="y2", relprob_subsets_matrix)
#' yfunc(charcell="y3", relprob_subsets_matrix)
#' yfunc(charcell="y4", relprob_subsets_matrix)
#' 
#' 
#' # Probabilities of different descendant rangesizes, for the smaller descendant, 
#' # under vicariant speciation
#' relprob_subsets_matrix = relative_probabilities_of_vicariants(max_numareas=6, 
#' maxent_constraint_01v=0.0001, NA_val=NA)
#' relprob_subsets_matrix
#' vfunc(charcell="v1_1", relprob_subsets_matrix)
#' vfunc(charcell="v1_2", relprob_subsets_matrix)
#' vfunc(charcell="v1_3", relprob_subsets_matrix)
#' vfunc(charcell="v1_4", relprob_subsets_matrix)
#' vfunc(charcell="v2_4", relprob_subsets_matrix)
#' vfunc(charcell="v2_2", relprob_subsets_matrix)
#' vfunc(charcell="v1_6", relprob_subsets_matrix)
#' vfunc(charcell="v2_6", relprob_subsets_matrix)
#' vfunc(charcell="v3_6", relprob_subsets_matrix)
#' 
#' relprob_subsets_matrix = relative_probabilities_of_vicariants(max_numareas=6, 
#' maxent_constraint_01v=0.5, NA_val=NA)
#' relprob_subsets_matrix
#' vfunc(charcell="v1_1", relprob_subsets_matrix)
#' vfunc(charcell="v1_2", relprob_subsets_matrix)
#' vfunc(charcell="v1_3", relprob_subsets_matrix)
#' vfunc(charcell="v1_4", relprob_subsets_matrix)
#' vfunc(charcell="v2_4", relprob_subsets_matrix)
#' vfunc(charcell="v2_2", relprob_subsets_matrix)
#' vfunc(charcell="v1_6", relprob_subsets_matrix)
#' vfunc(charcell="v2_6", relprob_subsets_matrix)
#' vfunc(charcell="v3_6", relprob_subsets_matrix)
#' 
#' relprob_subsets_matrix = relative_probabilities_of_vicariants(max_numareas=6, 
#' maxent_constraint_01v=0.9999, NA_val=NA)
#' relprob_subsets_matrix
#' vfunc(charcell="v1_1", relprob_subsets_matrix)
#' vfunc(charcell="v1_2", relprob_subsets_matrix)
#' vfunc(charcell="v1_3", relprob_subsets_matrix)
#' vfunc(charcell="v1_4", relprob_subsets_matrix)
#' vfunc(charcell="v2_4", relprob_subsets_matrix)
#' vfunc(charcell="v2_2", relprob_subsets_matrix)
#' vfunc(charcell="v1_6", relprob_subsets_matrix)
#' vfunc(charcell="v2_6", relprob_subsets_matrix)
#' vfunc(charcell="v3_6", relprob_subsets_matrix)
#' 
vfunc <- function(charcell, relprob_vicar_matrix)
	{
	defaults='
	charcell="v1_6"
	maxent_constraint_01v=0.00001
	relprob_vicar_matrix = relative_probabilities_of_subsets(max_numareas=6, maxent_constraint_01v)
	'
	# Remove the "v"
	#"v1_3"
	charcell2 = gsub(pattern="v", replacement="", x=charcell)
	
	words = strsplit(charcell2, split="_")[[1]]
	#"1" "3"
	
	smaller_daughter_range = as.numeric(words[1])
	widespread_ancestor_range = as.numeric(words[2])
	
	
	# The relative probability of the different options
	# subset daughter areas 1, 2, 3, 4...
	#rangesizes = seq(1, widespread_ancestor_range, by=1)
	
	prob_of_this_v = relprob_vicar_matrix[widespread_ancestor_range, smaller_daughter_range]
	prob_of_this_v
	
	return(prob_of_this_v)
	}




#######################################################
# nullsym_to_NA
#######################################################
#' Convert a specified null range code to NA
#' 
#' Takes a matrix \code{mat}, converts any instances of the \code{nullsym} symbol to \code{NA}.
#' 
#' @param mat A matrix.
#' @param nullsym A character specifying the null symbol.
#' @return \code{mat} The revised matrix
#' @export
#' @seealso \code{\link{remove_null_rowcols_from_mat}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' mat = matrix(c("-",1,1,1,"-",1,1,1,"-"), nrow=3)
#' mat
#' mat2 = nullsym_to_NA(mat, nullsym="-")
#' mat2
#'
nullsym_to_NA <- function(mat, nullsym="-")
	{
	mat[mat==nullsym] = NA
	return(mat)
	}





#######################################################
# Process optim results
#######################################################

#######################################################
# process_optim
#######################################################
#' Extract \code{optim} results to a row
#' 
#' After running an ML (maximum likelihood) search with \code{\link[stats]{optim}}, \code{\link[stats]{optim}} returns
#' a list with a variety of objects.  It is often handy to have the parameter values, log-likelihood, etc., extracted
#' to a table for comparison with other optimization runs.  \code{process_optim} does this.
#' 
#' @param optim_results A results object from \code{\link[stats]{optim}}
#' @param max_num_params Specify the number of parameters, if known. If \code{NULL}, the code will try to guess.
#' @return \code{tmprow3} A row holding the \code{\link[stats]{optim}} results, which can then be added to a table with \code{\link[base]{rbind}}.
#' @export
#' @seealso \code{\link[stats]{optim}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' testval=1
#' # Any optim() for a biogeography scenario would take too long to run for R CMD check.
#' 
process_optim <- function(optim_results, max_num_params=NULL)
	{
	# Get the number of params
	num_params = length(optim_results$par)


	if (is.null(max_num_params) || (num_params == max_num_params))
		{
		# Re-arrange optim results
		tmprow = unlist(optim_results)
		tmprow2_nonparam_indices = (num_params+1):length(tmprow)
		tmprow2 = tmprow[tmprow2_nonparam_indices]
		tmprow3 = c(tmprow2, tmprow[1:num_params])
		} else {
		
		# Error check
		if (max_num_params < num_params)
			{
			stop("ERROR: num_params in optim_results is GREATER than max_num_params")
			}

		# Fill out extra params
		num_zeros_to_add = max_num_params - length(optim_results$par)
		zeros_to_add = rep(0, num_zeros_to_add)
		
		parnums_to_add = (length(optim_results$par)+1) : max_num_params
		parnums_to_add_names = paste("par", parnums_to_add, sep="")
		
		# Re-arrange optim results
		tmprow = unlist(optim_results)
		tmprow2_nonparam_indices = (num_params+1):length(tmprow)
		tmprow2 = tmprow[tmprow2_nonparam_indices]
		tmprow2 = c(tmprow2, tmprow[1:num_params])
		names_tmprow2 = names(tmprow2) 
		names_tmprow3 = c(names_tmprow2, parnums_to_add_names)
		
		tmprow3 = c(tmprow2, zeros_to_add)
		names(tmprow3) = names_tmprow3
		}
	
	return(tmprow3)
	}







