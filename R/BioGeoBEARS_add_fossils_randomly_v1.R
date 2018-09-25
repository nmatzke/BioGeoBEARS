#######################################################
# Error checks for dumb things in fossils xls file
#######################################################

error_checks_for_fossils_xls <- function(xls, fix_duplicates=TRUE)
	{
	defaults='
	error_checks_for_fossils_xls(xls)
	'

	# Remove rows you don't want to use
	# Also remove rows where use is blank!
	TF = isblank_TF(xls$use)
	if (sum(TF) > 0)
		{
		txt = paste0("WARNING from error_checks_for_fossils_xls(): You have ", sum(TF), " blanks (or NAs etc.) in the 'use' column of your input Excel spreadsheet.  These are all being auto-set to 'no', meaning these rows will not be used. Edit the spreadsheet if you want a different result.")
		warning(txt)

		xls$use[TF] = "no"
		}
	
	xls = xls[xls$use=="yes",]
	
	# Blanks or NAs in max_age
	if (any(isblank_TF(xls$max_age)) == TRUE)
		{
		stoptxt = "STOP ERROR in error_checks_for_fossils_xls(): You have blanks (or NAs etc.) in some row(s) for max_age. Check your input Excel spreadsheet."
		cat("\n\n")
		cat(stoptxt)
		print(xls[TF,])
		cat("\n\n")
		} # END if (any(isblank_TF(xls$max_age)) == TRUE)
	
	TF = xls$max_age < xls$min_age
	if (any(TF) == TRUE)
		{
		stoptxt = "STOP ERROR in error_checks_for_fossils_xls(): You have max_age < min_age in some row(s). Printing them below:"
		cat("\n\n")
		cat(stoptxt)
		print(xls[TF,])
		cat("\n\n")
		} # END if (any(TF) == TRUE)
	
	# No blanks allowed
	TF = isblank_TF(xls$Fossil_Taxon)
	if (any(TF) == TRUE)
		{
		stoptxt = "STOP ERROR in error_checks_for_fossils_xls(): No cells in column xls$Fossil_Taxon can be blank. Go check your input spreadsheet!!"
		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
		stop(stoptxt)
		}

	# No spaces allowed
	fossil_taxa_names = xls$Fossil_Taxon
	TF = grepl(pattern=" ", x=fossil_taxa_names)
	if (any(TF) == TRUE)
		{
		stoptxt = "STOP ERROR in error_checks_for_fossils_xls(): No cells in xls$Fossil_Taxon can have spaces (' ') in the fossil taxon names. Go check your input spreadsheet!!"
		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
		stop(stoptxt)
		}
	
	# No commas allowed
	TF = grepl(pattern=",", x=fossil_taxa_names)
	if (any(TF) == TRUE)
		{
		stoptxt = "STOP ERROR in error_checks_for_fossils_xls(): No cells in xls$Fossil_Taxon can have commas (',') in the fossil taxon names. Go check your input spreadsheet!!"
		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
		stop(stoptxt)
		}


	# No ? allowed
	TF = grepl(pattern="\\?", x=fossil_taxa_names)
	if (any(TF) == TRUE)
		{
		stoptxt = "STOP ERROR in error_checks_for_fossils_xls(): No cells in xls$Fossil_Taxon can have commas (',') in the fossil taxon names. Go check your input spreadsheet!!"
		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
		stop(stoptxt)
		}

	# No slashes allowed
	TF = grepl(pattern="/", x=fossil_taxa_names)
	if (any(TF) == TRUE)
		{
		stoptxt = "STOP ERROR in error_checks_for_fossils_xls(): No cells in xls$Fossil_Taxon can have slashes ('/') in the fossil taxon names. Go check your input spreadsheet!!"
		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
		stop(stoptxt)
		}

	# No parenthesis allowed
	TF = grepl(pattern="\\(", x=fossil_taxa_names)
	if (any(TF) == TRUE)
		{
		stoptxt = "STOP ERROR in error_checks_for_fossils_xls(): No cells in xls$Fossil_Taxon can have parentheses ('(' or ')') in the fossil taxon names. Go check your input spreadsheet!!"
		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
		stop(stoptxt)
		}

	# No parenthesis allowed
	TF = grepl(pattern="\\)", x=fossil_taxa_names)
	if (any(TF) == TRUE)
		{
		stoptxt = "STOP ERROR in error_checks_for_fossils_xls(): No cells in xls$Fossil_Taxon can have parentheses ('(' or ')') in the fossil taxon names. Go check your input spreadsheet!!"
		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
		stop(stoptxt)
		}

	# No backslashes allowed
	TF = grepl(pattern="\\\\", x=fossil_taxa_names)
	if (any(TF) == TRUE)
		{
		stoptxt = "STOP ERROR in error_checks_for_fossils_xls(): No cells in xls$Fossil_Taxon can have backslashes ('\') in the fossil taxon names. Go check your input spreadsheet!!"
		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
		stop(stoptxt)
		}


	
	# No periods allowed, except 'indet.' and 'sp.'
	TF = grepl(pattern="\\.", x=fossil_taxa_names)
	TF2 = grepl(pattern="indet.", x=fossil_taxa_names)
	TF3 = grepl(pattern="sp.", x=fossil_taxa_names)
	TF[TF2] = FALSE
	TF[TF3] = FALSE
	fossil_taxa_names = gsub(pattern="indet.", replacement="indet", x=fossil_taxa_names)
	fossil_taxa_names = gsub(pattern="sp.", replacement="sp", x=fossil_taxa_names)
	xls$Fossil_Taxon = fossil_taxa_names
	
	if (any(TF) == TRUE)
		{
		stoptxt = "WARNING from error_checks_for_fossils_xls(): No cells in xls$Fossil_Taxon can have periods ('.') in the fossil taxon names (except for 'indet.' or 'sp.', where the periods will be auto-removed). Go check your input spreadsheet!!"
		cat("\n")
		cat(stoptxt)
		cat("\n")
		warning(stoptxt)
		}

	# No repeated names allowed
	uniq_fossil_taxa = unique(fossil_taxa_names)

	if ( (length(uniq_fossil_taxa) != length(fossil_taxa_names)) && (fix_duplicates == FALSE) )
		{
		counts = table(fossil_taxa_names)
		duplicate_names = names(counts)[counts > 1]
		duplicate_names
				
		stoptxt = "STOP ERROR in error_checks_for_fossils_xls(): There can be *NO* duplicate names in xls$Fossil_Taxon. Check your input spreadsheet!! The duplicated names are listed below..."
		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
		cat(duplicate_names, sep="\n")
		cat("\n\n")
		stop(stoptxt)	
		}
	
	# Fix repeated names
	if (length(uniq_fossil_taxa) != length(fossil_taxa_names))
		{
		counts = table(fossil_taxa_names)
		duplicate_names = names(counts)[counts > 1]
		duplicate_names
		
		# Fix the duplicates
		for (i in 1:length(duplicate_names))
			{
			name_to_add_numbers_to = duplicate_names[i]
			TF = fossil_taxa_names == name_to_add_numbers_to
			nums_to_add = 1:sum(TF)
			nums_as_strings = sprintf(fmt="%04.0f", (nums_to_add))
			nums_as_strings = paste("_specimen", nums_as_strings, sep="")
			new_names = paste(fossil_taxa_names[TF], nums_as_strings, sep="")
			
			# Insert new names
			fossil_taxa_names[TF] = new_names
			} # END for (i in 1:length(duplicate_names))
		cat("\n")
		txt = "WARNING from error_checks_for_fossils_xls(): added specimen numbers to OTU names in xls, in order to make all OTUs unique."
		cat(txt)
		cat("\n")
		warning(txt)
		xls$Fossil_Taxon = fossil_taxa_names
		} # END if (length(uniq_fossil_taxa) != length(fossil_taxa_names))
	
	# Blank cells
	TF = isblank_TF(xls$range_coded)
	if (any(TF) == TRUE)
		{
		stoptxt = "STOP ERROR in error_checks_for_fossils_xls(): No cells in column xls$range_coded can be blank. Go check your input spreadsheet!!"
		cat("\n\n")
		cat(stoptxt)
		cat("\n\n")
		stop(stoptxt)
		}
	
		
	cat("\nPassed error checks in error_checks_for_fossils_xls(), but this only checks for the problems Nick thought of. Returning xls (perhaps modified!)...\n\n")
	return(xls)
	} # END error_checks_for_fossils_xls <- function(xls)


# Add fossils to some or all of the trees in a NEXUS treefile,
# outputing a Newick file for each
add_fossils_to_many_trees <- function(nexfn, burnin, numtrees, xls, living_tipranges_fn, fn_prefix="tr", fn_suffix="", brlen_of_side_branch=1e-07, sample_trs_w_replacement=TRUE, plottree=TRUE)
	{
	defaults='	fn_suffix=""

	fn_prefix="tr"
	fn_suffix=""
	sample_trs_w_replacement=TRUE
	'
	
	# Make a list of the output tree filenames
	output_trfns = NULL
	
	txt = paste("\n\nadd_fossils_to_many_trees():\n\n", sep="")
	cat(txt)
		
	# Read the input geography file with the geography of living species 
	# (or whatever you have before the fossils are added)
	txt = paste("...is reading in the pre-fossils geographic ranges in '", living_tipranges_fn, "'.\n", sep="")
	cat(txt)
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=living_tipranges_fn)
	living_tipranges_df = tipranges@df

	
	# Process the input trees
	txt = paste("...is processing the trees in NEXUS file '", nexfn, "':\n\n", sep="")
	cat(txt)
	trs = read.nexus(nexfn)
	trs
	
	if (class(trs) == "multiPhylo")
		{
		num_nexus_trees = length(trs)
		} else {
		num_nexus_trees = 1
		} # END if (class(trs) == "multiPhylo")

	

	# Cut the burnin
	if (class(trs) == "multiPhylo")
		{
		trs = trs[(burnin+1):num_nexus_trees]
		numtrees_to_sample_from = length(trs)
		} else {
		trs = trs
		numtrees_to_sample_from = 1
		}
	
	
	# If the number of trees matches exactly, no sampling is needed, but I guess it doesn't matter
	treenums_to_use = sort(sample(x=1:numtrees_to_sample_from, size=numtrees, replace=sample_trs_w_replacement))

	# Cycle through the list of input trees
	for (i in 1:length(treenums_to_use))
		{
		txt = paste("\n\n\nAdding fossils to tree #", i, "/", length(treenums_to_use), "\n", sep="")
		cat(txt)
		
		# Re-load tipranges from original tipranges file
		tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=living_tipranges_fn)
		living_tipranges_df = tipranges@df
		
		
		# Select a tree from the sample
		treenum_to_use = treenums_to_use[i]
		
		if (class(trs) == "multiPhylo")
			{
			tr = trs[[ treenum_to_use ]]
			} else {
			tr = trs
			} # END if (class(trs) == "multiPhylo")
		trnum_txt = sprintf(fmt="%04.0f", (treenum_to_use+burnin))
		trnum_txt2 = sprintf(fmt="%04.0f", (i))


		# Are there any taxa in the tree that are not in the geography file?  That's OK, if these
		# taxa are added via the add_fossils file.
		treetips_missing_from_geog_TF = (tr$tip.label %in% row.names(tipranges@df)) == FALSE
		tips_to_drop = tr$tip.label[treetips_missing_from_geog_TF]
		tips_to_drop_txt = paste0(sort(tips_to_drop), collapse=", ")
		num_missing = sum(treetips_missing_from_geog_TF)
		if (num_missing > 0)
			{
			txt = paste0("Warning: add_fossils_to_many_trees() has noticed that your tree file '", nexfn, "', tree #", treenum_to_use, " has ", num_missing, " tipnames that are not in the geography file, '", living_tipranges_fn, "'. The tipnames that are in the tree, but missing from the geography file, are: ", tips_to_drop_txt, ".\n")
			cat(txt)
			warning(txt)
			#tr = drop.tip(phy=tr, tip=tips_to_drop)
			cat("\nTips removed from tree:\n")
			cat(sort(tips_to_drop), sep="\n")
			} # END if (num_missing > 0)

		# Are there any taxa in the geography file that are not in the tree?  Cut those.
		geogtips_missing_from_tree_TF = (row.names(tipranges@df) %in% tr$tip.label) == FALSE
		geogtips_to_drop = row.names(tipranges@df)[geogtips_missing_from_tree_TF]
		geogtips_to_drop_txt = paste0(sort(geogtips_to_drop), collapse=", ")
		num_missing = sum(geogtips_missing_from_tree_TF)
		if (num_missing > 0)
			{
			txt = paste0("Wait. add_fossils_to_many_trees() has noticed that your geography file '", living_tipranges_fn, "' has ", num_missing, " names that are not in the treefile, '", nexfn, "', tree #", treenum_to_use, ". Removing them from the geography file.  Really, you should probably get this kind of basic thing straightened out before running a complex analysis like this!  The tipnames that are in the geography file, but missing from the tree, are: ", geogtips_to_drop_txt, ".\n")

			cat(txt)
			warning(txt)
			
			tmpdf = tipranges@df
			tmpdf2 = tmpdf[geogtips_missing_from_tree_TF==FALSE, ]
			tipranges@df = tmpdf2

			cat("\nNames removed from geography file:\n")
			cat(sort(geogtips_to_drop), sep="\n")
			} # END if (num_missing > 0)



		
		# Add the list of fossils randomly, and save the tree
		tr_w_fossils_and_fossils_added_list = add_fossils_from_xls_randomly(tr=tr, xls=xls, brlen_of_side_branch=brlen_of_side_branch, plottree=plottree)
		tr_w_fossil_tips_added = tr_w_fossils_and_fossils_added_list$tr_w_fossil
		list_of_fossil_tips_added =	tr_w_fossils_and_fossils_added_list$list_of_fossil_tips_added
		outtr_fn = paste(fn_prefix, trnum_txt, "_wFossils", trnum_txt2, fn_suffix, ".newick", sep="")
		write.tree(phy=tr_w_fossil_tips_added, file=outtr_fn)
		
		cat("\n\n")
		txt= paste("Tree with fossils added saved to '", outtr_fn, "', which has ", length(tr_w_fossil_tips_added$tip.label), " tips.\n", sep="")
		cat(txt)
		
		# Add the list of fossils added to the tree, to the geography filename
		new_tipranges_df = add_list_of_fossil_tips_to_tipranges(tipranges_df=living_tipranges_df, list_of_fossil_tips_added=list_of_fossil_tips_added, xls=xls, numareas=ncol(living_tipranges_df))
		tipranges@df = new_tipranges_df
		# Save the new tipranges
		outgeog_fn = paste(fn_prefix, trnum_txt, "_wFossils", trnum_txt2, fn_suffix, ".geog", sep="")
		save_tipranges_to_LagrangePHYLIP(tipranges_object=tipranges, lgdata_fn=outgeog_fn)

		txt= paste("Geog with fossils added saved to '", outgeog_fn, "'.\n", sep="")
		cat(txt)
		#moref(outgeog_fn)
		
		# Add to list
		output_trfns = c(output_trfns, outtr_fn)
		} # end for (i in 1:length(treenums_to_use))
		# End cycle through input trees
	
	return(output_trfns)
	} # end add_fossils_to_many_trees()


# Make a geography file by adding fossils to a tipranges data frame
# tipranges_df: just a table, e.g. tipranges@df
# list_of_fossil_tips_added: a vector of names matching paste("fossil_", xls$Fossil_Taxon, sep="")
# xls: e.g. xls_full = read.xls("/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/z_help/Ruud_Scharn/RuudScharn_fossilsLiving_for_biogeography_v3.xlsx", sheet="fossils", skip=5)
add_list_of_fossil_tips_to_tipranges <- function(tipranges_df, list_of_fossil_tips_added, xls, numareas=ncol(tipranges_df))
	{
	defaults='
	tipranges_df=living_tipranges_df
	numareas=ncol(living_tipranges_df)
	
	# ...OR...
	tipranges_df = tipranges@df
	numareas=ncol(tipranges_df)
	xls = read.xls("/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/z_help/Ruud_Scharn/RuudScharn_fossilsLiving_for_biogeography_v3.xlsx", sheet="fossils", skip=5)

	list_of_fossil_tips_added
	'
	
	
	# The names of the fossils listed in xls
	fossils_in_xls = paste("fossil_", xls$Fossil_Taxon, sep="")
	
	
	# Check for an empty list of fossil tips to add; return original tipranges if so.
	if (length(list_of_fossil_tips_added) == 0)
		{
		txt = paste0("WARNING: add_list_of_fossil_tips_to_tipranges() found *0* new fossil ranges to add. Returning the original list of tipranges, with 0 fossils added.")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		
		warning(txt)
		
		new_tipranges_df = tipranges_df
		return(new_tipranges_df)
		} # END if (length(list_of_fossil_tips_added) == 0)
	
	
	# Add the fossil ranges to the list of ranges
	# Initialize
	tipranges_to_add = NULL
	for (i in 1:length(list_of_fossil_tips_added))
		{
		# Find the row in the xls table:
		TF = list_of_fossil_tips_added[i] == fossils_in_xls
		
		# STUPID Excel readin converts "'0000100" to "100" in R
		# Fix with sprintf
		#######################################################
		# THIS IS NO LONGER NEEDED, WITH XLConnect
		#######################################################
		
		if (numareas < 10)
			{
			sprintf_str = paste("%0", numareas, ".0f", sep="")
			}
		if (numareas >= 10)
			{
			sprintf_str = paste("%", numareas, ".0f", sep="")
			}
		sprintf_str = "%s"
		
		range_coded = sprintf(fmt=sprintf_str, xls$range_coded[TF])
		ranges = strsplit(range_coded, split="")[[1]]
		
		# What sort of geographic constraint?
		geog_constraint = xls$geog_constraint[TF]
		
		#print(geog_constraint)
		
		# If the fossil data is presence-only, then the absences don't mean anything,
		# just the presences
		if (geog_constraint == "presence-only")
			{
			ranges[ranges == "0"] = "?"
			}
		
		# If the fossil data is absence-only, then the presences don't mean anything,
		# just the absences
		# (which would be VERY WEIRD -- but MAYBE sometimes you could exclude a fossil
		#  from a well-sampled region)
		if (geog_constraint == "absence-only")
			{
			ranges[ranges == "1"] = "?"
			}

		# If the fossil data is "both", then the presences and absences
		# are both meaningful; i.e., you know the true range.
		# An alternative to the above strategies is to use a 
		# detection model, see Matzke Ph.D., chapter 4
		if (geog_constraint == "both")
			{
			# Make no change
			ranges = ranges
			}

		# detection model, see Matzke Ph.D., chapter 4
		if (geog_constraint == "detection_model")
			{
			# Not implemented yet -- you'd have to write an observation counts file
			# and a taphonomic control counts file
			error_txt = "\n\nWARNING: add_list_of_fossil_tips_to_tipranges() does not yet implement use of the detection model\n\n"
			cat(error_txt)
			
			# Make no change
			ranges = ranges
			}

		
		# OK, now add the modified ranges to the table of ranges to add
		tipranges_to_add = rbind(tipranges_to_add, ranges)
		}
	
	# Convert the tipranges to a data.frame
	tipranges_to_add = adf2(tipranges_to_add)
	tipranges_to_add
	
	names(tipranges_to_add) = names(tipranges_df)
	rownames(tipranges_to_add) = list_of_fossil_tips_added
	
	# Merge the tipranges
	new_tipranges_df = rbind(tipranges_df, tipranges_to_add)
	new_tipranges_df
	
	return(new_tipranges_df)
	}

# Add fossils from a data.frame, e.g. as read in from Excel with read.xls() from library(gdata)
# The data.frame MUST have the following headings:
# 
# Fossil_Taxon (no spaces, ".", "/", or any other punctuation except "_")
# stem_or_crown
# min_age
# max_age
# range_coded
# is_MRCA_of
# and_MRCA_of
# geog_constraint
# constraint_type
# taxon_specifiers
# specifier_type
#
# See example Excel file for what goes in each column.
# 
add_fossils_from_xls_randomly <- function(tr, xls, brlen_of_side_branch=1e-07, plottree=TRUE)
	{
	list_of_fossil_tips_added = NULL

	# Error check: no spaces in fossil names etc.
	colnames_to_check = c("Fossil_Taxon", "stem_or_crown", "range_coded", "is_MRCA_of", "and_MRCA_of", "geog_constraint", "constraint_type", "taxon_specifiers", "specifier_type")

	for (i in 1:length(colnames_to_check))
		{
		# Build and execute a TRUE/FALSE (TF) command:
		cmdstr = paste('TF1 = grepl(pattern=" ", x=xls$', colnames_to_check[i], ")", sep="")
		eval(parse(text=cmdstr))

		cmdstr = paste('TF2 = grepl(pattern="/", x=xls$', colnames_to_check[i], ")", sep="")
		eval(parse(text=cmdstr))

		cmdstr = paste('TF3 = grepl(pattern="\\\\.", x=xls$', colnames_to_check[i], ")", sep="")
		eval(parse(text=cmdstr))

		cmdstr = paste('TF4 = grepl(pattern="\\\\?", x=xls$', colnames_to_check[i], ")", sep="")
		eval(parse(text=cmdstr))
	
		cmdstr = paste('TF5 = grepl(pattern="\\\\(", x=xls$', colnames_to_check[i], ")", sep="")
		eval(parse(text=cmdstr))

		cmdstr = paste('TF6 = grepl(pattern="\\\\)", x=xls$', colnames_to_check[i], ")", sep="")
		eval(parse(text=cmdstr))
	
		TFsum = TF1 + TF2 + TF3 + TF4 + TF5 + TF6
		TF = TFsum >= 1
	
		# Report error if it occurs
		if (sum(TF) > 0)
			{
			error_txt = paste("\n\nWARNING in add_fossils_from_xls_randomly(): Spaces, punctuation, etc. are illegal in xls$", colnames_to_check[i], ".\n You might keep them if 'indet.' etc., though, so those will be removed. Printing the offending entries:\n\n", sep="")
			cat(error_txt)
		
			cmdstr = paste('offending_entries = xls$', colnames_to_check[i], "[TF]", sep="")
			eval(parse(text=cmdstr))
	
			for (j in 1:length(offending_entries))
				{
				txt = paste("'", offending_entries[j], "'\n", sep="")
				cat(txt)
				}
			
			# Fix "indet."
			offending_entries_new = gsub(pattern="indet.", replacement="indet", x=offending_entries)
			cmdstr = paste('xls$', colnames_to_check[i], "[TF] = offending_entries_new", sep="")
			eval(parse(text=cmdstr))

			# Fix "sp."
			offending_entries_new = gsub(pattern="sp.", replacement="sp", x=offending_entries)
			cmdstr = paste('xls$', colnames_to_check[i], "[TF] = offending_entries_new", sep="")
			eval(parse(text=cmdstr))
			
			warning("Warning from add_fossils_from_xls_randomly(): punctuation in OTU names is a bad idea, normally.")
			} # end: if (sum(TF) > 0)
		} # end: for (i in 1:length(colnames_to_check))



	# Add fossils (skipping fixnodes)
	# options are side_branch, ancestor, fixnode
	# skip fixnode cases here
	TF = xls$constraint_type != "fixnode"
	xls = xls[TF,]
	head(xls)

	numrows = nrow(xls)
	numrows
	# 33 fossils to add
	
	if (numrows < 1)
		{
		txt = paste0("WARNING: add_fossils_from_xls_randomly() found *0* active rows in 'xls'. Returning the original tree, with 0 fossils added.")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		
		tr_w_fossils_and_fossils_added_list = NULL
		tr_w_fossils_and_fossils_added_list$tr_w_fossil = tr
		tr_w_fossils_and_fossils_added_list$list_of_fossil_tips_added = list()
			
		warning(txt)
		return(tr_w_fossils_and_fossils_added_list)
		}
	
	# Loop through all fossils, adding the fossils one by one
	i=1
	list_of_fossil_tips_added = NULL
	#numrows = 1
	tr_orig = tr
	tr_w_fossil = tr
	
	for (i in 1:numrows)
		{
		counter_txt = paste("\n\nAttempting to add fossil #", i, "/", numrows, "...", sep="")
		cat(counter_txt)

		# Randomly draw the age of the fossil
		# Ages are times before present
		fossil_name = paste("fossil_", xls$Fossil_Taxon[i], sep="")
		max_age = as.numeric(xls$max_age[i])
		min_age = as.numeric(xls$min_age[i])
		constraint_type = as.character(xls$constraint_type[i])
		stem_or_crown = as.character(xls$stem_or_crown[i])



		# Error checks - 2017-03-14
		if (xls$specifier_type[i] == "specifiers")
			{
			# if ((xls$is_MRCA_of[i] %in% tr$tip.label) == FALSE)
			if ((xls$is_MRCA_of[i] %in% tr_w_fossil$tip.label) == FALSE)
				{
				txt = paste0("STOP ERROR in add_fossils_from_xls_randomly(): the tip specifier '", xls$is_MRCA_of[i], "' DOES NOT EXIST as a tip in the tree, even after adding any previous fossils in the table. Fix the Excel fossils file, and/or the tree!")
				cat("\n\n")
				cat(txt)
				cat("\n\n")
				stop(txt)
				}

			# if ((xls$and_MRCA_of[i] %in% tr$tip.label) == FALSE)
			if ((xls$and_MRCA_of[i] %in% tr_w_fossil$tip.label) == FALSE)
				{
				txt = paste0("STOP ERROR in add_fossils_from_xls_randomly(): the tip specifier '", xls$and_MRCA_of[i], "' DOES NOT EXIST as a tip in the tree, even after adding any previous fossils in the table. Fix the Excel fossils file, and/or the tree!")
				cat("\n\n")
				cat(txt)
				cat("\n\n")
				stop(txt)
				}
			} # END if (xls$specifier_type[i] == "specifiers")
		if (xls$specifier_type[i] == "species")
			{
			#if ((xls$taxon_specifiers[i] %in% tr$tip.label) == FALSE)
			if ((xls$taxon_specifiers[i] %in% tr_w_fossil$tip.label) == FALSE)
				{
				txt = paste0("STOP ERROR in add_fossils_from_xls_randomly(): the tip specifier '", xls$taxon_specifiers[i], "' DOES NOT EXIST as a tip in the tree, even after adding any previous fossils in the table. Fix the Excel fossils file, and/or the tree!")
				cat("\n\n")
				cat(txt)
				cat("\n\n")
				
				cat("Sorted tip labels in tree:\n")
				print(sort(tr$tip.label))
				stop(txt)
				}			
			} # END if (xls$specifier_type[i] == "species")



		# Get the tip_specifiers
		if (xls$specifier_type[i] == "specifiers")
			{
			# Get the two tip species that are tip specifiers
			tip_specifiers = c(xls$is_MRCA_of[i], xls$and_MRCA_of[i])
			} else {
			if ((xls$specifier_type[i] == "species") || (xls$specifier_type[i] == "specifiers"))
				{
				# Fossil is related to a terminal branch on the tree
				# Only ONE tip specified (stem of a single species is being used)
				# Input the species name
				# Be sure add_fossil_randomly() can accept a single-tip specifier
				tip_specifiers = xls$taxon_specifiers[i]
				} else {
				# Using Genus or Family as the specifier
				# Some ranked group; MUST be in the tip names
				# Note that this REQUIRES that names be NON-REDUNDANT; I do
				# this by searching on e.g. "FamilyName_"
				pattern_to_search_on = paste(xls$taxon_specifiers[i], "_", sep="")
				tips_in_clade_TF = grepl(pattern=pattern_to_search_on, x=tr$tip.label)
				tips_in_clade = tr$tip.label[tips_in_clade_TF]
				tip_specifiers = tips_in_clade
				} # end tip specifiers as ranks
			} # end tip specifiers

		#plottree=TRUE; fixnode_clades=NULL
		#tr=tr_w_fossil
		tr_w_fossil_added = add_fossil_randomly(tr=tr_w_fossil, tip_specifiers=tip_specifiers, fossil_name=fossil_name, max_age=max_age, min_age=min_age, constraint_type=constraint_type, stem_or_crown=stem_or_crown, plottree=plottree, fixnode_clades=NULL, brlen_of_side_branch=brlen_of_side_branch)

		# Check that the new tree is different
		if (length(tr_w_fossil_added$tip.label) == (length(tr_w_fossil$tip.label)+1) )
			{
			tr_w_fossil = tr_w_fossil_added
			list_of_fossil_tips_added = c(list_of_fossil_tips_added, fossil_name)
			} else {
			warning_txt = paste("\n\nWarning: fossil '", fossil_name, "' was NOT added. In the Excel fossils file, check constraints, age range, tip specifiers, etc.\n", sep="")
			cat(warning_txt)
			}
		} # END for (i in 1:numrows)
	
	tr_w_fossils_and_fossils_added_list = NULL
	tr_w_fossils_and_fossils_added_list$tr_w_fossil = tr_w_fossil
	tr_w_fossils_and_fossils_added_list$list_of_fossil_tips_added = list_of_fossil_tips_added
	
	return(tr_w_fossils_and_fossils_added_list)
	} # end add_fossils_from_xls_randomly()





add_fossil_randomly <- function(tr, tip_specifiers, fossil_name, max_age, min_age, constraint_type="side_branch", stem_or_crown="both", plottree=FALSE, fixnode_clades=NULL, fossil_age_within_constraints=TRUE, brlen_of_side_branch=1e-07)
	{
	defaults='
	tr=tr_w_fossil
	tip_specifiers
	fossil_name
	max_age
	min_age
	constraint_type
	stem_or_crown
	plottree=FALSE
	fixnode_clades=NULL
	'
	
	
	# 
	if (max_age < min_age)
		{
		errortxt = paste("\nERROR in add_fossil_randomly(): fossil '", fossil_name, "' max_age (", max_age, ") is less than fossil min_age (", min_age, ").\nmax_age must be greater than or equal to min_age.\n", sep="")
		stop(errortxt)
		}

	
	# Get the tree table
	trtable = prt(tr, printflag=FALSE)
	
	# Get the list of branches that are possible attachment points
	#print(head(trtable))
	#print(tail(trtable))
	
	
	#print("tip_specifiers:")
	#print(tip_specifiers)
	#print("stem_or_crown:")
	#print(stem_or_crown)
	branches = try(get_possible_branches_to_add_fossils_to(tr, trtable, tip_specifiers, stem_or_crown))
	
	if (class(branches) == "try-error")
		{
		txt = paste("\n\nWARNING in add_fossil_randomly(): In add_fossil_randomly(), 'get_possible_branches_to_add_fossils_to()' failed. The input tree is being returned. The fossil name is: '", fossil_name, " (min=", min_age, " / max=", max_age, "). tip_specifiers='", paste0(tip_specifiers, sep=","), "', stem_or_crown='", stem_or_crown, "'. Printing try-error...\n\n", sep="")
		cat(warning_txt)
		print("Printing 'branches' result (result was of class 'try-error'):")
		print(branches)
		return(tr)
		}	


	# Randomly draw the age of the fossil
	# Ages are times before present
	if (fossil_age_within_constraints == FALSE)
		{
		# Old strategy, 2014
		fossil_age = runif(n=1, min=min_age, max=max_age)
		} else {
		# New strategy, 2017-03-14
		# NOTE: INCLUSION OR NOT OF THE BRANCH BELOW THE CROWN NODE
		# IS ALL DECIDED in get_possible_branches_to_add_fossils_to()
		# Should we include the branch below the crown node??

		age_oldest_node = max(branches$time_bp)
# 		if (stem_or_crown == "crown")
# 			{
# 			# Don't include that branch!
# 			use_branch_below_this_node_TF = (branches$time_bp != age_oldest_node) == TRUE
# 			} else {
# 			# Do include that branch!
# 			use_branch_below_this_node_TF = (branches$time_bp <= age_oldest_node) == TRUE
# 			}
			
		# Filter the branches to include/exclude the branch below the crown node
# 		branches = branches[use_branch_below_this_node_TF, ]

		# Max age of each branch
		branch_min_ages = branches$time_bp
		
		temp_brlens = branches$edge.length
		temp_brlens[is.na(temp_brlens)] = 0
		branch_max_ages = branch_min_ages + temp_brlens
		oldest_branch_age = max(branch_max_ages)
		youngest_branch_age = min(branch_min_ages)
		
		# If the oldest branch age is younger than the oldest age of the fossil,
		# make max_age the oldest_branch_age
		original_max_age = max_age
		if (max_age > oldest_branch_age)
			{
			max_age = oldest_branch_age
			}
		
		# Check if the new max_age screws up the age range!
		if (max_age < min_age)
			{
			txt = paste0("WARNING in add_fossil_randomly(): For fossil '", fossil_name, "', your specified age range was min_age=", min_age, " Ma, max_age=", original_max_age, " Ma.  However, the date range of the available branches matching the constraints in your Excel file is only ", youngest_branch_age, "-", oldest_branch_age, " Ma.  This fossil will be skipped for producing this particular tree and geography file, but you should consider revising your tree and/or your Excel fossils file!")
			
			cat("\n")
			cat(txt)
			cat("\n")
			
			# Restore max_age for purposes of later script
			max_age = original_max_age
			} # END if (max_age < min_age)
		
		fossil_age = runif(n=1, min=min_age, max=max_age)
		} # END if (fossil_age_within_constraints == FALSE)


	# If you have NO branches that are possible attachment points for this tree,
	# then print WARNING, and return the input tree
	if ((nrow(branches) == 0) || (length(branches) == 0))
		{
		warning_txt = paste("\n\nWARNING: In add_fossil_randomly(), fossil has no possible branches to attach to\nin this tree, given the constraints and fossil date. The input tree is being returned.\nThe fossil name is: '", fossil_name, "', age=", fossil_age, " (min=", min_age, " / max=", max_age, ").", sep="")
		cat(warning_txt)
		return(tr)
		}


	# # Possible constraint types:
	# 1. side_branch
	# 2. direct_ancestor
	# 3. fixnode
	#constraint_type = "side_branch"
	
	
	
	txt = paste0("\nAttempting to add fossil '", fossil_name, "' as a constraint_type='", constraint_type, "', age=", fossil_age, " (min=", min_age, " / max=", max_age, ").")
	cat(txt)
	
	
	# Error check on constraint type
	constraint_type_options = c("side_branch", "direct_ancestor", "ancestor", "fixnode", "specifiers")
	TF = constraint_type %in% constraint_type_options
	
	if (TF == FALSE)
		{
		txt = paste0("STOP ERROR in add_fossil_randomly(): the input constraint_type for this row was '", constraint_type, "', but this is not an allowed constraint_type.")
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		}
	
	if (constraint_type == "side_branch")
		{
		# Get the branches that are possible attachment points
		# Any branch that extends below the fossil age, but is 
		# in the identified clade / branch list
		
		# Add random side branch
		tr_w_fossil = add_random_side_branch(tr, trtable, branches, fossil_name, fossil_age, plottree=plottree)
		}

	if ((constraint_type == "direct_ancestor") || (constraint_type == "ancestor") )
		{
		# Get the branches that are possible hosts of the ancestor hook
		# Any branch that extends below the fossil age, but has
		# min_age above fossil_age
		# in the identified clade / branch list
		tr_w_fossil = add_random_direct_ancestor_hook(tr, trtable, branches, fossil_name, fossil_age, plottree=plottree, brlen_of_side_branch=brlen_of_side_branch)
		}

	if (constraint_type == "fixnode")
		{
		# Then add to the list of fixnodes
		fixnode_clades = c(fixnode_clades, tip_specifiers)
		# Nothing to output here, really, do fixnodes in a different way
		# return(fixnode_clades)
		}
	
	return(tr_w_fossil)
	} # end add_fossil_randomly()


#' Get a list of branches that are possible 
#' attachment points for a fossil with a phylogenetic
#' position that is only approximately known
get_possible_branches_to_add_fossils_to <- function(tr, trtable, tip_specifiers, stem_or_crown="both")
	{
	#######################################################
	# Get the branches in the clade
	#######################################################
	#stem_or_crown = "both"	# could be stem, crown, or both
							# default (i.e. BEAST stem) is "both"

	# tip_specifiers must have a length of 1 or greater
	if (length(tip_specifiers) <= 0)
		{
		errortxt = "\n\nERROR in get_possible_branches_to_add_fossils_to(): length(tip_specifiers) cannot be 0 or less.\nPrinting tip_specifiers:\n\n"
		cat(errortxt)
		print(tip_specifiers)
		stop("ERROR caused stop.")
		} # end if (length(tip_specifiers) <= 0)

	# If there is just ONE tip_specifier, you have a terminal branch; getMRCA won't work
	# If there is just ONE tip_specifier, you have a terminal branch; the "stem_or_crown" cannot be "crown"
	if (length(tip_specifiers) == 1)
		{ 
		if (stem_or_crown == "crown")
			{
			errortxt = "\n\nERROR in get_possible_branches_to_add_fossils_to(): length(tip_specifiers)==1 and stem_or_crown=='crown'.\nCrown groups require two or more tips to be validly defined.\nValid 'stem_or_crown' options are 'stem' or 'both'.\nPrinting tip_specifiers:\n\n"
			cat(errortxt)
			print(tip_specifiers)
			stop("ERROR caused stop.")
			} # end if (stem_or_crown == "crown")
		
		# Otherwise, just return the branch for the tip
		TF = trtable$label == tip_specifiers
		nodenum_of_tip = trtable$node[TF]
		branches = trtable[nodenum_of_tip, ]
		return(branches)
		} # end if (length(tip_specifiers) == 1)
	
	nodenum_of_crown = getMRCA(phy=tr, tip=tip_specifiers)
	nodenum_of_crown

	# Get the branches in the crown group
	nodes = get_daughter_nodes(nodenum=nodenum_of_crown, tr=tr, nodes=NULL)
	nodes
	
	# Get the branches in the crown group, stem group, or both
	if (stem_or_crown == "crown")
		{
		branches = trtable[nodes,]
		}

	if (stem_or_crown == "both")
		{
		nodes = c(nodes, nodenum_of_crown)
		branches = trtable[nodes,]
		# And add the ancestor branch
		# 2017-05-16_BUG: this is getting the node BELOW the crown node;
		# we just need the branch ATTACHED to the crown node
		#stem_branch_nodenum = trtable$ancestor[nodenum_of_crown]
		#stem_branch = trtable[stem_branch_nodenum,]
		#branches = rbind(branches, stem_branch)
		
		# NOTE! We will want the nodenum_of_crown, so that we ONLY get the branch
		# that is the stem of 'nodenum_of_crown'
		}

	if (stem_or_crown == "stem")
		{
		# Use JUST the ancestor branch
		# 2017-05-16_BUG: this is getting the node BELOW the crown node;
		# we just need the branch ATTACHED to the crown node
		#stem_branch_nodenum = trtable$ancestor[nodenum_of_crown]
		#stem_branch = trtable[stem_branch_nodenum,]
		#branches = stem_branch
		branches = trtable[nodenum_of_crown,]
		
		# NOTE! We will want the nodenum_of_crown, so that we ONLY get the branch
		# that is the stem of 'nodenum_of_crown'
		}

		naTF = is.na(branches$node)
		if (sum(naTF) > 0)
			{
			cat("ERROR: get_possible_branches_to_add_fossils_to() produces NAs for \n", sep="")
			cat("trtable nodes: \n", sep="")
			print(nodes[naTF])
			
			cat("Causing this problem in the list of branches that are possible attachment points:\n", sep="")
			print(branches)
			stop()
			}
	
	
	return(branches)
	} # end get_possible_branches_to_add_fossils_to()



#' Given a fossil time and a clade, randomly add a 
#' side branch somewhere below the fossil to the clade
#' (uniformly on available branches; not perfect, but 
#'  workable)
add_random_side_branch <- function(tr, trtable, branches, fossil_name, fossil_age, plottree=FALSE)
	{
	# Get the branches that are possible attachment points
	# Any branch that extends below the fossil age, but is 
	# in the identified clade / branch list

	# The fossil cannot be the exact same age as any internal node
	internal_TF = trtable$node.type != "tip"
	matchTF = trtable$time_bp[internal_TF] == fossil_age
	if (sum(matchTF) > 0)
		{
		errormsg = "ERROR in add_random_side_branch(): the fossil tip you are adding is the exact same age\nas an internal node. This is not allowed.\n"
		cat("\n", errormsg, "\n", sep="")

		errormsg = "Returning the input tree...\n"
		cat("\n", errormsg, "\n", sep="")

		return(tr)
		}

	# Really, all that matters is that the max_age of the 
	# branch be older than the fossil
	
	# Max age of each branch
	branch_min_ages = branches$time_bp
	temp_brlens = branches$edge.length
	temp_brlens[is.na(temp_brlens)] = 0
	branch_max_ages = branch_min_ages + temp_brlens




	# Branches are possible attachment points for the side branch if
	# the fossil age is less than branch max_age
	fossil_age_LT_branch_max_age_TF = fossil_age < branch_max_ages
	branches_to_attach_to = branches[fossil_age_LT_branch_max_age_TF, ]

	# ERROR CHECK:
	if ((nrow(branches_to_attach_to) == 0) || length(branches_to_attach_to) == 0)
		{
		warning_txt = paste("\n\nWARNING: In add_random_side_branch(), fossil has no possible branches to attach to\nin this tree, given the constraints and fossil date. The input tree is being returned.\nThe fossil name is: '", fossil_name, "', age=", fossil_age, ".\n\n", sep="")
		cat(warning_txt)

		warning_txt = paste("The branches in the clade have these maximum ages:\n\n", sep="")
		cat(warning_txt)
		
		#print(branch_max_ages)

		warning_txt = paste("\n\nReturning the tree unchanged...\n", sep="")
		cat(warning_txt)
		
		return(tr)
		}

	
	# Make a new column, with available branch lengths (anything older than the fossil)
	branch_min_ages = branches_to_attach_to$time_bp
	temp_brlens = branches_to_attach_to$edge.length
	branch_max_ages = branch_min_ages + temp_brlens
	
	branch_lengths_available_to_attach_to = branches_to_attach_to$edge.length
	diff_between_fossil_and_branch_min_ages = branch_min_ages - fossil_age
	# Only when the result is negative, do you chop down the branch
	negative_TF = diff_between_fossil_and_branch_min_ages < 0

	branch_lengths_available_to_attach_to[negative_TF] = branch_lengths_available_to_attach_to[negative_TF] + diff_between_fossil_and_branch_min_ages[negative_TF]
	naTF = is.na(branch_lengths_available_to_attach_to)
	branch_lengths_available_to_attach_to[naTF] = 0
	branches_to_attach_to$attach_brlen = branch_lengths_available_to_attach_to
	
	# Proportions to sample from
	prop = branches_to_attach_to$attach_brlen / sum(branches_to_attach_to$attach_brlen)
	prop[is.na(prop)] = 0
	branchrow_to_attach_to = sample(1:length(prop), size=1, replace=FALSE, prob=prop)
	branchrow_to_attach_to
	
	# Uniformly sample a distance above the base of the branch
	height_above_base = runif(n=1, min=0, max=branches_to_attach_to$attach_brlen[branchrow_to_attach_to])

	# Pick the highest tip in the clade to count the age back from
	nodenum_at_top_of_branch = branches_to_attach_to$node[branchrow_to_attach_to]
	tipnums = get_daughter_tipnums(nodenum=nodenum_at_top_of_branch, tr=tr)
	tipnums
	highest_tip_in_clade_hts = trtable$node_ht[tipnums]
	highest_tip_in_clade_TF = highest_tip_in_clade_hts == max(highest_tip_in_clade_hts)
	tipnum_highest_tip_in_clade = tipnums[highest_tip_in_clade_TF][1]
	tipname_of_highest_tip_in_clade = trtable$label[tipnum_highest_tip_in_clade]
	
	# Time before present of the attachment point
	time_bp_attachment_point_absolute = trtable$time_bp[nodenum_at_top_of_branch] + trtable$edge.length[nodenum_at_top_of_branch] - height_above_base
	# Time before the highest tip above the attachment point
	time_bp_attachment_point_relative_to_tip = time_bp_attachment_point_absolute - trtable$time_bp[tipnum_highest_tip_in_clade]
	
	# Add the side branch
	length_of_side_branch = time_bp_attachment_point_absolute - fossil_age
	
	#print("length_of_side_branch:")
	#print(length_of_side_branch)
	
	# 
	tr_w_fossil = add_hook(t=tr, tipname=tipname_of_highest_tip_in_clade, brlen_of_side_branch=length_of_side_branch, depthtime=time_bp_attachment_point_relative_to_tip, plottree=FALSE, newtipname=fossil_name)
	
	if (plottree == TRUE)
		{
		if (length(tr_w_fossil$tip.label) > 100)
			{
			plot(tr_w_fossil, show.tip.label=FALSE)
			axisPhylo()
			titletxt = paste("Tree with ", fossil_name, " added", sep="")
			title(titletxt)
			} else {
			plot(tr_w_fossil, show.tip.label=TRUE)
			axisPhylo()			
			titletxt = paste("Tree with ", fossil_name, " added", sep="")
			title(titletxt)
			}
		}
	
	return(tr_w_fossil)
	} # end add_random_side_branch()


#' Given a time and a clade, randomly add a
#' direct ancestor (in the form of a hook)
add_random_direct_ancestor_hook <- function(tr, trtable, branches, fossil_name, fossil_age, plottree=FALSE, brlen_of_side_branch=1e-07)
	{
	# Get the branches that are possible hosts of the ancestor hook
	# Any branch that extends below the fossil age, but has
	# min_age above fossil_age
	# in the identified clade / branch list

	# The fossil cannot be the exact same age as any internal node
	internal_TF = trtable$node.type != "tip"
	matchTF = trtable$time_bp[internal_TF] == fossil_age
	if (sum(matchTF) > 0)
		{
		errormsg = "ERROR in add_random_direct_ancestor_hook(): the fossil tip you are adding is the exact same age\nas an internal node. This is not allowed.\n"
		cat("\n", errormsg, "\n", sep="")

		errormsg = "Returning the input tree...\n"
		cat("\n", errormsg, "\n", sep="")

		return(tr)
		}

	# Really, all that matters is that the max_age of the 
	# branch be older than the fossil

	# Max age of each branch
	branch_min_ages = branches$time_bp
	temp_brlens = branches$edge.length
	temp_brlens[is.na(temp_brlens)] = 0
	branch_max_ages = branch_min_ages + temp_brlens

	# Branches are possible attachment points for the ancestor if
	# the fossil age is less than branch max_age
	fossil_age_LT_branch_max_age_TF = fossil_age < branch_max_ages
	fossil_age_GT_branch_min_age_TF = fossil_age > branch_min_ages
	branches_containing_fossil_age_TF = (fossil_age_LT_branch_max_age_TF + fossil_age_GT_branch_min_age_TF) == 2
	branches_to_attach_to = branches[branches_containing_fossil_age_TF, ]

	# ERROR CHECK:
	if ((nrow(branches_to_attach_to) == 0) || length(branches_to_attach_to) == 0)
		{
		warning_txt = paste("\n\nWARNING: In add_random_direct_ancestor_hook(), fossil has no possible branches to attach to\nin this tree, given the constraints and fossil date. The input tree is being returned.\nThe fossil name is: '", fossil_name, "', age=", fossil_age, ".\n\n", sep="")
		cat(warning_txt)

		warning_txt = paste("The branches in the clade have these maximum/minimum ages:\n\n", sep="")
		cat(warning_txt)
		
		print(cbind(branch_max_ages,branch_min_ages))

		warning_txt = paste("\n\nReturning the tree unchanged...\n\n", sep="")
		cat(warning_txt)
		
		return(tr)
		}
	
	# Make a new column, with available branch lengths (anything older than the fossil)
	branch_min_ages = branches_to_attach_to$time_bp
	temp_brlens = branches_to_attach_to$edge.length
	branch_max_ages = branch_min_ages + temp_brlens
	
	branch_lengths_available_to_attach_to = branches_to_attach_to$edge.length
	naTF = is.na(branch_lengths_available_to_attach_to)
	branch_lengths_available_to_attach_to[naTF] = 0
	branches_to_attach_to$attach_brlen = branch_lengths_available_to_attach_to
	
	# Proportions to sample from -- EQUAL probability for each branch crossed
	prop = branches_to_attach_to$attach_brlen / sum(branches_to_attach_to$attach_brlen)
	prop[is.na(prop)] = 0
	#print(prop)
	branchrow_to_attach_to = sample(1:length(prop), size=1, replace=FALSE, prob=prop)
	branchrow_to_attach_to
	
	# The distance above the base of the branch is calculated
	height_above_base = branches_to_attach_to$time_bp[branchrow_to_attach_to] + branches_to_attach_to$edge.length[branchrow_to_attach_to] - fossil_age
	height_above_base
	fossil_age
	
	# Pick the highest tip in the clade to count the age back from
	nodenum_at_top_of_branch = branches_to_attach_to$node[branchrow_to_attach_to]
	tipnums = get_daughter_tipnums(nodenum=nodenum_at_top_of_branch, tr=tr)
	tipnums
	highest_tip_in_clade_hts = trtable$node_ht[tipnums]
	highest_tip_in_clade_TF = highest_tip_in_clade_hts == max(highest_tip_in_clade_hts)
	tipnum_highest_tip_in_clade = tipnums[highest_tip_in_clade_TF][1]
	tipname_of_highest_tip_in_clade = trtable$label[tipnum_highest_tip_in_clade]
	
	# Time before present of the attachment point
	time_bp_attachment_point_absolute = trtable$time_bp[nodenum_at_top_of_branch] + trtable$edge.length[nodenum_at_top_of_branch] - height_above_base
	# Time before the highest tip above the attachment point
	time_bp_attachment_point_relative_to_tip = time_bp_attachment_point_absolute - trtable$time_bp[tipnum_highest_tip_in_clade]
	
	# Add the side branch (as a hook, brlen=1e-07)
	tr_w_fossil = add_hook(t=tr, tipname=tipname_of_highest_tip_in_clade, brlen_of_side_branch=brlen_of_side_branch, depthtime=time_bp_attachment_point_relative_to_tip, plottree=FALSE, newtipname=fossil_name)
	
	if (plottree == TRUE)
		{
		if (length(tr_w_fossil$tip.label) > 100)
			{
			plot(tr_w_fossil, show.tip.label=FALSE)
			axisPhylo()
			titletxt = paste("Tree with ", fossil_name, " added", sep="")
			title(titletxt)
			} else {
			plot(tr_w_fossil, show.tip.label=TRUE)
			axisPhylo()			
			titletxt = paste("Tree with ", fossil_name, " added", sep="")
			title(titletxt)
			}
		}
	
	return(tr_w_fossil)	
	} # end add_random_direct_ancestor_hook()





get_daughter_tipnums <- function(nodenum, tr)
	{
	nodes = get_daughter_nodes(nodenum, tr, nodes=NULL)
	tips_TF = nodes <= length(tr$tip.label)
	tipnums = nodes[tips_TF]
	return(tipnums)
	}


get_daughter_tipnames <- function(nodenum, tr)
	{
	tipnums = get_daughter_tipnums(nodenum, tr)
	tipnames = tr$tip.label[tipnums]
	return(tipnames)
	}



get_daughter_nodes <- function(nodenum, tr, nodes=NULL)
	{
	if(is.null(nodes))
		{
		nodes = vector()
		}
	daughter_nodes = tr$edge[which(tr$edge[,1]==nodenum),2]
	
	# Error check, in case the starting nodenum is a tip
	if ((length(daughter_nodes) == 0) && (length(nodes)==0))
		{
		nodes = c(nodes, nodenum)
		} else {
		nodes = c(nodes, daughter_nodes)
		}
		
	internal_nodes_indices = which(daughter_nodes > length(tr$tip.label))
	if(length(internal_nodes_indices) > 0)
		{
		for (i in 1:length(internal_nodes_indices))
			{
			nodes = get_daughter_nodes(nodenum=daughter_nodes[internal_nodes_indices[i]], tr=tr, nodes=nodes)
			}
		}
	return(nodes)
	}




trace_parents_up <- function(nodenum, t, depthtime)
	{
	# Trace from a node up to its parents etc., a specified distance
	parent_node = get_parent_for_trace_parents_up(nodenum, t)
	
	# print nodenum
	#print(nodenum)
	#print(parent_node)
	#print(depthtime)
	
	length_to_parent = t$edge.length[t$edge[,2] == nodenum]
	#cat("length_to_parent: ", length_to_parent, ", length=", length(length_to_parent), "\n", sep="")
	#cat("depthtime: ", depthtime, "\n", sep="")
	if (length(length_to_parent) == 0)
		{
		print("ERROR: trace_parents_up() -- no length_to_parent returned, probably overshot bottom of tree")
		return(NA)
		}
	if (length_to_parent == depthtime)
		{
		print("ERROR: trace_parents_up() doesn't want to find an EXACT match between depthtime and a node; this will lead to problems in hook addition!")
		}
	if (length_to_parent > depthtime)
		{
		# you're done!
		return(nodenum)
		}
	else
		{
		# burrow up to parents
		depthtime = depthtime - length_to_parent
		parent_node = trace_parents_up(parent_node, t, depthtime)
		return(parent_node)
		}
	return(nodenum)
	}


get_parent_for_trace_parents_up <- function(nodenum, t, printflag=FALSE)
	{
	matching_edges = findall(nodenum, t$edge[,2])
	parent_nodenum = t$edge[,1][matching_edges][1]
	
	if (printflag)
		{
		print(paste("nodenum=", nodenum, " parent_nodenum=", parent_nodenum, sep=""))
		}
	if (is.na(parent_nodenum))
		{
		if (printflag)
			{
			print(paste("get_parent(): node ", nodenum, " has no parent, it's probably the root!\nAnd you missed whatever parent you were actually trying to find!", sep=""))
			}
		}
	return(parent_nodenum)
	}
	

dist_between_direct_ancestors <- function(ancestor_node, descendant_node, t, totaldist=0, printflag=FALSE)
	{
	# Recursive algorithm to get distance between descendent and ancestor
	
	# Error trap should operate before this
	
	if (ancestor_node == descendant_node)
		{
		if (printflag)
			{
			print("dist_between_direct_ancestors(): ancestor_node == descendant_node")
			}
		return(totaldist)
		}
	
	parent_node = get_parent(descendant_node, t)
	dist_to_parent = t$edge.length[t$edge[,2] == descendant_node]
	
	totaldist = dist_to_parent + totaldist
	
	#print(paste(parent_node, ancestor_node, sep=""))
	if (parent_node == ancestor_node)
		{
		return(totaldist)
		}
	else
		{
		totaldist = dist_between_direct_ancestors(ancestor_node, parent_node, t, totaldist)
		return(totaldist)
		}
	}

#' t			= tree
#' tipname		= tip to add hook below
#' depthtime	= time_bp_attachment_point_relative_to_tip  This is the time RELATIVE to the tip
add_hook <- function(t, tipname, depthtime, brlen_of_side_branch=0.0000001, plottree = FALSE, printflag=0, newtipname="default")
	{
	# Add a hook (a small side tip) to a phylogeny
	#
	# e.g.:
	# Do spatial variogram by doing points from many different species
	# add tips to tree
	#cat("owls(((Strix_aluco:4.2,Asio_otus:4.2):3.1,Athene_noctua:7.3):6.3,Tyto_alba:13.5);", file = "ex.tre", sep = "\n")
	#t <- read.tree("ex.tre")
	#prt(t)
	#newtree = add_hook(t, brlen_of_side_branch=0.0000001, plottree = TRUE)

	if (printflag >= 2)
		{
		cat("add_hook(): Adding below tipname ",  tipname, "\n", sep="")
		}

	newtree = t
	#daughter_nodenum_to_add_hook_below = 4
	#height_below_daughter_at_which_to_add_it = 1.0
	#brlen_of_side_branch=0.0000001
	
	# find the node, below this you will add the 
	tip_nodenum = which(t$tip.label == tipname)
	
	if (printflag >= 2)
		{
		print(paste("addhook(): tip_nodenum = ", tip_nodenum, sep=""))
		}

	daughter_nodenum_to_add_hook_below = trace_parents_up(tip_nodenum, t, depthtime)
	
	# Error trap in case of e.g. overrun of bottom of tree
	if (is.na(daughter_nodenum_to_add_hook_below))
		{
		cat("add_hook() error: daughter_nodenum_to_add_hook_below is NA, probably overshot bottom of tree (?); not adding hook\n")
		return(t)
		}
	
	if (printflag >= 2)
		{
		print(paste("addhook(): internal nodenum to add below = ", daughter_nodenum_to_add_hook_below, sep=""))
		}
	
	tip_to_ancestor_node_dist = dist_between_direct_ancestors(daughter_nodenum_to_add_hook_below, tip_nodenum, t, totaldist=0)
	
	
	
	height_below_daughter_at_which_to_add_it = depthtime - tip_to_ancestor_node_dist
	
	
	# add a new tip to the list of tips (this is the hook)
	new_tip_nodenum = get_nodenum_structural_root(t)
	
	# bump up all of the nodenums above the node tip by 1
	newtree$edge[t$edge >= new_tip_nodenum] = t$edge[t$edge >= new_tip_nodenum] + 1
	
	# add a new internal node at the end
	new_inNode = max(newtree$edge) + 1
	
	# add two new edges, and replace the old edge
	#print(t$edge[,2])
	#print(daughter_nodenum_to_add_hook_below)
	#print(t$edge[,2] == daughter_nodenum_to_add_hook_below)
	old_edge_num = which(t$edge[,2] == daughter_nodenum_to_add_hook_below)
	
	# extract the edgenums before and after this insertion (exceptions for in case
	# if the modified row is the first or last row)
	if (old_edge_num == 1)
		{
		first_old_edges_rownums = NULL
		} else {
		first_old_edges_rownums = 1:(old_edge_num-1) #newtree$edge[1:(old_edge_num-1), ]
		}
	if (old_edge_num == nrow(t$edge))
		{
		second_old_edges_rownums = NULL
		} else {
		second_old_edges_rownums = (old_edge_num+1):nrow(t$edge) # newtree$edge[, ]
		}
	
	
	
	# replace the edge, keeping the old parent (which may have increased by 1! use newtree!!), put the new internal node as the daughter)
	replacement_edge_row = newtree$edge[old_edge_num, ]
	replacement_edge_row[2] = new_inNode
	
	# subtract the distance below the daughter, from the top
	replacement_edge_length = t$edge.length[old_edge_num] - height_below_daughter_at_which_to_add_it
	
	
	# make the new edge, which goes below the old daughter node
	# you have to bump the daughter_nodenum_to_add_hook_below if it is
	# >= to the new_tip_nodenum
	if (daughter_nodenum_to_add_hook_below >= new_tip_nodenum)
		{
		daughter_nodenum_to_add_hook_below = daughter_nodenum_to_add_hook_below + 1
		}
	new_edge_below_old_daughter_node = c(new_inNode, daughter_nodenum_to_add_hook_below)
	new_edge_below_old_daughter_node_edge_length = height_below_daughter_at_which_to_add_it
	
	# make the new edge, which goes below the new tip: c(parent, daughter)
	new_edge_below_new_tip = c(new_inNode, new_tip_nodenum)
	new_edge_below_new_tip_edge_length = brlen_of_side_branch
	
	
	# add the edge rows before the one that is replaced, then the replaced edge, then the other old edges, then the 2 new edges
	new_edge_table = rbind(newtree$edge[first_old_edges_rownums, ], replacement_edge_row, newtree$edge[second_old_edges_rownums, ], new_edge_below_old_daughter_node, new_edge_below_new_tip)
	
	new_edgelength_list = c(t$edge.length[first_old_edges_rownums], replacement_edge_length, t$edge.length[second_old_edges_rownums], new_edge_below_old_daughter_node_edge_length, new_edge_below_new_tip_edge_length)
	
	# it MAY be important that the node numbers be INTEGER, not NUMERIC
	newtree$edge = matrix(as.integer(new_edge_table), ncol=2, byrow=FALSE)
	
	#row.names(newtree$edge) = NULL
	#newtree$edge[,1] = as.integer(newtree$edge[,1])
	#newtree$edge[,2] = as.integer(newtree$edge[,2])
	
	newtree$edge.length = new_edgelength_list
	#row.names(newtree$edge.length) = NULL
	
	# update number of internal nodes
	newtree$Nnode = t$Nnode + 1
	
	# add the new tip to the end of the list of tips
	if (newtipname == "default")
		{
		newtipname = paste("hook", new_tip_nodenum, sep="")
		} else {
		newtipname = newtipname
		}
	newtree$tip.label = c(t$tip.label, newtipname)
	
	if (printflag >= 2)
		{
		cat("\nAdding tip: ",  newtipname, sep="")
		}	
	
	# some crap to fix the tree formatting somehow
	# I mean, really, the tree was fucking logically correct, but
	# hanging plot and dist.nodes and probably anything
	# using reorder(phy, "pruningwise"), but I couldn't figure out why
	# I guess the order of the tips is important for some reason?
	# like maybe leftmost tip is leftmost in branching?
	# wtf kind of data architecture is this?
	# anyway, FUCK IT, writing to Newick and reading back fixes it.
	newtree = reorder(newtree)
	#tmpfn = "tmp_junktree.tree"
	#write.tree(newtree, tmpfn)
	newtree2 = read.tree(file="", text=write.tree(newtree, file=""))

	
	#cat("\nDone adding tips.\n")
	
	# plot, if desired:
	if (plottree == TRUE)
		{
		cat("add_hook(): plotting/printing the resulting tree...\n", sep="")
		prt(newtree)
		plot(newtree)
		}
	
	return(newtree2)
	} # END add_hook <- function(t, tipname, depthtime, brlen_of_side_branch=0.0000001, plottree = FALSE, printflag=0, newtipname="default")





################################################################################
# TREE MODIFICATION FUNCTIONS (e.g. adding hooks, choosing certain branches)
################################################################################
# newtipname="default" means new tips get "hook" and then the node number, e.g. hook414
# newtipname="tipnames" means new tips get either 
#    tipname_age      (for hooks on tip branches), or 
#    hook_tipname_age (for hooks on internal branches; just one of the tips is listed)
add_hooks <- function(tr, list_of_times_before_present, brlen_of_side_branch=0.0000001, plottree=FALSE, printflag=0, newtipnames="default")
	{
	# Take a list of ages, add hooks to any branch existing at that age

	# OK, RE-FUCKING DO
	# Gather a list of tip labels, and then record the distance below those tips that
	# you would go down to attach a hook
	hooktree = tr
	list_of_daughter_tipnames_to_add_hooks_below = c()
	list_of_ages_below_daughter = c()   # assumes ultrametric
	list_of_ages_absolute = c()
	list_of_new_tipnames = c()
	ntips = length(hooktree$tip.label)
	
	# Gather the list of tip labels, distance below each one, and new tip labels
	for (i in 1:length(list_of_times_before_present))
		{
		
		# Get the edges that exist at the time_slice in question
		time_slice = as.numeric(list_of_times_before_present[i])
		edge_times_bp = get_edge_times_before_present(hooktree)
		edges_that_exist_in_the_right_time = edges_existing_at_correct_time_bp_TF(time_slice, edge_times_bp)
		
		# get the nodes daughter to the branches that match
		nodenums_to_add_hooks_to = hooktree$edge[,2][edges_that_exist_in_the_right_time]

		if (printflag >= 1.5)
			{
			txt = paste0("Adding ", length(nodenums_to_add_hooks_to), " hooks for time #", i, "/", length(list_of_times_before_present), ": ", time_slice, " m.y.a.")
			cat("\n")
			cat(txt)
			}


		# calculate the times parent to these daughters at which to insert the hooks		
		#times_before_daughter_nodes = time_slice - edge_times_bp[edges_that_exist_in_the_right_time, 2]
		
		# trace these nodes to their tips in the (UNALTERED ORIGINAL) tree
# 		if (printflag >= 1.5)
# 			{
# 			cat("\n\tj(nodenum): ")
# 			}
		for (j in 1:length(nodenums_to_add_hooks_to))
			{
			nodenum = nodenums_to_add_hooks_to[j]

			if (printflag >= 2)
				{
				txt = paste0(j, ":", nodenum, " ")
				cat(txt)
				}

			
			# If the node is a tip
			if (nodenum <= ntips)
				{
				daughter_tipname_to_add_hooks_below = hooktree$tip.label[nodenum]
				
				edgenums = 1:nrow(edge_times_bp)
				edgenum = edgenums[tr$edge[,2] == nodenum]
				node_time_of_daughter = edge_times_bp[edgenum, 2]
				
				time_before_daughter_nodes = time_slice - node_time_of_daughter
				
				list_of_daughter_tipnames_to_add_hooks_below = c(list_of_daughter_tipnames_to_add_hooks_below, daughter_tipname_to_add_hooks_below)
				list_of_ages_below_daughter = c(list_of_ages_below_daughter, time_before_daughter_nodes)
				list_of_ages_absolute = c(list_of_ages_absolute, time_slice)
				list_of_new_tipnames = c(list_of_new_tipnames, daughter_tipname_to_add_hooks_below)
				} else {
				# Get *a* tip to compare time slice to
				temp_tipnames = get_all_daughter_tips_of_a_node(nodenum, hooktree)
				temp_tipname = temp_tipnames[1]
				namenums = 1:length(hooktree$tip.label)
				temp_tipnum = namenums[hooktree$tip.label == temp_tipname]
				list_of_daughter_tipnames_to_add_hooks_below = c(list_of_daughter_tipnames_to_add_hooks_below, hooktree$tip.label[temp_tipnum])
				
				#if (newtipnames == "tipnames")
				#	{
				hook_name_base = paste0("nodeBelow_", hooktree$tip.label[temp_tipnum])
				list_of_new_tipnames = c(list_of_new_tipnames, hook_name_base)
				#	list_of_daughter_tipnames_to_add_hooks_below = c(list_of_daughter_tipnames_to_add_hooks_below, hook_name_base)
				#	} # END if (newtipnames == "tipnames")


				#
				edgenums = 1:nrow(edge_times_bp)
				edgenum = edgenums[tr$edge[,2] == temp_tipnum]
				node_time_of_daughter = edge_times_bp[edgenum, 2]
				
				time_before_daughter_nodes = time_slice - node_time_of_daughter

				list_of_ages_below_daughter = c(list_of_ages_below_daughter, time_before_daughter_nodes)
				list_of_ages_absolute = c(list_of_ages_absolute, time_slice)
				} # END if (nodenum <= ntips)
			
			if (printflag >= 3)
				{
				print("print(list_of_ages_below_daughter):")
				print(list_of_ages_below_daughter)
				}
			}
		}
		

	if (printflag >= 1.5)
		{
		txt = paste0("\nAdding ", length(list_of_daughter_tipnames_to_add_hooks_below), " hooks to tree: ")
		cat(txt)
		}

	
	# Now, attach the hooks
	for (i in 1:length(list_of_daughter_tipnames_to_add_hooks_below))
		{
		#print(paste("i=", i, sep=""))
		tipname = list_of_daughter_tipnames_to_add_hooks_below[i]
		depthtime = as.numeric(list_of_ages_below_daughter[i])
		
		if (printflag >= 3)
			{
			print("print(depthtime):")
			print(depthtime)
			print("print(list_of_ages_below_daughter):")
			print(list_of_ages_below_daughter)
			}

		if (printflag >= 1.5)
			{
			txt = paste0(i, " ")
			cat(txt)
			}
			
		if (newtipnames == "default")
			{
			hooktree = add_hook(hooktree, tipname, depthtime, plottree=plottree, printflag=printflag, newtipname=newtipnames)
			}
		if (newtipnames == "tipnames")
			{
			new_name = paste0(list_of_new_tipnames[i], "_", list_of_ages_absolute[i])
			hooktree = add_hook(hooktree, tipname, depthtime, plottree=plottree, printflag=printflag, newtipname=new_name)
			}
		}
	
	return(hooktree)
	}


# Needed for add_hooks
edges_existing_at_correct_time_bp_TF <- function(time_slice, edge_times_bp, roundto=5)
	{
	# find the edges that exist in the right time
	
	# (note: round to a default of 6 digits (a single year) with roundto; this is 
	#  important for whether or not lineages exist at time=0 before present)

	# timepoint is younger or equal to the oldest end of the branch
	edges_that_start_below_time = round(edge_times_bp[, 1], digits=roundto) > time_slice
	
	# timepoint is older than the youngest end of the branch	
	edges_that_end_after_time = round(edge_times_bp[, 2], digits=roundto) <= time_slice
	edges_that_exist_in_the_right_time = edges_that_start_below_time + edges_that_end_after_time == 2
	return(edges_that_exist_in_the_right_time)
	}









