#######################################################
# Installing BioGeoBEARS from GitHub latest version
#######################################################
# CUT 2018: The old instructions to source() online upgrade .R files have been deleted,
#         all updates are now on the GitHub version of the package, version 1.1+
#######################################################
# Paste the stuff below, INSIDE the single-quote (') marks
# but NOT the single-quote marks themselves.
install_cmds_that_work_as_of_2023='

# Installation-from-scratch commands
install.packages("devtools")
install.packages("ape")
install.packages("optimx")
install.packages("GenSA")
install.packages("rexpokit")   
install.packages("cladoRcpp")
install.packages("snow")
install.packages("MultinomialCI")

library(devtools)
devtools::install_github(repo="nmatzke/BioGeoBEARS")

# NOTE: If you get a message like this
# * select "2. CRAN packages only" on "3. None"
# * If you get asked about "binaries" vs. "source", choose "binaries" 
#   (binaries are precompiled and easy to install; installing from source
#    requires that your computer have the correct compilers, which can be
#    challenging if you are not fairly expert)
' # END install_cmds_that_work_as_of_2023

#######################################################
#######################################################

#######################################################
# SETUP -- libraries/BioGeoBEARS updates
#######################################################

# Load the package (after installation, see above).
library(ape)
library(optimx)   # optimx seems better than R's default optim()
library(GenSA)    # GenSA seems better than optimx (but slower) on 5+ parameters, 
                  # seems to sometimes fail on simple problems (2-3 parameters)
library(rexpokit)
library(cladoRcpp)
#library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)
library(BioGeoBEARS)






#######################################################
# SETUP: YOUR WORKING DIRECTORY
#######################################################
# You will need to set your working directory to match your local system

# Note these very handy functions!
# Command "setwd(x)" sets your working directory
# Command "getwd()" gets your working directory and tells you what it is.
# Command "list.files()" lists the files in your working directory
# To get help on any command, use "?".  E.g., "?list.files"

# Set your working directory for output files
# default here is your home directory ("~")
# Change this as you like
# wd = "/drives/GDrive/z_help/Furnariidae_Bio_Geo/BSM_small/M3areas_allowed/"
# setwd(wd)

# Double-check your working directory with getwd()
getwd()

#######################################################
# SETUP: Extension data directory
#######################################################
# When R packages contain extra files, they are stored in the "extdata" directory 
# inside the installed package.
#
# BioGeoBEARS contains various example files and scripts in its extdata directory.
# 
# Each computer operating system might install BioGeoBEARS in a different place, 
# depending on your OS and settings. 
# 
# However, you can find the extdata directory like this:
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
extdata_dir
list.files(extdata_dir)


wd = slashslash(paste0(addslash(extdata_dir), "/examples/BSM_3taxa/M3areas_allowed/"))
setwd(wd)


# "system.file" looks in the directory of a specified package (in this case BioGeoBEARS)
# The function "np" is just a shortcut for normalizePath(), which converts the 
# path to the format appropriate for your system (e.g., Mac/Linux use "/", but 
# Windows uses "\\", if memory serves).

# Even when using your own data files, you should KEEP these commands in your 
# script, since the plot_BioGeoBEARS_results function needs a script from the 
# extdata directory to calculate the positions of "corners" on the plot. This cannot
# be made into a straight up BioGeoBEARS function because it uses C routines 
# from the package APE which do not pass R CMD check for some reason.

#######################################################
# SETUP: YOUR TREE FILE AND GEOGRAPHY FILE
#######################################################
# Example files are given below. To run your own data,
# make the below lines point to your own files, e.g.
# trfn = "/mydata/frogs/frogBGB/tree.newick"
# geogfn = "/mydata/frogs/frogBGB/geog.data"

#######################################################
# Phylogeny file
# Notes: 
# 1. Must be binary/bifurcating: no polytomies
# 2. No negative branchlengths (e.g. BEAST MCC consensus trees sometimes have negative branchlengths)
# 3. Be careful of very short branches, as BioGeoBEARS will interpret ultrashort branches as direct ancestors
# 4. You can use non-ultrametric trees, but BioGeoBEARS will interpret any tips significantly below the 
#    top of the tree as fossils!  This is only a good idea if you actually do have fossils in your tree,
#    as in e.g. Wood, Matzke et al. (2013), Systematic Biology.
# 5. The default settings of BioGeoBEARS make sense for trees where the branchlengths are in units of 
#    millions of years, and the tree is 1-1000 units tall. If you have a tree with a total height of
#    e.g. 0.00001, you will need to adjust e.g. the max values of d and e, or (simpler) multiply all
#    your branchlengths to get them into reasonable units.
# 6. DON'T USE SPACES IN SPECIES NAMES, USE E.G. "_"
#######################################################
# This is the example Newick file for Hawaiian 3taxa
# (from Ree & Smith 2008)
# "trfn" = "tree file name"
trfn = "tree.newick"

# Look at the raw Newick file:
moref(trfn)

dev.off()
dev.off()
dev.off()

# Look at your phylogeny:
pdffn = "tree.pdf"
pdf(file=pdffn, width=6, height=6)

tr = read.tree(trfn)
tr
plot(tr)
title("Example 3-taxon tree")
axisPhylo() # plots timescale
dev.off()

#######################################################
# Geography file
# Notes:
# 1. This is a PHYLIP-formatted file. This means that in the 
#    first line, 
#    - the 1st number equals the number of rows (species)
#    - the 2nd number equals the number of columns (number of areas)
#    - after a tab, put the areas in parentheses, with spaces: (A B C D)
#
# 1.5. Example first line:
#    10    4    (A B C D)
# 
# 2. The second line, and subsequent lines:
#    speciesA    0110
#    speciesB    0111
#    speciesC    0001
#         ...
# 
# 2.5a. This means a TAB between the species name and the area 0/1s
# 2.5b. This also means NO SPACE AND NO TAB between the area 0/1s.
# 
# 3. See example files at:
#    http://phylo.wikidot.com/biogeobears#files
# 
# 4. Make you understand what a PLAIN-TEXT EDITOR is:
#    http://phylo.wikidot.com/biogeobears#texteditors
#
# 3. The PHYLIP format is the same format used for C++ LAGRANGE geography files.
#
# 4. All names in the geography file must match names in the phylogeny file.
#
# 5. DON'T USE SPACES IN SPECIES NAMES, USE E.G. "_"
#
# 6. Operational taxonomic units (OTUs) should ideally be phylogenetic lineages, 
#    i.e. genetically isolated populations.  These may or may not be identical 
#    with species.  You would NOT want to just use specimens, as each specimen 
#    automatically can only live in 1 area, which will typically favor DEC+J 
#    models.  This is fine if the species/lineages really do live in single areas,
#    but you wouldn't want to assume this without thinking about it at least. 
#    In summary, you should collapse multiple specimens into species/lineages if 
#    data indicates they are the same genetic population.
######################################################

# This is the example geography file for Hawaiian 3taxa
# (from Ree & Smith 2008)
geogfn = "geog.data"

# Look at the raw geography text file:
moref(geogfn)

# Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

# Maximum range size observed:
max(rowSums(dfnums_to_numeric(tipranges@df)))

# Set the maximum number of areas any species may occupy; this cannot be larger 
# than the number of areas you set up, but it can be smaller.
max_range_size = 3

####################################################
####################################################
# KEY HINT: The number of states (= number of different possible geographic ranges)
# depends on (a) the number of areas and (b) max_range_size.
# If you have more than about 500-600 states, the calculations will get REALLY slow,
# since the program has to exponentiate a matrix of e.g. 600x600.  Often the computer
# will just sit there and crunch, and never get through the calculation of the first
# likelihood.
# 
# (this is also what is usually happening when LAGRANGE hangs: you have too many states!)
#
# To check the number of states for a given number of ranges, try:
numstates_from_numareas(numareas=4, maxareas=4, include_null_range=TRUE)
numstates_from_numareas(numareas=4, maxareas=4, include_null_range=FALSE)
numstates_from_numareas(numareas=4, maxareas=3, include_null_range=TRUE)
numstates_from_numareas(numareas=4, maxareas=2, include_null_range=TRUE)

# Large numbers of areas have problems:
numstates_from_numareas(numareas=10, maxareas=10, include_null_range=TRUE)

# ...unless you limit the max_range_size:
numstates_from_numareas(numareas=10, maxareas=2, include_null_range=TRUE)
####################################################
####################################################

#######################################################
#######################################################
# DEC AND DEC+J ANALYSIS
#######################################################
#######################################################
# NOTE: The BioGeoBEARS "DEC" model is identical with 
# the Lagrange DEC model, and should return identical
# ML estimates of parameters, and the same 
# log-likelihoods, for the same datasets.
#
# Ancestral state probabilities at nodes will be slightly 
# different, since BioGeoBEARS is reporting the 
# ancestral state probabilities under the global ML
# model, and Lagrange is reporting ancestral state
# probabilities after re-optimizing the likelihood
# after fixing the state at each node. These will 
# be similar, but not identical. See Matzke (2014),
# Systematic Biology, for discussion.
#
# Also see Matzke (2014) for presentation of the 
# DEC+J model.
#######################################################
#######################################################

#######################################################
#######################################################

#######################################################
# Run DEC
#######################################################

# Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
# 1. Here, un-comment ONLY the files you want to use.
# 2. Also un-comment "BioGeoBEARS_run_object = section_the_tree(...", below.
# 3. For example files see (a) extdata_dir, 
#  or (b) http://phylo.wikidot.com/biogeobears#files
#  and BioGeoBEARS Google Group posts for further hints)
#
# Uncomment files you wish to use in time-stratified analyses:
BioGeoBEARS_run_object$timesfn = "timeperiods_C2.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed_noC2.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1
# (use more cores to speed it up; this requires
# library(parallel) and/or library(snow). The package "parallel" 
# is now default on Macs in R 3.0+, but apparently still 
# has to be typed on some Windows machines. Note: apparently 
# parallel works on Mac command-line R, but not R.app.
# BioGeoBEARS checks for this and resets to 1
# core with R.app)

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
# I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
# but the results are imprecise and so I haven't explored it further.
# In a Bayesian analysis, it might work OK, but the ML point estimates are
# not identical.
# Also, I have not implemented all functions to work with force_sparse=TRUE.
# Volunteers are welcome to work on it!!
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE, fossils_older_than=0.001, cut_fossils=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DEC model
# (nothing to do; defaults)

# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

# This contains the model object
BioGeoBEARS_run_object$BioGeoBEARS_model_object

# This table contains the parameters of the model 
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

# Run this to check inputs. Read the error messages if you get them!
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "3taxa_DEC_M3_areasallowed_v1.Rdata"
if (runslow)
    {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    

    save(res, file=resfn)
    resDEC = res
    } else {
    # Loads to "res"
    load(resfn)
    resDEC = res
    }

#######################################################
# Run DEC+J
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
BioGeoBEARS_run_object$timesfn = "timeperiods_C2.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed_noC2.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE, fossils_older_than=0.001, cut_fossils=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DEC+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "3taxa_DEC+J_M3_areasallowed_v1.Rdata"
runslow = TRUE
if (runslow)
    {
    #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")

    res = bears_optim_run(BioGeoBEARS_run_object)
    res    

    save(res, file=resfn)

    resDECj = res
    } else {
    # Loads to "res"
    load(resfn)
    resDECj = res
    }

#######################################################
# PDF plots
#######################################################
pdffn = "3taxa_DEC_vs_DEC+J_M3_areasallowed_v1.pdf"
pdf(pdffn, width=6, height=6)

#######################################################
# Plot ancestral states - DEC
#######################################################
analysis_titletxt ="BioGeoBEARS DEC on 3taxa M3_areasallowed"

# Setup
results_object = resDEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

#######################################################
# Plot ancestral states - DECJ
#######################################################
analysis_titletxt ="BioGeoBEARS DEC+J on 3taxa M3_areasallowed"

# Setup
results_object = resDECj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it











#######################################################
# Time-stratified Biogeographic Stochastic Mapping (BSM)
#######################################################
model_name = "DEC_M3_timestrat"
res = resDEC

pdffn = paste0("Example_", model_name, "_v1.pdf")
pdf(pdffn, height=6, width=6)

analysis_titletxt = paste0(model_name, " on Example")

# Setup
results_object = res
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#######################################################
# Stochastic mapping on DEC M3b stratified with islands coming up
#######################################################
clado_events_tables = NULL
ana_events_tables = NULL
lnum = 0

#######################################################
# Get the inputs for Biogeographical Stochastic Mapping
# Note: this can be slow for large state spaces and trees, since 
# the independent likelihoods for each branch are being pre-calculated
# E.g., for 10 areas, this requires calculation of a 1024x1024 matrix
# for each branch.  On a tree with ~800 tips and thus ~1600 branches, this was about 1.6 gigs
# for storage of "BSM_inputs_file.Rdata".
# Update: 2015-09-23 -- now, if you used multicore functionality for the ML analysis,
# the same settings will be used for get_inputs_for_stochastic_mapping().
#######################################################
BSM_inputs_fn = "BSM_inputs_file.Rdata"
runInputsSlow = TRUE
if (runInputsSlow)
	{
	# debug:
	# cluster_already_open=FALSE; rootedge=FALSE; statenum_bottom_root_branch_1based=NULL; printlevel=1; min_branchlength=0.000001
	stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=res)
	save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
	} else {
	# Loads to "stochastic_mapping_inputs_list"
	load(BSM_inputs_fn)
	} # END if (runInputsSlow)

# Check inputs (doesn't work the same on unconstr)
names(stochastic_mapping_inputs_list)
stochastic_mapping_inputs_list$phy2
stochastic_mapping_inputs_list$COO_weights_columnar
stochastic_mapping_inputs_list$unconstr
set.seed(seed=as.numeric(Sys.time()))

runBSMslow = TRUE
if (runBSMslow == TRUE)
    {
    # Saves to: RES_clado_events_tables.Rdata
    # Saves to: RES_ana_events_tables.Rdata
		# Bug check:
		# stochastic_mapping_inputs_list=stochastic_mapping_inputs_list; maxnum_maps_to_try=100; nummaps_goal=50; maxtries_per_branch=40000; save_after_every_try=TRUE; savedir=getwd(); seedval=12345; wait_before_save=0.01; master_nodenum_toPrint=0
		
    BSM_output = runBSM(res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=1, nummaps_goal=50, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=12345, wait_before_save=0.01, master_nodenum_toPrint=0)

    RES_clado_events_tables = BSM_output$RES_clado_events_tables
    RES_ana_events_tables = BSM_output$RES_ana_events_tables
    } else {
    # Load previously saved...

    # Loads to: RES_clado_events_tables
    load(file="RES_clado_events_tables.Rdata")
    # Loads to: RES_ana_events_tables
    load(file="RES_ana_events_tables.Rdata")
    BSM_output = NULL
    BSM_output$RES_clado_events_tables = RES_clado_events_tables
    BSM_output$RES_ana_events_tables = RES_ana_events_tables
    } # END if (runBSMslow == TRUE)

# Extract BSM output
clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables
head(clado_events_tables[[1]])
head(ana_events_tables[[1]])
length(clado_events_tables)
length(ana_events_tables)

include_null_range = TRUE
areanames = names(tipranges@df)
areas = areanames
max_range_size = 3

# Note: If you did something to change the states_list from the default given the number of areas, you would
# have to manually make that change here as well! (e.g., areas_allowed matrix, or manual reduction of the states_list)
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)

colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range=TRUE)

############################################
# Setup for painting a single stochastic map
############################################
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
stratified = TRUE
clado_events_table = clado_events_tables[[1]]
ana_events_table = ana_events_tables[[1]]

# cols_to_get = names(clado_events_table[,-ncol(clado_events_table)])
# colnums = match(cols_to_get, names(ana_events_table))
# ana_events_table_cols_to_add = ana_events_table[,colnums]
# anagenetic_events_txt_below_node = rep("none", nrow(ana_events_table_cols_to_add))
# ana_events_table_cols_to_add = cbind(ana_events_table_cols_to_add, anagenetic_events_txt_below_node)
# rows_to_get_TF = ana_events_table_cols_to_add$node <= length(tr$tip.label)
# master_table_cladogenetic_events = rbind(ana_events_table_cols_to_add[rows_to_get_TF,], clado_events_table)

############################################
# Open a PDF
############################################
pdffn = paste0(model_name, "_single_stochastic_map_n1.pdf")
pdf(file=pdffn, height=6, width=6)

# Convert the BSM into a modified res object
master_table_cladogenetic_events = clado_events_tables[[1]]
resmod = stochastic_map_states_into_res(res=res, master_table_cladogenetic_events=master_table_cladogenetic_events, stratified=stratified)

plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=FALSE, show.tip.label=TRUE)

# Paint on the branch states
paint_stochastic_map_branches(res=resmod, master_table_cladogenetic_events=master_table_cladogenetic_events, colors_list_for_states=colors_list_for_states, lwd=5, lty=par("lty"), root.edge=TRUE, stratified=stratified)

plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=TRUE, show.tip.label=TRUE)

############################################
# Close PDF
############################################
dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

#######################################################
# Plot all 50 stochastic maps to PDF
#######################################################
# Setup
include_null_range = include_null_range
areanames = areanames
areas = areanames
max_range_size = max_range_size
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)
colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range=TRUE)
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
stratified = stratified

# Loop through the maps and plot to PDF
pdffn = paste0(model_name, "_", length(clado_events_tables), "BSMs_v1.pdf")
pdf(file=pdffn, height=6, width=6)

nummaps_goal = 50
for (i in 1:nummaps_goal)
    {
    clado_events_table = clado_events_tables[[i]]
    analysis_titletxt = paste0(model_name, " - Stochastic Map #", i, "/", nummaps_goal)
    plot_BSM(results_object=res, clado_events_table=clado_events_table, stratified=stratified, analysis_titletxt=analysis_titletxt, addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, show.tip.label=TRUE, include_null_range=include_null_range)
    } # END for (i in 1:nummaps_goal)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

#######################################################
# Summarize stochastic map tables
#######################################################
length(clado_events_tables)
length(ana_events_tables)

head(clado_events_tables[[1]][,-20])
tail(clado_events_tables[[1]][,-20])

head(ana_events_tables[[1]])
tail(ana_events_tables[[1]])

areanames = names(tipranges@df)
actual_names = areanames
actual_names

# Get the dmat and times (if any)
dmat_times = get_dmat_times_from_res(res=res, numstates=NULL)
dmat_times

# Extract BSM output
clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables

# Simulate the source areas
BSMs_w_sourceAreas = simulate_source_areas_ana_clado(res, clado_events_tables, ana_events_tables, areanames)
clado_events_tables = BSMs_w_sourceAreas$clado_events_tables
ana_events_tables = BSMs_w_sourceAreas$ana_events_tables

# Count all anagenetic and cladogenetic events
counts_list = count_ana_clado_events(clado_events_tables, ana_events_tables, areanames, actual_names)

summary_counts_BSMs = counts_list$summary_counts_BSMs
print(conditional_format_table(summary_counts_BSMs))

# Histogram of event counts
hist_event_counts(counts_list, pdffn=paste0(model_name, "_histograms_of_event_counts.pdf"))

#######################################################
# Print counts to files
#######################################################
tmpnames = names(counts_list)
cat("\n\nWriting tables* of counts to tab-delimited text files:\n(* = Tables have dimension=2 (rows and columns). Cubes (dimension 3) and lists (dimension 1) will not be printed to text files.) \n\n")
for (i in 1:length(tmpnames))
    {
    cmdtxt = paste0("item = counts_list$", tmpnames[i])
    eval(parse(text=cmdtxt))

    # Skip cubes
    if (length(dim(item)) != 2)
        {
        next()
        }

    outfn = paste0(tmpnames[i], ".txt")
    if (length(item) == 0)
        {
        cat(outfn, " -- NOT written, *NO* events recorded of this type", sep="")
        cat("\n")
        } else {
        cat(outfn)
        cat("\n")
        write.table(conditional_format_table(item), file=outfn, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
        } # END if (length(item) == 0)
    } # END for (i in 1:length(tmpnames))
cat("...done.\n")

#######################################################
# Check that ML ancestral state/range probabilities and
# the mean of the BSMs approximately line up
#######################################################
library(MultinomialCI)    # For 95% CIs on BSM counts
check_ML_vs_BSM(res, clado_events_tables, model_name, tr=NULL, plot_each_node=FALSE, linreg_plot=TRUE, MultinomialCI=TRUE)


