#######################################################
# This is an introductory example script for the 
# R package "BioGeoBEARS" by Nick Matzke
#
# The package is designed for ML and Bayesian inference
# of 
# 
# (a) ancestral geographic ranges, and 
# 
# (b) perhaps more importantly, models for the 
#     evolution of geographic range across a phylogeny.
#
# The example below implements and compares:
# 
# (1) The standard 2-parameter DEC model implemented in 
#     the program LAGRANGE (Ree & Smith 2008); users will
#     notice that the ML parameter inference and log-
#     likelihoods are identical
#
# (2) A DEC+J model implemented in BioGeoBEARS, wherein
#     a third parameter, j, is added, representing the 
#     relative per-event weight of founder-event / jump
#     speciation events at cladogenesis events.  The 
#     higher j is, the more probability these events have,
#     and the less probability the standard LAGRANGE
#     cladogenesis events have.
#
# (3) Some standard model-testing (LRT and AIC) is 
#     implemented at the end so that users may compare models
#
# (4) The script does similar tests of a DIVA-like model (Ronquist 1997)
#     and a BAYAREA-like model (Landis, Matzke, Moore, & Huelsenbeck, 2013)
# 
#######################################################

# Load the package (do install.packages("BioGeoBEARS") to install it).
library(BioGeoBEARS)

#######################################################
# UNIVERSAL SETUP
#######################################################

# You will need to set your working directory to match your local system
# You can find the input files at:
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
extdata_dir
list.files(extdata_dir)

# Set your working directory for output files
# default here is your home direcotry
wd = np("~")
setwd(wd)

max_range_size = 4

# Results table
restable = NULL
teststable = NULL



#######################################################
# DEC AND DECJ ANALYSIS
#######################################################

#######################################################
# Run DEC
#######################################################

BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$speedup=TRUE		# seems to work OK
BioGeoBEARS_run_object$calc_ancprobs=TRUE		# get ancestral states from optim run


# Set up the stratified part
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"

BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$num_cores_to_use=4

BioGeoBEARS_run_object$force_sparse=FALSE
BioGeoBEARS_run_object

# Input geography text file
geogfn = np(paste(addslash(extdata_dir), "Psychotria_geog.data", sep=""))
BioGeoBEARS_run_object$geogfn = geogfn

tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=np(BioGeoBEARS_run_object$geogfn))
tipranges

# Input tree
trfn = np(paste(addslash(extdata_dir), "Psychotria_5.2.newick", sep=""))
BioGeoBEARS_run_object$trfn = trfn

tr = read.tree(BioGeoBEARS_run_object$trfn)
# tipnames1 = sort(tr$tip.label)
# length(tipnames1)
# tipnames2 = sort(rownames(tipranges@df))
# length(tipnames2)
# tipnames1 == tipnames2
# 
# tipnames1
# tipnames2



# Run to check inputs
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
# Divide the tree up by strata
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$master_table



# Set up DEC model
# (nothing to do; defaults)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "Psychotria_DEC_M0_unconstrained_v1.Rdata"
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
# Run DECj
#######################################################

BioGeoBEARS_run_object = define_BioGeoBEARS_run()

BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$num_cores_to_use=4
BioGeoBEARS_run_object$geogfn = geogfn

# Set up the stratified part
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"

BioGeoBEARS_run_object$force_sparse=FALSE	# sparse=FALSE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$speedup=TRUE		# seems to work OK
BioGeoBEARS_run_object$calc_ancprobs=TRUE		# get ancestral states from optim run


# Run to check inputs
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
# Divide the tree up by strata
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$master_table



# Set up DEC+J model
# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "Psychotria_DECJ_M0_unconstrained_v1.Rdata"
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














pdffn = "Psychotria_DEC_vs_DECj_M0_unconstrained_v1.pdf"
pdf(pdffn, width=8.5, height=11)


#######################################################
# Plot ancestral states - DEC
#######################################################
analysis_titletxt ="BioGeoBEARS DEC on Psychotria M0_unconstrained"

# Setup
results_object = resDEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)





#######################################################
# Plot ancestral states - DECJ
#######################################################
analysis_titletxt ="BioGeoBEARS DEC+J on Psychotria M0_unconstrained"

# Setup
results_object = resDECj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

junk='
results_object; analysis_titletxt; addl_params=list("j"); plotwhat="text"; label.offset=NULL; tipcex=0.45; statecex=0.7; splitcex=0.6; titlecex=0.8; plotsplits=TRUE; cornercoords_loc=scriptdir; include_null_range=TRUE; tr=tr; tipranges=tipranges
'

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)


dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)










#######################################################
# Stats
#######################################################
LnL_2 = as.numeric(resDEC$optim_result$fvalues)
LnL_1 = as.numeric(resDECj$optim_result$fvalues)
numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats


res1
res2

rbind(res1, res2)
tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res1, res2)
teststable = rbind(teststable, tmp_tests)






#######################################################
# DIVA AND DIVAJ ANALYSIS
#######################################################



#######################################################
# Run DIVA
#######################################################

BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object
BioGeoBEARS_run_object$speedup=TRUE		# seems to work OK
BioGeoBEARS_run_object$calc_ancprobs=TRUE		# get ancestral states from optim run


# Set up the stratified part
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"

BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$num_cores_to_use=4

BioGeoBEARS_run_object$force_sparse=FALSE
BioGeoBEARS_run_object

# Input geography text file
geogfn = np(paste(addslash(extdata_dir), "Psychotria_geog.data", sep=""))
BioGeoBEARS_run_object$geogfn = geogfn

tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=np(BioGeoBEARS_run_object$geogfn))
tipranges

# Input tree
trfn = np(paste(addslash(extdata_dir), "Psychotria_5.2.newick", sep=""))
BioGeoBEARS_run_object$trfn = trfn





# Run to check inputs
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
# Divide the tree up by strata
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$master_table


# Set up DIVA model
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "Psychotria_DIVA_M0_unconstrained_v1.Rdata"
if (runslow)
	{
	res = bears_optim_run(BioGeoBEARS_run_object)
	res	
	
	save(res, file=resfn)
	resDIVA = res
	} else {
	# Loads to "res"
	load(resfn)
	resDIVA = res
	}






#######################################################
# Run DIVAj
#######################################################

BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object

BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$num_cores_to_use=4
BioGeoBEARS_run_object$geogfn = geogfn

# Set up the stratified part
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"

BioGeoBEARS_run_object$force_sparse=FALSE	# sparse=FALSE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$speedup=TRUE		# seems to work OK
BioGeoBEARS_run_object$calc_ancprobs=TRUE		# get ancestral states from optim run


# Run to check inputs
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
# Divide the tree up by strata
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$master_table



# Set up DIVA model
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "Psychotria_DIVAJ_M0_unconstrained_v1.Rdata"
runslow = TRUE
if (runslow)
	{
	#sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")

	
	res = bears_optim_run(BioGeoBEARS_run_object)
	res	
	
	save(res, file=resfn)
	
	resDIVAj = res
	} else {
	# Loads to "res"
	load(resfn)
	resDIVAj = res
	}







pdffn = "Psychotria_DIVA_vs_DIVAj_M0_unconstrained_v1.pdf"
pdf(pdffn, width=8.5, height=11)


#######################################################
# Plot ancestral states - DIVA
#######################################################
analysis_titletxt ="BioGeoBEARS DIVA on Psychotria M0_unconstrained"

# Setup
results_object = resDIVA
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)








#######################################################
# Plot ancestral states - DIVAJ
#######################################################
analysis_titletxt ="BioGeoBEARS DIVA+J on Psychotria M0_unconstrained"

# Setup
results_object = resDIVAj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

junk='
results_object; analysis_titletxt; addl_params=list("j"); plotwhat="text"; label.offset=NULL; tipcex=0.45; statecex=0.7; splitcex=0.6; titlecex=0.8; plotsplits=TRUE; cornercoords_loc=scriptdir; include_null_range=TRUE; tr=tr; tipranges=tipranges
'

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)


dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)










#######################################################
# Stats
#######################################################
LnL_2 = as.numeric(resDIVA$optim_result$fvalues)
LnL_1 = as.numeric(resDIVAj$optim_result$fvalues)
numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

res1
res2

rbind(res1, res2)
conditional_format_table(stats)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res1, res2)
teststable = rbind(teststable, tmp_tests)








#######################################################
# BAYAREA AND BAYAREAJ ANALYSIS
#######################################################

#######################################################
# Run BAYAREA
#######################################################

BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object
BioGeoBEARS_run_object$speedup=TRUE		# seems to work OK
BioGeoBEARS_run_object$calc_ancprobs=TRUE		# get ancestral states from optim run


# Set up the stratified part
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"

BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$num_cores_to_use=4

BioGeoBEARS_run_object$force_sparse=FALSE
BioGeoBEARS_run_object

# Input geography text file
geogfn = np(paste(addslash(extdata_dir), "Psychotria_geog.data", sep=""))
BioGeoBEARS_run_object$geogfn = geogfn

tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=np(BioGeoBEARS_run_object$geogfn))
tipranges

# Input tree
trfn = np(paste(addslash(extdata_dir), "Psychotria_5.2.newick", sep=""))
BioGeoBEARS_run_object$trfn = trfn


# Run to check inputs
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
# Divide the tree up by strata
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$master_table



# Set up BAYAREA+J model
# Only sympatric/range-copying (y) events allowed
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Exact copying (smaller descendant always the same size as the larger descendant)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "Psychotria_BAYAREA_M0_unconstrained_v1.Rdata"
if (runslow)
	{
	res = bears_optim_run(BioGeoBEARS_run_object)
	res	
	
	save(res, file=resfn)
	resBAYAREA = res
	} else {
	# Loads to "res"
	load(resfn)
	resBAYAREA = res
	}






#######################################################
# Run BAYAREAj
#######################################################

BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object

BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$num_cores_to_use=4
BioGeoBEARS_run_object$geogfn = geogfn

# Set up the stratified part
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"

BioGeoBEARS_run_object$force_sparse=FALSE	# sparse=FALSE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$speedup=TRUE		# seems to work OK
BioGeoBEARS_run_object$calc_ancprobs=TRUE		# get ancestral states from optim run

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01


# Run to check inputs
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
# Divide the tree up by strata
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$master_table



# Set up BAYAREA+J model
# Only sympatric/range-copying (y) events allowed
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Exact copying (smaller descendant always the same size as the larger descendant)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "Psychotria_BAYAREAJ_M0_unconstrained_v1.Rdata"
runslow = TRUE
if (runslow)
	{
	res = bears_optim_run(BioGeoBEARS_run_object)
	res	
	
	save(res, file=resfn)
	
	resBAYAREAj = res
	} else {
	# Loads to "res"
	load(resfn)
	resBAYAREAj = res
	}









pdffn = "Psychotria_BAYAREA_vs_BAYAREAj_M0_unconstrained_v1.pdf"
pdf(pdffn, width=8.5, height=11)


#######################################################
# Plot ancestral states - BAYAREA
#######################################################
analysis_titletxt ="BioGeoBEARS BAYAREA on Psychotria M0_unconstrained"

# Setup
results_object = resBAYAREA
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)









#######################################################
# Plot ancestral states - BAYAREAJ
#######################################################
analysis_titletxt ="BioGeoBEARS BAYAREA+J on Psychotria M0_unconstrained"

# Setup
results_object = resBAYAREAj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

junk='
results_object; analysis_titletxt; addl_params=list("j"); plotwhat="text"; label.offset=NULL; tipcex=0.45; statecex=0.7; splitcex=0.6; titlecex=0.8; plotsplits=TRUE; cornercoords_loc=scriptdir; include_null_range=TRUE; tr=tr; tipranges=tipranges
'

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)


dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)










#######################################################
# Stats
#######################################################
LnL_2 = as.numeric(resBAYAREA$optim_result$fvalues)
LnL_1 = as.numeric(resBAYAREAj$optim_result$fvalues)
numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

res1
res2

rbind(res1, res2)
conditional_format_table(stats)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res1, res2)
teststable = rbind(teststable, tmp_tests)


# RESULTS: DEC, DECJ, DIVA, DIVAJ, BAYAREA, BAYAREAJ

restable
teststable



#######################################################
# Your results should be:
#######################################################

# > restable
# 
# LnL numparams            d            e         j
# -34.5         2 3.506337e-02 2.847676e-02 0.0000000
# -20.9         3 1.000000e-15 1.000000e-15 0.1142643
# -33.1         2 4.474911e-02 1.000000e-15 0.0000000
# -21.1         3 1.000000e-15 1.000000e-15 0.1157198
# -40.3         2 1.875408e-02 3.057997e-01 0.0000000
# -21.6         3 1.000000e-15 2.000001e-09 0.1080933


#######################################################
# The p-value of the LRT (Likelihood Ratio Test) tells you whether or not you can reject the
# null hypothesis that DEC and DEC+J confer equal likelihoods on the data
#
# AIC and AIC model weights are also shown, giving a sense of the relative probability of the two models.
#
# (One could easily do model weights between all 6 models, but this is not done here.)
#######################################################

# > teststable
# 
# LnLalt LnLnull DFalt DFnull DF Dstatistic    pval        test       tail  AIC1  AIC2 AICwt1  AICwt2 AICweight_ratio_model1 AICweight_ratio_model2
# 
# -20.95  -34.54     3      2  1      27.19 1.8e-07 chi-squared one-tailed  47.9 73.08   1.00 3.4e-06                 294895                3.4e-06
# 
# -21.09  -33.15     3      2  1      24.13 9.0e-07 chi-squared one-tailed 48.17  70.3   1.00 1.6e-05                  63797                1.6e-05
# 
# -21.55  -40.33     3      2  1      37.56 8.8e-10 chi-squared one-tailed 49.11 84.67   1.00 1.9e-08               5.28e+07                1.9e-08
# 


