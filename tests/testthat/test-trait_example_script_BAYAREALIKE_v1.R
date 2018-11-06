

test_that(desc="Check that cladoRcpp version number is >= 0.15", code={

version_number = packageVersion("cladoRcpp")
TF = version_number >= 0.15

if (TF == FALSE)
	{
	txt = 'STOP ERROR inside test_that(desc="Check that cladoRcpp version number is >= 0.15"): the BioGeoBEARS "testthat" tests, located in BioGeoBEARS/tests, require that cladoRcpp have version 0.15 or higher to work. To get the new version, try "devtools::install_github(repo="nmatzke/cladoRcpp", quick=TRUE, dependencies=FALSE, build_vignettes=FALSE, keep_source=TRUE, local=FALSE, force=TRUE)".'
	
	cat("\n\n")
	cat(txt)
	cat("\n\n")
	}

expect_equal(object=TF, expected=TRUE)

}) # END test_that









test_that(desc="Check that GenSA is installed", code={

TF = is.element("GenSA", installed.packages()[,1])

if (TF == FALSE)
	{
	txt = 'STOP ERROR inside test_that(desc="Check that GenSA is installed"): the BioGeoBEARS "testthat" tests of trait-based models, located in BioGeoBEARS/tests, require that GenSA be installed. To get it, try \n\ninstall.packages("GenSA")\n.\n'
	
	cat("\n\n")
	cat(txt)
	cat("\n\n")
	}

expect_equal(object=TF, expected=TRUE)

}) # END test_that





test_that(desc="Check trait-dependent dispersal inference, for base models derived from 'BAYAREALIKE'", code={

# Skip the slow tests in online checks
testthat::skip_on_cran()
testthat::skip_on_travis()



setup='
library(parallel)
library(cladoRcpp)
library(BioGeoBEARS)
library(GenSA)
sourceall("/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)
'


# library(parallel)
# library(cladoRcpp)
# library(BioGeoBEARS)
# library(GenSA)

library(cladoRcpp)
#source("/drives/Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp.R")
library(BioGeoBEARS)
library(parallel) # for detectCores
library(GenSA)


###############
# Specify the input files
###############
set.seed(54321)

#wd = "/drives/GDrive/__GDrive_projects/2016-12-07_Kristina_Klaus_Podocarpaceae/__doc3/Supplementary_Material/example/BAYAREALIKE_wTraits/"
#setwd(wd)

# Example directory with files for traits-based analysis
extdata_dir = np(system.file("extdata/examples/trait_examples/BAYAREALIKE_wTraits", package="BioGeoBEARS"))
extdata_dir
list.files(extdata_dir)


trfn = np(paste(addslash(extdata_dir), "simtr_observed.newick", sep=""))
tr = read.tree(trfn)

geog_1area_fn = np(paste(addslash(extdata_dir), "geog_1area.data", sep=""))
geogfn = np(paste(addslash(extdata_dir), "geog_sim_observed.txt", sep=""))
traitsfn = np(paste(addslash(extdata_dir), "traits_sim_observed.txt", sep=""))

#######################################################
# Inference
#######################################################
max_range_size = 4


#######################################################
# Traits-only model -- 1 rate
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$print_optim = TRUE
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
BioGeoBEARS_run_object$max_range_size = 1
BioGeoBEARS_run_object$num_cores_to_use=23
BioGeoBEARS_run_object$use_optimx="GenSA"
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = geog_1area_fn
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1000000
BioGeoBEARS_run_object$force_sparse = FALSE  # works with kexpmv, but compare to dense,
# time-stratify to break up long branches if you see major differences in lnL

# Set up BAYAREALIKE model, but set all rates to 0 (data are 1 invariant area)
# (nothing to do; defaults)

# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

# This contains the model object
BioGeoBEARS_run_object$BioGeoBEARS_model_object

# This table contains the parameters of the model 
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = 0.0

# Set up BAYAREALIKE model
# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# No jump dispersal/founder-event speciation
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999


tr = read.tree(BioGeoBEARS_run_object$trfn)
#plot(tr); axisPhylo()

geog_values = getranges_from_LagrangePHYLIP(traitsfn)

trait_fn = traitsfn
trait_values = getranges_from_LagrangePHYLIP(lgdata_fn=trait_fn)
trait_values

# Add the traits data and model
BioGeoBEARS_run_object = add_trait_to_BioGeoBEARS_run_object(BioGeoBEARS_run_object, trait_fn=trait_fn)


# Look at the params table
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

#######################################################
# Manual modifications of trait-based model
#######################################################
# Edit t12 and t21 rates

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "init"] = 0.001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "est"] = 0.001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "max"] = 1

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "type"] = "t12"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "init"] = 0.001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "est"] = 0.001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "max"] = 1

# No multipliers on geog (set m1 and m2 to 1)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "init"] = 1.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "init"] = 1.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "est"] = 1.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "est"] = 1.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "desc"] = "trait-based dispersal rate multipliers m1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "desc"] = "trait-based dispersal rate multipliers m2"

# Run this to check inputs. Read the error messages if you get them!
BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# No printing during testthat checks
BioGeoBEARS_run_object$print_optim = FALSE
BioGeoBEARS_run_object$printlevel = 0

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "sim_traitsOnly_1rate_v1.Rdata"
if (runslow)
		{
		res = bears_optim_run(BioGeoBEARS_run_object)
		res    

		save(res, file=resfn)
		resTrait_1rate = res
		} else {
		# Loads to "res"
		load(resfn)
		resTrait_1rate = res
		} # END if (runslow)



#######################################################
# Traits-only model -- 2 rates
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$print_optim = TRUE
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
BioGeoBEARS_run_object$max_range_size = 1
BioGeoBEARS_run_object$num_cores_to_use=23
BioGeoBEARS_run_object$use_optimx="GenSA"
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = geog_1area_fn
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1000000
BioGeoBEARS_run_object$force_sparse = FALSE  # works with kexpmv, but compare to dense,
# time-stratify to break up long branches if you see major differences in lnL

# Set up BAYAREALIKE model, but set all rates to 0 (data are 1 invariant area)
# (nothing to do; defaults)

# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

# This contains the model object
BioGeoBEARS_run_object$BioGeoBEARS_model_object

# This table contains the parameters of the model 
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = 0.0

# Set up BAYAREALIKE model
# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# No jump dispersal/founder-event speciation
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999


tr = read.tree(BioGeoBEARS_run_object$trfn)
#plot(tr); axisPhylo()

trait_fn = traitsfn
trait_values = getranges_from_LagrangePHYLIP(lgdata_fn=trait_fn)
trait_values

# Add the traits data and model
BioGeoBEARS_run_object = add_trait_to_BioGeoBEARS_run_object(BioGeoBEARS_run_object, trait_fn=trait_fn)


# Look at the params table
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

#######################################################
# Manual modifications of trait-based model
#######################################################
# Edit t12 and t21 rates

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "init"] = 0.001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "est"] = 0.001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "max"] = 1

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "init"] = 0.001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "est"] = 0.001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "max"] = 1

# No multipliers on geog (set m1 and m2 to 1)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "init"] = 1.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "init"] = 1.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "est"] = 1.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "est"] = 1.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "desc"] = "trait-based dispersal rate multipliers m1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "desc"] = "trait-based dispersal rate multipliers m2"

# Run this to check inputs. Read the error messages if you get them!
BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)


# No printing during testthat checks
BioGeoBEARS_run_object$print_optim = FALSE
BioGeoBEARS_run_object$printlevel = 0

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "sim_traitsOnly_2rates_v1.Rdata"
if (runslow)
		{
		res = bears_optim_run(BioGeoBEARS_run_object)
		res    

		save(res, file=resfn)
		resTrait_2rates = res
		} else {
		# Loads to "res"
		load(resfn)
		resTrait_2rates = res
		} # END if (runslow)



#######################################################
# Run BAYAREALIKE (on geography only)
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$print_optim = TRUE
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$num_cores_to_use=23
BioGeoBEARS_run_object$use_optimx="GenSA"
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1000000
BioGeoBEARS_run_object$force_sparse = FALSE  # works with kexpmv, but compare to dense,
# time-stratify to break up long branches if you see major differences in lnL

#tr = read.tree(BioGeoBEARS_run_object$trfn)

# Set up BAYAREALIKE model
# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# No jump dispersal/founder-event speciation
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999


BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)


# No printing during testthat checks
BioGeoBEARS_run_object$print_optim = FALSE
BioGeoBEARS_run_object$printlevel = 0

runslow = TRUE
resfn = "BAYAREALIKE_inf.Rdata"
if (runslow)
	{
	res = bears_optim_run(BioGeoBEARS_run_object)
	res    

	save(res, file=resfn)
	resBAYAREALIKE = res
	} else {
	# Loads to "res"
	load(resfn)
	resBAYAREALIKE = res
	}


#######################################################
# Run BAYAREALIKE+J (on geography only)
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$print_optim = TRUE
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$num_cores_to_use=23
BioGeoBEARS_run_object$use_optimx="GenSA"
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1000000
BioGeoBEARS_run_object$force_sparse = FALSE  # works with kexpmv, but compare to dense,
# time-stratify to break up long branches if you see major differences in lnL

#tr = read.tree(BioGeoBEARS_run_object$trfn)

# Set up BAYAREALIKE+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resBAYAREALIKE$outputs@params_table["d","est"]
estart = resBAYAREALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in BAYAREALIKE+J) or 2 (as in DIVALIKE+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows 
# machines. I can't replicate this on my Mac machines, but it is almost certainly
# just some precision under-run issue, when optim/optimx tries some parameter value 
# just below zero.  The "min" and "max" options on each parameter are supposed to
# prevent this, but apparently optim/optimx sometimes go slightly beyond 
# these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" 
# slightly for each parameter:
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999




# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart


BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)


# No printing during testthat checks
BioGeoBEARS_run_object$print_optim = FALSE
BioGeoBEARS_run_object$printlevel = 0

print("Printing warnings: 'warnings()':")
print(warnings())

runslow = TRUE
resfn = "BAYAREALIKEj_inf.Rdata"
if (runslow)
	{
	res = bears_optim_run(BioGeoBEARS_run_object)
	res    

	save(res, file=resfn)
	resBAYAREALIKEj = res
	} else {
	# Loads to "res"
	load(resfn)
	resBAYAREALIKEj = res
	}






#######################################################
# Run BAYAREALIKE + t12 + t21 + m2, starting from BAYAREALIKE-geog and 2rates
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$print_optim = TRUE
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$num_cores_to_use=23
BioGeoBEARS_run_object$use_optimx="GenSA"
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$trfn = trfn
#BioGeoBEARS_run_object$timesfn = "times_v2.txt"
#BioGeoBEARS_run_object$distsfn = "geological_distances_v3_div100_stay_same.txt"
#BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1000000
BioGeoBEARS_run_object$force_sparse = FALSE  # works with kexpmv, but compare to dense,
# time-stratify to break up long branches if you see major differences in lnL

tr = read.tree(BioGeoBEARS_run_object$trfn)
#plot(tr); axisPhylo()

# Set up BAYAREALIKE model
# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# No jump dispersal/founder-event speciation
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999



trait_fn = traitsfn
trait_values = getranges_from_LagrangePHYLIP(lgdata_fn=trait_fn)
trait_values

# Add the traits data and model
BioGeoBEARS_run_object = add_trait_to_BioGeoBEARS_run_object(BioGeoBEARS_run_object, trait_fn=trait_fn)


# Look at the params table
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table


# Starting values from ML results of simpler run
t12_start = resTrait_2rates$outputs@params_table["t12","est"]
t21_start = resTrait_2rates$outputs@params_table["t21","est"]
m2_start = 1
dstart = resBAYAREALIKE$outputs@params_table["d","est"]
estart = max(c(resBAYAREALIKE$outputs@params_table["e","est"], 1.1*BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"]))
jstart = 0.0001

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "est"] = estart


#######################################################
# Manual modifications of trait-based model
#######################################################
# Edit t12 and t21 rates

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "init"] = t12_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "est"] = t12_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "max"] = round(max(t12_start, t21_start)*5, 3)


BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "init"] = t21_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "est"] = t21_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "max"] = round(max(t12_start, t21_start)*5, 3)


# Set 0/1 multipliers on dispersal rate
# For flightlessness (m2), max multiplier is 1, and
# fix to a small value, or estimate
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "init"] = 1
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "est"] = 1
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "min"] = 0.01
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "max"] = 1

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "init"] = m2_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "est"] = m2_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "min"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "max"] = 10

BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)


# No printing during testthat checks
BioGeoBEARS_run_object$print_optim = FALSE
BioGeoBEARS_run_object$printlevel = 0


runslow = TRUE
resfn = "BAYAREALIKE+t12+t21+m2_inf.Rdata"
if (runslow)
	{
	res = bears_optim_run(BioGeoBEARS_run_object, skip_optim=FALSE)
	
	save(res, file=resfn)
	resBAYAREALIKE_t12_t21_m2 = res
	} else {
	# Loads to "res"
	load(resfn)
	resBAYAREALIKE_t12_t21_m2 = res
	}




#######################################################
# Run BAYAREALIKEj + t12 + t21 + m2, starting from BAYAREALIKEj-geog and 2rates
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$print_optim = TRUE
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$num_cores_to_use=23
BioGeoBEARS_run_object$use_optimx="GenSA"
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$trfn = trfn
#BioGeoBEARS_run_object$timesfn = "times_v2.txt"
#BioGeoBEARS_run_object$distsfn = "geological_distances_v3_div100_stay_same.txt"
#BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1000000
BioGeoBEARS_run_object$force_sparse = FALSE  # works with kexpmv, but compare to dense,
# time-stratify to break up long branches if you see major differences in lnL

tr = read.tree(BioGeoBEARS_run_object$trfn)
#plot(tr); axisPhylo()


# Set up BAYAREALIKE+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resBAYAREALIKEj$outputs@params_table["d","est"]
estart = resBAYAREALIKEj$outputs@params_table["e","est"]
jstart = resBAYAREALIKEj$outputs@params_table["j","est"]

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in BAYAREALIKE+J) or 2 (as in DIVALIKE+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows 
# machines. I can't replicate this on my Mac machines, but it is almost certainly
# just some precision under-run issue, when optim/optimx tries some parameter value 
# just below zero.  The "min" and "max" options on each parameter are supposed to
# prevent this, but apparently optim/optimx sometimes go slightly beyond 
# these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" 
# slightly for each parameter:
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999


trait_fn = traitsfn
trait_values = getranges_from_LagrangePHYLIP(lgdata_fn=trait_fn)
trait_values

# Add the traits data and model
BioGeoBEARS_run_object = add_trait_to_BioGeoBEARS_run_object(BioGeoBEARS_run_object, trait_fn=trait_fn)



# Starting values from ML results of simpler run
t12_start = resTrait_2rates$outputs@params_table["t12","est"]
t21_start = resTrait_2rates$outputs@params_table["t21","est"]
m2_start = 1
dstart = resBAYAREALIKEj$outputs@params_table["d","est"]
estart = max(c(resBAYAREALIKEj$outputs@params_table["e","est"], 1.1*BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"]))
jstart = resBAYAREALIKEj$outputs@params_table["j","est"]


# Set up BAYAREALIKE+J model
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "est"] = estart

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Crash fix
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.0001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 1e-13
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 1e-13

#######################################################
# Manual modifications of trait-based model
#######################################################
# Edit t12 and t21 rates
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "init"] = t12_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "est"] = t12_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "max"] = round(max(t12_start, t21_start)*5, 3)


BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "init"] = t21_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "est"] = t21_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "max"] = round(max(t12_start, t21_start)*5, 3)


# Set 0/1 multipliers on dispersal rate
# For flightlessness (m2), max multiplier is 1, and
# fix to a small value, or estimate
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "init"] = 1
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "est"] = 1
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "min"] = 0.01
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "max"] = 1

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "init"] = m2_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "est"] = m2_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "min"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "max"] = 10

BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)


# No printing during testthat checks
BioGeoBEARS_run_object$print_optim = FALSE
BioGeoBEARS_run_object$printlevel = 0



resfn = "BAYAREALIKEJ+t12+t21+m2_inf.Rdata"
runslow = TRUE
if (runslow)
	{
	res = bears_optim_run(BioGeoBEARS_run_object, skip_optim=FALSE)

	save(res, file=resfn)
	resBAYAREALIKEj_t12_t21_m2 = res
	} else {
	# Loads to "res"
	load(resfn)
	resBAYAREALIKEj_t12_t21_m2 = res
	}







#######################################################
# Run BAYAREALIKEj + t12 + t21 + m2, starting from BAYAREALIKE + t12 + t21 + m2
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$print_optim = TRUE
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$num_cores_to_use=23
BioGeoBEARS_run_object$use_optimx="GenSA"
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$trfn = trfn
#BioGeoBEARS_run_object$timesfn = "times_v2.txt"
#BioGeoBEARS_run_object$distsfn = "geological_distances_v3_div100_stay_same.txt"
#BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1000000
BioGeoBEARS_run_object$force_sparse = FALSE  # works with kexpmv, but compare to dense,
# time-stratify to break up long branches if you see major differences in lnL

tr = read.tree(BioGeoBEARS_run_object$trfn)
#plot(tr); axisPhylo()



# Set up BAYAREALIKE+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resBAYAREALIKE$outputs@params_table["d","est"]
estart = resBAYAREALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in BAYAREALIKE+J) or 2 (as in DIVALIKE+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows 
# machines. I can't replicate this on my Mac machines, but it is almost certainly
# just some precision under-run issue, when optim/optimx tries some parameter value 
# just below zero.  The "min" and "max" options on each parameter are supposed to
# prevent this, but apparently optim/optimx sometimes go slightly beyond 
# these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" 
# slightly for each parameter:
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999



trait_fn = traitsfn
trait_values = getranges_from_LagrangePHYLIP(lgdata_fn=trait_fn)
trait_values

# Add the traits data and model
BioGeoBEARS_run_object = add_trait_to_BioGeoBEARS_run_object(BioGeoBEARS_run_object, trait_fn=trait_fn)



# Starting values from ML results of simpler run
t12_start = resBAYAREALIKE_t12_t21_m2$outputs@params_table["t12","est"]
t21_start = resBAYAREALIKE_t12_t21_m2$outputs@params_table["t21","est"]
m2_start = resBAYAREALIKE_t12_t21_m2$outputs@params_table["m2","est"]
dstart = resBAYAREALIKE_t12_t21_m2$outputs@params_table["d","est"]
estart = max(c(resBAYAREALIKE_t12_t21_m2$outputs@params_table["e","est"], 1.1*BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"]))
jstart = 0.0001


# Set up BAYAREALIKE+J model
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "est"] = estart

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Crash fix
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.0001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 1e-13
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 1e-13

#######################################################
# Manual modifications of trait-based model
#######################################################
# Edit t12 and t21 rates
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "init"] = t12_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "est"] = t12_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "max"] = round(max(t12_start, t21_start)*5, 3)

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "init"] = t21_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "est"] = t21_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "max"] = round(max(t12_start, t21_start)*5, 3)


# Set 0/1 multipliers on dispersal rate
# For flightlessness (m2), max multiplier is 1, and
# fix to a small value, or estimate
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "init"] = 1
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "est"] = 1
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "min"] = 0.01
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "max"] = 1

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "init"] = m2_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "est"] = m2_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "min"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "max"] = 10



BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)


# No printing during testthat checks
BioGeoBEARS_run_object$print_optim = FALSE
BioGeoBEARS_run_object$printlevel = 0

resfn = "BAYAREALIKEJ+t12+t21+m2_rep2_inf.Rdata"
runslow = TRUE
if (runslow)
	{
	res = bears_optim_run(BioGeoBEARS_run_object)
	res    

	save(res, file=resfn)

	resBAYAREALIKEj_t12_t21_m2_rep2 = res
	} else {
	# Loads to "res"
	load(resfn)
	resBAYAREALIKEj_t12_t21_m2_rep2 = res
	}








param_names = c("lnL", "d", "e", "j", "t12", "t21", "m1", "m2")

Trait_1rate_results = c(
resTrait_1rate$total_loglikelihood,
resTrait_1rate$output@params_table["d", "est"], 
resTrait_1rate$output@params_table["e", "est"], 
resTrait_1rate$output@params_table["j", "est"], 
resTrait_1rate$output@params_table["t12", "est"], 
resTrait_1rate$output@params_table["t21", "est"], 
resTrait_1rate$output@params_table["m1", "est"], 
resTrait_1rate$output@params_table["m2", "est"]
)
names(Trait_1rate_results) = paste("Trait_1rate_", param_names, sep="")

Trait_2rates_results = c(
resTrait_2rates$total_loglikelihood,
resTrait_2rates$output@params_table["d", "est"], 
resTrait_2rates$output@params_table["e", "est"], 
resTrait_2rates$output@params_table["j", "est"], 
resTrait_2rates$output@params_table["t12", "est"], 
resTrait_2rates$output@params_table["t21", "est"], 
resTrait_2rates$output@params_table["m1", "est"], 
resTrait_2rates$output@params_table["m2", "est"]
)
names(Trait_2rates_results) = paste("Trait_2rates_", param_names, sep="")

BAYAREALIKE_results = c(
resBAYAREALIKE$total_loglikelihood,
resBAYAREALIKE$output@params_table["d", "est"], 
resBAYAREALIKE$output@params_table["e", "est"], 
resBAYAREALIKE$output@params_table["j", "est"], 
resBAYAREALIKE$output@params_table["t12", "est"], 
resBAYAREALIKE$output@params_table["t21", "est"], 
resBAYAREALIKE$output@params_table["m1", "est"], 
resBAYAREALIKE$output@params_table["m2", "est"]
)
names(BAYAREALIKE_results) = paste("BAYAREALIKE_", param_names, sep="")


BAYAREALIKEj_results = c(
resBAYAREALIKEj$total_loglikelihood,
resBAYAREALIKEj$output@params_table["d", "est"], 
resBAYAREALIKEj$output@params_table["e", "est"], 
resBAYAREALIKEj$output@params_table["j", "est"], 
resBAYAREALIKEj$output@params_table["t12", "est"], 
resBAYAREALIKEj$output@params_table["t21", "est"], 
resBAYAREALIKEj$output@params_table["m1", "est"], 
resBAYAREALIKEj$output@params_table["m2", "est"]
)
names(BAYAREALIKEj_results) = paste("BAYAREALIKEj_", param_names, sep="")

BAYAREALIKE_t12_t21_m2_results = c(
resBAYAREALIKE_t12_t21_m2$total_loglikelihood,
resBAYAREALIKE_t12_t21_m2$output@params_table["d", "est"], 
resBAYAREALIKE_t12_t21_m2$output@params_table["e", "est"], 
resBAYAREALIKE_t12_t21_m2$output@params_table["j", "est"], 
resBAYAREALIKE_t12_t21_m2$output@params_table["t12", "est"], 
resBAYAREALIKE_t12_t21_m2$output@params_table["t21", "est"], 
resBAYAREALIKE_t12_t21_m2$output@params_table["m1", "est"], 
resBAYAREALIKE_t12_t21_m2$output@params_table["m2", "est"]
)
names(BAYAREALIKE_t12_t21_m2_results) = paste("BAYAREALIKE_t12_t21_m2_", param_names, sep="")

BAYAREALIKEj_t12_t21_m2_results = c(
resBAYAREALIKEj_t12_t21_m2$total_loglikelihood,
resBAYAREALIKEj_t12_t21_m2$output@params_table["d", "est"], 
resBAYAREALIKEj_t12_t21_m2$output@params_table["e", "est"], 
resBAYAREALIKEj_t12_t21_m2$output@params_table["j", "est"], 
resBAYAREALIKEj_t12_t21_m2$output@params_table["t12", "est"], 
resBAYAREALIKEj_t12_t21_m2$output@params_table["t21", "est"], 
resBAYAREALIKEj_t12_t21_m2$output@params_table["m1", "est"], 
resBAYAREALIKEj_t12_t21_m2$output@params_table["m2", "est"]
)
names(BAYAREALIKEj_t12_t21_m2_results) = paste("BAYAREALIKEj_t12_t21_m2_", param_names, sep="")

BAYAREALIKEj_t12_t21_m2_rep2_results = c(
resBAYAREALIKEj_t12_t21_m2_rep2$total_loglikelihood,
resBAYAREALIKEj_t12_t21_m2_rep2$output@params_table["d", "est"], 
resBAYAREALIKEj_t12_t21_m2_rep2$output@params_table["e", "est"], 
resBAYAREALIKEj_t12_t21_m2_rep2$output@params_table["j", "est"], 
resBAYAREALIKEj_t12_t21_m2_rep2$output@params_table["t12", "est"], 
resBAYAREALIKEj_t12_t21_m2_rep2$output@params_table["t21", "est"], 
resBAYAREALIKEj_t12_t21_m2_rep2$output@params_table["m1", "est"], 
resBAYAREALIKEj_t12_t21_m2_rep2$output@params_table["m2", "est"]
)
names(BAYAREALIKEj_t12_t21_m2_rep2_results) = paste("BAYAREALIKEj_t12_t21_m2_rep2_", param_names, sep="")


tmp_results = c(Trait_1rate_results, Trait_2rates_results, BAYAREALIKE_results, BAYAREALIKEj_results, BAYAREALIKE_t12_t21_m2_results, BAYAREALIKEj_t12_t21_m2_results, BAYAREALIKEj_t12_t21_m2_rep2_results)
tmp_results_mat = as.matrix(tmp_results, nrow=7, byrow=TRUE)
tmp_results_mat = as.data.frame(tmp_results_mat, stringsAsFactors=FALSE)

outfn = slashslash("params_inferred.txt")
write.table(x=tmp_results, file=outfn, append=FALSE, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
#moref(outfn)



# Rdata_fn = "BAYAREALIKE+t12+t21+m2_inf.Rdata"
# rerun_optim_table_BAYAREALIKE = rerun_optimization_w_HiLow(res=NULL, Rdata_fn=Rdata_fn, runslow=TRUE)
# 
# Rdata_fn = "BAYAREALIKEJ+t12+t21+m2_inf.Rdata"
# rerun_optim_table_BAYAREALIKEj = rerun_optimization_w_HiLow(res=NULL, Rdata_fn=Rdata_fn, runslow=TRUE)
# 
# Rdata_fn = "BAYAREALIKEJ+t12+t21+m2_rep2_inf.Rdata"
# rerun_optim_table_BAYAREALIKEj_rep2 = rerun_optimization_w_HiLow(res=NULL, Rdata_fn=Rdata_fn, runslow=TRUE)


#####################################
# Archived answers
#####################################
trait_1rate_lnL = -22.2601116898632
trait_1rate_t12 = 0.0119494911453998
trait_1rate_t21 = 0.0119494911453998

trait_2rate_lnL = -22.2516171160283
trait_2rate_t12 = 0.0124468030348559
trait_2rate_t21 = 0.0114012761732464

BAYAREALIKE_lnL = -156.081201427971
BAYAREALIKE_d = 0.0151784236784355
BAYAREALIKE_e = 0.0395752033623117
BAYAREALIKE_j = 0.000000e+00

BAYAREALIKEj_lnL = -106.063972271043
BAYAREALIKEj_d = 0.00536930031696263
BAYAREALIKEj_e = 9.99999999994061e-08
BAYAREALIKEj_j = 0.161074344802673

BAYAREALIKE_t12_t21_m2_lnL = -177.623522128533
BAYAREALIKE_t12_t21_m2_d = 0.0216429924891068
BAYAREALIKE_t12_t21_m2_e = 0.0419960196268601
BAYAREALIKE_t12_t21_m2_j = 0.000000e+00
BAYAREALIKE_t12_t21_m2_t12 = 0.0129710407591378
BAYAREALIKE_t12_t21_m2_t21 = 0.0105856903746298
BAYAREALIKE_t12_t21_m2_m2 = 0.52743168240104

BAYAREALIKEj_t12_t21_m2_rep1_lnL = -126.087073551735
BAYAREALIKEj_t12_t21_m2_rep1_d = 0.00788955637983855
BAYAREALIKEj_t12_t21_m2_rep1_e = 1e-13
BAYAREALIKEj_t12_t21_m2_rep1_j = 0.28452725365985
BAYAREALIKEj_t12_t21_m2_rep1_t12 = 0.0128706537541464
BAYAREALIKEj_t12_t21_m2_rep1_t21 = 0.0118195373271578
BAYAREALIKEj_t12_t21_m2_rep1_m2 = 0.32380017884888

BAYAREALIKEj_t12_t21_m2_rep2_lnL = -126.087152218089
BAYAREALIKEj_t12_t21_m2_rep2_d = 0.00787935823645791
BAYAREALIKEj_t12_t21_m2_rep2_e = 1e-13
BAYAREALIKEj_t12_t21_m2_rep2_j = 0.283726786050578
BAYAREALIKEj_t12_t21_m2_rep2_t12 = 0.0128542227522938
BAYAREALIKEj_t12_t21_m2_rep2_t21 = 0.0118162460770523
BAYAREALIKEj_t12_t21_m2_rep2_m2 = 0.325945808223546


txt = "Checking the maximum lnL inferences for 7 models in example BAYAREALIKE_wTraits trait-dependent dispersal script..."
cat("\n")
cat(txt)
archive_lnLs = c(trait_1rate_lnL, trait_2rate_lnL, BAYAREALIKE_lnL, BAYAREALIKEj_lnL, BAYAREALIKE_t12_t21_m2_lnL, BAYAREALIKEj_t12_t21_m2_rep1_lnL, BAYAREALIKEj_t12_t21_m2_rep2_lnL)
test_lnLs = c(resTrait_1rate$total_loglikelihood, resTrait_2rates$total_loglikelihood, resBAYAREALIKE$total_loglikelihood, resBAYAREALIKEj$total_loglikelihood, resBAYAREALIKE_t12_t21_m2$total_loglikelihood, resBAYAREALIKEj_t12_t21_m2$total_loglikelihood, resBAYAREALIKEj_t12_t21_m2_rep2$total_loglikelihood)

expect_equal(object=round(archive_lnLs, digits=4), expected=round(test_lnLs, digits=4))
cat("...PASSED  ")



txt = "Checking the ML inferences of 'd' for 5 models in example BAYAREALIKE_wTraits trait-dependent dispersal script..."
cat("\n")
cat(txt)
archive_vals = c(BAYAREALIKE_d, BAYAREALIKEj_d, BAYAREALIKE_t12_t21_m2_d, BAYAREALIKEj_t12_t21_m2_rep1_d, BAYAREALIKEj_t12_t21_m2_rep2_d)
test_vals = c(BAYAREALIKE_results[2], BAYAREALIKEj_results[2], BAYAREALIKE_t12_t21_m2_results[2], BAYAREALIKEj_t12_t21_m2_results[2], BAYAREALIKEj_t12_t21_m2_rep2_results[2])

expect_equal(object=round(archive_vals, digits=4), expected=unname(round(test_vals, digits=4)))
cat("...PASSED  ")


txt = "Checking the ML inferences of 'e' for 5 models in example BAYAREALIKE_wTraits trait-dependent dispersal script..."
cat("\n")
cat(txt)
archive_vals = c(BAYAREALIKE_e, BAYAREALIKEj_e, BAYAREALIKE_t12_t21_m2_e, BAYAREALIKEj_t12_t21_m2_rep1_e, BAYAREALIKEj_t12_t21_m2_rep2_e)
test_vals = c(BAYAREALIKE_results[3], BAYAREALIKEj_results[3], BAYAREALIKE_t12_t21_m2_results[3], BAYAREALIKEj_t12_t21_m2_results[3], BAYAREALIKEj_t12_t21_m2_rep2_results[3])

expect_equal(object=round(archive_vals, digits=4), expected=unname(round(test_vals, digits=4)))
cat("...PASSED  ")


txt = "Checking the ML inferences of 'j' for 5 models in example BAYAREALIKE_wTraits trait-dependent dispersal script..."
cat("\n")
cat(txt)
archive_vals = c(BAYAREALIKE_j, BAYAREALIKEj_j, BAYAREALIKE_t12_t21_m2_j, BAYAREALIKEj_t12_t21_m2_rep1_j, BAYAREALIKEj_t12_t21_m2_rep2_j)
test_vals = c(BAYAREALIKE_results[4], BAYAREALIKEj_results[4], BAYAREALIKE_t12_t21_m2_results[4], BAYAREALIKEj_t12_t21_m2_results[4], BAYAREALIKEj_t12_t21_m2_rep2_results[4])

expect_equal(object=round(archive_vals, digits=4), expected=unname(round(test_vals, digits=4)))
cat("...PASSED  ")





txt = "Checking the ML inferences of 't12' for 5 models in example BAYAREALIKE_wTraits trait-dependent dispersal script..."
cat("\n")
cat(txt)
archive_vals = c(trait_1rate_t12, trait_2rate_t12, BAYAREALIKE_t12_t21_m2_t12, BAYAREALIKEj_t12_t21_m2_rep1_t12, BAYAREALIKEj_t12_t21_m2_rep2_t12)
test_vals = c(Trait_1rate_results[5], Trait_2rates_results[5], BAYAREALIKE_t12_t21_m2_results[5], BAYAREALIKEj_t12_t21_m2_results[5], BAYAREALIKEj_t12_t21_m2_rep2_results[5])

expect_equal(object=round(archive_vals, digits=3), expected=unname(round(test_vals, digits=3)))
cat("...PASSED  ")


txt = "Checking the ML inferences of 't21' for 5 models in example BAYAREALIKE_wTraits trait-dependent dispersal script..."
cat("\n")
cat(txt)
archive_vals = c(trait_1rate_t21, trait_2rate_t21, BAYAREALIKE_t12_t21_m2_t21, BAYAREALIKEj_t12_t21_m2_rep1_t21, BAYAREALIKEj_t12_t21_m2_rep2_t21)
test_vals = c(Trait_1rate_results[6], Trait_2rates_results[6], BAYAREALIKE_t12_t21_m2_results[6], BAYAREALIKEj_t12_t21_m2_results[6], BAYAREALIKEj_t12_t21_m2_rep2_results[6])

expect_equal(object=round(archive_vals, digits=3), expected=unname(round(test_vals, digits=3)))
cat("...PASSED  ")



txt = "Checking the ML inferences of 'm2' for 3 models in example BAYAREALIKE_wTraits trait-dependent dispersal script..."
cat("\n")
cat(txt)
archive_vals = c(BAYAREALIKE_t12_t21_m2_m2, BAYAREALIKEj_t12_t21_m2_rep1_m2, BAYAREALIKEj_t12_t21_m2_rep2_m2)
test_vals = c(BAYAREALIKE_t12_t21_m2_results[8], BAYAREALIKEj_t12_t21_m2_results[8], BAYAREALIKEj_t12_t21_m2_rep2_results[8])

expect_equal(object=round(archive_vals, digits=3), expected=unname(round(test_vals, digits=3)))
cat("...PASSED  ")


#######################################################
# Other checks -- "sanity checks" on ML inference with trait-based models
#######################################################
txt = "Checking that the 2rates traits lnL >= 1rates traits lnL (which it must be, because the 1rates model is nested inside 2rates) (archive values)..."
cat("\n")
cat(txt)
expect_gt(object=trait_2rate_lnL, expected=trait_1rate_lnL)
cat("...PASSED  ")

txt = "Checking that the 2rates traits lnL >= 1rates traits lnL (which it must be, because the 1rates model is nested inside 2rates) (test values)..."
cat("\n")
cat(txt)
expect_gt(object=Trait_2rates_results[1], expected=Trait_1rate_results[1])
cat("...PASSED  ")


txt = "Checking that BAYAREALIKE+J lnL >= BAYAREALIKE lnL (which it must be, as BAYAREALIKE is nested inside BAYAREALIKE+J) (archive values)..."
cat("\n")
cat(txt)
expect_gt(object=BAYAREALIKEj_lnL, expected=BAYAREALIKE_lnL)
cat("...PASSED  ")

txt = "Checking that BAYAREALIKE+J lnL >= BAYAREALIKE lnL (which it must be, as BAYAREALIKE is nested inside BAYAREALIKE+J) (test values)..."
cat("\n")
cat(txt)
expect_gt(object=BAYAREALIKEj_results[1], expected=BAYAREALIKE_results[1])
cat("...PASSED  ")






txt = "Checking that BAYAREALIKE+J+t12+t21+m2 lnL (on rep. 2) == BAYAREALIKE+J+t12+t21+m2 lnL (on rep. 1) (which it must be, because these are the same model, with optimization started from two different starting points) (archive values)..."
cat("\n")
cat(txt)
expect_equal(object=round(BAYAREALIKEj_t12_t21_m2_rep2_lnL, digits=3), expected=round(BAYAREALIKEj_t12_t21_m2_rep1_lnL, digits=3))
cat("...PASSED  ")

# lnL
txt = "Checking that BAYAREALIKE+J+t12+t21+m2 lnL (on rep. 2) == BAYAREALIKE+J+t12+t21+m2 lnL (on rep. 1) (which it must be, because these are the same model, with optimization started from two different starting points) (test values)..."
cat("\n")
cat(txt)
expect_equal(object=unname(round(BAYAREALIKEj_t12_t21_m2_rep2_results[1], digits=3)), expected=unname(round(BAYAREALIKEj_t12_t21_m2_results[1], digits=3)))
cat("...PASSED  ")



txt = "Checking that BAYAREALIKE+J+t12+t21+m2 lnL (on rep. 1) >= BAYAREALIKE+t12+t21+m2 lnL (which it must be, as BAYAREALIKE+t12+t21+m2 is nested inside BAYAREALIKE+J+t12+t21+m2) (archive values)..."
cat("\n")
cat(txt)
expect_gt(object=BAYAREALIKEj_t12_t21_m2_rep1_lnL, expected=BAYAREALIKE_t12_t21_m2_lnL)
cat("...PASSED  ")


txt = "Checking that BAYAREALIKE+J+t12+t21+m2 lnL (on rep. 1) >= BAYAREALIKE+t12+t21+m2 lnL (which it must be, as BAYAREALIKE+t12+t21+m2 is nested inside BAYAREALIKE+J+t12+t21+m2) (test values)..."
cat("\n")
cat(txt)
expect_gt(object=BAYAREALIKEj_t12_t21_m2_results[1], expected=BAYAREALIKE_t12_t21_m2_results[1])
cat("...PASSED  ")


# d
txt = "Checking that BAYAREALIKE+J+t12+t21+m2 d (on rep. 2) == BAYAREALIKE+J+t12+t21+m2 d (on rep. 1) (which it must be, because these are the same model, with optimization started from two different starting points)..."
cat("\n")
cat(txt)
expect_equal(object=unname(round(BAYAREALIKEj_t12_t21_m2_rep2_results[2], digits=4)), expected=unname(round(BAYAREALIKEj_t12_t21_m2_results[2], digits=4)))
cat("...PASSED  ")

# e
txt = "Checking that BAYAREALIKE+J+t12+t21+m2 e (on rep. 2) == BAYAREALIKE+J+t12+t21+m2 e (on rep. 1) (which it must be, because these are the same model, with optimization started from two different starting points)..."
cat("\n")
cat(txt)
expect_equal(object=unname(round(BAYAREALIKEj_t12_t21_m2_rep2_results[3], digits=4)), expected=unname(round(BAYAREALIKEj_t12_t21_m2_results[3], digits=4)))
cat("...PASSED  ")

# j
txt = "Checking that BAYAREALIKE+J+t12+t21+m2 j (on rep. 2) == BAYAREALIKE+J+t12+t21+m2 j (on rep. 1) (which it must be, because these are the same model, with optimization started from two different starting points)..."
cat("\n")
cat(txt)
expect_equal(object=unname(round(BAYAREALIKEj_t12_t21_m2_rep2_results[4], digits=1)), expected=unname(round(BAYAREALIKEj_t12_t21_m2_results[4], digits=1)))
cat("...PASSED  ")

# t12
txt = "Checking that BAYAREALIKE+J+t12+t21+m2 t12 (on rep. 2) == BAYAREALIKE+J+t12+t21+m2 t12 (on rep. 1) (which it must be, because these are the same model, with optimization started from two different starting points)..."
cat("\n")
cat(txt)
expect_equal(object=unname(round(BAYAREALIKEj_t12_t21_m2_rep2_results[5], digits=4)), expected=unname(round(BAYAREALIKEj_t12_t21_m2_results[5], digits=4)))
cat("...PASSED  ")

# t21
txt = "Checking that BAYAREALIKE+J+t12+t21+m2 t21 (on rep. 2) == BAYAREALIKE+J+t12+t21+m2 t21 (on rep. 1) (which it must be, because these are the same model, with optimization started from two different starting points)..."
cat("\n")
cat(txt)
expect_equal(object=unname(round(BAYAREALIKEj_t12_t21_m2_rep2_results[6], digits=4)), expected=unname(round(BAYAREALIKEj_t12_t21_m2_results[6], digits=4)))
cat("...PASSED  ")

# m2
txt = "Checking that BAYAREALIKE+J+t12+t21+m2 m2 (on rep. 2) == BAYAREALIKE+J+t12+t21+m2 m2 (on rep. 1) (which it must be, because these are the same model, with optimization started from two different starting points)..."
cat("\n")
cat(txt)
expect_equal(object=unname(round(BAYAREALIKEj_t12_t21_m2_rep2_results[8], digits=1)), expected=unname(round(BAYAREALIKEj_t12_t21_m2_results[8], digits=1)))
cat("...PASSED  ")


cat("\n")
}) # END test_that



