

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





test_that(desc="Check trait-dependent dispersal inference, for base models derived from 'DIVALIKE'", code={


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

#wd = "/drives/GDrive/__GDrive_projects/2016-12-07_Kristina_Klaus_Podocarpaceae/__doc3/Supplementary_Material/example/DIVALIKE_wTraits/"
#setwd(wd)

# Example directory with files for traits-based analysis
extdata_dir = np(system.file("extdata/examples/trait_examples/DIVALIKE_wTraits", package="BioGeoBEARS"))
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

# Set up DIVALIKE model, but set all rates to 0 (data are 1 invariant area)
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

BioGeoBEARS_run_object$printlevel = 0
BioGeoBEARS_run_object$print_optim = FALSE

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

# Set up DIVALIKE model, but set all rates to 0 (data are 1 invariant area)
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

BioGeoBEARS_run_object$printlevel = 0
BioGeoBEARS_run_object$print_optim = FALSE

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
# Run DIVALIKE (on geography only)
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

# Set up DIVALIKE model
# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# No jump dispersal/founder-event speciation
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01


BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

BioGeoBEARS_run_object$printlevel = 0
BioGeoBEARS_run_object$print_optim = FALSE

runslow = TRUE
resfn = "DIVALIKE_inf.Rdata"
if (runslow)
	{
	res = bears_optim_run(BioGeoBEARS_run_object)
	res    

	save(res, file=resfn)
	resDIVALIKE = res
	} else {
	# Loads to "res"
	load(resfn)
	resDIVALIKE = res
	}


#######################################################
# Run DIVALIKE+J (on geography only)
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

# Set up DIVALIKE+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDIVALIKE$outputs@params_table["d","est"]
estart = resDIVALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# Add jump dispersal/founder-event speciation
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DIVALIKE+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999


BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

BioGeoBEARS_run_object$printlevel = 0
BioGeoBEARS_run_object$print_optim = FALSE

print("Printing warnings: 'warnings()':")
print(warnings())

runslow = TRUE
resfn = "DIVALIKEj_inf.Rdata"
if (runslow)
	{
	res = bears_optim_run(BioGeoBEARS_run_object)
	res    

	save(res, file=resfn)
	resDIVALIKEj = res
	} else {
	# Loads to "res"
	load(resfn)
	resDIVALIKEj = res
	}






#######################################################
# Run DIVALIKE + t12 + t21 + m2, starting from DIVALIKE-geog and 2rates
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

# Set up DIVALIKE model
# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# No jump dispersal/founder-event speciation
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01


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
dstart = resDIVALIKE$outputs@params_table["d","est"]
estart = max(c(resDIVALIKE$outputs@params_table["e","est"], 1.1*BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"]))
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

BioGeoBEARS_run_object$printlevel = 0
BioGeoBEARS_run_object$print_optim = FALSE


runslow = TRUE
resfn = "DIVALIKE+t12+t21+m2_inf.Rdata"
if (runslow)
	{
	res = bears_optim_run(BioGeoBEARS_run_object, skip_optim=FALSE)
	
	save(res, file=resfn)
	resDIVALIKE_t12_t21_m2 = res
	} else {
	# Loads to "res"
	load(resfn)
	resDIVALIKE_t12_t21_m2 = res
	}




#######################################################
# Run DIVALIKEj + t12 + t21 + m2, starting from DIVALIKEj-geog and 2rates
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


# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# Add jump dispersal/founder-event speciation
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DIVALIKE+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999



trait_fn = traitsfn
trait_values = getranges_from_LagrangePHYLIP(lgdata_fn=trait_fn)
trait_values

# Add the traits data and model
BioGeoBEARS_run_object = add_trait_to_BioGeoBEARS_run_object(BioGeoBEARS_run_object, trait_fn=trait_fn)



# Starting values from ML results of simpler run
t12_start = resTrait_2rates$outputs@params_table["t12","est"]
t21_start = resTrait_2rates$outputs@params_table["t21","est"]
m2_start = 1
dstart = resDIVALIKEj$outputs@params_table["d","est"]
estart = max(c(resDIVALIKEj$outputs@params_table["e","est"], 1.1*BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"]))
jstart = resDIVALIKEj$outputs@params_table["j","est"]


# Set up DIVALIKE+J model
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

BioGeoBEARS_run_object$printlevel = 0
BioGeoBEARS_run_object$print_optim = FALSE



resfn = "DIVALIKEJ+t12+t21+m2_inf.Rdata"
runslow = TRUE
if (runslow)
	{
	res = bears_optim_run(BioGeoBEARS_run_object, skip_optim=FALSE)

	save(res, file=resfn)
	resDIVALIKEj_t12_t21_m2 = res
	} else {
	# Loads to "res"
	load(resfn)
	resDIVALIKEj_t12_t21_m2 = res
	}







#######################################################
# Run DIVALIKEj + t12 + t21 + m2, starting from DIVALIKE + t12 + t21 + m2
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


# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# Add jump dispersal/founder-event speciation
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DIVALIKE+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999



trait_fn = traitsfn
trait_values = getranges_from_LagrangePHYLIP(lgdata_fn=trait_fn)
trait_values

# Add the traits data and model
BioGeoBEARS_run_object = add_trait_to_BioGeoBEARS_run_object(BioGeoBEARS_run_object, trait_fn=trait_fn)



# Starting values from ML results of simpler run
t12_start = resDIVALIKE_t12_t21_m2$outputs@params_table["t12","est"]
t21_start = resDIVALIKE_t12_t21_m2$outputs@params_table["t21","est"]
m2_start = resDIVALIKE_t12_t21_m2$outputs@params_table["m2","est"]
dstart = resDIVALIKE_t12_t21_m2$outputs@params_table["d","est"]
estart = max(c(resDIVALIKE_t12_t21_m2$outputs@params_table["e","est"], 1.1*BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"]))
jstart = 0.0001


# Set up DIVALIKE+J model
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

BioGeoBEARS_run_object$printlevel = 0
BioGeoBEARS_run_object$print_optim = FALSE

resfn = "DIVALIKEJ+t12+t21+m2_rep2_inf.Rdata"
runslow = TRUE
if (runslow)
	{
	res = bears_optim_run(BioGeoBEARS_run_object)
	res    

	save(res, file=resfn)

	resDIVALIKEj_t12_t21_m2_rep2 = res
	} else {
	# Loads to "res"
	load(resfn)
	resDIVALIKEj_t12_t21_m2_rep2 = res
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

DIVALIKE_results = c(
resDIVALIKE$total_loglikelihood,
resDIVALIKE$output@params_table["d", "est"], 
resDIVALIKE$output@params_table["e", "est"], 
resDIVALIKE$output@params_table["j", "est"], 
resDIVALIKE$output@params_table["t12", "est"], 
resDIVALIKE$output@params_table["t21", "est"], 
resDIVALIKE$output@params_table["m1", "est"], 
resDIVALIKE$output@params_table["m2", "est"]
)
names(DIVALIKE_results) = paste("DIVALIKE_", param_names, sep="")


DIVALIKEj_results = c(
resDIVALIKEj$total_loglikelihood,
resDIVALIKEj$output@params_table["d", "est"], 
resDIVALIKEj$output@params_table["e", "est"], 
resDIVALIKEj$output@params_table["j", "est"], 
resDIVALIKEj$output@params_table["t12", "est"], 
resDIVALIKEj$output@params_table["t21", "est"], 
resDIVALIKEj$output@params_table["m1", "est"], 
resDIVALIKEj$output@params_table["m2", "est"]
)
names(DIVALIKEj_results) = paste("DIVALIKEj_", param_names, sep="")

DIVALIKE_t12_t21_m2_results = c(
resDIVALIKE_t12_t21_m2$total_loglikelihood,
resDIVALIKE_t12_t21_m2$output@params_table["d", "est"], 
resDIVALIKE_t12_t21_m2$output@params_table["e", "est"], 
resDIVALIKE_t12_t21_m2$output@params_table["j", "est"], 
resDIVALIKE_t12_t21_m2$output@params_table["t12", "est"], 
resDIVALIKE_t12_t21_m2$output@params_table["t21", "est"], 
resDIVALIKE_t12_t21_m2$output@params_table["m1", "est"], 
resDIVALIKE_t12_t21_m2$output@params_table["m2", "est"]
)
names(DIVALIKE_t12_t21_m2_results) = paste("DIVALIKE_t12_t21_m2_", param_names, sep="")

DIVALIKEj_t12_t21_m2_results = c(
resDIVALIKEj_t12_t21_m2$total_loglikelihood,
resDIVALIKEj_t12_t21_m2$output@params_table["d", "est"], 
resDIVALIKEj_t12_t21_m2$output@params_table["e", "est"], 
resDIVALIKEj_t12_t21_m2$output@params_table["j", "est"], 
resDIVALIKEj_t12_t21_m2$output@params_table["t12", "est"], 
resDIVALIKEj_t12_t21_m2$output@params_table["t21", "est"], 
resDIVALIKEj_t12_t21_m2$output@params_table["m1", "est"], 
resDIVALIKEj_t12_t21_m2$output@params_table["m2", "est"]
)
names(DIVALIKEj_t12_t21_m2_results) = paste("DIVALIKEj_t12_t21_m2_", param_names, sep="")

DIVALIKEj_t12_t21_m2_rep2_results = c(
resDIVALIKEj_t12_t21_m2_rep2$total_loglikelihood,
resDIVALIKEj_t12_t21_m2_rep2$output@params_table["d", "est"], 
resDIVALIKEj_t12_t21_m2_rep2$output@params_table["e", "est"], 
resDIVALIKEj_t12_t21_m2_rep2$output@params_table["j", "est"], 
resDIVALIKEj_t12_t21_m2_rep2$output@params_table["t12", "est"], 
resDIVALIKEj_t12_t21_m2_rep2$output@params_table["t21", "est"], 
resDIVALIKEj_t12_t21_m2_rep2$output@params_table["m1", "est"], 
resDIVALIKEj_t12_t21_m2_rep2$output@params_table["m2", "est"]
)
names(DIVALIKEj_t12_t21_m2_rep2_results) = paste("DIVALIKEj_t12_t21_m2_rep2_", param_names, sep="")


tmp_results = c(Trait_1rate_results, Trait_2rates_results, DIVALIKE_results, DIVALIKEj_results, DIVALIKE_t12_t21_m2_results, DIVALIKEj_t12_t21_m2_results, DIVALIKEj_t12_t21_m2_rep2_results)
tmp_results_mat = as.matrix(tmp_results, nrow=7, byrow=TRUE)
tmp_results_mat = as.data.frame(tmp_results_mat, stringsAsFactors=FALSE)

outfn = slashslash("params_inferred.txt")
write.table(x=tmp_results, file=outfn, append=FALSE, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
#moref(outfn)



# Rdata_fn = "DIVALIKE+t12+t21+m2_inf.Rdata"
# rerun_optim_table_DIVALIKE = rerun_optimization_w_HiLow(res=NULL, Rdata_fn=Rdata_fn, runslow=TRUE)
# 
# Rdata_fn = "DIVALIKEJ+t12+t21+m2_inf.Rdata"
# rerun_optim_table_DIVALIKEj = rerun_optimization_w_HiLow(res=NULL, Rdata_fn=Rdata_fn, runslow=TRUE)
# 
# Rdata_fn = "DIVALIKEJ+t12+t21+m2_rep2_inf.Rdata"
# rerun_optim_table_DIVALIKEj_rep2 = rerun_optimization_w_HiLow(res=NULL, Rdata_fn=Rdata_fn, runslow=TRUE)


#####################################
# Archived answers
#####################################
trait_1rate_lnL = -22.26011
trait_1rate_t12 = 0.01194949
trait_1rate_t21 = 0.01194949

trait_2rate_lnL = -22.25162
trait_2rate_t12 = 0.01244680
trait_2rate_t21 = 0.01140128

DIVALIKE_lnL = -98.3675464013305
DIVALIKE_d = 0.0116422984205537
DIVALIKE_e = 1e-12
DIVALIKE_j = 0.000000e+00

DIVALIKEj_lnL = -96.2306138973692
DIVALIKEj_d = 0.00941669853856458
DIVALIKEj_e = 1e-12
DIVALIKEj_j = 0.0491188224855976

DIVALIKE_t12_t21_m2_lnL = -119.725602339355
DIVALIKE_t12_t21_m2_d = 0.0147351137021721
DIVALIKE_t12_t21_m2_e = 1.000000e-12
DIVALIKE_t12_t21_m2_j = 0.000000e+00
DIVALIKE_t12_t21_m2_t12 = 0.0122855084043226
DIVALIKE_t12_t21_m2_t21 = 0.0111865810098732
DIVALIKE_t12_t21_m2_m2 = 0.588197547605291

DIVALIKEj_t12_t21_m2_rep1_lnL = -117.09978122649
DIVALIKEj_t12_t21_m2_rep1_d = 0.0123404589203886
DIVALIKEj_t12_t21_m2_rep1_e = 1.000000e-13
DIVALIKEj_t12_t21_m2_rep1_j = 0.0841168179750863
DIVALIKEj_t12_t21_m2_rep1_t12 = 0.0124688975806053
DIVALIKEj_t12_t21_m2_rep1_t21 = 0.0112368654697434
DIVALIKEj_t12_t21_m2_rep1_m2 = 0.48893052313498

DIVALIKEj_t12_t21_m2_rep2_lnL = -117.099777296057
DIVALIKEj_t12_t21_m2_rep2_d = 0.0123454457815759
DIVALIKEj_t12_t21_m2_rep2_e = 1.000000e-13
DIVALIKEj_t12_t21_m2_rep2_j = 0.084232734376256
DIVALIKEj_t12_t21_m2_rep2_t12 = 0.0124742258672918
DIVALIKEj_t12_t21_m2_rep2_t21 = 0.0112314563505419
DIVALIKEj_t12_t21_m2_rep2_m2 = 0.488081895137516


txt = "Checking the maximum lnL inferences for 7 models in example DIVALIKE_wTraits trait-dependent dispersal script..."
cat("\n")
cat(txt)
archive_lnLs = c(trait_1rate_lnL, trait_2rate_lnL, DIVALIKE_lnL, DIVALIKEj_lnL, DIVALIKE_t12_t21_m2_lnL, DIVALIKEj_t12_t21_m2_rep1_lnL, DIVALIKEj_t12_t21_m2_rep2_lnL)
test_lnLs = c(resTrait_1rate$total_loglikelihood, resTrait_2rates$total_loglikelihood, resDIVALIKE$total_loglikelihood, resDIVALIKEj$total_loglikelihood, resDIVALIKE_t12_t21_m2$total_loglikelihood, resDIVALIKEj_t12_t21_m2$total_loglikelihood, resDIVALIKEj_t12_t21_m2_rep2$total_loglikelihood)

expect_equal(object=round(archive_lnLs, digits=4), expected=round(test_lnLs, digits=4))
cat("...PASSED  ")



txt = "Checking the ML inferences of 'd' for 5 models in example DIVALIKE_wTraits trait-dependent dispersal script..."
cat("\n")
cat(txt)
archive_vals = c(DIVALIKE_d, DIVALIKEj_d, DIVALIKE_t12_t21_m2_d, DIVALIKEj_t12_t21_m2_rep1_d, DIVALIKEj_t12_t21_m2_rep2_d)
test_vals = c(DIVALIKE_results[2], DIVALIKEj_results[2], DIVALIKE_t12_t21_m2_results[2], DIVALIKEj_t12_t21_m2_results[2], DIVALIKEj_t12_t21_m2_rep2_results[2])

expect_equal(object=round(archive_vals, digits=4), expected=unname(round(test_vals, digits=4)))
cat("...PASSED  ")


txt = "Checking the ML inferences of 'e' for 5 models in example DIVALIKE_wTraits trait-dependent dispersal script..."
cat("\n")
cat(txt)
archive_vals = c(DIVALIKE_e, DIVALIKEj_e, DIVALIKE_t12_t21_m2_e, DIVALIKEj_t12_t21_m2_rep1_e, DIVALIKEj_t12_t21_m2_rep2_e)
test_vals = c(DIVALIKE_results[3], DIVALIKEj_results[3], DIVALIKE_t12_t21_m2_results[3], DIVALIKEj_t12_t21_m2_results[3], DIVALIKEj_t12_t21_m2_rep2_results[3])

expect_equal(object=round(archive_vals, digits=4), expected=unname(round(test_vals, digits=4)))
cat("...PASSED  ")


txt = "Checking the ML inferences of 'j' for 5 models in example DIVALIKE_wTraits trait-dependent dispersal script..."
cat("\n")
cat(txt)
archive_vals = c(DIVALIKE_j, DIVALIKEj_j, DIVALIKE_t12_t21_m2_j, DIVALIKEj_t12_t21_m2_rep1_j, DIVALIKEj_t12_t21_m2_rep2_j)
test_vals = c(DIVALIKE_results[4], DIVALIKEj_results[4], DIVALIKE_t12_t21_m2_results[4], DIVALIKEj_t12_t21_m2_results[4], DIVALIKEj_t12_t21_m2_rep2_results[4])

expect_equal(object=round(archive_vals, digits=4), expected=unname(round(test_vals, digits=4)))
cat("...PASSED  ")





txt = "Checking the ML inferences of 't12' for 5 models in example DIVALIKE_wTraits trait-dependent dispersal script..."
cat("\n")
cat(txt)
archive_vals = c(trait_1rate_t12, trait_2rate_t12, DIVALIKE_t12_t21_m2_t12, DIVALIKEj_t12_t21_m2_rep1_t12, DIVALIKEj_t12_t21_m2_rep2_t12)
test_vals = c(Trait_1rate_results[5], Trait_2rates_results[5], DIVALIKE_t12_t21_m2_results[5], DIVALIKEj_t12_t21_m2_results[5], DIVALIKEj_t12_t21_m2_rep2_results[5])

expect_equal(object=round(archive_vals, digits=3), expected=unname(round(test_vals, digits=3)))
cat("...PASSED  ")


txt = "Checking the ML inferences of 't21' for 5 models in example DIVALIKE_wTraits trait-dependent dispersal script..."
cat("\n")
cat(txt)
archive_vals = c(trait_1rate_t21, trait_2rate_t21, DIVALIKE_t12_t21_m2_t21, DIVALIKEj_t12_t21_m2_rep1_t21, DIVALIKEj_t12_t21_m2_rep2_t21)
test_vals = c(Trait_1rate_results[6], Trait_2rates_results[6], DIVALIKE_t12_t21_m2_results[6], DIVALIKEj_t12_t21_m2_results[6], DIVALIKEj_t12_t21_m2_rep2_results[6])

expect_equal(object=round(archive_vals, digits=3), expected=unname(round(test_vals, digits=3)))
cat("...PASSED  ")



txt = "Checking the ML inferences of 'm2' for 3 models in example DIVALIKE_wTraits trait-dependent dispersal script..."
cat("\n")
cat(txt)
archive_vals = c(DIVALIKE_t12_t21_m2_m2, DIVALIKEj_t12_t21_m2_rep1_m2, DIVALIKEj_t12_t21_m2_rep2_m2)
test_vals = c(DIVALIKE_t12_t21_m2_results[8], DIVALIKEj_t12_t21_m2_results[8], DIVALIKEj_t12_t21_m2_rep2_results[8])

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


txt = "Checking that DIVALIKE+J lnL >= DIVALIKE lnL (which it must be, as DIVALIKE is nested inside DIVALIKE+J) (archive values)..."
cat("\n")
cat(txt)
expect_gt(object=DIVALIKEj_lnL, expected=DIVALIKE_lnL)
cat("...PASSED  ")

txt = "Checking that DIVALIKE+J lnL >= DIVALIKE lnL (which it must be, as DIVALIKE is nested inside DIVALIKE+J) (test values)..."
cat("\n")
cat(txt)
expect_gt(object=DIVALIKEj_results[1], expected=DIVALIKE_results[1])
cat("...PASSED  ")






txt = "Checking that DIVALIKE+J+t12+t21+m2 lnL (on rep. 2) == DIVALIKE+J+t12+t21+m2 lnL (on rep. 1) (which it must be, because these are the same model, with optimization started from two different starting points) (archive values)..."
cat("\n")
cat(txt)
expect_equal(object=round(DIVALIKEj_t12_t21_m2_rep2_lnL, digits=3), expected=round(DIVALIKEj_t12_t21_m2_rep1_lnL, digits=3))
cat("...PASSED  ")

# lnL
txt = "Checking that DIVALIKE+J+t12+t21+m2 lnL (on rep. 2) == DIVALIKE+J+t12+t21+m2 lnL (on rep. 1) (which it must be, because these are the same model, with optimization started from two different starting points) (test values)..."
cat("\n")
cat(txt)
expect_equal(object=unname(round(DIVALIKEj_t12_t21_m2_rep2_results[1], digits=4)), expected=unname(round(DIVALIKEj_t12_t21_m2_results[1], digits=4)))
cat("...PASSED  ")



txt = "Checking that DIVALIKE+J+t12+t21+m2 lnL (on rep. 1) >= DIVALIKE+t12+t21+m2 lnL (which it must be, as DIVALIKE+t12+t21+m2 is nested inside DIVALIKE+J+t12+t21+m2) (archive values)..."
cat("\n")
cat(txt)
expect_gt(object=DIVALIKEj_t12_t21_m2_rep1_lnL, expected=DIVALIKE_t12_t21_m2_lnL)
cat("...PASSED  ")


txt = "Checking that DIVALIKE+J+t12+t21+m2 lnL (on rep. 1) >= DIVALIKE+t12+t21+m2 lnL (which it must be, as DIVALIKE+t12+t21+m2 is nested inside DIVALIKE+J+t12+t21+m2) (test values)..."
cat("\n")
cat(txt)
expect_gt(object=DIVALIKEj_t12_t21_m2_results[1], expected=DIVALIKE_t12_t21_m2_results[1])
cat("...PASSED  ")



# d
txt = "Checking that DIVALIKE+J+t12+t21+m2 d (on rep. 2) == DIVALIKE+J+t12+t21+m2 d (on rep. 1) (which it must be, because these are the same model, with optimization started from two different starting points)..."
cat("\n")
cat(txt)
expect_equal(object=unname(round(DIVALIKEj_t12_t21_m2_rep2_results[2], digits=4)), expected=unname(round(DIVALIKEj_t12_t21_m2_results[2], digits=4)))
cat("...PASSED  ")

# e
txt = "Checking that DIVALIKE+J+t12+t21+m2 e (on rep. 2) == DIVALIKE+J+t12+t21+m2 e (on rep. 1) (which it must be, because these are the same model, with optimization started from two different starting points)..."
cat("\n")
cat(txt)
expect_equal(object=unname(round(DIVALIKEj_t12_t21_m2_rep2_results[3], digits=4)), expected=unname(round(DIVALIKEj_t12_t21_m2_results[3], digits=4)))
cat("...PASSED  ")

# j
txt = "Checking that DIVALIKE+J+t12+t21+m2 j (on rep. 2) == DIVALIKE+J+t12+t21+m2 j (on rep. 1) (which it must be, because these are the same model, with optimization started from two different starting points)..."
cat("\n")
cat(txt)
expect_equal(object=unname(round(DIVALIKEj_t12_t21_m2_rep2_results[4], digits=3)), expected=unname(round(DIVALIKEj_t12_t21_m2_results[4], digits=3)))
cat("...PASSED  ")

# t12
txt = "Checking that DIVALIKE+J+t12+t21+m2 t12 (on rep. 2) == DIVALIKE+J+t12+t21+m2 t12 (on rep. 1) (which it must be, because these are the same model, with optimization started from two different starting points)..."
cat("\n")
cat(txt)
expect_equal(object=unname(round(DIVALIKEj_t12_t21_m2_rep2_results[5], digits=4)), expected=unname(round(DIVALIKEj_t12_t21_m2_results[5], digits=4)))
cat("...PASSED  ")

# t21
txt = "Checking that DIVALIKE+J+t12+t21+m2 t21 (on rep. 2) == DIVALIKE+J+t12+t21+m2 t21 (on rep. 1) (which it must be, because these are the same model, with optimization started from two different starting points)..."
cat("\n")
cat(txt)
expect_equal(object=unname(round(DIVALIKEj_t12_t21_m2_rep2_results[6], digits=4)), expected=unname(round(DIVALIKEj_t12_t21_m2_results[6], digits=4)))
cat("...PASSED  ")

# m2
txt = "Checking that DIVALIKE+J+t12+t21+m2 m2 (on rep. 2) == DIVALIKE+J+t12+t21+m2 m2 (on rep. 1) (which it must be, because these are the same model, with optimization started from two different starting points)..."
cat("\n")
cat(txt)
expect_equal(object=unname(round(DIVALIKEj_t12_t21_m2_rep2_results[8], digits=2)), expected=unname(round(DIVALIKEj_t12_t21_m2_results[8], digits=2)))
cat("...PASSED  ")



cat("\n")
}) # END test_that



