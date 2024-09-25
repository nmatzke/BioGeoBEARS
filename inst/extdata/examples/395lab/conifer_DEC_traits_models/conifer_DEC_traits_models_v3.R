#######################################################
# Setup - beginning
#######################################################

# Install some necessary R packages from CRAN (these are precompiled,
# so easy to install on all systems)

# Run the "install.packages" commands inside the '' marks
# (but WITHOUT the '' marks ONCE.
run_text_inside_quote_once = '
install.packages("ape")
install.packages("devtools")
install.packages("Rcpp")
install.packages("ape")
install.packages("FD")
install.packages("snow")
install.packages("rexpokit")
install.packages("cladoRcpp")

install.packages("phytools")
install.packa
ges("phangorn")
install.packages("phylobase")
install.packages("optimx")
install.packages("GenSA")


################################################
# OLD:
# Install BioGeoBEARS from GitHub
################################################
# (BioGeoBEARS is pure R, so installation is easy *if* the above 
#  packages have been installed)
#library(devtools)
#install_github(repo="nmatzke/BioGeoBEARS", upgrade="never")


################################################
# NEW:
# To avoid GitHub overload, download BioGeoBEARS from Canvas, then install locally:
################################################
# Download from:
# Canvas -> BIOSCI 395 -> Files -> Southern_Conifer_Biogeog
# https://canvas.auckland.ac.nz/courses/106014/files/folder/Southern_conifer_biogeog
# BioGeoBEARS_1.1.3.tar.gz
install.packages("BioGeoBEARS_1.1.3.tar.gz", repos=NULL, type="source")

' # END installation commands


# Check that your BioGeoBEARS installation loads
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)


# Get the file locations for 395 lab files, inside the BioGeoBEARS installation
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
extdata_dir
list.files(extdata_dir)

# Get the locations of the 395 lab files from GitHub install
labdir = paste(extdata_dir, "examples/395lab/", sep="/")
labpt1a = paste(extdata_dir, "examples/395lab/Psychotria_M0_equalRates/", sep="/")
labpt1b = paste(extdata_dir, "examples/395lab/Psychotria_M2_oneWayDispersal/", sep="/")
labpt1c = paste(extdata_dir, "examples/395lab/Psychotria_M4_DistanceDispersal/", sep="/")
labpt2a = paste(extdata_dir, "examples/395lab/conifer_DEC_traits_models/", sep="/")
labpt2b = paste(extdata_dir, "examples/395lab/conifer_DEC+x_traits_models/", sep="/")

labpt1a_script = paste(extdata_dir, "examples/395lab/Psychotria_M0_equalRates/Psychotria_M0_v3a.R", sep="/")
labpt1b_script = paste(extdata_dir, "examples/395lab/Psychotria_M2_oneWayDispersal/Psychotria_M2_oneWayDispersal_v3a.R", sep="/")
labpt1c_script = paste(extdata_dir, "examples/395lab/Psychotria_M4_DistanceDispersal/Psychotria_M4_DistanceDispersal_v3a.R/", sep="/")
labpt2a_script = paste(extdata_dir, "examples/395lab/conifer_DEC_traits_models/conifer_DEC_traits_models_v3a.R", sep="/")
labpt2b_script = paste(extdata_dir, "examples/395lab/conifer_DEC+x_traits_models/conifer_DEC+x_traits_models_v3a.R", sep="/")

#######################################################
# Setup - ENDing
#######################################################


#######################################################
# Loading packages
#######################################################
# Close any open graphics
dev.off(); dev.off(); dev.off(); dev.off(); dev.off(); 

# Load the package (after installation, see above).
library(ape)						# for read.tree
library(GenSA)
library(optimx)         # You need to have some version of optimx available
                        # as it is a BioGeoBEARS dependency; however, if you
                        # don't want to use optimx, and use optim() (from R core) 
                        # you can set:
                        # BioGeoBEARS_run_object$use_optimx = FALSE
                        # ...everything should work either way -- NJM 2014-01-08
library(FD)       # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)
library(BioGeoBEARS)



# Set working directory
#wd = "/drives/GDrive/__classes/BIOSCI395/lab/BGBlab/conifer_DEC_traits_models/"
wd = "~/Downloads/"
setwd(wd)

# Get 395 locations in GitHub install
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
extdata_dir
list.files(extdata_dir)

trfn = slashslash(paste(labpt2a, "tree.newick", sep="/"))
tr = read.tree(trfn)
tr
plot(tr, show.tip.label=FALSE)
axisPhylo()
title("Phylogeny of 197 southern conifers")
mtext(text="Millions of years ago (mega-annum; Ma)", side=1, line=3, cex=1.0)

# Get the locations of the 395 lab files from GitHub install
labdir = paste(extdata_dir, "examples/395lab/", sep="/")
labpt2a = paste(extdata_dir, "examples/395lab/conifer_DEC_traits_models/", sep="/")
labpt2b = paste(extdata_dir, "examples/395lab/conifer_DEC+x_traits_models/", sep="/")



#######################################################
# Inference
#######################################################
max_range_size = 3
geogfn = slashslash(paste(labpt2b, "geog.data", sep="/"))
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
range_sizes = rowSums(dfnums_to_numeric(tipranges@df))
range_sizes[range_sizes > 2]

areas = getareas_from_tipranges_object(tipranges)
#areas = c("A", "B", "C", "D")

# This is the list of states/ranges, where each state/range
# is a list of areas, counting from 0
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)

# Remove all 3-area ranges except DFG / 0-based 123
states_list_0based_NEW = states_list_0based[c(1:22,33)]



#######################################################
# Traits-only model -- 1 rate
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$print_optim = TRUE
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
BioGeoBEARS_run_object$max_range_size = 1
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$use_optimx=TRUE
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = slashslash(paste(labpt2a, "geog_1area.data", sep="/"))
BioGeoBEARS_run_object$trfn = slashslash(paste(labpt2a, "tree.newick", sep="/"))
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1000000
BioGeoBEARS_run_object$force_sparse = FALSE  # works with kexpmv, but compare to dense,
# time-stratify to break up long branches if you see major differences in lnL

# Set up DEC model, but set all rates to 0 (data are 1 invariant area)
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

trait_fn = slashslash(paste(labpt2a, "trait.data", sep="/"))
geog_values = getranges_from_LagrangePHYLIP(trait_fn)

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

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "traitsOnly_1rate_v3a.Rdata"
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
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$use_optimx=TRUE
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = slashslash(paste(labpt2a, "geog_1area.data", sep="/"))
BioGeoBEARS_run_object$trfn = slashslash(paste(labpt2a, "tree.newick", sep="/"))
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1000000
BioGeoBEARS_run_object$force_sparse = FALSE  # works with kexpmv, but compare to dense,
# time-stratify to break up long branches if you see major differences in lnL

# Set up DEC model, but set all rates to 0 (data are 1 invariant area)
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

trait_fn = slashslash(paste(labpt2a, "trait.data", sep="/"))
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

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "traitsOnly_2rates_v3a.Rdata"
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
# Run DEC (on geography only)
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$print_optim = TRUE
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$use_optimx=TRUE
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = slashslash(paste(labpt2a, "geog.data", sep="/"))
BioGeoBEARS_run_object$trfn = slashslash(paste(labpt2a, "tree.newick", sep="/"))
#BioGeoBEARS_run_object$timesfn = "times_v2.txt"
#BioGeoBEARS_run_object$distsfn = "modern_distances_subset.txt"
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1000000
BioGeoBEARS_run_object$force_sparse = FALSE
BioGeoBEARS_run_object$states_list = states_list_0based_NEW
#tr = read.tree(BioGeoBEARS_run_object$trfn)


# Add x as a free parameter
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = 0
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = 0
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","max"] = 0


BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "DEC_inf_v3a.Rdata"
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
# Run DEC+J (on geography only)
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$print_optim = TRUE
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$use_optimx=TRUE
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = slashslash(paste(labpt2a, "geog.data", sep="/"))
BioGeoBEARS_run_object$trfn = slashslash(paste(labpt2a, "tree.newick", sep="/"))
#BioGeoBEARS_run_object$timesfn = "times_v2.txt"
#BioGeoBEARS_run_object$distsfn = "modern_distances_subset.txt"
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1000000
BioGeoBEARS_run_object$force_sparse = FALSE

#tr = read.tree(BioGeoBEARS_run_object$trfn)

dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001
#xstart = resDEC$outputs@params_table["x","est"]

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
BioGeoBEARS_run_object$states_list = states_list_0based_NEW

# Add x as a free parameter
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = xstart
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = xstart
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","max"] = 0
BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

print("Printing warnings: 'warnings()':")
print(warnings())

runslow = TRUE
resfn = "DECj_inf_v3a.Rdata"
if (runslow)
	{
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
# Run DEC + t12 + t21 + m2, starting from DEC-geog and 2rates
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$print_optim = TRUE
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$use_optimx=TRUE
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = slashslash(paste(labpt2a, "geog.data", sep="/"))
BioGeoBEARS_run_object$trfn = slashslash(paste(labpt2a, "tree.newick", sep="/"))
#BioGeoBEARS_run_object$timesfn = "times_v2.txt"
#BioGeoBEARS_run_object$distsfn = "modern_distances_subset.txt"
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1000000
BioGeoBEARS_run_object$force_sparse = FALSE
BioGeoBEARS_run_object$states_list = states_list_0based_NEW

tr = read.tree(BioGeoBEARS_run_object$trfn)
#plot(tr); axisPhylo()

trait_fn = slashslash(paste(labpt2a, "trait.data", sep="/"))
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
dstart = resDEC$outputs@params_table["d","est"]
estart = max(c(resDEC$outputs@params_table["e","est"], 1.1*BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"]))
jstart = 0.0001
#xstart = resDEC$outputs@params_table["x","est"]

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "est"] = estart

# Add x as a free parameter
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = xstart
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = xstart
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","max"] = 0
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","min"] = -2


#######################################################
# Manual modifications of trait-based model
#######################################################
# Edit t12 and t21 rates

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "init"] = t12_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "est"] = t12_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "max"] = round(max(t12_start, t21_start)* 10, 3) 

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "init"] = t21_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "est"] = t21_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "max"] = round(max(t12_start, t21_start)* 10, 3) 


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
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "max"] = 1

BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)





runslow = TRUE
resfn = "DEC+t12+t21+m2_inf_v3a.Rdata"
if (runslow)
	{
	# Calculate the lnL for the parameters, and store in text file
	res = bears_optim_run(BioGeoBEARS_run_object, skip_optim=FALSE, skip_optim_option="return_all")
	save(res, file=resfn)
	resDEC_t12_t21_m2 = res
	} else {
	# Loads to "res"
	load(resfn)
	resDEC_t12_t21_m2 = res
	}




#######################################################
# Run DECj + t12 + t21 + m2, starting from DECj-geog and 2rates
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$print_optim = TRUE
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$use_optimx=TRUE
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = slashslash(paste(labpt2a, "geog.data", sep="/"))
BioGeoBEARS_run_object$trfn = slashslash(paste(labpt2a, "tree.newick", sep="/"))
#BioGeoBEARS_run_object$timesfn = "times_v2.txt"
#BioGeoBEARS_run_object$distsfn = "modern_distances_subset.txt"
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1000000
BioGeoBEARS_run_object$force_sparse = FALSE
BioGeoBEARS_run_object$states_list = states_list_0based_NEW

tr = read.tree(BioGeoBEARS_run_object$trfn)
#plot(tr); axisPhylo()

trait_fn = slashslash(paste(labpt2a, "trait.data", sep="/"))
trait_values = getranges_from_LagrangePHYLIP(lgdata_fn=trait_fn)
trait_values

# Add the traits data and model
BioGeoBEARS_run_object = add_trait_to_BioGeoBEARS_run_object(BioGeoBEARS_run_object, trait_fn=trait_fn)



# Starting values from ML results of simpler run
t12_start = resTrait_2rates$outputs@params_table["t12","est"]
t21_start = resTrait_2rates$outputs@params_table["t21","est"]
m2_start = 1
dstart = resDECj$outputs@params_table["d","est"]
estart = max(c(resDECj$outputs@params_table["e","est"], 1.1*BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"]))
jstart = resDECj$outputs@params_table["j","est"]
#xstart = resDECj$outputs@params_table["x","est"]


# Set up DEC+J model
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "est"] = estart

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Add x as a free parameter
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = xstart
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = xstart
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","max"] = 0
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","min"] = -2

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
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "max"] = round(max(t12_start, t21_start)* 10, 3) 

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "init"] = t21_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "est"] = t21_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "max"] = round(max(t12_start, t21_start)* 10, 3) 


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
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "max"] = 1




BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)



resfn = "DECJ+t12+t21+m2_inf_v3a.Rdata"
runslow = TRUE
if (runslow)
	{
	# Calculate the lnL for the parameters, and store in text file
	res = bears_optim_run(BioGeoBEARS_run_object, skip_optim=FALSE, skip_optim_option="return_all")
	save(res, file=resfn)

	resDECj_t12_t21_m2 = res
	} else {
	# Loads to "res"
	load(resfn)
	resDECj_t12_t21_m2 = res
	}










#######################################################
# Plot: best model to screen (may look squashed)
#######################################################

pdffn = "southern_conifers_DEC+J+trait_v3a.pdf"
pdf(file=pdffn, width=10, height=30)

#######################################################
# Extract just geography ancestral states from geog+trait ancestral states
#######################################################
geog_res = get_geog_from_traitgeog_results(res=resDECj_t12_t21_m2, num_trait_states=2)

#######################################################
# Extract just trait ancestral states from geog+trait ancestral states
#######################################################
trait_res = get_trait_from_traitgeog_results(res=resDECj_t12_t21_m2, num_trait_states=2)
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=slashslash(paste(labpt2a, "geog.data", sep="/")))

#######################################################
# Plot the geographic ancestral states
#######################################################
analysis_titletxt = "Southern conifers: Geog reconstruction under DEC+J+trait"

scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
results_object = geog_res
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)


#######################################################
# Plot the trait ancestral states
#######################################################
analysis_titletxt = "Southern conifers: Trait reconstruction under DEC+J+trait"
results_object = trait_res
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("t12","t21"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=FALSE, tr=tr, tipranges=trait_values)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=FALSE, tr=tr, tipranges=trait_values)


dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)












#######################################################
# Plot: best model to PDF (should look better; but you may have to open the PDF
#       manually, if system(cmdstr) does not successfully open PDF; if you are
#       Rstudio Cloud in a browser, you will have to download the PDF to your 
#       hard drive to view it.)
#######################################################
dev.off(); dev.off(); # (close previous graphics devices to make way for the PDF)

pdffn = "southern_conifers_DEC+J+trait_v3a.pdf"
pdf(file=pdffn, width=10, height=30)

#######################################################
# Extract just geography ancestral states from geog+trait ancestral states
#######################################################
geog_res = get_geog_from_traitgeog_results(res=resDECj_t12_t21_m2, num_trait_states=2)

#######################################################
# Extract just trait ancestral states from geog+trait ancestral states
#######################################################
trait_res = get_trait_from_traitgeog_results(res=resDECj_t12_t21_m2, num_trait_states=2)
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=slashslash(paste(labpt2a, "geog.data", sep="/")))

#######################################################
# Plot the geographic ancestral states
#######################################################
analysis_titletxt = "Southern conifers: Geog reconstruction under DEC+J+trait"

scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
results_object = geog_res
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)


#######################################################
# Plot the trait ancestral states
#######################################################
analysis_titletxt = "Southern conifers: Trait reconstruction under DEC+J+trait"
results_object = trait_res
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("t12","t21"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=FALSE, tr=tr, tipranges=trait_values)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=FALSE, tr=tr, tipranges=trait_values)

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)







#######################################################
# Extract parameters, lnLs, etc. like you did for
# Psychotria, and calculate AICs & AIC weights
#######################################################

resTrait_1rate$total_loglikelihood
resTrait_2rates$total_loglikelihood
resDEC$total_loglikelihood
resDECj$total_loglikelihood
resDEC_t12_t21_m2$total_loglikelihood
resDECj_t12_t21_m2$total_loglikelihood



#######################################################
# Extract M.L. parameter values
#######################################################
params_to_get = c("d", "e", "j", "x", "t12", "t21", "m2")
# Column "est" has the estimates
resTrait_1rate$output@params_table[params_to_get,]

print_param_ests <- function(res, params_to_get)
	{
	paramvals_to_print = res$output@params_table[params_to_get,]$est
	names(paramvals_to_print) = params_to_get
	print(paramvals_to_print)
	}


# Print to screen for pasting into table
param_ests_table = rbind(print_param_ests(resTrait_1rate, params_to_get=params_to_get), 
print_param_ests(resTrait_2rates, params_to_get=params_to_get), 
print_param_ests(resDEC, params_to_get=params_to_get), 
print_param_ests(resDECj, params_to_get=params_to_get), 
print_param_ests(resDEC_t12_t21_m2, params_to_get=params_to_get), 
print_param_ests(resDECj_t12_t21_m2, params_to_get=params_to_get))

param_ests_table_df = adf(param_ests_table)
row.names(param_ests_table_df) = NULL
param_ests_table_df



cat("\n\n")
cat("...end of script, printing any warnings() to screen.")
print(warnings)
cat("\n")

