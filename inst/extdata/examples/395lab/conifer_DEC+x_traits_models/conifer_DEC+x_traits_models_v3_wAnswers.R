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

labpt1a_script = paste(extdata_dir, "examples/395lab/Psychotria_M0_equalRates/Psychotria_M0_v3b.R", sep="/")
labpt1b_script = paste(extdata_dir, "examples/395lab/Psychotria_M2_oneWayDispersal/Psychotria_M2_oneWayDispersal_v3b.R", sep="/")
labpt1c_script = paste(extdata_dir, "examples/395lab/Psychotria_M4_DistanceDispersal/Psychotria_M4_DistanceDispersal_v3b.R/", sep="/")
labpt2a_script = paste(extdata_dir, "examples/395lab/conifer_DEC_traits_models/conifer_DEC_traits_models_v3b.R", sep="/")
labpt2b_script = paste(extdata_dir, "examples/395lab/conifer_DEC+x_traits_models/conifer_DEC+x_traits_models_v3b.R", sep="/")

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
wd = "~/Downloads/395lab23/"
setwd(wd)

# Get 395 locations in GitHub install
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
extdata_dir
list.files(extdata_dir)

trfn = slashslash(paste(labpt2b, "tree.newick", sep="/"))
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
BioGeoBEARS_run_object$num_cores_to_use = 4
BioGeoBEARS_run_object$use_optimx=TRUE
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = slashslash(paste(labpt2b, "geog_1area.data", sep="/"))
BioGeoBEARS_run_object$trfn = slashslash(paste(labpt2b, "tree.newick", sep="/"))
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

trait_fn = slashslash(paste(labpt2b, "trait.data", sep="/"))
geog_values = getranges_from_LagrangePHYLIP(trait_fn)

trait_fn = slashslash(paste(labpt2b, "trait.data", sep="/"))
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
resfn = "traitsOnly_1rate_v3b.Rdata"
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
BioGeoBEARS_run_object$num_cores_to_use = 4
BioGeoBEARS_run_object$use_optimx=TRUE
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = slashslash(paste(labpt2b, "geog_1area.data", sep="/"))
BioGeoBEARS_run_object$trfn = slashslash(paste(labpt2b, "tree.newick", sep="/"))
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

trait_fn = slashslash(paste(labpt2b, "trait.data", sep="/"))
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
resfn = "traitsOnly_2rates_v3b.Rdata"
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
BioGeoBEARS_run_object$num_cores_to_use = 4
BioGeoBEARS_run_object$use_optimx=TRUE
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = slashslash(paste(labpt2b, "geog.data", sep="/"))
BioGeoBEARS_run_object$trfn = slashslash(paste(labpt2b, "tree.newick", sep="/"))
#BioGeoBEARS_run_object$timesfn = "times_v2.txt"
BioGeoBEARS_run_object$distsfn = slashslash(paste(labpt2b, "modern_distances_subset.txt", sep="/"))
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1000000
BioGeoBEARS_run_object$force_sparse = FALSE

#tr = read.tree(BioGeoBEARS_run_object$trfn)


# Add x as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","max"] = 0


BioGeoBEARS_run_object$states_list = states_list_0based_NEW

BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "DECx_inf_v3b.Rdata"
if (runslow)
	{
	res = bears_optim_run(BioGeoBEARS_run_object)
	res    

	save(res, file=resfn)
	resDECx = res
	} else {
	# Loads to "res"
	load(resfn)
	resDECx = res
	}


#######################################################
# Run DEC+J (on geography only)
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$print_optim = TRUE
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$num_cores_to_use = 4
BioGeoBEARS_run_object$use_optimx=TRUE
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = slashslash(paste(labpt2b, "geog.data", sep="/"))
BioGeoBEARS_run_object$trfn = slashslash(paste(labpt2b, "tree.newick", sep="/"))
#BioGeoBEARS_run_object$timesfn = "times_v2.txt"
BioGeoBEARS_run_object$distsfn = slashslash(paste(labpt2b, "modern_distances_subset.txt", sep="/"))
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1000000
BioGeoBEARS_run_object$force_sparse = FALSE

#tr = read.tree(BioGeoBEARS_run_object$trfn)

dstart = resDECx$outputs@params_table["d","est"]
estart = resDECx$outputs@params_table["e","est"]
jstart = 0.0001
xstart = resDECx$outputs@params_table["x","est"]

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Add x as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = xstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = xstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","max"] = 0

BioGeoBEARS_run_object$states_list = states_list_0based_NEW


BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

print("Printing warnings: 'warnings()':")
print(warnings())

runslow = TRUE
resfn = "DECxj_inf_v3b.Rdata"
if (runslow)
	{
	res = bears_optim_run(BioGeoBEARS_run_object)
	res    

	save(res, file=resfn)
	resDECxj = res
	} else {
	# Loads to "res"
	load(resfn)
	resDECxj = res
	}






#######################################################
# Run DEC + x + t12 + t21 + m2, starting from DEC-geog and 2rates
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$print_optim = TRUE
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$num_cores_to_use = 4
BioGeoBEARS_run_object$use_optimx=TRUE
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = slashslash(paste(labpt2b, "geog.data", sep="/"))
BioGeoBEARS_run_object$trfn = slashslash(paste(labpt2b, "tree.newick", sep="/"))
#BioGeoBEARS_run_object$timesfn = "times_v2.txt"
BioGeoBEARS_run_object$distsfn = slashslash(paste(labpt2b, "modern_distances_subset.txt", sep="/"))
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



trait_fn = slashslash(paste(labpt2b, "trait.data", sep="/"))
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
#dstart = resDECx$outputs@params_table["d","est"]
#estart = max(c(resDECx$outputs@params_table["e","est"], 1.1*BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"]))
#jstart = 0.0001
#xstart = resDECx$outputs@params_table["x","est"]
dstart = 0.2
estart = 0.004
jstart = 0.0001
xstart = -1.2

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "est"] = estart

# Add x as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = xstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = xstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","max"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","min"] = -2


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



# Calculate the lnL for the parameters, and store in text file


runslow = TRUE
resfn = "DECx+t12+t21+m2_inf_v3b.Rdata"
if (runslow)
	{
	# Calculate the lnL for the parameters, and store in text file
	res = bears_optim_run(BioGeoBEARS_run_object, skip_optim=FALSE, skip_optim_option="return_all")

	save(res, file=resfn)
	resDECx_t12_t21_m2 = res
	} else {
	# Loads to "res"
	load(resfn)
	resDECx_t12_t21_m2 = res
	}




#######################################################
# Run DECj + x + t12 + t21 + m2, starting from DECj-geog and 2rates
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$print_optim = TRUE
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$num_cores_to_use = 4
BioGeoBEARS_run_object$use_optimx=TRUE
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = slashslash(paste(labpt2b, "geog.data", sep="/"))
BioGeoBEARS_run_object$trfn = slashslash(paste(labpt2b, "tree.newick", sep="/"))
#BioGeoBEARS_run_object$timesfn = "times_v2.txt"
BioGeoBEARS_run_object$distsfn = slashslash(paste(labpt2b, "modern_distances_subset.txt", sep="/"))
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

trait_fn = slashslash(paste(labpt2b, "trait.data", sep="/"))
trait_values = getranges_from_LagrangePHYLIP(lgdata_fn=trait_fn)
trait_values

# Add the traits data and model
BioGeoBEARS_run_object = add_trait_to_BioGeoBEARS_run_object(BioGeoBEARS_run_object, trait_fn=trait_fn)



# Starting values from ML results of simpler run
t12_start = resTrait_2rates$outputs@params_table["t12","est"]
t21_start = resTrait_2rates$outputs@params_table["t21","est"]
m2_start = 0.5
#dstart = resDECxj$outputs@params_table["d","est"]
#estart = max(c(resDECxj$outputs@params_table["e","est"], 1.1*BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"]))
#jstart = resDECxj$outputs@params_table["j","est"]
#xstart = resDECxj$outputs@params_table["x","est"]
dstart = 0.08
estart = 0.00001
jstart = 0.6
xstart = -1

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
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = xstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = xstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","max"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","min"] = -2

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



resfn = "DECJx+t12+t21+m2_inf_v3b.Rdata"
runslow = TRUE
if (runslow)
	{
	# Calculate the lnL for the parameters, and store in text file
	res = bears_optim_run(BioGeoBEARS_run_object, skip_optim=FALSE, skip_optim_option="return_all")
	save(res, file=resfn)

	resDECxj_t12_t21_m2 = res
	} else {
	# Loads to "res"
	load(resfn)
	resDECxj_t12_t21_m2 = res
	}



#######################################################
# Plot: best model to screen (may look squashed)
#######################################################

pdffn = "southern_conifers_DEC+J+x+trait_v3b.pdf"
#pdf(file=pdffn, width=10, height=30)


#######################################################
# Extract just geography ancestral states from geog+trait ancestral states
#######################################################
geog_res = get_geog_from_traitgeog_results(res=resDECxj_t12_t21_m2, num_trait_states=2)

#######################################################
# Extract just trait ancestral states from geog+trait ancestral states
#######################################################
trait_res = get_trait_from_traitgeog_results(res=resDECxj_t12_t21_m2, num_trait_states=2)
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=slashslash(paste(labpt2a, "geog.data", sep="/")))

#######################################################
# Plot the geographic ancestral states
#######################################################
analysis_titletxt = "Southern conifers: Geog reconstruction under DEC+J+x+trait"

scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
results_object = geog_res
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
# plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)


#######################################################
# Plot the trait ancestral states
#######################################################
analysis_titletxt = "Southern conifers: Trait reconstruction under DEC+J+x+trait"
results_object = trait_res
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("t12","t21"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=FALSE, tr=tr, tipranges=trait_values)

# Pie chart
# plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=FALSE, tr=tr, tipranges=trait_values)


#dev.off()
# cmdstr = paste0("open ", pdffn)
# system(cmdstr)





#######################################################
# Plot: best model to PDF (should look better; but you may have to open the PDF
#       manually, if system(cmdstr) does not successfully open PDF; if you are
#       Rstudio Cloud in a browser, you will have to download the PDF to your 
#       hard drive to view it.)
#######################################################
dev.off(); dev.off(); # (close previous graphics devices to make way for the PDF)

pdffn = "southern_conifers_DEC+J+x+trait_v3b.pdf"
pdf(file=pdffn, width=10, height=30)


#######################################################
# Extract just geography ancestral states from geog+trait ancestral states
#######################################################
geog_res = get_geog_from_traitgeog_results(res=resDECxj_t12_t21_m2, num_trait_states=2)

#######################################################
# Extract just trait ancestral states from geog+trait ancestral states
#######################################################
trait_res = get_trait_from_traitgeog_results(res=resDECxj_t12_t21_m2, num_trait_states=2)
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=slashslash(paste(labpt2a, "geog.data", sep="/")))

#######################################################
# Plot the geographic ancestral states
#######################################################
analysis_titletxt = "Southern conifers: Geog reconstruction under DEC+J+x+trait"

scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
results_object = geog_res
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)


#######################################################
# Plot the trait ancestral states
#######################################################
analysis_titletxt = "Southern conifers: Trait reconstruction under DEC+J+x+trait"
results_object = trait_res
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("t12","t21"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=FALSE, tr=tr, tipranges=trait_values)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=FALSE, tr=tr, tipranges=trait_values)


dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)











#######################################################
# Extract maximized log-likelihoods
#######################################################

# t12=t21 - Binary trait only (1=Fleshy cone, 2=non-fleshy cone), 1-rate model
resTrait_1rate$total_loglikelihood
# t12+t21 -- Binary trait only, 2-rate model
resTrait_2rates$total_loglikelihood
# DEC+x (distance-dependent dispersal) model on geography only
resDECx$total_loglikelihood
# DEC+x+j (distance-dependent, with jump dispersal) model on geography only
resDECxj$total_loglikelihood
# Trait+geography joint model, DEC+x+t12+t21+m2
resDECx_t12_t21_m2$total_loglikelihood
# Trait+geography joint model, DEC+x+j+t12+t21+m2
resDECxj_t12_t21_m2$total_loglikelihood

# 2022-2023 numbers:
-25.67922
-25.08393
-299.1735
-288.5104
-323.5684
-308.0184

# 2024 numbers:
-33.0481
-33.04532
-299.1735
-288.5104
-331.4065
-314.308


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
print_param_ests(resDECx, params_to_get=params_to_get), 
print_param_ests(resDECxj, params_to_get=params_to_get), 
print_param_ests(resDECx_t12_t21_m2, params_to_get=params_to_get), 
print_param_ests(resDECxj_t12_t21_m2, params_to_get=params_to_get))

param_ests_table_df = adf(param_ests_table)
row.names(param_ests_table_df) = NULL
param_ests_table_df

# 2022-2023 numbers:
#          d           e         j          x         t12         t21        m2
# 0.00000000 0.000000000 0.0000000  0.0000000 0.004071375 0.004071375 1.0000000
# 0.00000000 0.000000000 0.0000000  0.0000000 0.005386846 0.001570499 1.0000000
# 0.11476838 0.005355614 0.0000000 -0.8698849          NA          NA        NA
# 0.12347404 0.003303965 0.2059576 -0.9803483          NA          NA        NA
# 0.26682580 0.005663280 0.0000000 -1.1378666 0.005389503 0.001572181 0.9995106
# 0.07674597 0.000297265 0.6148225 -0.9965451 0.005051829 0.002708328 0.5074576

# 2024 numbers:
         d           e           j           x         t12         t21          m2 
0.000000000 0.000000000 0.000000000 0.000000000 0.005377814 0.005377814 1.000000000 
0.000000000 0.000000000 0.000000000 0.000000000 0.005249353 0.005535828 1.000000000 
0.114768383  0.005355614  0.000000000 -0.869884936           NA           NA           NA 
0.123474044  0.003303965  0.205957648 -0.980348339           NA           NA           NA 
0.287038079  0.005691870  0.000000000 -1.158972552  0.005132391  0.005497680  0.976246081 
0.0791837954  0.0001134913  0.6192036837 -0.9968830615  0.0043901973  0.0059324178  0.4732876151 


# Several independent re-runs to check optimizations
# (not needed for lab exercise)
# 
Rdata_fn = "DECx+t12+t21+m2_inf_v3b.Rdata"
rerun_optim_table_DECx = rerun_optimization_w_HiLow(res=NULL, Rdata_fn=Rdata_fn, runslow=TRUE)
# 
Rdata_fn = "DECJx+t12+t21+m2_inf_v3b.Rdata"
rerun_optim_table_DECjx = rerun_optimization_w_HiLow(res=NULL, Rdata_fn=Rdata_fn, runslow=TRUE)
# 
#Rdata_fn = "DECJx+t12+t21+m2_rep2_inf.Rdata"
#rerun_optim_table_DECj_rep2 = rerun_optimization_w_HiLow(res=NULL, #Rdata_fn=Rdata_fn, runslow=FALSE)
# 




# 2022-2023 results

# Printing 'DECx+t12+t21+m2_inf_rerun_optim_table_v3b.txt':
#                 lnL nparam         d           e
# orig_inf  -321.5659      6 0.2505165 0.005633778
# orig_redo -321.5432      6 0.2662884 0.005792895
# lowstart  -321.6107      6 0.2352627 0.005680150
# histart   -321.5240      6 0.3191225 0.005749231
#                   x         t12         t21
# orig_inf  -1.018379 0.005436032 0.001560916
# orig_redo -1.027870 0.005441661 0.001569041
# lowstart  -1.002477 0.005404404 0.001732090
# histart   -1.092792 0.005449204 0.001589304
#                  m2
# orig_inf  0.5146364
# orig_redo 0.5139577
# lowstart  0.5329633
# histart   0.4926514
# 
# Printing 'DECJx+t12+t21+m2_inf_rerun_optim_table_v3b.txt':
#                 lnL nparam          d            e
# orig_inf  -308.4762      7 0.09887329 1.338152e-03
# orig_redo -307.7471      7 0.07422314 1.000000e-13
# lowstart  -308.0513      7 0.09810287 1.073554e-03
# histart   -307.6879      7 0.08231410 1.000000e-13
#                   x         j         t12
# orig_inf  -1.083895 0.5657454 0.005655982
# orig_redo -1.065512 0.7005708 0.005332781
# lowstart  -1.060149 0.5789154 0.005592259
# histart   -1.078181 0.7197148 0.005333396
#                   t21        m2
# orig_inf  0.002020766 0.8011699
# orig_redo 0.001836048 0.6797034
# lowstart  0.001699287 0.6733107
# histart   0.002006096 0.6247004


# 2024-09-17 results

orig_inf  -314.3080      7 0.07918380 0.0001134913 -0.9968831 0.6192037 0.004390197
orig_redo -314.2975      7 0.08410638 0.0002665114 -1.0027131 0.6228974 0.004376951
lowstart  -316.5862      7 0.04761638 0.0001158004 -0.7665762 0.4350984 0.004644678
histart   -314.4221      7 0.13800569 0.0006551540 -1.1667794 0.8131723 0.004416693
                  t21        m2
orig_inf  0.005932418 0.4732876
orig_redo 0.005997088 0.4700550
lowstart  0.005751463 0.3623133
histart   0.005820861 0.5111106


cat("\n\n")
cat("...end of script, printing any warnings() to screen.")
print(warnings)
cat("\n")

