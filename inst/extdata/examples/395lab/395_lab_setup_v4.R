

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
install.packages("FD")
install.packages("snow")
install.packages("rexpokit")
install.packages("cladoRcpp")

install.packages("phytools")
install.packages("phangorn")
install.packages("phylobase")
install.packages("optimx")
install.packages("GenSA")

# Install additional dependencies
install.packages(c("plotrix","gdata","minqa","fdrtool","statmod","SparseM","spam","MultinomialCI"))

################################################
# Small numbers of students:
# Install BioGeoBEARS from GitHub
################################################
# (BioGeoBEARS is pure R, so installation is easy *if* the above 
#  packages have been installed)
library(devtools)
install_github(repo="nmatzke/BioGeoBEARS", dependencies=TRUE, upgrade="never")

################################################
# Large numbers of students (eg a university class)
# To avoid GitHub overload, download BioGeoBEARS from Canvas, then install locally:
################################################
# Download from:
# Canvas -> BIOSCI 395 -> Files -> Southern_Conifer_Biogeog
# https://canvas.auckland.ac.nz/courses/106014/files/folder/Southern_conifer_biogeog
# BioGeoBEARS_1.1.3.tar.gz

# Install BioGeoBEARS from Canvas-downloaded zipfile
install.packages("BioGeoBEARS_1.1.3.tar.gz", repos=NULL, type="source", dependencies=TRUE)

' # END run_text_inside_quote_once

# Check that your BioGeoBEARS installation loads
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

labpt1a_script = paste(extdata_dir, "examples/395lab/Psychotria_M0_equalRates/Psychotria_M0_v1.R", sep="/")
labpt1b_script = paste(extdata_dir, "examples/395lab/Psychotria_M2_oneWayDispersal/Psychotria_M2_oneWayDispersal_v1.R", sep="/")
labpt1c_script = paste(extdata_dir, "examples/395lab/Psychotria_M4_DistanceDispersal/Psychotria_M4_DistanceDispersal_v1.R/", sep="/")
labpt2a_script = paste(extdata_dir, "examples/395lab/conifer_DEC_traits_models/conifer_DEC_traits_models_v3.R", sep="/")
labpt2b_script = paste(extdata_dir, "examples/395lab/conifer_DEC+x_traits_models/conifer_DEC+x_traits_models_v3.R", sep="/")


# Open an R script in R Studio
file.edit(labpt1a_script)
# ONCE LINE 70 WORKS, TAKE A BREAK, WE WILL CONTINUE ONCE EVERYONE IS HERE

