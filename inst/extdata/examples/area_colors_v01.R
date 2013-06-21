
#######################################################
# Source code
#######################################################
library(LaplacesDemon)
library(BioGeoBEARS)
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_basics_v1.R', chdir = TRUE)
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_generics_v1.R', chdir = TRUE)
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_classes_v1.R', chdir = TRUE)
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_univ_model_v1.R', chdir = TRUE)
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_readwrite_v1.R', chdir = TRUE)
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_models_v1.R', chdir = TRUE)
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/calc_loglike_sp_v01.R', chdir = TRUE)
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_LaplacesDemon_v1.R', chdir = TRUE)
#source('/Dropbox/_njm/_biogeog_sim_utils_v1.R', chdir = TRUE)


#######################################################
# File locations
#######################################################
# In package
extdata_dir = system.file("extdata", package="BioGeoBEARS")
# development
extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"


#######################################################
# Input files (tree and tipdata)
#######################################################
trfn = paste(addslash(extdata_dir), "Psychotria_5.2.newick", sep="")
tr = read.tree(trfn)
geogfn = paste(addslash(extdata_dir), "Psychotria_geog.data", sep="")


# Get ranges from geogfn
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges






#######################################################
# Do colors
#######################################################
# library(grDevices)
col2rgb(col="red", alpha=TRUE)	# alpha = opacity
col2rgb(col="red", alpha=FALSE)

# Get colors, from blue to red
numareas = 3



# Get the default areas
areas = getareas_from_tipranges_object(tipranges)
areas

# Name the areas (Kauai, Oahu, Maui-Nui, Hawaii)
#areanames = c("K","O","M","H")
areanames = areas

# Make states and state indices
states = areas_list_to_states_list_old(areas=areas)
states
states2 = rcpp_areas_list_to_states_list(areas=areas)
states2

colors_matrix = get_colors_for_numareas(numareas=length(areas), use_rainbow=FALSE)
colors_matrix

states_list_0based_index = rcpp_areas_list_to_states_list(areas=areas, include_null_range=FALSE)


colors_list_for_states = mix_colors_for_states(colors_matrix, states_list_0based_index, exclude_null=TRUE)
colors_list_for_states

possible_ranges_list_txt = states_list_indexes_to_areastxt(states_list=states_list_0based_index, areanames, counting_base=0, sep="")
possible_ranges_list_txt


# Make a legend:
#plot(1:4,1:4) 
#lines(1:4,4:1, col="blue") 
#legend("top",leg=c("a","b"),col=c("black","blue"), fill=TRUE) 


# Plot, no borders (bty="n"), no labels (xlab, ylab), no tick marks (xaxt, yaxt)
plot(1:10, 1:10, pch=".", col="white", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
#lines(1:4,4:1, col="blue") 
#legend("top", leg=c("a","b"),col=c("black","blue"), fill=TRUE) 
legend("top", leg=possible_ranges_list_txt, fill=colors_list_for_states, ncol=2, title="Legend", cex=2.5)#, fill=TRUE) 


