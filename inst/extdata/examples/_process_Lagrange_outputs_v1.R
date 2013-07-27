#######################################################
# Process example LAGRANGE runs
#######################################################

library(BioGeoBEARS)
source('/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_readwrite_v1.R', chdir = TRUE)

wd = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/"
setwd()

fns = c(
"Psychotria_M0/LGpy_2012/psychotria_M0.results.txt",
"Psychotria_M0/LGpy_2013/psychotria_M0.results.txt",
"Psychotria_M0/LGcpp/Psychotria_M0_lgcpp_out.txt",
"Psychotria_M0strat/LGpy_2012/psychotria_M0_strat.lagrange.results.txt",
"Psychotria_M0strat/LGpy_2013/psychotria_M0_strat.lagrange.results.txt",
"Psychotria_M0strat/LGcpp/Psychotria_M0strat_lgcpp_out.txt",
"Psychotria_M1/LGpy_2012/psychotria_M1.results.txt",
"Psychotria_M1/LGpy_2013/psychotria_M1.results.txt",
"Psychotria_M1/LGcpp/Psychotria_M1_lgcpp_out.txt",
"Psychotria_M1strat/LGpy_2012/psychotria_M1_strat.lagrange.results.txt",
"Psychotria_M1strat/LGpy_2013/psychotria_M1_strat.lagrange.results.txt",
"Psychotria_M1strat/LGcpp/Psychotria_M1strat_lgcpp_out.txt",
"Psychotria_M2/LGpy_2012/psychotria_M2.results.txt",
"Psychotria_M2/LGpy_2013/psychotria_M2.results.txt",
"Psychotria_M2/LGcpp/Psychotria_M2_lgcpp_out.txt",
"Psychotria_M2strat/LGpy_2012/psychotria_M2_strat.lagrange.results.txt",
"Psychotria_M2strat/LGpy_2013/psychotria_M2_strat.lagrange.results.txt",
"Psychotria_M2strat/LGcpp/Psychotria_M2strat_lgcpp_out.txt",
"Psychotria_M3strat/LGpy_2012/psychotria_M3_strat.lagrange.results.txt",
"Psychotria_M3strat/LGpy_2013/psychotria_M3_strat.lagrange.results.txt",
"Psychotria_M3strat/LGcpp/Psychotria_M3strat_lgcpp_out.txt",
"Psychotria_M3strat/LGcpp_v2/Psychotria_M3strat_lgcpp_out.txt",
"Lonicera_M0/LGpy_2012/Lonicera_M0.results.txt",
"Lonicera_M0/LGpy_2013/Lonicera_M0.results.txt",
"Lonicera_M0/LGcpp/Lonicera_M0_lgcpp_out.txt",
"Cyrtandra_M0/LGpy_2012/Cyrtandra_M0.results.txt",
"Cyrtandra_M0/LGpy_2013/Cyrtandra_M0.results.txt",
"Cyrtandra_M0/LGcpp/Cyrtandra_M0_lgcpp_out.txt"
)


# Go through the files and extract the summary statistics for each run.
# This is useful to see that BioGeoBEARS ML searches produce the 
# same log-likelihoods (LnLs) and parameters inferences as LAGRANGE.
# (The exception seems to be that C++ LAGRANGE has a different LnL and
#  slightly different parameter inferences in stratified analysis WITH
#  a varying dispersal probability matrix, i.e. dispersal mask, i.e. "Dmask"
#  in the code. It appears that perhaps the Dmask is being given a dual use in the
#  code -- dispersal probability multiplier, and ranges allowed, and the latter is
#  used to zero out the likelihoods of prohibited ancestral states.  However, the 
#  rate matrix should do this automatically, if EXPOKIT matrix exponentiation is being run in 
#  the correct direction; if it's not, you may have to manually zero out impossible ancestral
#  states.  This code zeroing out impossible ancestral states is blanked out in the C++
#  code for dense matrix exponentiation.  Perhaps (?) this allows more ancestral states,
#  thus higher log-likelihoods.
sumstats = NULL
for (i in 1:length(fns))
	{
	tmpfn = np(paste(addslash(wd), fns[i], sep=""))
	if (grepl(pattern="LGpy", x=tmpfn) == TRUE)
		{
		tmprow = parse_lagrange_python_output(outfn=tmpfn, outputfiles=FALSE, new_splits_fn=FALSE, new_states_fn=FALSE)
		print(sumstats)
		}
	if (grepl(pattern="LGcpp", x=tmpfn) == TRUE)
		{
		tmprow = parse_lagrange_output(outfn=tmpfn, outputfiles=FALSE, new_splits_fn=FALSE, new_states_fn=FALSE)
		print(sumstats)
		}
	path = get_path_first(inpath=fns[i])
	filename = get_path_last(path=fns[i])
	program = get_path_last(path=path)
	model = get_path_first(inpath=path)
	tmprow = cbind(model, program, path, filename, tmprow)
	sumstats = rbind(sumstats, tmprow)
	}

sumstats = adf2(sumstats)
sumstats



#######################################################
# Write the sumstats file to a table
#######################################################
sumstats_fn = paste(wd, "LAGRANGE_example_param_LnL_results.txt", sep="")
write.table(x=sumstats, file=sumstats_fn, append=FALSE, quote=FALSE, sep="	", row.names=FALSE, col.names=TRUE)

# preview the file
cat(sumstats_fn)
moref(sumstats_fn)

# read the file
sumstats_fn = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples/LAGRANGE_example_param_LnL_results.txt"
LAGRANGE_example_sumstats = read.table(sumstats_fn, header=TRUE, stringsAsFactors=FALSE)
LAGRANGE_example_sumstats
class(LAGRANGE_example_sumstats)


