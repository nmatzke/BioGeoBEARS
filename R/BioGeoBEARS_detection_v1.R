#######################################################
# DETECTION PROBABILITIES
#######################################################

# source("/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_detection_v1.R")

require("ape")
require("rexpokit")
require("cladoRcpp")


# Qmat contains the d & e probs
# COO_probs_columnar contains the simulation probs


#######################################################
# read_detections
#######################################################
#' Read a file with detection counts per area
#' 
#' This function reads in a tab-delimited text file containing counts of detections of 
#' each OTU in each region.  These could be from database records or some other source.
#' 
#' @param detects_fn The filename of the detections file.
#' @param OTUnames Default \code{NULL}, in which case the first column of the text file is used as row names/OTU names.
#' @param areanames Default \code{NULL}, in which case the text file column headings are used.
#' @param tmpskip How many lines should be skipped before reading the text file?  Default \code{0}.
#' @param phy An ape phylo object. If included, the rows will be sorted to match the order of tree tip labels.
#' @return \code{dtf} 
#' @export
#' @seealso \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Bottjer_Jablonski_1988
#' @examples
#' testval=1
#' 
read_detections <- function(detects_fn, OTUnames=NULL, areanames=NULL, tmpskip=0, phy=NULL)
	{
	runjunk='
	detects_fn = "/Library/Frameworks/R.framework/Versions/2.15/Resources/library/BioGeoBEARS/extdata/Psychotria_detections_v1.txt"
	detects_fn = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/Psychotria_detections_v1.txt"
	OTUnames=NULL
	areanames=NULL
	tmpskip=0
	'
	if (is.null(OTUnames) == TRUE)
		{
		dtf = read.table(file=detects_fn, header=TRUE, skip=tmpskip, sep="	", quote="", stringsAsFactors = FALSE, strip.white=TRUE, fill=TRUE, row.names=1)
		dtf
		} else {
		# Assumes no OTU names
		dtf = read.table(file=detects_fn, header=TRUE, skip=tmpskip, sep="	", quote="", stringsAsFactors = FALSE, strip.white=TRUE, fill=TRUE, row.names=OTUnames)
		dtf
		}

	if (is.null(areanames) == TRUE)
		{
		pass=1
		} else {
		names(dtf) = areanames
		}

	# Order tip labels
	if (!is.null(phy))
		{
		tipranges_df_order = match(phy$tip.label, rownames(dtf))
		dtf = dtf[tipranges_df_order, ]
		}
		
	# Standardize and return
	dtf = adf2(dtf)
	dtf = dfnums_to_numeric(dtf)
	
	return(dtf)
	}


#######################################################
# read_controls
#######################################################
#' Read a file with the total number of detections in a taphonomic control
#' 
#' This function reads in a tab-delimited text file containing counts of detections of 
#' the taphonomic controls in each region.  These numbers should always be equal to or
#' larger than the counts in the detections file.
#' 
#' The idea of taphonomic controls dates back at least to work of Bottjer & Jablonski (1988).  The basic
#' idea is that if you have taxa of roughly similar detectability, then detections of other 
#' taxa give some idea of overall detection effort.  Obviously this is a very simple 
#' model that can be criticized in any number of ways (different alpha diversity in each 
#' region, different detectability of individual taxa, etc.), but it is a useful starting
#' point as there has been no implementation of any detection model in historical/phylogenetic 
#' biogeography to date.
#' 
#' One could imagine (a) every OTU and area has a different count of detections and taphonomic control 
#' detections, or (b) the taphonomic control detections are specified by area, and shared across all OTUs.
#' Situation (b) is likely more common, but this function implements (a).  Behavior (b) 
#' could be reproduced by summing each column, and/or copying this sum to all cells for a 
#' particular area.
#' 
#' @param controls_fn The filename of the file containing the counts of taphonomic control detections.
#' @param OTUnames Default \code{NULL}, in which case the first column of the text file is used as row names/OTU names.
#' @param areanames Default \code{NULL}, in which case the text file column headings are used.
#' @param tmpskip How many lines should be skipped before reading the text file?  Default \code{0}.
#' @param phy An ape phylo object. If included, the rows will be sorted to match the order of tree tip labels.
#' @return \code{dtf} 
#' @export
#' @seealso \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Bottjer_Jablonski_1988
#' @examples
#' testval=1
#' 
read_controls <- function(controls_fn, OTUnames=NULL, areanames=NULL, tmpskip=0, phy=NULL)
	{
	runjunk='
	controls_fn = "/Library/Frameworks/R.framework/Versions/2.15/Resources/library/BioGeoBEARS/extdata/Psychotria_controls_v1.txt"
	controls_fn = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/Psychotria_controls_v1.txt"
	OTUnames=NULL
	areanames=NULL
	tmpskip=0
	'
	if (is.null(OTUnames) == TRUE)
		{
		dtf = read.table(file=controls_fn, header=TRUE, skip=tmpskip, sep="	", quote="", stringsAsFactors = FALSE, strip.white=TRUE, fill=TRUE, row.names=1)
		dtf
		} else {
		# Assumes no OTU names
		dtf = read.table(file=controls_fn, header=TRUE, skip=tmpskip, sep="	", quote="", stringsAsFactors = FALSE, strip.white=TRUE, fill=TRUE, row.names=OTUnames)
		dtf
		}

	if (is.null(areanames) == TRUE)
		{
		pass=1
		} else {
		names(dtf) = areanames
		}
	
	# Order tip labels
	if (!is.null(phy))
		{
		tipranges_df_order = match(phy$tip.label, rownames(dtf))
		dtf = dtf[tipranges_df_order, ]
		}
	
	# Standardize and return
	dtf = adf2(dtf)
	dtf = dfnums_to_numeric(dtf)
	
	return(dtf)
	}






# Likelihood model
# Calculate P(obs | presence)
# Source: 
# "/Users/nickm/Desktop/__projects/___phylokriging_example/__scripts/_ocean_data_predict_v6.R', chdir = TRUE)

#######################################################
# calc_obs_like
#######################################################
#' Calculate likelihood of count data given true presence/absence and parameters
#' 
#' This function calculates P(data|presence,parameters), i.e. the probability of some detection and 
#' taphonomic control counts, given the true geographic range/state, and parameters such as \code{dp}, a 
#' detection probability (and, optionally, a false detection probability, \code{fdp}).
#' 
#' The idea of taphonomic controls dates back at least to work of Bottjer & Jablonski (1988).  The basic
#' idea is that if you have taxa of roughly similar detectability, then detections of other 
#' taxa give some idea of overall detection effort.  Obviously this is a very simple 
#' model that can be criticized in any number of ways (different alpha diversity in each 
#' region, different detectability of individual taxa, etc.), but it is a useful starting
#' point as there has been no implementation of any detection model in historical/phylogenetic 
#' biogeography to date.
#' 
#' One could imagine (a) every OTU and area has a different count of detections and taphonomic control 
#' detections, or (b) the taphonomic control detections are specified by area, and shared across all OTUs.
#' Situation (b) is likely more common, but this function assumes (a) as this is the more thorough case.
#' Behavior (b) could be reproduced by summing each column, and/or copying this sum to all cells for a 
#' particular area.
#' 
#' @param truly_present Is the OTU of interest known/conditionally assumed to be truly present (\code{TRUE}) or truly absent (\code{FALSE})?
#' @param obs_target_species A count of detections of your OTU of interest, e.g. as produced from a cell of the matrix output from \code{\link{read_detections}}.
#' @param obs_all_species A count of detections of your taphonomic controls, e.g. as produced from a cell of the output from \code{\link{read_controls}}.
#' @param mean_frequency This is the proportion of samples from the taphonomic control group that will truly be from this OTU, GIVEN that the OTU is present.
#' This could be estimated, but a decent first guess is (total # samples of OTU of interest / total # of samples in the taphonomic control group where
#' the OTU is known to be present).  All that is really needed is some reasonable value, such that more sampling without detection lowers the 
#' likelihood of the data on the hypothesis of true presence, and vice versa.  This value can only be 1 when the number of detections = the number 
#' of taphonomic control detections, for every OTU and area.  This is the implicit assumption in e.g. standard historical biogeography analyses in 
#' LAGRANGE or BioGeoBEARS.
#' @param dp The detection probability.  This is the per-sample probability that you will correctly detect the OTU in question, 
#' when you are looking at it.  Default is 1, which is the implicit assumption in standard analyses.
#' @param fdp The false detection probability.  This is probability of falsely concluding a detection of the OTU of interest occurred, when in 
#' fact the specimen was of something else.  The default is 0, which assumes zero error rate, 
#' i.e. the assumption being made in all historical biogeography analyses that do not take into account detection 
#' probability.  These options are being included for completeness, but it may not be wise to try to infer \code{mean_frequency},  
#' \code{dp} and \code{fdp} all at once due to identifiability issues (and estimation of fdp may take a very large amount of data).  However, 
#' fixing some of these parameters to reasonable values can allow the user to effectively include beliefs about the uncertainty of the 
#' input data into the analysis, if desired.
#' @return \code{lnlike_allobs_given_absence} The natural log-likelihood of the data, given the model & assumption of true presence or absence.
#' @export
#' @seealso \code{\link{mapply_calc_post_prob_presence}}, \code{\link{calc_post_prob_presence}}, \code{\link{mapply_calc_obs_like}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Bottjer_Jablonski_1988
#' @examples
#' 
#' # Example: 10 observations of the species mean dramatically higher likelihood of the
#' # data on the hypothesis that it is truly present.
#' 
#' # With zero error rate
#' obs_target_species = 10
#' obs_all_species = 100
#' mean_frequency=0.1
#' dp=1
#' fdp=0
#' LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_presence
#' LnL_under_absence
#' # Note that the probability of getting detections, under the hypothesis of 
#' # true absence, is -Inf
#' 
#' 
#' # With a small error rate, there is some small but positive probability of
#' # falsely getting 10 detections
#' obs_target_species = 10
#' obs_all_species = 100
#' mean_frequency=0.1
#' dp=0.99
#' fdp=0.001
#' LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_presence
#' LnL_under_absence
#' # i.e. the prob. of the data is 1 under the hypothesis of presence, and 0 
#' # under the hypothesis of absence (ln(prob) = 0 & -Inf, respectively)
#' 
#' 
#' # Note that with very high error rates, your conclusion could reverse
#' obs_target_species = 10
#' obs_all_species = 100
#' mean_frequency=0.1
#' dp=0.5
#' fdp=0.3
#' LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_presence
#' LnL_under_absence
#' 
#' 
#' # Example #2 -- what if you have ZERO detections, but lots of detections 
#' # of your taphonomic control?
#' obs_target_species = 0
#' obs_all_species = 1
#' mean_frequency=0.1
#' dp=1
#' fdp=0
#' LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_presence
#' LnL_under_absence
#' 
#' # With a slight error rate
#' obs_target_species = 0
#' obs_all_species = 1
#' mean_frequency=0.1
#' dp=0.99
#' fdp=0.001
#' LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_presence
#' LnL_under_absence
#' 
#' 
#' obs_target_species = 0
#' obs_all_species = 2
#' mean_frequency=0.1
#' dp=1
#' fdp=0
#' LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_presence
#' LnL_under_absence
#' 
#' # With a slight error rate
#' obs_target_species = 0
#' obs_all_species = 2
#' mean_frequency=0.1
#' dp=0.99
#' fdp=0.001
#' LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_presence
#' LnL_under_absence
#' 
#' 
#' 
#' 
#' 
#' # Example #3 -- what if you have ZERO detections, but only a few 
#' # detections of your taphonomic control?
#' obs_target_species = 0
#' obs_all_species = 100
#' mean_frequency=0.1
#' dp=1
#' fdp=0
#' LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_presence
#' LnL_under_absence
#' 
#' # With a slight error rate
#' obs_target_species = 0
#' obs_all_species = 100
#' mean_frequency=0.1
#' dp=0.99
#' fdp=0.001
#' LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_presence
#' LnL_under_absence
#' 
#' 
#' 
#' # Special cases -- e.g., no data
#' # Prob(data)=1, ln(prob)=0
#' obs_target_species = 0
#' obs_all_species = 0
#' mean_frequency=0.1
#' dp=0.99
#' fdp=0.001
#' LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_presence
#' LnL_under_absence
#' 
#' obs_target_species = 0
#' obs_all_species = 0
#' mean_frequency=0.1
#' dp=1
#' fdp=0
#' LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_presence
#' LnL_under_absence
#' 
#' 
#' # What if, for some reason, you put in identical detections and taphonomic control
#' # counts? (e.g., you load in a standard tipranges file)
#' obs_target_species = 1
#' obs_all_species = 1
#' mean_frequency=1
#' dp=1
#' fdp=0
#' LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_presence
#' LnL_under_absence
#' 
#' # What if, for some reason, you put in identical detections and taphonomic control
#' # counts? (e.g., you load in a standard tipranges file)
#' obs_target_species = 1
#' obs_all_species = 1
#' mean_frequency=1
#' dp=0.99
#' fdp=0.001
#' LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' LnL_under_presence
#' LnL_under_absence
#' 
#' 
calc_obs_like <- function(truly_present=TRUE, obs_target_species, obs_all_species, mean_frequency=0.1, dp=1, fdp=0)
	{
	# Error check -- if you have NO taphonomic controls, you have no information to calculate the likelihood.
	# Therefore just return 1 for any hypothesis (implicit flat prior)
	if (obs_all_species == 0)
		{
		# Probability of these "data" is 1 under any hypothesis; ln(prob=1) is 0
		return(0)
		}
	
	# Hard-coded model parameters
	#########################################
	# Sampling probability
	# Given a single sample, what is the probability that it is the target species?
	Psample_is_the_species_given_presence = mean_frequency
	Psample_is_NOT_species_given_presence = 1-mean_frequency
	Psample_is_the_species_given_absence = 0
	Psample_is_NOT_species_given_absence = 1-Psample_is_the_species_given_absence

	# Identification probability
	# Given an ID, what is the probability that it is accurate?
	P_ID_is_accurate_given_species_is_target = dp
	P_ID_is_NOT_accurate_given_species_is_target = 1-P_ID_is_accurate_given_species_is_target

	P_ID_is_accurate_given_NONtarget_species = 1-fdp
	P_ID_is_NOT_accurate_given_NONtarget_species = 1-P_ID_is_accurate_given_NONtarget_species	# false detection prob
	#########################################
	
	
	# Get observations
	
	# count of the observations of the target species
	obs_target_species = obs_target_species
	
	# count of the observations of the other species (non-target)
	obs_NONtarget_species = obs_all_species-obs_target_species
	
	
	# Calculate P(obs | true_presence)
	if (truly_present == TRUE)
		{
		# prob of presence data (given presence in pixel)
		# prob of ways to observe species correctly (true positive)
		p1 = Psample_is_the_species_given_presence * P_ID_is_accurate_given_species_is_target
		# prob of ways to observe species incorrectly (false positive)
		p2 = Psample_is_NOT_species_given_presence * P_ID_is_NOT_accurate_given_NONtarget_species
		
		lnlike_presencedata_given_presence = log(p1 + p2)*obs_target_species
		exp(lnlike_presencedata_given_presence)
	
	
		# prob of absence data (given presence in pixel)
		# prob of ways to observe NOT species correctly (true negative)
		# still given presence
		p1 = Psample_is_NOT_species_given_presence * P_ID_is_accurate_given_NONtarget_species
		# prob of ways to observe NOT species incorrectly (false negative)
		p2 = Psample_is_the_species_given_presence * P_ID_is_NOT_accurate_given_species_is_target
		
		lnlike_absensedata_given_presence = log(p1 + p2)*obs_NONtarget_species
		# -Inf * 0 trap; prob is 1, log prob is 0
		if ( (log(p1 + p2) == -Inf) && obs_NONtarget_species == 0)
			{
			lnlike_absensedata_given_presence = 0
			}
		exp(lnlike_absensedata_given_presence)
		
		lnlike_allobs_given_presence = lnlike_presencedata_given_presence + lnlike_absensedata_given_presence
		exp(lnlike_allobs_given_presence)
		
		return(lnlike_allobs_given_presence)
		}
		
		
	# Calculate P(obs | true_absence)
	if (truly_present == FALSE)
		{
		# prob of ways to observe species correctly (true result, which can't be positive)
		p1 = Psample_is_the_species_given_absence * P_ID_is_accurate_given_species_is_target
		# prob of ways to observe species incorrectly (false result, which might positive)
		p2 = Psample_is_NOT_species_given_absence * P_ID_is_NOT_accurate_given_NONtarget_species
		
		lnlike_presencedata_given_absence = log(p1 + p2)*obs_target_species
		# -Inf * 0 trap; prob is 1, log prob is 0
		if ( (log(p1 + p2) == -Inf) && obs_target_species == 0)
			{
			lnlike_presencedata_given_absence = 0
			}
		exp(lnlike_presencedata_given_absence)
		
	
		# prob of ways to observe NOT species correctly (true negative)
		# still given absence
		p1 = Psample_is_NOT_species_given_absence * P_ID_is_accurate_given_NONtarget_species
		# prob of ways to observe NOT species incorrectly (false negative)
		p2 = Psample_is_the_species_given_absence * P_ID_is_NOT_accurate_given_species_is_target
		
		lnlike_absensedata_given_absence = log(p1 + p2)*obs_NONtarget_species
		exp(lnlike_absensedata_given_absence)

		
		lnlike_allobs_given_absence = lnlike_presencedata_given_absence + lnlike_absensedata_given_absence
		exp(lnlike_allobs_given_absence)
		
		return(lnlike_allobs_given_absence)
		}
	}




#######################################################
# calc_post_prob_presence
#######################################################
#' Calculate posterior probability of presence, given count data and parameters
#' 
#' This function calculates P(presence|count data,parameters), i.e. the posterior 
#' probability of presence in an area, given data on detection counts and 
#' taphonomic control counts, and a detection model with the parameters mean_frequency, \code{dp}, a 
#' detection probability (and, optionally, a false detection probability, \code{fdp}).
#' 
#' Essentially, this function combines a prior probability, with the likelihood function 
#' (coded in \code{\link{calc_obs_like}}) to produce a posterior probability of presence
#' given Bayes' Theorem (Bayes & Price, 1763).
#' 
#' The idea of taphonomic controls dates back at least to work of Bottjer & Jablonski (1988).  The basic
#' idea is that if you have taxa of roughly similar detectability, then detections of other 
#' taxa give some idea of overall detection effort.  Obviously this is a very simple 
#' model that can be criticized in any number of ways (different alpha diversity in each 
#' region, different detectability of individual taxa, etc.), but it is a useful starting
#' point as there has been no implementation of any detection model in historical/phylogenetic 
#' biogeography to date.
#' 
#' One could imagine (a) every OTU and area has a different count of detections and taphonomic control 
#' detections, or (b) the taphonomic control detections are specified by area, and shared across all OTUs.
#' Situation (b) is likely more common, but this function assumes (a) as this is the more thorough case.
#' Behavior (b) could be reproduced by summing each column, and/or copying this sum to all cells for a 
#' particular area.
#' 
#' @param prior_prob_presence The prior probability of presence, i.e. when no detection or taphonomic control data 
#' whatsoever is available.  Default is set to 0.01 which expresses my totally uninformed bias that 
#' in whatever your data is, your species of interest probably doesn't live in the typical area you are 
#' looking at.
#' @param obs_target_species A count of detections of your OTU of interest, e.g. as produced from a cell of the matrix output from \code{\link{read_detections}}.
#' @param obs_all_species A count of detections of your taphonomic controls, e.g. as produced from a cell of the output from \code{\link{read_controls}}.
#' @param mean_frequency This is the proportion of samples from the taphonomic control group that will truly be from this OTU, GIVEN that the OTU is present.
#' This could be estimated, but a decent first guess is (total # samples of OTU of interest / total # of samples in the taphonomic control group where
#' the OTU is known to be present).  All that is really needed is some reasonable value, such that more sampling without detection lowers the 
#' likelihood of the data on the hypothesis of true presence, and vice versa.  This value can only be 1 when the number of detections = the number 
#' of taphonomic control detections, for every OTU and area.  This is the implicit assumption in e.g. standard historical biogeography analyses in 
#' LAGRANGE or BioGeoBEARS.
#' @param dp The detection probability.  This is the per-sample probability that you will correctly detect the OTU in question, 
#' when you are looking at it.  Default is 1, which is the implicit assumption in standard analyses.
#' @param fdp The false detection probability.  This is probability of falsely concluding a detection of the OTU of interest occurred, when in 
#' fact the specimen was of something else.  The default is 0, which assumes zero error rate, 
#' i.e. the assumption being made in all historical biogeography analyses that do not take into account detection 
#' probability.  These options are being included for completeness, but it may not be wise to try to infer \code{mean_frequency},  
#' \code{dp} and \code{fdp} all at once due to identifiability issues (and estimation of fdp may take a very large amount of data).  However, 
#' fixing some of these parameters to reasonable values can allow the user to effectively include beliefs about the uncertainty of the 
#' input data into the analysis, if desired.
#' @param print_progress If not the default (\code{""}), print whatever is in print_progress, followed by a space (for error checking/surveying results).
#' @return \code{post_prob} The posterior probability of presence, given the prior probability, the model parameters, and the data.
#' @export
#' @seealso \code{\link{calc_obs_like}},  \code{\link{mapply_calc_post_prob_presence}}, \code{\link{mapply_calc_obs_like}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://en.wikipedia.org/wiki/Bayes'_theorem}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Bottjer_Jablonski_1988
#' 	 @cite Bayes_1763
#' @examples
#'
#' # Calculate posterior probability of presence in an area, 
#' # given a dp (detection probability) and detection model. 
#' 
#' # With zero error rate
#' obs_target_species = 10
#' obs_all_species = 100
#' mean_frequency=0.1
#' dp=1
#' fdp=0
#' prior_prob_presence = 0.01
#' post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' post_prob
#' # i.e., with perfect detection, the prob. of presence is 1 under the 
#' # hypothesis of presence, and 0 under the hypothesis of 
#' # (This is because the likelihood of the data under 
#' # presence and absence are ln(prob) = 0 & -Inf, respectively.)
#' 
#' 
#' # Note that with very high error rates, your conclusion could reverse
#' obs_target_species = 10
#' obs_all_species = 100
#' mean_frequency=0.1
#' dp=0.5
#' fdp=0.3
#' prior_prob_presence = 0.01
#' post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' post_prob
#' 
#' 
#' # With 0 error rate, even 1 observation makes P(presence) = 1
#' obs_target_species = 1
#' obs_all_species = 100
#' mean_frequency=0.1
#' dp=1
#' fdp=0
#' prior_prob_presence = 0.01
#' post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' post_prob
#' 
#' 
#' # With a small error rate, there is some small but positive probability of
#' # falsely getting 10 detections; but it may be effectively 0
#' obs_target_species = 10
#' obs_all_species = 100
#' mean_frequency=0.1
#' dp=0.99
#' fdp=0.001
#' prior_prob_presence = 0.01
#' post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' post_prob
#' 
#' 
#' # If you have only 1 detection, and you have 100 taphonomic controls and
#' # a mean_frequency of sampling the OTU of interest of 0.1, then there is 
#' # still a very low probability of presence (since, under your model, 
#' # you should expect to see about 10 detections, not 1)
#' obs_all_species = 100
#' mean_frequency=0.1
#' dp=0.99
#' fdp=0.001
#' prior_prob_presence = 0.01
#' 
#' obs_target_species = 0
#' calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' 
#' obs_target_species = 1
#' calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' 
#' obs_target_species = 2
#' calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' 
#' obs_target_species = 3
#' calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' 
#' 
#' 
#' # Note how quickly this chances if you drop the mean_frequency from 0.1 
#' # to 0.01. This means that if you want single detections to count for 
#' # a lot, you need either a low mean_frequency which matches the observed 
#' # frequency, or an extremely high/perfect detection probability (dp).
#' obs_all_species = 100
#' mean_frequency=0.01
#' dp=0.99
#' fdp=0.001
#' prior_prob_presence = 0.01
#' 
#' obs_target_species = 0
#' calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' 
#' obs_target_species = 1
#' calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' 
#' obs_target_species = 2
#' calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' 
#' obs_target_species = 3
#' calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' 
#' 
#' # Changing mean_frequency from 0.01 to 0.001 actually LOWERS the posterior
#' # probability of presence based on 1 detection, as we have a somewhat 
#' # significant false detection rate:
#' obs_target_species = 1
#' obs_all_species = 100
#' mean_frequency=0.001
#' dp=0.99
#' fdp=0.001
#' prior_prob_presence = 0.01
#' 
#' obs_target_species = 0
#' calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' 
#' obs_target_species = 1
#' calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' 
#' obs_target_species = 2
#' calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' 
#' obs_target_species = 3
#' calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' 
#' 
#' # Change false detection probability to a much lower value
#' obs_all_species = 100
#' mean_frequency=0.001
#' dp=0.99
#' fdp=0.00001
#' prior_prob_presence = 0.01
#' 
#' obs_target_species = 0
#' calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' 
#' obs_target_species = 1
#' calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' 
#' obs_target_species = 2
#' calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' 
#' obs_target_species = 3
#' calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' 
#' 
#' 
#' 
#' # Change false detection probability to 0
#' obs_all_species = 100
#' mean_frequency=0.001
#' dp=0.99
#' fdp=0.0
#' prior_prob_presence = 0.01
#' 
#' 
#' obs_target_species = 0
#' calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' 
#' obs_target_species = 1
#' calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' 
#' obs_target_species = 2
#' calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' 
#' obs_target_species = 3
#' calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' # Change mean_frequency to 0.001
#' obs_all_species = 100
#' mean_frequency=0.001
#' dp=0.99
#' fdp=0.0
#' prior_prob_presence = 0.01
#' 
#' obs_target_species = 0
#' calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' 
#' obs_target_species = 1
#' calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' 
#' obs_target_species = 2
#' calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' 
#' obs_target_species = 3
#' calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' # Example #2 -- what if you have ZERO detections, but lots of detections 
#' # of your taphonomic control?
#' obs_target_species = 0
#' obs_all_species = 100
#' mean_frequency=0.1
#' dp=1
#' fdp=0
#' prior_prob_presence = 0.01
#' post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' post_prob
#' 
#' # With a slight error rate
#' obs_target_species = 0
#' obs_all_species = 100
#' mean_frequency=0.1
#' dp=0.99
#' fdp=0.001
#' prior_prob_presence = 0.01
#' post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' post_prob
#' 
#' 
#' obs_target_species = 0
#' obs_all_species = 2
#' mean_frequency=0.1
#' dp=1
#' fdp=0
#' prior_prob_presence = 0.01
#' post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' post_prob
#' 
#' # With a slight error rate
#' obs_target_species = 0
#' obs_all_species = 2
#' mean_frequency=0.1
#' dp=0.99
#' fdp=0.001
#' prior_prob_presence = 0.01
#' post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' post_prob
#' 
#' 
#' 
#' 
#' 
#' # Example #3 -- what if you have ZERO detections, but only a few 
#' # detections of your taphonomic control?
#' obs_target_species = 0
#' obs_all_species = 1
#' mean_frequency=0.1
#' dp=1
#' fdp=0
#' prior_prob_presence = 0.01
#' post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' post_prob
#' 
#' # With a slight error rate
#' obs_target_species = 0
#' obs_all_species = 1
#' mean_frequency=0.1
#' dp=0.99
#' fdp=0.001
#' prior_prob_presence = 0.01
#' post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' post_prob
#' 
#' 
#' 
#' # Special cases -- e.g., no data
#' # Prob(data)=1, ln(prob)=0
#' obs_target_species = 0
#' obs_all_species = 0
#' mean_frequency=0.1
#' dp=0.99
#' fdp=0.001
#' prior_prob_presence = 0.01
#' post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' post_prob
#' 
#' obs_target_species = 0
#' obs_all_species = 0
#' mean_frequency=0.1
#' dp=1
#' fdp=0
#' prior_prob_presence = 0.01
#' post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' post_prob
#' 
#' 
#' # What if, for some reason, you put in identical detections and taphonomic control
#' # counts? (e.g., you load in a standard tipranges file)
#' obs_target_species = 1
#' obs_all_species = 1
#' mean_frequency=1
#' dp=1
#' fdp=0
#' prior_prob_presence = 0.01
#' post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' post_prob
#' 
#' # What if, for some reason, you put in identical detections and taphonomic control
#' # counts? (e.g., you load in a standard tipranges file)
#' obs_target_species = 1
#' obs_all_species = 1
#' mean_frequency=1
#' dp=0.99
#' fdp=0.001
#' prior_prob_presence = 0.01
#' post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
#' obs_all_species, mean_frequency, dp, fdp)
#' post_prob
#' 
#' 
calc_post_prob_presence <- function(prior_prob_presence=0.01, obs_target_species, obs_all_species, mean_frequency=0.1, dp=1, fdp=0, print_progress="")
	{
	# Priors
	prior_prob_presence = prior_prob_presence
	prior_prob_absence = 1-prior_prob_presence

	# likelihood of data (assuming truly present)
	lnlike_allobs_given_presence = calc_obs_like(truly_present=TRUE, obs_target_species, obs_all_species, mean_frequency=mean_frequency, dp=dp, fdp=fdp)
	
	numerator = lnlike_allobs_given_presence + log(prior_prob_presence)
	
	
	# prob(obs)
	# denominator
	d1 = lnlike_allobs_given_presence + log(prior_prob_presence)
	
	lnlike_allobs_given_absence = calc_obs_like(truly_present=FALSE, obs_target_species, obs_all_species, mean_frequency=mean_frequency, dp=dp, fdp=fdp)
	
	d2 = lnlike_allobs_given_absence + log(prior_prob_absence)
	
	denominator = log(exp(d1) + exp(d2))
	
	ln_post_prob = numerator - denominator
	post_prob = exp(ln_post_prob)
	
	if (print_progress == "")
		{
		nothing_happens = 0
		} else {
		cat(print_progress, " ", sep="")
		}
	
	return(post_prob)
	}










#######################################################
# mapply_calc_post_prob_presence
#######################################################
#' Mapply version of calc_post_prob_presence()
#' 
#' This function applies \code{\link{calc_post_prob_presence}} to all cells of 
#' the input matrices \code{obs_target_species} and \code{obs_all_species}. These matrices
#' obviously must have the same dimensions.
#' 
#' The inputs are the same as for \code{\link{calc_post_prob_presence}}, except that
#' \code{obs_target_species} and \code{obs_all_species} can be matrices.
#' 
#' @param prior_prob_presence The prior probability of presence, i.e. when no detection or taphonomic control data 
#' whatsoever is available.  Default is set to 0.01 which expresses my totally uninformed bias that 
#' in whatever your data is, your species of interest probably doesn't live in the typical area you are 
#' looking at.
#' @param obs_target_species A scalar or column/vector/matrix of detection counts, e.g. as produced from the output from \code{\link{read_detections}}.
#' @param obs_all_species A scalar or column/vector/matrix of detection counts, e.g. as produced from the output from \code{\link{read_controls}}.
#' @param mean_frequency This is the proportion of samples from the taphonomic control group that will truly be from this OTU, GIVEN that the OTU is present.
#' This could be estimated, but a decent first guess is (total # samples of OTU of interest / total # of samples in the taphonomic control group where
#' the OTU is known to be present).  All that is really needed is some reasonable value, such that more sampling without detection lowers the 
#' likelihood of the data on the hypothesis of true presence, and vice versa.  This value can only be 1 when the number of detections = the number 
#' of taphonomic control detections, for every OTU and area.  This is the implicit assumption in e.g. standard historical biogeography analyses in 
#' LAGRANGE or BioGeoBEARS.
#' @param dp The detection probability.  This is the per-sample probability that you will correctly detect the OTU in question, 
#' when you are looking at it.  Default is 1, which is the implicit assumption in standard analyses.
#' @param fdp The false detection probability.  This is probability of falsely concluding a detection of the OTU of interest occurred, when in 
#' fact the specimen was of something else.  The default is 0, which assumes zero error rate, 
#' i.e. the assumption being made in all historical biogeography analyses that do not take into account detection 
#' probability.  These options are being included for completeness, but it may not be wise to try to infer \code{mean_frequency},  
#' \code{dp} and \code{fdp} all at once due to identifiability issues (and estimation of fdp may take a very large amount of data).  However, 
#' fixing some of these parameters to reasonable values can allow the user to effectively include beliefs about the uncertainty of the 
#' input data into the analysis, if desired.
#' @param print_progress If not the default (\code{""}), print whatever is in print_progress, followed by a space (for error checking/surveying results).
#' @return \code{pp_df} A matrix of the posterior probability of presence, given the prior probability, the model parameters, and the data.
#' @export
#' @seealso \code{\link{calc_obs_like}}, \code{\link{calc_post_prob_presence}}, \code{\link{mapply_calc_obs_like}}
#' \code{\link{Pdata_given_rangerow}}, \code{\link[base]{mapply}}, \code{\link{tiplikes_wDetectionModel}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://en.wikipedia.org/wiki/Bayes'_theorem}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Bottjer_Jablonski_1988
#' 	 @cite Bayes_1763
#' @examples
#' 
#' # Calculate posterior probability of presence in an area, 
#' # given a dp (detection probability) and detection model. 
#' 
#' # soft-coded input files
#' extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
#' detects_fn = np(paste(extdata_dir, "/Psychotria_detections_v1.txt", sep=""))
#' controls_fn = np(paste(extdata_dir, "/Psychotria_controls_v1.txt", sep=""))
#' 
#' OTUnames=NULL
#' areanames=NULL
#' tmpskip=0
#' 
#' detects_df = read_detections(detects_fn, OTUnames=NULL, areanames=NULL, tmpskip=0)
#' controls_df = read_controls(controls_fn, OTUnames=NULL, areanames=NULL, tmpskip=0)
#' 
#' detects_df
#' controls_df
#' detects_df / controls_df
#' 
#' 
#' # Calculate data likelihoods, and posterior probability of presence=TRUE
#' mean_frequency=0.1
#' dp=1
#' fdp=0
#' 
#' mapply_calc_obs_like(truly_present=TRUE, obs_target_species=detects_df,
#' obs_all_species=controls_df, mean_frequency, dp, fdp)
#' 
#' mapply_calc_obs_like(truly_present=FALSE, obs_target_species=detects_df,
#' obs_all_species=controls_df, mean_frequency, dp, fdp)
#' 
#' mapply_calc_post_prob_presence(prior_prob_presence=0.01, 
#' obs_target_species=detects_df,
#' obs_all_species=controls_df, mean_frequency, dp, fdp)
#' 
mapply_calc_post_prob_presence <- function(prior_prob_presence=0.01, obs_target_species, obs_all_species, mean_frequency=0.1, dp=1, fdp=0, print_progress="")
	{
	# Calculate the posterior probability of presence in each area in a matrix / data.frame
	tmpdf = mapply(FUN=calc_post_prob_presence, obs_target_species=unlist(obs_target_species),
	obs_all_species=unlist(obs_all_species), MoreArgs=list(prior_prob_presence=prior_prob_presence, mean_frequency=mean_frequency, dp=dp, fdp=fdp), USE.NAMES=TRUE)

	pp_df = matrix(data=tmpdf, ncol=ncol(obs_target_species), nrow=nrow(obs_target_species), byrow=FALSE)
	pp_df
	
	pp_df = adf2(pp_df)
	names(pp_df) = names(obs_target_species)
	rownames(pp_df) = rownames(obs_target_species)
	
	return(pp_df)
	}


#######################################################
# mapply_calc_obs_like
#######################################################
#' Mapply version of calc_obs_like()
#' 
#' This function applies \code{\link{calc_obs_like}} to all cells of 
#' the input matrices \code{obs_target_species} and \code{obs_all_species}. These matrices
#' obviously must have the same dimensions.
#' 
#' The inputs are the same as for \code{\link{calc_obs_like}}, except that
#' \code{obs_target_species} and \code{obs_all_species} can be matrices.
#' 
#' @param truly_present Is the OTU of interest known/conditionally assumed to be truly present (\code{TRUE}) or truly absent (\code{FALSE})?
#' @param obs_target_species A scalar or column/vector/matrix of detection counts, e.g. as produced from the output from \code{\link{read_detections}}.
#' @param obs_all_species A scalar or column/vector/matrix of detection counts, e.g. as produced from the output from \code{\link{read_controls}}.
#' @param mean_frequency This is the proportion of samples from the taphonomic control group that will truly be from this OTU, GIVEN that the OTU is present.
#' This could be estimated, but a decent first guess is (total # samples of OTU of interest / total # of samples in the taphonomic control group where
#' the OTU is known to be present).  All that is really needed is some reasonable value, such that more sampling without detection lowers the 
#' likelihood of the data on the hypothesis of true presence, and vice versa.  This value can only be 1 when the number of detections = the number 
#' of taphonomic control detections, for every OTU and area.  This is the implicit assumption in e.g. standard historical biogeography analyses in 
#' LAGRANGE or BioGeoBEARS.
#' @param dp The detection probability.  This is the per-sample probability that you will correctly detect the OTU in question, 
#' when you are looking at it.  Default is 1, which is the implicit assumption in standard analyses.
#' @param fdp The false detection probability.  This is probability of falsely concluding a detection of the OTU of interest occurred, when in 
#' fact the specimen was of something else.  The default is 0, which assumes zero error rate, 
#' i.e. the assumption being made in all historical biogeography analyses that do not take into account detection 
#' probability.  These options are being included for completeness, but it may not be wise to try to infer \code{mean_frequency},  
#' \code{dp} and \code{fdp} all at once due to identifiability issues (and estimation of fdp may take a very large amount of data).  However, 
#' fixing some of these parameters to reasonable values can allow the user to effectively include beliefs about the uncertainty of the 
#' input data into the analysis, if desired.
#' @return \code{pp_df} A matrix of the natural log-likelihood of the data, given the model & assumption of true presence or absence.
#' @export
#' @seealso \code{\link{calc_obs_like}}, \code{\link{calc_post_prob_presence}}, \code{\link{mapply_calc_post_prob_presence}}, 
#' \code{\link{Pdata_given_rangerow}}, \code{\link[base]{mapply}}, \code{\link{tiplikes_wDetectionModel}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://en.wikipedia.org/wiki/Bayes'_theorem}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Bottjer_Jablonski_1988
#' 	 @cite Bayes_1763
#' @examples
#' test=1
#' # Calculate likelihood of data, given presence in an area, 
#' # given a dp (detection probability) and detection model. 
#' 
#' # soft-coded input files
#' extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
#' detects_fn = np(paste(extdata_dir, "/Psychotria_detections_v1.txt", sep=""))
#' controls_fn = np(paste(extdata_dir, "/Psychotria_controls_v1.txt", sep=""))
#' 
#' OTUnames=NULL
#' areanames=NULL
#' tmpskip=0
#' 
#' detects_df = read_detections(detects_fn, OTUnames=NULL, areanames=NULL, tmpskip=0)
#' controls_df = read_controls(controls_fn, OTUnames=NULL, areanames=NULL, tmpskip=0)
#' 
#' detects_df
#' controls_df
#' detects_df / controls_df
#' 
#' 
#' # Calculate data likelihoods, and posterior probability of presence=TRUE
#' mean_frequency=0.1
#' dp=1
#' fdp=0
#' 
#' mapply_calc_obs_like(truly_present=TRUE, obs_target_species=detects_df,
#' obs_all_species=controls_df, mean_frequency, dp, fdp)
#' 
#' mapply_calc_obs_like(truly_present=FALSE, obs_target_species=detects_df,
#' obs_all_species=controls_df, mean_frequency, dp, fdp)
#' 
#' mapply_calc_post_prob_presence(prior_prob_presence=0.01, 
#' obs_target_species=detects_df,
#' obs_all_species=controls_df, mean_frequency, dp, fdp)
#' 
mapply_calc_obs_like <- function(truly_present=TRUE, obs_target_species, obs_all_species, mean_frequency=0.1, dp=1, fdp=0)
	{
	# Calculate the posterior probability of presence in each area in a matrix / data.frame
	tmpdf = mapply(FUN=calc_obs_like, obs_target_species=unlist(obs_target_species),
	obs_all_species=unlist(obs_all_species), MoreArgs=list(truly_present=truly_present, mean_frequency=mean_frequency, dp=dp, fdp=fdp), USE.NAMES=TRUE)

	like_df = matrix(data=tmpdf, ncol=ncol(obs_target_species), nrow=nrow(obs_target_species), byrow=FALSE)
	like_df
	
	like_df = adf2(like_df)
	names(like_df) = names(obs_target_species)
	rownames(like_df) = rownames(obs_target_species)
	
	return(like_df)
	}



#######################################################
# Pdata_given_range_dp
#######################################################
#' Calculate probability of detection data given a true geographic range and a detection probability
#' 
#' This function calculates P(data|range,dp), i.e. the probability of some detection and 
#' taphonomic control counts, given the true geographic range/state, and \code{dp}, a 
#' detection probability (and, optionally, a false detection probability, \code{fdp}).
#' 
#' The idea of taphonomic controls dates back at least to work of Bottjer & Jablonski (1988).  The basic
#' idea is that if you have taxa of roughly similar detectability, then detections of other 
#' taxa give some idea of overall detection effort.  Obviously this is a very simple 
#' model that can be criticized in any number of ways (different alpha diversity in each 
#' region, different detectability of individual taxa, etc.), but it is a useful starting
#' point as there has been no implementation of any detection model in historical/phylogenetic 
#' biogeography to date.
#' 
#' One could imagine (a) every OTU and area has a different count of detections and taphonomic control 
#' detections, or (b) the taphonomic control detections are specified by area, and shared across all OTUs.
#' Situation (b) is likely more common, but this function assumes (a) as this is the more thorough case.
#' Behavior (b) could be reproduced by summing each column, and/or copying this sum to all cells for a 
#' particular area.
#' 
#' @param truerange_areas The list of areas (as 0-based numbers) in this geographic range/state.
#' @param numareas The function needs to know the total number of areas in the analysis.
#' @param detects_df A column/vector of detection counts, as produced from a column of the output from \code{\link{read_detections}}.
#' @param controls_df A column/vector of detection counts, as produced from a column of the output from \code{\link{read_controls}}.
#' @param dp The detection probability.  This is the per-sample probability that you will detect the OTU in question.
#' In other words, the model assumes that each specimen from the taphonomic control group has a chance of being
#' a representative of the OTU you are looking for.  The default is 1, which assumes perfect detection, i.e. the assumption 
#' being made in all historical biogeography analyses that do not take into account detection probability.  A value of 1
#' will only work when the taphonomic control count equals the detection count; any other data would have likelihood=0.
#' @param fdp The false detection probability.  This is probability of falsely concluding a detection occurred, when in 
#' fact the specimen was of something else.  The default is 0, which assumes zero error rate, 
#' i.e. the assumption being made in all historical biogeography analyses that do not take into account detection 
#' probability.  This option is being included for completeness, but it may not be wise to try to infer both 
#' \code{dp} and \code{fdp} at once due to identifiability issues (and estimation of fdp may take a very large amount of data).
#' @param mean_frequency This is the proportion of samples from the taphonomic control group that will truly be from this OTU, GIVEN that the OTU is present.
#' This could be estimated, but a decent first guess is (total # samples of OTU of interest / total # of samples in the taphonomic control group where
#' the OTU is known to be present).  All that is really needed is some reasonable value, such that more sampling without detection lowers the 
#' likelihood of the data on the hypothesis of true presence, and vice versa.  This value can only be 1 when the number of detections = the number 
#' of taphonomic control detections, for every OTU and area.  This is the implicit assumption in e.g. standard historical biogeography analyses in 
#' LAGRANGE or BioGeoBEARS.
#' @return \code{dtf} 
#' @export
#' @seealso \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Bottjer_Jablonski_1988
#' @examples
#' testval=1
#' 
Pdata_given_rangerow_dp <- function(truerange_areas, numareas, detects_df, controls_df, mean_frequency=0.1, dp=1, fdp=0)
	{
	runjunk='
	controls_fn = "/Library/Frameworks/R.framework/Versions/2.15/Resources/library/BioGeoBEARS/extdata/Psychotria_controls_v1.txt"
	controls_fn = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/Psychotria_controls_v1.txt"
	OTUnames=NULL
	areanames=NULL
	tmpskip=0
	' # end runjunk

	# Likelihood of a state = the probability of data at ALL areas (occupied, and not, in the state/geographic range) given the
	# conditionally assumed true presences/absences in the specified range
	
	# Build a TRUE/FALSE matrix specifying the ranges in this assumed true state/geographic range
	# (copy for every row, this does 1 state/geographic range across all OTUs)
	range_as_areas_TF = matrix(data=FALSE, nrow=nrow(detects_df), ncol=numareas)
	range_as_areas_TF[truerange_areas+1] = TRUE
	
	# OK, calculate likelihood OF THE DATA ON THAT RANGE
	LnLs_of_data_in_each_area = mapply(FUN=calc_obs_like, truly_present=range_as_areas_TF, obs_target_species=detects_df,
obs_all_species=controls_df, MoreArgs=list(mean_frequency=mean_frequency, dp=dp, fdp=fdp), USE.NAMES=TRUE)
	
	# LnL of the data at each tip, for this state
	LnLs_of_this_range_for_each_OTU = rowSums(LnLs_of_data_in_each_area)
	LnLs_of_this_range_for_each_OTU	
	
	return(LnLs_of_this_range_for_each_OTU)
	}





#######################################################
# Pdata_given_rangerow
#######################################################
#' Calculate probability of detection data given a true geographic range and a detection probability
#' 
#' This function calculates P(data|range,dp), i.e. the probability of some detection and 
#' taphonomic control counts, given the true geographic range/state, and \code{dp}, a 
#' detection probability (and, optionally, a false detection probability, \code{fdp}).
#' 
#' The idea of taphonomic controls dates back at least to work of Bottjer & Jablonski (1988).  The basic
#' idea is that if you have taxa of roughly similar detectability, then detections of other 
#' taxa give some idea of overall detection effort.  Obviously this is a very simple 
#' model that can be criticized in any number of ways (different alpha diversity in each 
#' region, different detectability of individual taxa, etc.), but it is a useful starting
#' point as there has been no implementation of any detection model in historical/phylogenetic 
#' biogeography to date.
#' 
#' One could imagine (a) every OTU and area has a different count of detections and taphonomic control 
#' detections, or (b) the taphonomic control detections are specified by area, and shared across all OTUs.
#' Situation (b) is likely more common, but this function assumes (a) as this is the more thorough case.
#' Behavior (b) could be reproduced by summing each column, and/or copying this sum to all cells for a 
#' particular area.
#' 
#' @param range_as_areas_TF The list of areas (as TRUE/FALSE) in this geographic range/state.
#' @param detects_df_row A column/vector of detection counts, as produced from a row of the output from \code{\link{read_detections}}.
#' @param controls_df_row A column/vector of detection counts, as produced from a row of the output from \code{\link{read_controls}}.
#' @param mean_frequency This is the proportion of samples from the taphonomic control group that will truly be from this OTU, GIVEN that the OTU is present.
#' This could be estimated, but a decent first guess is (total # samples of OTU of interest / total # of samples in the taphonomic control group where
#' the OTU is known to be present).  All that is really needed is some reasonable value, such that more sampling without detection lowers the 
#' likelihood of the data on the hypothesis of true presence, and vice versa.  This value can only be 1 when the number of detections = the number 
#' of taphonomic control detections, for every OTU and area.  This is the implicit assumption in e.g. standard historical biogeography analyses in 
#' LAGRANGE or BioGeoBEARS.
#' @param dp The detection probability.  This is the per-sample probability that you will correctly detect the OTU in question, 
#' when you are looking at it.  Default is 1, which is the implicit assumption in standard analyses.
#' @param fdp The false detection probability.  This is probability of falsely concluding a detection of the OTU of interest occurred, when in 
#' fact the specimen was of something else.  The default is 0, which assumes zero error rate, 
#' i.e. the assumption being made in all historical biogeography analyses that do not take into account detection 
#' probability.  These options are being included for completeness, but it may not be wise to try to infer \code{mean_frequency},  
#' \code{dp} and \code{fdp} all at once due to identifiability issues (and estimation of fdp may take a very large amount of data).  However, 
#' fixing some of these parameters to reasonable values can allow the user to effectively include beliefs about the uncertainty of the 
#' input data into the analysis, if desired.
#' @param return_LnLs If \code{FALSE} (default), return exp(sum(LnLs of data in each area)), i.e. the likelihood of the data, 
#' non-logged. If \code{TRUE}, return the LnLs of the data in each area.
#' @return \code{likelihood_of_data_given_range} The (non-logged!) likelihood of the data given the input range, and the 
#' detection model parameters. If \code{return_LnLs=TRUE}, returns \code{LnLs_of_data_in_each_area}, the LnLs of the data in each area.
#' @export
#' @seealso \code{\link{calc_obs_like}}, \code{\link[base]{mapply}}, \code{\link{tiplikes_wDetectionModel}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Bottjer_Jablonski_1988
#' @examples
#' testval=1
#' 
#' # soft-coded input files
#' extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
#' detects_fn = np(paste(extdata_dir, "/Psychotria_detections_v1.txt", sep=""))
#' controls_fn = np(paste(extdata_dir, "/Psychotria_controls_v1.txt", sep=""))
#' 
#' OTUnames=NULL
#' areanames=NULL
#' tmpskip=0
#' 
#' detects_df = read_detections(detects_fn, OTUnames=NULL, areanames=NULL, tmpskip=0)
#' controls_df = read_controls(controls_fn, OTUnames=NULL, areanames=NULL, tmpskip=0)
#' 
#' detects_df
#' controls_df
#' detects_df / controls_df
#' 
#' mean_frequency=0.1
#' dp=1
#' fdp=0
#' 
#' mapply_calc_obs_like(truly_present=TRUE, obs_target_species=detects_df,
#' obs_all_species=controls_df, mean_frequency, dp, fdp)
#' 
#' mapply_calc_obs_like(truly_present=FALSE, obs_target_species=detects_df,
#' obs_all_species=controls_df, mean_frequency, dp, fdp)
#' 
#' mapply_calc_post_prob_presence(prior_prob_presence=0.01, 
#' obs_target_species=detects_df,
#' obs_all_species=controls_df, mean_frequency, dp, fdp)
#' 
#' 
#' 
#' # Now, calculate the likelihood of the data given a geographic range
#' numareas = 4
#' tmpranges = list(c(0), c(1), c(0,1))
#' truerange_areas = tmpranges[[3]]
#' truerange_areas
#' 
#' 	
#' # Build a TRUE/FALSE row specifying the ranges in this assumed true 
#' # state/geographic range
#' range_as_areas_TF = matrix(data=FALSE, nrow=1, ncol=numareas)
#' range_as_areas_TF[truerange_areas+1] = TRUE
#' range_as_areas_TF
#' 
#' detects_df_row = detects_df[1,]
#' controls_df_row = controls_df[1,]
#' 
#' # Manual method, superceded by Pdata_given_rangerow():
#' # LnLs_of_data_in_each_area = mapply(FUN=calc_obs_like, 
#' # obs_target_species=detects_df_row,
#' # obs_all_species=controls_df_row, truly_present=range_as_areas_TF, 
#' # MoreArgs=list(mean_frequency=mean_frequency, dp=dp, fdp=fdp), 
#' # USE.NAMES=TRUE)
#' 
#' 
#' # Calculate data likelihoods on for this geographic range
#' mean_frequency=0.1
#' dp=1
#' fdp=0
#' 

#' # Get the likelihood (the probability of the data, given this range)
#' likelihood_of_data_given_range = Pdata_given_rangerow(
#' range_as_areas_TF=range_as_areas_TF,
#' detects_df_row=detects_df_row, 
#' controls_df_row=controls_df_row, mean_frequency=mean_frequency, dp=dp, fdp=fdp)
#' likelihood_of_data_given_range
#' 
#' # Return the raw log-likelihoods:
#' LnLs_of_data_in_each_area = Pdata_given_rangerow(range_as_areas_TF=range_as_areas_TF,
#' detects_df_row=detects_df_row, 
#' controls_df_row=controls_df_row, mean_frequency=mean_frequency, dp=dp, fdp=fdp, 
#' return_LnLs=TRUE)
#' 
#' detects_df_row
#' controls_df_row
#' LnLs_of_data_in_each_area
#' 
#' # The likelihood: the probability of the data in each area:
#' exp(LnLs_of_data_in_each_area)
#' 
Pdata_given_rangerow <- function(range_as_areas_TF, detects_df_row, controls_df_row, mean_frequency=0.1, dp=1, fdp=0, return_LnLs=FALSE)
	{
	runjunk='
	controls_fn = "/Library/Frameworks/R.framework/Versions/2.15/Resources/library/BioGeoBEARS/extdata/Psychotria_controls_v1.txt"
	controls_fn = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/Psychotria_controls_v1.txt"
	OTUnames=NULL
	areanames=NULL
	tmpskip=0
	' # end runjunk

	# Likelihood of a state = the probability of data at ALL areas (occupied, and not, in the state/geographic range) given the
	# conditionally assumed true presences/absences in the specified range
	
	# Build a TRUE/FALSE matrix specifying the ranges in this assumed true state/geographic range
	# (copy for every row, this does 1 state/geographic range across all OTUs)
	#ranges_as_TF = matrix(data=FALSE, nrow=nrow(detects_df), ncol=numareas)
	#ranges_as_TF[truerange_areas+1] = TRUE

	# OK, calculate likelihood OF THE DATA ON THAT RANGE
	LnLs_of_data_in_each_area = mapply(FUN=calc_obs_like, truly_present=range_as_areas_TF, obs_target_species=detects_df_row,
obs_all_species=controls_df_row, MoreArgs=list(mean_frequency=mean_frequency, dp=0.99, fdp=0.001), USE.NAMES=TRUE)
	
	if (return_LnLs == FALSE)
		{
		LnL_range = sum(LnLs_of_data_in_each_area)
		likelihood_of_data_given_range = exp(LnL_range)
		return(likelihood_of_data_given_range)
		} else {
		return(LnLs_of_data_in_each_area)
		}
	}







#######################################################
# tiplikes_wDetectionModel
#######################################################
#' Calculate probability of detection data for each OTU at each range in a list of states/geographic ranges
#' 
#' This function calculates P(data|range,dp), i.e. the probability of some detection and 
#' taphonomic control counts, given the true geographic range/state, and \code{dp}, a 
#' detection probability (and, optionally, a false detection probability, \code{fdp}).
#' 
#' This function performs the operation for all states/ranges for all tips.
#' 
#' The idea of taphonomic controls dates back at least to work of Bottjer & Jablonski (1988).  The basic
#' idea is that if you have taxa of roughly similar detectability, then detections of other 
#' taxa give some idea of overall detection effort.  Obviously this is a very simple 
#' model that can be criticized in any number of ways (different alpha diversity in each 
#' region, different detectability of individual taxa, etc.), but it is a useful starting
#' point as there has been no implementation of any detection model in historical/phylogenetic 
#' biogeography to date.
#' 
#' One could imagine (a) every OTU and area has a different count of detections and taphonomic control 
#' detections, or (b) the taphonomic control detections are specified by area, and shared across all OTUs.
#' Situation (b) is likely more common, but this function assumes (a) as this is the more thorough case.
#' Behavior (b) could be reproduced by summing each column, and/or copying this sum to all cells for a 
#' particular area.
#' 
#' @param states_list_0based_index A states_list, 0-based, e.g. from \code{\link[cladoRcpp]{rcpp_areas_list_to_states_list}}.
#' @param numareas The number of areas being considered in the analysis. If \code{NULL} (default), this is calculated to be the maximum range length, or 
#' one plus the maximum 0-based index in any of the ranges.
#' @param detects_df A matrix/data.frame of detection counts, as produced from the output from \code{\link{read_detections}}.
#' @param controls_df A matrix/data.frame of detection counts, as produced from the output from \code{\link{read_controls}}.
#' @param mean_frequency This is the proportion of samples from the taphonomic control group that will truly be from this OTU, GIVEN that the OTU is present.
#' This could be estimated, but a decent first guess is (total # samples of OTU of interest / total # of samples in the taphonomic control group where
#' the OTU is known to be present).  All that is really needed is some reasonable value, such that more sampling without detection lowers the 
#' likelihood of the data on the hypothesis of true presence, and vice versa.  This value can only be 1 when the number of detections = the number 
#' of taphonomic control detections, for every OTU and area.  This is the implicit assumption in e.g. standard historical biogeography analyses in 
#' LAGRANGE or BioGeoBEARS.
#' @param dp The detection probability.  This is the per-sample probability that you will correctly detect the OTU in question, 
#' when you are looking at it.  Default is 1, which is the implicit assumption in standard analyses.
#' @param fdp The false detection probability.  This is probability of falsely concluding a detection of the OTU of interest occurred, when in 
#' fact the specimen was of something else.  The default is 0, which assumes zero error rate, 
#' i.e. the assumption being made in all historical biogeography analyses that do not take into account detection 
#' probability.  These options are being included for completeness, but it may not be wise to try to infer \code{mean_frequency},  
#' \code{dp} and \code{fdp} all at once due to identifiability issues (and estimation of fdp may take a very large amount of data).  However, 
#' fixing some of these parameters to reasonable values can allow the user to effectively include beliefs about the uncertainty of the 
#' input data into the analysis, if desired.
#' @param null_range_gets_0_like If \code{TRUE} (default), then the data is given zero probability on the hypothesis that the 
#' range is a null range (i.e., no areas occupied).  This is equivalent to saying that you are sure/are willing to assume that 
#' the OTU exists somewhere in your study area, at the timepoint being considered.  Null ranges are identified by length=1, 
#' containing NULL, NA, "", "_", etc.
#' @return \code{tip_condlikes_of_data_on_each_state} The (non-logged!) likelihood of the data for each tip, given each possible range, and the 
#' detection model parameters.
#' @export
#' @seealso \code{\link{Pdata_given_rangerow}}, \code{\link{calc_obs_like}}, \code{\link[base]{mapply}}, \code{\link{read_detections}}, \code{\link{read_controls}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Bottjer_Jablonski_1988
#' @examples
#' testval=1
#' 
#' # soft-coded input files
#' extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
#' detects_fn = np(paste(extdata_dir, "/Psychotria_detections_v1.txt", sep=""))
#' controls_fn = np(paste(extdata_dir, "/Psychotria_controls_v1.txt", sep=""))
#' 
#' detects_df = read_detections(detects_fn, OTUnames=NULL, areanames=NULL, tmpskip=0)
#' controls_df = read_controls(controls_fn, OTUnames=NULL, areanames=NULL, tmpskip=0)
#' 
#' # Calculate the likelihood of the data at each tip, for each possible geographic range
#' numareas = 4
#' tmpranges = list(c(0), c(1), c(0,1))
#' 
#' mean_frequency=0.1
#' dp=1
#' fdp=0
#' 
#' tip_condlikes_of_data_on_each_state = 
#' tiplikes_wDetectionModel(states_list_0based_index=tmpranges, numareas=numareas, 
#' detects_df, controls_df, mean_frequency=mean_frequency, dp=dp, fdp=fdp, 
#' null_range_gets_0_like=TRUE)
#' 
#' tip_condlikes_of_data_on_each_state
#' 
tiplikes_wDetectionModel <- function(states_list_0based_index, numareas=NULL, detects_df, controls_df, mean_frequency=0.1, dp=1, fdp=0, null_range_gets_0_like=TRUE)
	{
	runjunk='
	controls_fn = "/Library/Frameworks/R.framework/Versions/2.15/Resources/library/BioGeoBEARS/extdata/Psychotria_controls_v1.txt"
	controls_fn = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/Psychotria_controls_v1.txt"
	OTUnames=NULL
	areanames=NULL
	tmpskip=0
	' # end runjunk
	
	
	# Get the number of areas, if needed
	if (is.null(numareas))
		{
		maxlength = max(unlist(lapply(X=states_list_0based_index,FUN=length)), na.rm=TRUE)
		maxnum = 1+max(unlist(lapply(X=states_list_0based_index,FUN=max, na.rm=TRUE)), na.rm=TRUE)
		numareas = max(c(maxnum, maxlength), na.rm=TRUE)
		}
	
	# Likelihood of a state = the probability of data at ALL areas (occupied, and not, in the state/geographic range) given the
	# conditionally assumed true presences/absences in the specified range
	
	# Go through the states
	#states_list_0based_index = tmpstates
	numstates = length(states_list_0based_index)
	numtips = nrow(detects_df)
	tip_condlikes_of_data_on_each_state = matrix(data=0, nrow=numtips, ncol=numstates)
	for (statenum in 1:numstates)
		{
		# Get the current state/geographic range
		truerange_areas = states_list_0based_index[[statenum]]
		
		# Check for null range, if you're supposed to, and the length is 1
		if ((null_range_gets_0_like==TRUE) && (length(truerange_areas)==1))
			{
			if ( (truerange_areas=="_") || (is.null(truerange_areas)==TRUE) || (is.na(truerange_areas)==TRUE) || (truerange_areas=="") )
				{
				# Then you have a null range, and you have decided to give the data 0 likelihood on all states
				likes_vector = rep(0, times=numtips)
				tip_condlikes_of_data_on_each_state[, statenum] = likes_vector
				
				# skip to next loop
				next()
				}
			}
		
		# Build a TRUE/FALSE row specifying the ranges in this assumed true 
		# state/geographic range
		range_as_areas_TF = matrix(data=FALSE, nrow=nrow(detects_df), ncol=numareas)
		range_as_areas_TF[, truerange_areas+1] = TRUE	# MUST BE OVER ALL ROWS
		range_as_areas_TF

		# iterate over the rows
		likes_vector = rep(0, times=numtips)
		for (i in 1:nrow(range_as_areas_TF))
			{
			likes_vector[i] = Pdata_given_rangerow(range_as_areas_TF=range_as_areas_TF[i,], detects_df_row=detects_df[i,], 
		controls_df_row=controls_df[i,], mean_frequency=mean_frequency, dp=dp, fdp=fdp, return_LnLs=FALSE)
			}
	
		tip_condlikes_of_data_on_each_state[, statenum] = likes_vector
		}

	tip_condlikes_of_data_on_each_state



	return(tip_condlikes_of_data_on_each_state)
	}







#######################################################
# prob_of_states_from_prior_prob_areas
#######################################################
#' Calculate probability of detection data for each OTU at each range in a list of states/geographic ranges
#' 
#' This function calculates P(data|range,dp), i.e. the probability of some detection and 
#' taphonomic control counts, given the true geographic range/state, and \code{dp}, a 
#' detection probability (and, optionally, a false detection probability, \code{fdp}).
#' 
#' This function performs the operation for all states/ranges for all tips.
#' 
#' The idea of taphonomic controls dates back at least to work of Bottjer & Jablonski (1988).  The basic
#' idea is that if you have taxa of roughly similar detectability, then detections of other 
#' taxa give some idea of overall detection effort.  Obviously this is a very simple 
#' model that can be criticized in any number of ways (different alpha diversity in each 
#' region, different detectability of individual taxa, etc.), but it is a useful starting
#' point as there has been no implementation of any detection model in historical/phylogenetic 
#' biogeography to date.
#' 
#' One could imagine (a) every OTU and area has a different count of detections and taphonomic control 
#' detections, or (b) the taphonomic control detections are specified by area, and shared across all OTUs.
#' Situation (b) is likely more common, but this function assumes (a) as this is the more thorough case.
#' Behavior (b) could be reproduced by summing each column, and/or copying this sum to all cells for a 
#' particular area.
#' 
#' @param states_list_0based_index A states_list, 0-based, e.g. from \code{\link[cladoRcpp]{rcpp_areas_list_to_states_list}}.
#' @param numareas The number of areas being considered in the analysis. If \code{NULL} (default), this is calculated to be the maximum range length, or 
#' one plus the maximum 0-based index in any of the ranges.
#' @param prior_prob_presence The prior probability of presence, i.e. when no detection or taphonomic control data 
#' whatsoever is available.  Default is set to 0.01 which expresses my totally uninformed bias that 
#' in whatever your data is, your species of interest probably doesn't live in the typical area you are 
#' looking at.
#' @param null_range_gets_0_prob If \code{TRUE} (default), then the null range is given zero probability. 
#' A null range has no areas occupied.  This is equivalent to saying that you are sure/are willing to assume that 
#' the OTU exists somewhere in your study area, at the timepoint being considered.  Null ranges are identified by length=1, 
#' containing NULL, NA, "", "_", etc.
#' @param normalize_probs If \code{TRUE}, the probabilities of each range will be normalized so that they sum to 1. Otherwise, they 
#' won't.
#' @return \code{prob_of_each_range} The probability of each range, given the prior probability of presence in each area.
#' @export
#' @seealso \code{\link[cladoRcpp]{rcpp_areas_list_to_states_list}}, \code{\link{Pdata_given_rangerow}}, 
#' \code{\link{calc_obs_like}}, \code{\link[base]{mapply}}, \code{\link{read_detections}}, \code{\link{read_controls}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Bottjer_Jablonski_1988
#' @examples
#' testval=1
#' 
#' prior_prob_presence = 0.01
#' 
#' areas = c("K", "O", "M", "H")
#' numareas = length(areas)
#' states_list_0based_index = 
#' rcpp_areas_list_to_states_list(areas=areas, maxareas=4, include_null_range=TRUE)
#' states_list_0based_index
#' numareas = 4
#' 
#' mean_frequency=0.1
#' dp=1
#' fdp=0
#' 
#' prob_of_states_from_prior_prob_areas(states_list_0based_index, numareas=numareas,
#' prior_prob_presence=prior_prob_presence, null_range_gets_0_prob=TRUE, 
#' normalize_probs=TRUE)
#' 
#' prob_of_states_from_prior_prob_areas(states_list_0based_index, numareas=numareas,
#' prior_prob_presence=prior_prob_presence, null_range_gets_0_prob=TRUE, 
#' normalize_probs=FALSE)
#' 
#' prob_of_states_from_prior_prob_areas(states_list_0based_index, numareas=numareas,
#' prior_prob_presence=prior_prob_presence, null_range_gets_0_prob=FALSE, 
#' normalize_probs=TRUE)
#' 
#' prob_of_states_from_prior_prob_areas(states_list_0based_index, numareas=numareas,
#' prior_prob_presence=prior_prob_presence, null_range_gets_0_prob=FALSE, 
#' normalize_probs=FALSE)
#' 
prob_of_states_from_prior_prob_areas <- function(states_list_0based_index, numareas=NULL, prior_prob_presence=0.01, null_range_gets_0_prob=TRUE, normalize_probs=TRUE)
	{
	runjunk='
	' # end runjunk
	
	
	# Get the number of areas, if needed
	if (is.null(numareas))
		{
		maxlength = max(unlist(lapply(X=states_list_0based_index,FUN=length)), na.rm=TRUE)
		maxnum = 1+max(unlist(lapply(X=states_list_0based_index,FUN=max, na.rm=TRUE)), na.rm=TRUE)
		numareas = max(c(maxnum, maxlength), na.rm=TRUE)
		}
	
	numstates = length(states_list_0based_index)

	range_as_areas_TF = matrix(data=FALSE, nrow=numstates, ncol=numareas)
	prob_of_TF = matrix(data=0, nrow=numstates, ncol=numareas)
	for (statenum in 1:numstates)
		{
		# Get the current state/geographic range
		truerange_areas = states_list_0based_index[[statenum]]

		# Check for null range, if you're supposed to, and the length is 1
		if ((null_range_gets_0_prob==TRUE) && (length(truerange_areas)==1))
			{
			if ( (truerange_areas=="_") || (is.null(truerange_areas)==TRUE) || (is.na(truerange_areas)==TRUE) || (truerange_areas=="") )
				{
				# Then you have a null range, and you have decided to give the data 0 likelihood on all states
				range_as_areas_TF[statenum, ] = NA
			
				# skip to next loop
				next()
				}
			}

	
		# Build a TRUE/FALSE row specifying the ranges in this assumed true 
		# state/geographic range
		range_as_areas_TF[statenum, truerange_areas+1] = TRUE	# MUST BE OVER ALL ROWS
		range_as_areas_TF
		}
	range_as_areas_TF

	# Prob of each range
	prob_of_TF[range_as_areas_TF == TRUE] = prior_prob_presence
	prob_of_TF[range_as_areas_TF == FALSE] = 1-prior_prob_presence

	# Set anything in the NULL range to 0
	range_as_areas_TF[is.na(range_as_areas_TF)] = 0

	prob_of_TF
	prob_of_each_range = exp(rowSums(log(prob_of_TF)))
	prob_of_each_range

	if (normalize_probs == TRUE)
		{
		prob_of_each_range = prob_of_each_range / sum(prob_of_each_range)
		prob_of_each_range
		}



	return(prob_of_each_range)
	}




#######################################################
# post_prob_states
#######################################################
#' Calculate posterior probability of each states/geographic ranges, given prior probabilities and data likelihoods
#' 
#' This function calculates P(range|data,detection model), i.e. the probability of each possible 
#' range, given a prior probability of each range, and the likelihood of each range.
#' 
#' The prior probability of each range should be considered by the user.  Note that putting the same prior on 
#' the probability of occurrence in each individual range does NOT mean a flat prior on each 
#' state/geographic range.  This fact is demonstrated in the function \code{\link{prob_of_states_from_prior_prob_areas}}.
#' 
#' @param prob_of_each_range The probability of each range, given the prior probability of presence in each area.
#' @param condlikes_of_data_on_each_range The probability of the data, conditional on each range (i.e., the likelihood), as 
#' found in e.g. a row of the output from \code{\link{tiplikes_wDetectionModel}}.
#' 
#' @return \code{posterior_probs} The posterior probability of each range.
#' @export
#' @seealso \code{\link{prob_of_states_from_prior_prob_areas}}, \code{\link{tiplikes_wDetectionModel}}, 
#' \code{\link[cladoRcpp]{rcpp_areas_list_to_states_list}}, \code{\link{Pdata_given_rangerow}}, 
#' \code{\link{calc_obs_like}}, \code{\link[base]{mapply}}, \code{\link{read_detections}}, \code{\link{read_controls}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://en.wikipedia.org/wiki/Log_probability}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Bottjer_Jablonski_1988
#' @examples
#' testval=1
#' 
#' # soft-coded input files
#' extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
#' detects_fn = np(paste(extdata_dir, "/Psychotria_detections_v1.txt", sep=""))
#' controls_fn = np(paste(extdata_dir, "/Psychotria_controls_v1.txt", sep=""))
#' 
#' detects_df = read_detections(detects_fn, OTUnames=NULL, areanames=NULL, tmpskip=0)
#' controls_df = read_controls(controls_fn, OTUnames=NULL, areanames=NULL, tmpskip=0)
#' 
#' 
#' # Setup
#' prior_prob_presence = 0.01
#' areas = c("K", "O", "M", "H")
#' numareas = length(areas)
#' maxareas = length(areas)
#' states_list_0based_index = 
#' rcpp_areas_list_to_states_list(areas=areas, maxareas=maxareas, 
#'                                include_null_range=TRUE)
#' states_list_0based_index
#' 
#' mean_frequency=0.1
#' dp=1
#' fdp=0
#' 
#' tip_condlikes_of_data_on_each_state = 
#' tiplikes_wDetectionModel(states_list_0based_index, numareas=numareas, 
#' detects_df, controls_df, mean_frequency=mean_frequency, dp=dp, fdp=fdp, 
#' null_range_gets_0_like=TRUE)
#' 
#' tip_condlikes_of_data_on_each_state
#' 
#' 
#' 
#' # To get denominator, just iterate over all the states
#' # Prior probability
#' prob_of_each_range = prob_of_states_from_prior_prob_areas(states_list_0based_index, 
#' numareas=numareas,
#' prior_prob_presence=prior_prob_presence, null_range_gets_0_prob=TRUE, 
#' normalize_probs=TRUE)
#' 
#' # Likelihoods of the data on each range
#' condlikes_of_data_on_each_range = tip_condlikes_of_data_on_each_state[1,]
#' 
#' posterior_probs = post_prob_states(prob_of_each_range, 
#'                   condlikes_of_data_on_each_range)
#' posterior_probs
#' 
#' # Should sum to 1
#' sum(posterior_probs)
#' 
post_prob_states <- function(prob_of_each_range, condlikes_of_data_on_each_range)
	{
	# The denominator in Bayes's Theorem is the P(data), i.e. the probability of the 
	# data conditional on all possible hypotheses, i.e. the sum of the data likelihood times the 
	# prior for all states.\

	# Sum probabilities in log space
	# http://en.wikipedia.org/wiki/Log_probability
	log_probs_in_denominator = log(prob_of_each_range) + log(condlikes_of_data_on_each_range)
	log_probs_in_denominator

	# Simplest way to add probabilities in log formet
	denominator = log(sum(exp(log_probs_in_denominator)))
	denominator


	numstates = length(prob_of_each_range)
	numerators = rep(NA, times=numstates)
	for (i in 1:numstates)
		{
		prior_prob_of_this_range = prob_of_each_range[i]
		like_of_data_on_this_range = condlikes_of_data_on_each_range[i]
	
		# Sum probabilities in log space
		# http://en.wikipedia.org/wiki/Log_probability
		numerators[i] = log(prior_prob_of_this_range * like_of_data_on_this_range)
		}

	# Calculate the posterior probabilities
	posterior_probs = exp(numerators - denominator)
	posterior_probs
	
	# They should sum to 1
	sum(posterior_probs)
	
	return(posterior_probs)
	}








#######################################################
# post_prob_states_matrix
#######################################################
#' Calculate posterior probability of each states/geographic ranges, given prior probabilities and data likelihoods
#' 
#' This function calculates P(range|data,detection model), i.e. the probability of each possible 
#' range, given a prior probability of each range, and the likelihood of each range.
#' 
#' The prior probability of each range should be considered by the user.  Note that putting the same prior on 
#' the probability of occurrence in each individual range does NOT mean a flat prior on each 
#' state/geographic range.  This fact is demonstrated in the function \code{\link{prob_of_states_from_prior_prob_areas}}.
#' 
#' @param prob_of_each_range The probability of each range, given the prior probability of presence in each area.
#' @param tip_condlikes_of_data_on_each_state The probability of the data, conditional on each range (i.e., the likelihood), as 
#' found in e.g. a row of the output from \code{\link{tiplikes_wDetectionModel}}.
#' @return \code{posterior_probs} The posterior probability of each range.
#' @export
#' @seealso \code{\link{prob_of_states_from_prior_prob_areas}}, \code{\link{tiplikes_wDetectionModel}}, 
#' \code{\link[cladoRcpp]{rcpp_areas_list_to_states_list}}, \code{\link{Pdata_given_rangerow}}, 
#' \code{\link{calc_obs_like}}, \code{\link[base]{mapply}}, \code{\link{read_detections}}, \code{\link{read_controls}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{http://en.wikipedia.org/wiki/Log_probability}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite Bottjer_Jablonski_1988
#' @examples
#' testval=1
#' 
#' # soft-coded input files
#' extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
#' #extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"
#' detects_fn = np(paste(extdata_dir, "/Psychotria_detections_v1.txt", sep=""))
#' controls_fn = np(paste(extdata_dir, "/Psychotria_controls_v1.txt", sep=""))
#' 
#' detects_df = read_detections(detects_fn, OTUnames=NULL, areanames=NULL, tmpskip=0)
#' controls_df = read_controls(controls_fn, OTUnames=NULL, areanames=NULL, tmpskip=0)
#' 
#' 
#' # Setup
#' prior_prob_presence = 0.01
#' areas = c("K", "O", "M", "H")
#' numareas = length(areas)
#' maxareas = length(areas)
#' states_list_0based_index = 
#' rcpp_areas_list_to_states_list(areas=areas, maxareas=maxareas, 
#'                                include_null_range=TRUE)
#' states_list_0based_index
#' 
#' mean_frequency=0.1
#' dp=1
#' fdp=0
#' 
#' tip_condlikes_of_data_on_each_state = 
#' tiplikes_wDetectionModel(states_list_0based_index, numareas=numareas, 
#' detects_df, controls_df, mean_frequency=mean_frequency, dp=dp, fdp=fdp, 
#' null_range_gets_0_like=TRUE)
#' 
#' tip_condlikes_of_data_on_each_state
#' 
#' 
#' 
#' # To get denominator, just iterate over all the states
#' # Prior probability
#' prob_of_each_range = prob_of_states_from_prior_prob_areas(states_list_0based_index, 
#' numareas=numareas,
#' prior_prob_presence=prior_prob_presence, null_range_gets_0_prob=TRUE, 
#' normalize_probs=TRUE)
#' 
#' posterior_probs_matrix = post_prob_states_matrix(prob_of_each_range, 
#'                   tip_condlikes_of_data_on_each_state)
#' posterior_probs_matrix
#' 
#' # Should sum to 1
#' rowSums(posterior_probs_matrix)
#' 
#' # How does posterior probability correlate with likelihood and prior probability?
#' par(mfrow=c(1,2))
#' plot(x=jitter(log(tip_condlikes_of_data_on_each_state)), 
#' y=jitter(log(posterior_probs_matrix)))
#' title("Correlation of data likelihoods\nand posterior probabilities")
#' 
#' prob_of_each_range_matrix = matrix(data=prob_of_each_range, 
#' nrow=nrow(posterior_probs_matrix), ncol=length(prob_of_each_range))
#' plot(x=jitter(log(prob_of_each_range_matrix)), 
#' y=jitter(log(posterior_probs_matrix)))
#' title("Correlation of prior probability\nand posterior probabilities")
#' 
post_prob_states_matrix <- function(prob_of_each_range, tip_condlikes_of_data_on_each_state)
	{
	# The denominator in Bayes's Theorem is the P(data), i.e. the probability of the 
	# data conditional on all possible hypotheses, i.e. the sum of the data likelihood times the 
	# prior for all states.\
	
	numtips = nrow(tip_condlikes_of_data_on_each_state)
	posterior_probs_matrix = matrix(NA, nrow=numtips, ncol=ncol(tip_condlikes_of_data_on_each_state))
	
	# Iterate through each tip
	for (tipnum in 1:numtips)
		{
		condlikes_of_data_on_each_range = tip_condlikes_of_data_on_each_state[tipnum, ]
	
		# Sum probabilities in log space
		# http://en.wikipedia.org/wiki/Log_probability
		log_probs_in_denominator = log(prob_of_each_range) + log(condlikes_of_data_on_each_range)
		log_probs_in_denominator

		# Simplest way to add probabilities in log formet
		denominator = log(sum(exp(log_probs_in_denominator)))
		denominator


		numstates = length(prob_of_each_range)
		numerators = rep(NA, times=numstates)
		for (i in 1:numstates)
			{
			prior_prob_of_this_range = prob_of_each_range[i]
			like_of_data_on_this_range = condlikes_of_data_on_each_range[i]
	
			# Sum probabilities in log space
			# http://en.wikipedia.org/wiki/Log_probability
			numerators[i] = log(prior_prob_of_this_range * like_of_data_on_this_range)
			}

		# Calculate the posterior probabilities
		posterior_probs = exp(numerators - denominator)
		posterior_probs
		
		posterior_probs_matrix[tipnum, ] = posterior_probs
		}
	
	# They should sum to 1
	rowSums(posterior_probs_matrix)
	
	return(posterior_probs_matrix)
	}





# add_probs_as_logs <- function(lnP_list)
# 	{
# 	
# 	# Complex:
# #	num_lnPs_to_add = length(lnP_list)
# # 	sum_lnPs = lnP_list[1]
# # 	for (i in 2:num_lnPs_to_add)
# # 		{
# # 		lnP = lnP_list[i]
# # 		
# # 		# Sum probabilities in log space
# # 		# http://en.wikipedia.org/wiki/Log_probability
# # 		#sum_lnPs = -1 * log( exp(-sum_lnPs) + exp(-lnP) )
# # 		
# # 		# Simpler:
# # 		sum_lnPs = 1 * log( exp(sum_lnPs) + exp(lnP) )
# # 		}
# 
# 	# Simple:
# 	sum_lnPs = log(sum(exp(lnP_list)))
# 
# 	return(sum_lnPs)
# 	}

