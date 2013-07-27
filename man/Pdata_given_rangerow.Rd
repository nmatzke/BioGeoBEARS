\name{Pdata_given_rangerow}
\alias{Pdata_given_rangerow}
\title{Calculate probability of detection data given a true geographic range and a detection probability}
\usage{
  Pdata_given_rangerow(range_as_areas_TF, detects_df_row,
    controls_df_row, mean_frequency = 0.1, dp = 1, fdp = 0,
    return_LnLs = FALSE)
}
\arguments{
  \item{range_as_areas_TF}{The list of areas (as
  TRUE/FALSE) in this geographic range/state.}

  \item{detects_df_row}{A column/vector of detection
  counts, as produced from a row of the output from
  \code{\link{read_detections}}.}

  \item{controls_df_row}{A column/vector of detection
  counts, as produced from a row of the output from
  \code{\link{read_controls}}.}

  \item{mean_frequency}{This is the proportion of samples
  from the taphonomic control group that will truly be from
  this OTU, GIVEN that the OTU is present. This could be
  estimated, but a decent first guess is (total # samples
  of OTU of interest / total # of samples in the taphonomic
  control group where the OTU is known to be present).  All
  that is really needed is some reasonable value, such that
  more sampling without detection lowers the likelihood of
  the data on the hypothesis of true presence, and vice
  versa.  This value can only be 1 when the number of
  detections = the number of taphonomic control detections,
  for every OTU and area.  This is the implicit assumption
  in e.g. standard historical biogeography analyses in
  LAGRANGE or BioGeoBEARS.}

  \item{dp}{The detection probability.  This is the
  per-sample probability that you will correctly detect the
  OTU in question, when you are looking at it.  Default is
  1, which is the implicit assumption in standard
  analyses.}

  \item{fdp}{The false detection probability.  This is
  probability of falsely concluding a detection of the OTU
  of interest occurred, when in fact the specimen was of
  something else.  The default is 0, which assumes zero
  error rate, i.e. the assumption being made in all
  historical biogeography analyses that do not take into
  account detection probability.  These options are being
  included for completeness, but it may not be wise to try
  to infer \code{mean_frequency}, \code{dp} and \code{fdp}
  all at once due to identifiability issues (and estimation
  of fdp may take a very large amount of data).  However,
  fixing some of these parameters to reasonable values can
  allow the user to effectively include beliefs about the
  uncertainty of the input data into the analysis, if
  desired.}

  \item{return_LnLs}{If \code{FALSE} (default), return
  exp(sum(LnLs of data in each area)), i.e. the likelihood
  of the data, non-logged. If \code{TRUE}, return the LnLs
  of the data in each area.}
}
\value{
  \code{likelihood_of_data_given_range} The (non-logged!)
  likelihood of the data given the input range, and the
  detection model parameters. If \code{return_LnLs=TRUE},
  returns \code{LnLs_of_data_in_each_area}, the LnLs of the
  data in each area.
}
\description{
  This function calculates P(data|range,dp), i.e. the
  probability of some detection and taphonomic control
  counts, given the true geographic range/state, and
  \code{dp}, a detection probability (and, optionally, a
  false detection probability, \code{fdp}).
}
\details{
  The idea of taphonomic controls dates back at least to
  work of Bottjer & Jablonski (1988).  The basic idea is
  that if you have taxa of roughly similar detectability,
  then detections of other taxa give some idea of overall
  detection effort.  Obviously this is a very simple model
  that can be criticized in any number of ways (different
  alpha diversity in each region, different detectability
  of individual taxa, etc.), but it is a useful starting
  point as there has been no implementation of any
  detection model in historical/phylogenetic biogeography
  to date.

  One could imagine (a) every OTU and area has a different
  count of detections and taphonomic control detections, or
  (b) the taphonomic control detections are specified by
  area, and shared across all OTUs. Situation (b) is likely
  more common, but this function assumes (a) as this is the
  more thorough case. Behavior (b) could be reproduced by
  summing each column, and/or copying this sum to all cells
  for a particular area.
}
\note{
  Go BEARS!
}
\examples{
testval=1

# soft-coded input files
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
detects_fn = np(paste(extdata_dir, "/Psychotria_detections_v1.txt", sep=""))
controls_fn = np(paste(extdata_dir, "/Psychotria_controls_v1.txt", sep=""))

OTUnames=NULL
areanames=NULL
tmpskip=0

detects_df = read_detections(detects_fn, OTUnames=NULL, areanames=NULL, tmpskip=0)
controls_df = read_controls(controls_fn, OTUnames=NULL, areanames=NULL, tmpskip=0)

detects_df
controls_df
detects_df / controls_df

mean_frequency=0.1
dp=1
fdp=0

mapply_calc_obs_like(truly_present=TRUE, obs_target_species=detects_df,
obs_all_species=controls_df, mean_frequency, dp, fdp)

mapply_calc_obs_like(truly_present=FALSE, obs_target_species=detects_df,
obs_all_species=controls_df, mean_frequency, dp, fdp)

mapply_calc_post_prob_presence(prior_prob_presence=0.01,
obs_target_species=detects_df,
obs_all_species=controls_df, mean_frequency, dp, fdp)



# Now, calculate the likelihood of the data given a geographic range
numareas = 4
tmpranges = list(c(0), c(1), c(0,1))
truerange_areas = tmpranges[[3]]
truerange_areas


# Build a TRUE/FALSE row specifying the ranges in this assumed true
# state/geographic range
range_as_areas_TF = matrix(data=FALSE, nrow=1, ncol=numareas)
range_as_areas_TF[truerange_areas+1] = TRUE
range_as_areas_TF

detects_df_row = detects_df[1,]
controls_df_row = controls_df[1,]

# Manual method, superceded by Pdata_given_rangerow():
# LnLs_of_data_in_each_area = mapply(FUN=calc_obs_like,
# obs_target_species=detects_df_row,
# obs_all_species=controls_df_row, truly_present=range_as_areas_TF,
# MoreArgs=list(mean_frequency=mean_frequency, dp=dp, fdp=fdp),
# USE.NAMES=TRUE)


# Calculate data likelihoods on for this geographic range
mean_frequency=0.1
dp=1
fdp=0

# Get the likelihood (the probability of the data, given this range)
likelihood_of_data_given_range = Pdata_given_rangerow(
range_as_areas_TF=range_as_areas_TF,
detects_df_row=detects_df_row,
controls_df_row=controls_df_row, mean_frequency=mean_frequency, dp=dp, fdp=fdp)
likelihood_of_data_given_range

# Return the raw log-likelihoods:
LnLs_of_data_in_each_area = Pdata_given_rangerow(range_as_areas_TF=range_as_areas_TF,
detects_df_row=detects_df_row,
controls_df_row=controls_df_row, mean_frequency=mean_frequency, dp=dp, fdp=fdp,
return_LnLs=TRUE)

detects_df_row
controls_df_row
LnLs_of_data_in_each_area

# The likelihood: the probability of the data in each area:
exp(LnLs_of_data_in_each_area)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS

  Bottjer_Jablonski_1988
}
\seealso{
  \code{\link{calc_obs_like}}, \code{\link[base]{mapply}},
  \code{\link{tiplikes_wDetectionModel}}
}

