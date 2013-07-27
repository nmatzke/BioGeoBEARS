\name{calc_obs_like}
\alias{calc_obs_like}
\title{Calculate likelihood of count data given true presence/absence and parameters}
\usage{
  calc_obs_like(truly_present = TRUE, obs_target_species,
    obs_all_species, mean_frequency = 0.1, dp = 1, fdp = 0)
}
\arguments{
  \item{truly_present}{Is the OTU of interest
  known/conditionally assumed to be truly present
  (\code{TRUE}) or truly absent (\code{FALSE})?}

  \item{obs_target_species}{A count of detections of your
  OTU of interest, e.g. as produced from a cell of the
  matrix output from \code{\link{read_detections}}.}

  \item{obs_all_species}{A count of detections of your
  taphonomic controls, e.g. as produced from a cell of the
  output from \code{\link{read_controls}}.}

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
}
\value{
  \code{lnlike_allobs_given_absence} The natural
  log-likelihood of the data, given the model & assumption
  of true presence or absence.
}
\description{
  This function calculates P(data|presence,parameters),
  i.e. the probability of some detection and taphonomic
  control counts, given the true geographic range/state,
  and parameters such as \code{dp}, a detection probability
  (and, optionally, a false detection probability,
  \code{fdp}).
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
# Example: 10 observations of the species mean dramatically higher likelihood of the
# data on the hypothesis that it is truly present.

# With zero error rate
obs_target_species = 10
obs_all_species = 100
mean_frequency=0.1
dp=1
fdp=0
LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_presence
LnL_under_absence
# Note that the probability of getting detections, under the hypothesis of
# true absence, is -Inf


# With a small error rate, there is some small but positive probability of
# falsely getting 10 detections
obs_target_species = 10
obs_all_species = 100
mean_frequency=0.1
dp=0.99
fdp=0.001
LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_presence
LnL_under_absence
# i.e. the prob. of the data is 1 under the hypothesis of presence, and 0
# under the hypothesis of absence (ln(prob) = 0 & -Inf, respectively)


# Note that with very high error rates, your conclusion could reverse
obs_target_species = 10
obs_all_species = 100
mean_frequency=0.1
dp=0.5
fdp=0.3
LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_presence
LnL_under_absence


# Example #2 -- what if you have ZERO detections, but lots of detections
# of your taphonomic control?
obs_target_species = 0
obs_all_species = 1
mean_frequency=0.1
dp=1
fdp=0
LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_presence
LnL_under_absence

# With a slight error rate
obs_target_species = 0
obs_all_species = 1
mean_frequency=0.1
dp=0.99
fdp=0.001
LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_presence
LnL_under_absence


obs_target_species = 0
obs_all_species = 2
mean_frequency=0.1
dp=1
fdp=0
LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_presence
LnL_under_absence

# With a slight error rate
obs_target_species = 0
obs_all_species = 2
mean_frequency=0.1
dp=0.99
fdp=0.001
LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_presence
LnL_under_absence





# Example #3 -- what if you have ZERO detections, but only a few
# detections of your taphonomic control?
obs_target_species = 0
obs_all_species = 100
mean_frequency=0.1
dp=1
fdp=0
LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_presence
LnL_under_absence

# With a slight error rate
obs_target_species = 0
obs_all_species = 100
mean_frequency=0.1
dp=0.99
fdp=0.001
LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_presence
LnL_under_absence



# Special cases -- e.g., no data
# Prob(data)=1, ln(prob)=0
obs_target_species = 0
obs_all_species = 0
mean_frequency=0.1
dp=0.99
fdp=0.001
LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_presence
LnL_under_absence

obs_target_species = 0
obs_all_species = 0
mean_frequency=0.1
dp=1
fdp=0
LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_presence
LnL_under_absence


# What if, for some reason, you put in identical detections and taphonomic control
# counts? (e.g., you load in a standard tipranges file)
obs_target_species = 1
obs_all_species = 1
mean_frequency=1
dp=1
fdp=0
LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_presence
LnL_under_absence

# What if, for some reason, you put in identical detections and taphonomic control
# counts? (e.g., you load in a standard tipranges file)
obs_target_species = 1
obs_all_species = 1
mean_frequency=1
dp=0.99
fdp=0.001
LnL_under_presence = calc_obs_like(truly_present=TRUE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_absence = calc_obs_like(truly_present=FALSE, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
LnL_under_presence
LnL_under_absence
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
  \code{\link{mapply_calc_post_prob_presence}},
  \code{\link{calc_post_prob_presence}},
  \code{\link{mapply_calc_obs_like}}
}

