\name{calc_post_prob_presence}
\alias{calc_post_prob_presence}
\title{Calculate posterior probability of presence, given count data and parameters}
\usage{
  calc_post_prob_presence(prior_prob_presence = 0.01,
    obs_target_species, obs_all_species,
    mean_frequency = 0.1, dp = 1, fdp = 0,
    print_progress = "")
}
\arguments{
  \item{prior_prob_presence}{The prior probability of
  presence, i.e. when no detection or taphonomic control
  data whatsoever is available.  Default is set to 0.01
  which expresses my totally uninformed bias that in
  whatever your data is, your species of interest probably
  doesn't live in the typical area you are looking at.}

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

  \item{print_progress}{If not the default (\code{""}),
  print whatever is in print_progress, followed by a space
  (for error checking/surveying results).}
}
\value{
  \code{post_prob} The posterior probability of presence,
  given the prior probability, the model parameters, and
  the data.
}
\description{
  This function calculates P(presence|count
  data,parameters), i.e. the posterior probability of
  presence in an area, given data on detection counts and
  taphonomic control counts, and a detection model with the
  parameters mean_frequency, \code{dp}, a detection
  probability (and, optionally, a false detection
  probability, \code{fdp}).
}
\details{
  Essentially, this function combines a prior probability,
  with the likelihood function (coded in
  \code{\link{calc_obs_like}}) to produce a posterior
  probability of presence given Bayes' Theorem (Bayes &
  Price, 1763).

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
# Calculate posterior probability of presence in an area,
# given a dp (detection probability) and detection model.

# With zero error rate
obs_target_species = 10
obs_all_species = 100
mean_frequency=0.1
dp=1
fdp=0
prior_prob_presence = 0.01
post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
post_prob
# i.e., with perfect detection, the prob. of presence is 1 under the
# hypothesis of presence, and 0 under the hypothesis of
# (This is because the likelihood of the data under
# presence and absence are ln(prob) = 0 & -Inf, respectively.)


# Note that with very high error rates, your conclusion could reverse
obs_target_species = 10
obs_all_species = 100
mean_frequency=0.1
dp=0.5
fdp=0.3
prior_prob_presence = 0.01
post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
post_prob


# With 0 error rate, even 1 observation makes P(presence) = 1
obs_target_species = 1
obs_all_species = 100
mean_frequency=0.1
dp=1
fdp=0
prior_prob_presence = 0.01
post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
post_prob


# With a small error rate, there is some small but positive probability of
# falsely getting 10 detections; but it may be effectively 0
obs_target_species = 10
obs_all_species = 100
mean_frequency=0.1
dp=0.99
fdp=0.001
prior_prob_presence = 0.01
post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
post_prob


# If you have only 1 detection, and you have 100 taphonomic controls and
# a mean_frequency of sampling the OTU of interest of 0.1, then there is
# still a very low probability of presence (since, under your model,
# you should expect to see about 10 detections, not 1)
obs_all_species = 100
mean_frequency=0.1
dp=0.99
fdp=0.001
prior_prob_presence = 0.01

obs_target_species = 0
calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)

obs_target_species = 1
calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)

obs_target_species = 2
calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)

obs_target_species = 3
calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)



# Note how quickly this chances if you drop the mean_frequency from 0.1
# to 0.01. This means that if you want single detections to count for
# a lot, you need either a low mean_frequency which matches the observed
# frequency, or an extremely high/perfect detection probability (dp).
obs_all_species = 100
mean_frequency=0.01
dp=0.99
fdp=0.001
prior_prob_presence = 0.01

obs_target_species = 0
calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)

obs_target_species = 1
calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)

obs_target_species = 2
calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)

obs_target_species = 3
calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)


# Changing mean_frequency from 0.01 to 0.001 actually LOWERS the posterior
# probability of presence based on 1 detection, as we have a somewhat
# significant false detection rate:
obs_target_species = 1
obs_all_species = 100
mean_frequency=0.001
dp=0.99
fdp=0.001
prior_prob_presence = 0.01

obs_target_species = 0
calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)

obs_target_species = 1
calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)

obs_target_species = 2
calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)

obs_target_species = 3
calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)


# Change false detection probability to a much lower value
obs_all_species = 100
mean_frequency=0.001
dp=0.99
fdp=0.00001
prior_prob_presence = 0.01

obs_target_species = 0
calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)

obs_target_species = 1
calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)

obs_target_species = 2
calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)

obs_target_species = 3
calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)




# Change false detection probability to 0
obs_all_species = 100
mean_frequency=0.001
dp=0.99
fdp=0.0
prior_prob_presence = 0.01


obs_target_species = 0
calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)

obs_target_species = 1
calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)

obs_target_species = 2
calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)

obs_target_species = 3
calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)









# Change mean_frequency to 0.001
obs_all_species = 100
mean_frequency=0.001
dp=0.99
fdp=0.0
prior_prob_presence = 0.01

obs_target_species = 0
calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)

obs_target_species = 1
calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)

obs_target_species = 2
calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)

obs_target_species = 3
calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)







# Example #2 -- what if you have ZERO detections, but lots of detections
# of your taphonomic control?
obs_target_species = 0
obs_all_species = 100
mean_frequency=0.1
dp=1
fdp=0
prior_prob_presence = 0.01
post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
post_prob

# With a slight error rate
obs_target_species = 0
obs_all_species = 100
mean_frequency=0.1
dp=0.99
fdp=0.001
prior_prob_presence = 0.01
post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
post_prob


obs_target_species = 0
obs_all_species = 2
mean_frequency=0.1
dp=1
fdp=0
prior_prob_presence = 0.01
post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
post_prob

# With a slight error rate
obs_target_species = 0
obs_all_species = 2
mean_frequency=0.1
dp=0.99
fdp=0.001
prior_prob_presence = 0.01
post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
post_prob





# Example #3 -- what if you have ZERO detections, but only a few
# detections of your taphonomic control?
obs_target_species = 0
obs_all_species = 1
mean_frequency=0.1
dp=1
fdp=0
prior_prob_presence = 0.01
post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
post_prob

# With a slight error rate
obs_target_species = 0
obs_all_species = 1
mean_frequency=0.1
dp=0.99
fdp=0.001
prior_prob_presence = 0.01
post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
post_prob



# Special cases -- e.g., no data
# Prob(data)=1, ln(prob)=0
obs_target_species = 0
obs_all_species = 0
mean_frequency=0.1
dp=0.99
fdp=0.001
prior_prob_presence = 0.01
post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
post_prob

obs_target_species = 0
obs_all_species = 0
mean_frequency=0.1
dp=1
fdp=0
prior_prob_presence = 0.01
post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
post_prob


# What if, for some reason, you put in identical detections and taphonomic control
# counts? (e.g., you load in a standard tipranges file)
obs_target_species = 1
obs_all_species = 1
mean_frequency=1
dp=1
fdp=0
prior_prob_presence = 0.01
post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
post_prob

# What if, for some reason, you put in identical detections and taphonomic control
# counts? (e.g., you load in a standard tipranges file)
obs_target_species = 1
obs_all_species = 1
mean_frequency=1
dp=0.99
fdp=0.001
prior_prob_presence = 0.01
post_prob = calc_post_prob_presence(prior_prob_presence, obs_target_species,
obs_all_species, mean_frequency, dp, fdp)
post_prob
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
  \url{http://en.wikipedia.org/wiki/Bayes'_theorem}

  Matzke_2012_IBS

  Bottjer_Jablonski_1988

  Bayes_1763
}
\seealso{
  \code{\link{calc_obs_like}},
  \code{\link{mapply_calc_post_prob_presence}},
  \code{\link{mapply_calc_obs_like}}
}

