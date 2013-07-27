\name{post_prob_states}
\alias{post_prob_states}
\title{Calculate posterior probability of each states/geographic ranges, given prior probabilities and data likelihoods}
\usage{
  post_prob_states(prob_of_each_range,
    condlikes_of_data_on_each_range)
}
\arguments{
  \item{prob_of_each_range}{The probability of each range,
  given the prior probability of presence in each area.}

  \item{condlikes_of_data_on_each_range}{The probability of
  the data, conditional on each range (i.e., the
  likelihood), as found in e.g. a row of the output from
  \code{\link{tiplikes_wDetectionModel}}.}
}
\value{
  \code{posterior_probs} The posterior probability of each
  range.
}
\description{
  This function calculates P(range|data,detection model),
  i.e. the probability of each possible range, given a
  prior probability of each range, and the likelihood of
  each range.
}
\details{
  The prior probability of each range should be considered
  by the user.  Note that putting the same prior on the
  probability of occurrence in each individual range does
  NOT mean a flat prior on each state/geographic range.
  This fact is demonstrated in the function
  \code{\link{prob_of_states_from_prior_prob_areas}}.
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

detects_df = read_detections(detects_fn, OTUnames=NULL, areanames=NULL, tmpskip=0)
controls_df = read_controls(controls_fn, OTUnames=NULL, areanames=NULL, tmpskip=0)


# Setup
prior_prob_presence = 0.01
areas = c("K", "O", "M", "H")
numareas = length(areas)
maxareas = length(areas)
states_list_0based_index =
rcpp_areas_list_to_states_list(areas=areas, maxareas=maxareas,
                               include_null_range=TRUE)
states_list_0based_index

mean_frequency=0.1
dp=1
fdp=0

tip_condlikes_of_data_on_each_state =
tiplikes_wDetectionModel(states_list_0based_index, numareas=numareas,
detects_df, controls_df, mean_frequency=mean_frequency, dp=dp, fdp=fdp,
null_range_gets_0_like=TRUE)

tip_condlikes_of_data_on_each_state



# To get denominator, just iterate over all the states
# Prior probability
prob_of_each_range = prob_of_states_from_prior_prob_areas(states_list_0based_index,
numareas=numareas,
prior_prob_presence=prior_prob_presence, null_range_gets_0_prob=TRUE,
normalize_probs=TRUE)

# Likelihoods of the data on each range
condlikes_of_data_on_each_range = tip_condlikes_of_data_on_each_state[1,]

posterior_probs = post_prob_states(prob_of_each_range,
                  condlikes_of_data_on_each_range)
posterior_probs

# Should sum to 1
sum(posterior_probs)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
  \url{http://en.wikipedia.org/wiki/Log_probability}

  Matzke_2012_IBS

  Bottjer_Jablonski_1988
}
\seealso{
  \code{\link{prob_of_states_from_prior_prob_areas}},
  \code{\link{tiplikes_wDetectionModel}},
  \code{\link[cladoRcpp]{rcpp_areas_list_to_states_list}},
  \code{\link{Pdata_given_rangerow}},
  \code{\link{calc_obs_like}}, \code{\link[base]{mapply}},
  \code{\link{read_detections}},
  \code{\link{read_controls}}
}

