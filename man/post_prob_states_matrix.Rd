\name{post_prob_states_matrix}
\alias{post_prob_states_matrix}
\title{Calculate posterior probability of each states/geographic ranges, given prior probabilities and data likelihoods}
\usage{
  post_prob_states_matrix(prob_of_each_range,
    tip_condlikes_of_data_on_each_state)
}
\arguments{
  \item{prob_of_each_range}{The probability of each range,
  given the prior probability of presence in each area.}

  \item{tip_condlikes_of_data_on_each_state}{The
  probability of the data, conditional on each range (i.e.,
  the likelihood), as found in e.g. a row of the output
  from \code{\link{tiplikes_wDetectionModel}}.}
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
#extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"
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

posterior_probs_matrix = post_prob_states_matrix(prob_of_each_range,
                  tip_condlikes_of_data_on_each_state)
posterior_probs_matrix

# Should sum to 1
rowSums(posterior_probs_matrix)

# How does posterior probability correlate with likelihood and prior probability?
par(mfrow=c(1,2))
plot(x=jitter(log(tip_condlikes_of_data_on_each_state)),
y=jitter(log(posterior_probs_matrix)))
title("Correlation of data likelihoods\\nand posterior probabilities")

prob_of_each_range_matrix = matrix(data=prob_of_each_range,
nrow=nrow(posterior_probs_matrix), ncol=length(prob_of_each_range))
plot(x=jitter(log(prob_of_each_range_matrix)),
y=jitter(log(posterior_probs_matrix)))
title("Correlation of prior probability\\nand posterior probabilities")
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

