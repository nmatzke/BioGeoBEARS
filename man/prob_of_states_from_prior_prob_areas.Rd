\name{prob_of_states_from_prior_prob_areas}
\alias{prob_of_states_from_prior_prob_areas}
\title{Calculate probability of detection data for each OTU at each range in a list of states/geographic ranges}
\usage{
  prob_of_states_from_prior_prob_areas(states_list_0based_index,
    numareas = NULL, prior_prob_presence = 0.01,
    null_range_gets_0_prob = TRUE, normalize_probs = TRUE)
}
\arguments{
  \item{states_list_0based_index}{A states_list, 0-based,
  e.g. from
  \code{\link[cladoRcpp]{rcpp_areas_list_to_states_list}}.}

  \item{numareas}{The number of areas being considered in
  the analysis. If \code{NULL} (default), this is
  calculated to be the maximum range length, or one plus
  the maximum 0-based index in any of the ranges.}

  \item{prior_prob_presence}{The prior probability of
  presence, i.e. when no detection or taphonomic control
  data whatsoever is available.  Default is set to 0.01
  which expresses my totally uninformed bias that in
  whatever your data is, your species of interest probably
  doesn't live in the typical area you are looking at.}

  \item{null_range_gets_0_prob}{If \code{TRUE} (default),
  then the null range is given zero probability. A null
  range has no areas occupied.  This is equivalent to
  saying that you are sure/are willing to assume that the
  OTU exists somewhere in your study area, at the timepoint
  being considered.  Null ranges are identified by
  length=1, containing NULL, NA, "", "_", etc.}

  \item{normalize_probs}{If \code{TRUE}, the probabilities
  of each range will be normalized so that they sum to 1.
  Otherwise, they won't.}
}
\value{
  \code{prob_of_each_range} The probability of each range,
  given the prior probability of presence in each area.
}
\description{
  This function calculates P(data|range,dp), i.e. the
  probability of some detection and taphonomic control
  counts, given the true geographic range/state, and
  \code{dp}, a detection probability (and, optionally, a
  false detection probability, \code{fdp}).
}
\details{
  This function performs the operation for all
  states/ranges for all tips.

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

prior_prob_presence = 0.01

areas = c("K", "O", "M", "H")
numareas = length(areas)
states_list_0based_index =
rcpp_areas_list_to_states_list(areas=areas, maxareas=4, include_null_range=TRUE)
states_list_0based_index
numareas = 4

mean_frequency=0.1
dp=1
fdp=0

prob_of_states_from_prior_prob_areas(states_list_0based_index, numareas=numareas,
prior_prob_presence=prior_prob_presence, null_range_gets_0_prob=TRUE,
normalize_probs=TRUE)

prob_of_states_from_prior_prob_areas(states_list_0based_index, numareas=numareas,
prior_prob_presence=prior_prob_presence, null_range_gets_0_prob=TRUE,
normalize_probs=FALSE)

prob_of_states_from_prior_prob_areas(states_list_0based_index, numareas=numareas,
prior_prob_presence=prior_prob_presence, null_range_gets_0_prob=FALSE,
normalize_probs=TRUE)

prob_of_states_from_prior_prob_areas(states_list_0based_index, numareas=numareas,
prior_prob_presence=prior_prob_presence, null_range_gets_0_prob=FALSE,
normalize_probs=FALSE)
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
  \code{\link[cladoRcpp]{rcpp_areas_list_to_states_list}},
  \code{\link{Pdata_given_rangerow}},
  \code{\link{calc_obs_like}}, \code{\link[base]{mapply}},
  \code{\link{read_detections}},
  \code{\link{read_controls}}
}

