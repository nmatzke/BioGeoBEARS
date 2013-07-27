\name{Pdata_given_rangerow_dp}
\alias{Pdata_given_rangerow_dp}
\title{Calculate probability of detection data given a true geographic range and a detection probability}
\usage{
  Pdata_given_rangerow_dp(truerange_areas, numareas,
    detects_df, controls_df, mean_frequency = 0.1, dp = 1,
    fdp = 0)
}
\arguments{
  \item{truerange_areas}{The list of areas (as 0-based
  numbers) in this geographic range/state.}

  \item{numareas}{The function needs to know the total
  number of areas in the analysis.}

  \item{detects_df}{A column/vector of detection counts, as
  produced from a column of the output from
  \code{\link{read_detections}}.}

  \item{controls_df}{A column/vector of detection counts,
  as produced from a column of the output from
  \code{\link{read_controls}}.}

  \item{dp}{The detection probability.  This is the
  per-sample probability that you will detect the OTU in
  question. In other words, the model assumes that each
  specimen from the taphonomic control group has a chance
  of being a representative of the OTU you are looking for.
  The default is 1, which assumes perfect detection, i.e.
  the assumption being made in all historical biogeography
  analyses that do not take into account detection
  probability.  A value of 1 will only work when the
  taphonomic control count equals the detection count; any
  other data would have likelihood=0.}

  \item{fdp}{The false detection probability.  This is
  probability of falsely concluding a detection occurred,
  when in fact the specimen was of something else.  The
  default is 0, which assumes zero error rate, i.e. the
  assumption being made in all historical biogeography
  analyses that do not take into account detection
  probability.  This option is being included for
  completeness, but it may not be wise to try to infer both
  \code{dp} and \code{fdp} at once due to identifiability
  issues (and estimation of fdp may take a very large
  amount of data).}

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
}
\value{
  \code{dtf}
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
  \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}}
}

