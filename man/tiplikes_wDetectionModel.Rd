\name{tiplikes_wDetectionModel}
\alias{tiplikes_wDetectionModel}
\title{Calculate probability of detection data for each OTU at each range in a list of states/geographic ranges}
\usage{
  tiplikes_wDetectionModel(states_list_0based_index,
    numareas = NULL, detects_df, controls_df,
    mean_frequency = 0.1, dp = 1, fdp = 0,
    null_range_gets_0_like = TRUE)
}
\arguments{
  \item{states_list_0based_index}{A states_list, 0-based,
  e.g. from
  \code{\link[cladoRcpp]{rcpp_areas_list_to_states_list}}.}

  \item{numareas}{The number of areas being considered in
  the analysis. If \code{NULL} (default), this is
  calculated to be the maximum range length, or one plus
  the maximum 0-based index in any of the ranges.}

  \item{detects_df}{A matrix/data.frame of detection
  counts, as produced from the output from
  \code{\link{read_detections}}.}

  \item{controls_df}{A matrix/data.frame of detection
  counts, as produced from the output from
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

  \item{null_range_gets_0_like}{If \code{TRUE} (default),
  then the data is given zero probability on the hypothesis
  that the range is a null range (i.e., no areas occupied).
  This is equivalent to saying that you are sure/are
  willing to assume that the OTU exists somewhere in your
  study area, at the timepoint being considered.  Null
  ranges are identified by length=1, containing NULL, NA,
  "", "_", etc.}
}
\value{
  \code{tip_condlikes_of_data_on_each_state} The
  (non-logged!) likelihood of the data for each tip, given
  each possible range, and the detection model parameters.
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

# soft-coded input files
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
detects_fn = np(paste(extdata_dir, "/Psychotria_detections_v1.txt", sep=""))
controls_fn = np(paste(extdata_dir, "/Psychotria_controls_v1.txt", sep=""))

detects_df = read_detections(detects_fn, OTUnames=NULL, areanames=NULL, tmpskip=0)
controls_df = read_controls(controls_fn, OTUnames=NULL, areanames=NULL, tmpskip=0)

# Calculate the likelihood of the data at each tip, for each possible geographic range
numareas = 4
tmpranges = list(c(0), c(1), c(0,1))

mean_frequency=0.1
dp=1
fdp=0

tip_condlikes_of_data_on_each_state =
tiplikes_wDetectionModel(states_list_0based_index=tmpranges, numareas=numareas,
detects_df, controls_df, mean_frequency=mean_frequency, dp=dp, fdp=fdp,
null_range_gets_0_like=TRUE)

tip_condlikes_of_data_on_each_state
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
  \code{\link{Pdata_given_rangerow}},
  \code{\link{calc_obs_like}}, \code{\link[base]{mapply}},
  \code{\link{read_detections}},
  \code{\link{read_controls}}
}

