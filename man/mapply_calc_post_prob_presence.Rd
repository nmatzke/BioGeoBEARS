\name{mapply_calc_post_prob_presence}
\alias{mapply_calc_post_prob_presence}
\title{Mapply version of calc_post_prob_presence()}
\usage{
  mapply_calc_post_prob_presence(prior_prob_presence = 0.01,
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

  \item{obs_target_species}{A scalar or
  column/vector/matrix of detection counts, e.g. as
  produced from the output from
  \code{\link{read_detections}}.}

  \item{obs_all_species}{A scalar or column/vector/matrix
  of detection counts, e.g. as produced from the output
  from \code{\link{read_controls}}.}

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
  \code{pp_df} A matrix of the posterior probability of
  presence, given the prior probability, the model
  parameters, and the data.
}
\description{
  This function applies
  \code{\link{calc_post_prob_presence}} to all cells of the
  input matrices \code{obs_target_species} and
  \code{obs_all_species}. These matrices obviously must
  have the same dimensions.
}
\details{
  The inputs are the same as for
  \code{\link{calc_post_prob_presence}}, except that
  \code{obs_target_species} and \code{obs_all_species} can
  be matrices.
}
\note{
  Go BEARS!
}
\examples{
# Calculate posterior probability of presence in an area,
# given a dp (detection probability) and detection model.

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


# Calculate data likelihoods, and posterior probability of presence=TRUE
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
  \code{\link{calc_post_prob_presence}},
  \code{\link{mapply_calc_obs_like}}
  \code{\link{Pdata_given_rangerow}},
  \code{\link[base]{mapply}},
  \code{\link{tiplikes_wDetectionModel}}
}

