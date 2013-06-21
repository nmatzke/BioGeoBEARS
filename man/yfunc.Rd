\name{yfunc}
\alias{yfunc}
\title{Extract the appropriate probability for a sympatric/range-copying speciation event, given text code for rangesize of smaller descendant, and ancestor}
\usage{
  yfunc(charcell, relprob_subsets_matrix)
}
\arguments{
  \item{charcell}{The text in the cell, indicating the type
  of speciation/cladogenesis range inheritance event.}

  \item{relprob_subsets_matrix}{A numeric matrix describing
  the relative probability of each smaller daughter range,
  conditional on the ancestral rangesize.}
}
\value{
  \code{prob_of_this_s}, a numeric value giving the
  relative probability of that descendent-ancestor
  rangesize pair.
}
\description{
  Extract the appropriate probability for a
  sympatric/range-copying speciation event, given text code
  for rangesize of smaller descendant, and ancestor
}
\note{
  Go BEARS!
}
\examples{
testval=1
# Examples

# Probabilities of different descendant rangesizes, for the smaller
# descendant, under sympatric/subset speciation
# (plus sympatric/range-copying, which is folded in):
relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6,
maxent_constraint_01=0.0001, NA_val=NA)
relprob_subsets_matrix
sfunc(charcell="s1_1", relprob_subsets_matrix)
sfunc(charcell="s1_2", relprob_subsets_matrix)
sfunc(charcell="s1_3", relprob_subsets_matrix)
sfunc(charcell="s2_3", relprob_subsets_matrix)

relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6,
maxent_constraint_01=0.5, NA_val=NA)
relprob_subsets_matrix
sfunc(charcell="s1_1", relprob_subsets_matrix)
sfunc(charcell="s1_2", relprob_subsets_matrix)
sfunc(charcell="s1_3", relprob_subsets_matrix)
sfunc(charcell="s2_3", relprob_subsets_matrix)

relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6,
maxent_constraint_01=0.9999, NA_val=NA)
relprob_subsets_matrix
sfunc(charcell="s1_1", relprob_subsets_matrix)
sfunc(charcell="s1_2", relprob_subsets_matrix)
sfunc(charcell="s1_3", relprob_subsets_matrix)
sfunc(charcell="s2_3", relprob_subsets_matrix)

relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6,
maxent_constraint_01=0.0001, NA_val=NA)
relprob_subsets_matrix
yfunc(charcell="y1", relprob_subsets_matrix)
yfunc(charcell="y2", relprob_subsets_matrix)
yfunc(charcell="y3", relprob_subsets_matrix)
yfunc(charcell="y4", relprob_subsets_matrix)

relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6,
maxent_constraint_01=0.5, NA_val=NA)
relprob_subsets_matrix
yfunc(charcell="y1", relprob_subsets_matrix)
yfunc(charcell="y2", relprob_subsets_matrix)
yfunc(charcell="y3", relprob_subsets_matrix)
yfunc(charcell="y4", relprob_subsets_matrix)

relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6,
maxent_constraint_01=0.9999, NA_val=NA)
relprob_subsets_matrix
yfunc(charcell="y1", relprob_subsets_matrix)
yfunc(charcell="y2", relprob_subsets_matrix)
yfunc(charcell="y3", relprob_subsets_matrix)
yfunc(charcell="y4", relprob_subsets_matrix)


# Probabilities of different descendant rangesizes, for the smaller descendant,
# under vicariant speciation
relprob_subsets_matrix = relative_probabilities_of_vicariants(max_numareas=6,
maxent_constraint_01v=0.0001, NA_val=NA)
relprob_subsets_matrix
vfunc(charcell="v1_1", relprob_subsets_matrix)
vfunc(charcell="v1_2", relprob_subsets_matrix)
vfunc(charcell="v1_3", relprob_subsets_matrix)
vfunc(charcell="v1_4", relprob_subsets_matrix)
vfunc(charcell="v2_4", relprob_subsets_matrix)
vfunc(charcell="v2_2", relprob_subsets_matrix)
vfunc(charcell="v1_6", relprob_subsets_matrix)
vfunc(charcell="v2_6", relprob_subsets_matrix)
vfunc(charcell="v3_6", relprob_subsets_matrix)

relprob_subsets_matrix = relative_probabilities_of_vicariants(max_numareas=6,
maxent_constraint_01v=0.5, NA_val=NA)
relprob_subsets_matrix
vfunc(charcell="v1_1", relprob_subsets_matrix)
vfunc(charcell="v1_2", relprob_subsets_matrix)
vfunc(charcell="v1_3", relprob_subsets_matrix)
vfunc(charcell="v1_4", relprob_subsets_matrix)
vfunc(charcell="v2_4", relprob_subsets_matrix)
vfunc(charcell="v2_2", relprob_subsets_matrix)
vfunc(charcell="v1_6", relprob_subsets_matrix)
vfunc(charcell="v2_6", relprob_subsets_matrix)
vfunc(charcell="v3_6", relprob_subsets_matrix)

relprob_subsets_matrix = relative_probabilities_of_vicariants(max_numareas=6,
maxent_constraint_01v=0.9999, NA_val=NA)
relprob_subsets_matrix
vfunc(charcell="v1_1", relprob_subsets_matrix)
vfunc(charcell="v1_2", relprob_subsets_matrix)
vfunc(charcell="v1_3", relprob_subsets_matrix)
vfunc(charcell="v1_4", relprob_subsets_matrix)
vfunc(charcell="v2_4", relprob_subsets_matrix)
vfunc(charcell="v2_2", relprob_subsets_matrix)
vfunc(charcell="v1_6", relprob_subsets_matrix)
vfunc(charcell="v2_6", relprob_subsets_matrix)
vfunc(charcell="v3_6", relprob_subsets_matrix)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
  \url{http://en.wikipedia.org/wiki/Maximum_entropy_probability_distribution}

  Matzke_2012_IBS

  Harte2011

  ReeSmith2008

  Ronquist1996_DIVA

  Ronquist_1997_DIVA

  Ronquist_Sanmartin_2011

  Landis_Matzke_etal_2013_BayArea
}
\seealso{
  \code{\link{sfunc}}, \code{\link{vfunc}},
  \code{\link{relative_probabilities_of_subsets}},
  \code{\link{symbolic_to_relprob_matrix_sp}},
  \code{\link{get_probvals}}, \code{\link[FD]{maxent}},
  \code{\link{calcZ_part}}, \code{\link{calcP_n}}
}

