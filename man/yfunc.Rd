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

# Probabilities of different descendant rangesizes, for the smaller descendant, under sympatric/subset speciation
# (plus sympatric/range-copying, which is folded in):
relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, maxent_constraint_01=0.0001, NA_val=NA)
relprob_subsets_matrix
sfunc(charcell="s1_1", relprob_subsets_matrix)
sfunc(charcell="s1_2", relprob_subsets_matrix)
sfunc(charcell="s1_3", relprob_subsets_matrix)
sfunc(charcell="s2_3", relprob_subsets_matrix)

relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, maxent_constraint_01=0.5, NA_val=NA)
relprob_subsets_matrix
sfunc(charcell="s1_1", relprob_subsets_matrix)
sfunc(charcell="s1_2", relprob_subsets_matrix)
sfunc(charcell="s1_3", relprob_subsets_matrix)
sfunc(charcell="s2_3", relprob_subsets_matrix)

relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, maxent_constraint_01=0.9999, NA_val=NA)
relprob_subsets_matrix
sfunc(charcell="s1_1", relprob_subsets_matrix)
sfunc(charcell="s1_2", relprob_subsets_matrix)
sfunc(charcell="s1_3", relprob_subsets_matrix)
sfunc(charcell="s2_3", relprob_subsets_matrix)

relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, maxent_constraint_01=0.0001, NA_val=NA)
relprob_subsets_matrix
yfunc(charcell="y1", relprob_subsets_matrix)
yfunc(charcell="y2", relprob_subsets_matrix)
yfunc(charcell="y3", relprob_subsets_matrix)
yfunc(charcell="y4", relprob_subsets_matrix)

relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, maxent_constraint_01=0.5, NA_val=NA)
relprob_subsets_matrix
yfunc(charcell="y1", relprob_subsets_matrix)
yfunc(charcell="y2", relprob_subsets_matrix)
yfunc(charcell="y3", relprob_subsets_matrix)
yfunc(charcell="y4", relprob_subsets_matrix)

relprob_subsets_matrix = relative_probabilities_of_subsets(max_numareas=6, maxent_constraint_01=0.9999, NA_val=NA)
relprob_subsets_matrix
yfunc(charcell="y1", relprob_subsets_matrix)
yfunc(charcell="y2", relprob_subsets_matrix)
yfunc(charcell="y3", relprob_subsets_matrix)
yfunc(charcell="y4", relprob_subsets_matrix)


# Probabilities of different descendant rangesizes, for the smaller descendant, under vicariant speciation
relprob_subsets_matrix = relative_probabilities_of_vicariants(max_numareas=6, maxent_constraint_01v=0.0001, NA_val=NA)
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

relprob_subsets_matrix = relative_probabilities_of_vicariants(max_numareas=6, maxent_constraint_01v=0.5, NA_val=NA)
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

relprob_subsets_matrix = relative_probabilities_of_vicariants(max_numareas=6, maxent_constraint_01v=0.9999, NA_val=NA)
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

  Matzke N (2012). "Founder-event speciation in BioGeoBEARS
  package dramatically improves likelihoods and alters
  parameter inference in Dispersal-Extinction-Cladogenesis
  (DEC) analyses." _Frontiers of Biogeography_, *4*(suppl.
  1), pp. 210. ISSN 1948-6596, Poster abstract published in
  the Conference Program and Abstracts of the International
  Biogeography Society 6th Biannual Meeting, Miami,
  Florida. Poster Session P10: Historical and
  Paleo-Biogeography. Poster 129B. January 11, 2013, <URL:
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}>.

  Harte J (2011). _Maximum entropy and ecology : a theory
  of abundance, distribution, and energetics_, series
  Oxford series in ecology and evolution. Oxford University
  Press, Oxford ; New York. ISBN 9780199593415 0199593418
  9780199593422 (pbk.) 0199593426 (pbk.)/ISSN, 2011284493
  John Harte.. ill. ; 24 cm. Includes bibliographical
  references (p. [244]-252) and index. Oxford series in
  ecology and evolution., <URL:
  http://socrates.berkeley.edu/~hartelab/MaxEnt.html>.

  Ree RH and Smith SA (2008). "Maximum likelihood inference
  of geographic range evolution by dispersal, local
  extinction, and cladogenesis." _Systematic Biology_,
  *57*(1), pp. 4-14. <URL:
  http://dx.doi.org/10.1080/10635150701883881>, <URL:
  http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&dopt=Citation&list_uids=18253896>.

  Ronquist F (1996). _DIVA version 1.1. Computer program
  and manual (online)._. Uppsala University. <URL:
  http://www.ebc.uu.se/systzoo/research/diva/manual/dmanual.html>.

  Ronquist F (1997). "Dispersal-Vicariance Analysis: A New
  Approach to the Quantification of Historical
  Biogeography." _Systematic Biology_, *46*(1), pp.
  195-203. <URL: http://www.jstor.org/stable/2413643>.

  Ronquist F and Sanmartin I (2011). "Phylogenetic Methods
  in Biogeography." _Annual Review of Ecology, Evolution,
  and Systematics_, *42*(1), pp. 441-464. <URL:
  http://dx.doi.org/10.1146/annurev-ecolsys-102209-144710>,
  <URL:
  http://www.annualreviews.org/doi/abs/10.1146/annurev-ecolsys-102209-144710>.

  Landis M, Matzke NJ, Moore BR and Huelsenbeck JP (2013).
  "Statistical Analysis of Biogeography when the Number of
  Areas is Large." _Systematic Biology_, *submitted*.
}
\seealso{
  \code{\link{sfunc}}, \code{\link{vfunc}},
  \code{\link{relative_probabilities_of_subsets}},
  \code{\link{symbolic_to_relprob_matrix_sp}},
  \code{\link{get_probvals}}, \code{\link[FD]{maxent}},
  \code{\link{calcZ_part}}, \code{\link{calcP_n}}
}

