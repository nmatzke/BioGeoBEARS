\name{relative_probabilities_of_subsets}
\alias{relative_probabilities_of_subsets}
\title{Calculate probability of different descendant rangesizes, for the smaller descendant, in subset speciation}
\usage{
  relative_probabilities_of_subsets(max_numareas = 6,
    maxent_constraint_01 = 0.5, NA_val = NA)
}
\arguments{
  \item{max_numareas}{The maximum number of areas possible
  allowed for the smaller-ranged-daughter in this type of
  cladogenesis/speciation.}

  \item{maxent_constraint_01}{The parameter describing the
  probability distribution on descendant rangesizes for the
  smaller descendant. See above.}

  \item{NA_val}{The output matrix consists of ancestral
  rangesizes and rangesizes of the smaller descendant.
  Some values are disallowed -- e.g. descendant ranges
  larger than the ancestor; or, in subset speciation,
  descendant ranges the same size as the ancestor are
  disallowed.  All disallowed descendant rangesizes get
  \code{NA_val}.}
}
\value{
  \code{\link{relative_probabilities_of_vicariants}},
  \code{relprob_subsets_matrix}, a numeric matrix giving
  the relative probability of each rangesize for the
  smaller descendant of an ancestral range, conditional on
  the ancestral rangesize.
}
\description{
  "Rangesize" here means "number of areas in a geographic
  range".  The \code{LAGRANGE} cladogenesis model requires
  that, during cladogenesis events, one daughter lineage
  will ALWAYS have a geographic range of size 1.  This is
  argued for in \cite{Ree et al. (2008)} on the grounds
  that new species usually get isolated and start in a new
  area.  This is a reasonable proposition, but still, it
  would be nice to test the assumption.  In addition, it
  could be that some speciation modes, especially
  vicariance, obey different rules.  E.g., \code{DIVA}
  (\cite{Ronquist (1996)}, \cite{Ronquist (1997)}) allows
  vicariant speciation to divide up the ancestral range in
  every possible way (e.g., ABCD-->AB|CD, or AC|BD, or
  A|BCD, or D|ABC, etc.), but \code{LAGRANGE} would only
  allow vicariance to split off areas of size 1:
  (ABCD-->A|BCD, B|ACD, etc.) (Ronquist_Sanmartin_2011).
}
\details{
  To test different models, the user has to have control of
  the relative probability of different descendant
  rangesizes.  The probability of each descendant rangesize
  could be parameterized individually, but we have a
  limited amount of observational data (essentially one
  character), so efficient parameterizations should be
  sought.

  One way to do this is with the Maximum Entropy
  (\cite{Harte (2011)}) discrete probability distribution
  of a number of ordered states.  Normally this is applied
  (in examples) to the problem of estimation of the
  relative probability of the different faces of a 6-sided
  die.  The input "knowledge" is the true mean of the dice
  rolls.  If the mean value is 3.5, then each face of the
  die will have probability 1/6.  If the mean value is
  close to 1, then the die is severely skewed such that the
  probability of rolling 1 is 99%+ and the probability of
  other die rolls is very small.  If the mean value is
  close to 6, then the probability distribution is skewed
  towards higher numbers.

  Here in \code{BioGeoBEARS}, we use the same Maximum
  Entropy function to specify the relative probability of
  geographic ranges of a number of different rangesizes.
  This is merely used so that a single parameter can
  control the probability distribution -- there is no
  MaxEnt estimation going on here.  The user specifies a
  value for the parameter \code{maxent_constraint_01}
  between 0.0001 and 0.9999.  This can then be applied to
  all of the different ancestor-descendant range
  combinations in the cladogenesis/speciation matrix.

  Example values of \code{maxent_constraint_01} would give
  the following results:

  \code{maxent_constraint_01 = 0.0001} -- The smaller
  descendant has rangesize 1 with 100% probability (as in
  \code{LAGRANGE})\cr \code{maxent_constraint_01 = 0.5} --
  The smaller descendant can be any rangesize equal
  probability. This is effectively what happens in
  \code{DIVA}'s version of vicariance speciation\cr
  \code{maxent_constraint_01 = 0.9999} -- The smaller
  descendant will take the largest possible rangesize for a
  given type of speciation, and a given ancestral
  rangesize.  E.g., for sympatric/range-copying speciation
  (the ancestor is simply copied to both descendants, as in
  a continuous-time model with no cladogenesis effect), an
  ancestor of size 3 would product two descendant lineages
  of size 3. Such a model is implemented in the program
  \code{BayArea} (\cite{Landis et al. (2013)}).
  \code{LAGRANGE}, on the other hand, would only allow
  range-copying for ancestral ranges of size 1.\cr

  \bold{Note:} In \code{LAGRANGE}-type models, at
  speciation/cladogenesis events, one descendant daughter
  branch ALWAYS has size 1, whereas the other descendant
  daughter branch either (a) is the same (in
  sympatric/range-copying speciation), (b) inherits the
  complete ancestral range (in sympatric/subset speciation)
  or (c) inherits the remainder of the range (in
  vicariant/range-division speciation).  LAGRANGE-type
  behavior (the smaller descendant has rangesize 1 with
  100% probability, and 0% probability for any other
  rangesize) can be achieved by setting the
  \code{maxent_constraint_01} parameter to 0.0001.

  \code{See also:} Maximum Entropy probability distribution
  for discrete variable with given mean (and discrete
  uniform flat prior)
  \url{http://en.wikipedia.org/wiki/Maximum_entropy_probability_distribution}

  Currently, the function \code{\link[FD]{maxent}} from the
  \code{\link[FD]{FD}} package is used to get the discrete
  probability distribution, given the number of states and
  the \code{maxent_constraint_01} parameter.  This could
  also be done with \code{\link{get_probvals}}, which uses
  \code{\link{calcZ_part}}, \code{\link{calcP_n}},
  following equations 6.3-6.4 of \cite{Harte (2011)},
  although this is not yet implemented.
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
relative_probabilities_of_subsets(max_numareas=6, maxent_constraint_01=0.0001,
NA_val=NA)
relative_probabilities_of_subsets(max_numareas=6, maxent_constraint_01=0.5,
NA_val=NA)
relative_probabilities_of_subsets(max_numareas=6, maxent_constraint_01=0.9999,
NA_val=NA)

# Probabilities of different descendant rangesizes, for the smaller descendant,
# under vicariant speciation
relative_probabilities_of_vicariants(max_numareas=6, maxent_constraint_01v=0.0001,
NA_val=NA)
relative_probabilities_of_vicariants(max_numareas=6, maxent_constraint_01v=0.5,
NA_val=NA)
relative_probabilities_of_vicariants(max_numareas=6, maxent_constraint_01v=0.9999,
NA_val=NA)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
  \url{http://en.wikipedia.org/wiki/Maximum_entropy_probability_distribution}

  ReeSmith2008

  Ronquist1996_DIVA

  Ronquist_1997_DIVA

  Harte2011

  Landis_Matzke_etal_2013_BayArea

  Matzke_2012_IBS

  Ronquist_Sanmartin_2011
}
\seealso{
  \code{\link{symbolic_to_relprob_matrix_sp}},
  \code{\link{get_probvals}}, \code{\link[FD]{maxent}},
  \code{\link{calcZ_part}}, \code{\link{calcP_n}}
}

