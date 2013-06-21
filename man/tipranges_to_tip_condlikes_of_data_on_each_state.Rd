\name{tipranges_to_tip_condlikes_of_data_on_each_state}
\alias{tipranges_to_tip_condlikes_of_data_on_each_state}
\title{Convert a tipranges object to the tip likelihoods}
\usage{
  tipranges_to_tip_condlikes_of_data_on_each_state(tipranges,
    phy, states_list = NULL,
    maxareas = length(getareas_from_tipranges_object(tipranges)))
}
\arguments{
  \item{tipranges}{An object of class \code{tipranges}.}

  \item{phy}{A phylogenetic tree (\code{\link[ape]{ape}}
  object of class \code{\link[ape]{phylo}})}

  \item{states_list}{A complete list of the different
  states, of class \code{\link{list}} form}

  \item{maxareas}{The maximum number of areas in a
  geographic range, if the user does}
}
\value{
  \code{tip_condlikes_of_data_on_each_state} For each
  tip/row, likelihood of that tip's data under each
  possible true geographic range (columns)
}
\description{
  This function takes a tipranges object, and converts it
  to tip likelihoods for input into the likelihood
  calculations of \code{\link{calc_loglike_sp}}.
}
\details{
  This (like LAGRANGE (\cite{Ree et al. (2008)}) and every
  other available program) assumes that the geographic
  ranges at the tips are known with certainty.  Reality may
  be different, particularly for sparsely-studied, scarce,
  or fossil taxa.  In such a case, a detection model is
  needed to specify the likelihood of the observation data
  under each possible geographic range at the tips.

  Note that data likelihoods under this or that hypothesis
  are not the same thing as probabilities.  E.g., with DNA,
  if sequencing machine says that the base could be either
  A or C, but not G or T, then the likelihood of the data
  for that nucleotide position for that species would be 1
  1 0 0, not 0.5 0.5 0 0.  See Felsenstein (2004), p. 255,
  for more.
}
\note{
  Go BEARS!
}
\examples{
testval=1
# Define a tipranges object
tipranges_object = define_tipranges_object()
tipranges_object

areanames = getareas_from_tipranges_object(tipranges_object)
areanames

# Specify phylogeny to go with default tipranges object
newick_str = "((tip1:1,tip2:1):1,tip3:2):1;"
phy = read.tree(file="", text=newick_str)

# Here, we will assume the maximum range size is all areas, but it could be smaller
maxareas = length(areanames)
\dontrun{
states_list = areas_list_to_states_list_old(areas=areanames, include_null_range=TRUE,
maxareas=maxareas)
states_list
}

states_list = areas_list_to_states_list_new(areas=areanames, include_null_range=TRUE,
maxareas=maxareas)
states_list

tip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(
tipranges=tipranges_object, phy=phy, states_list=states_list, maxareas=maxareas )
tip_condlikes_of_data_on_each_state
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  ReeSmith2008

  Matzke_2012_IBS

  Felsenstein2004
}
\seealso{
  \code{\link{define_tipranges_object}},
  \code{\link{getareas_from_tipranges_object}},
  \code{\link{areas_list_to_states_list_new}},
  \code{\link[cladoRcpp]{areas_list_to_states_list_old}},
  \code{\link{binary_ranges_to_letter_codes}}
}

