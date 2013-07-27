\name{calc_loglike_sp}
\alias{calc_loglike_sp}
\title{Calculate log-likelihood with a transition matrix and speciation events -- byte-compiled}
\usage{
  calc_loglike_sp(tip_condlikes_of_data_on_each_state, phy,
    Qmat, spPmat = NULL, min_branchlength = 1e-21,
    return_what = "loglike",
    probs_of_states_at_root = NULL, rootedge = FALSE,
    sparse = FALSE, printlevel = 1, use_cpp = TRUE,
    input_is_COO = FALSE, spPmat_inputs = NULL,
    cppSpMethod = 3, cluster_already_open = NULL,
    calc_ancprobs = FALSE, null_range_allowed = TRUE,
    fixnode = NULL, fixlikes = NULL, stratified = FALSE,
    states_allowed_TF = NULL)
}
\arguments{
  \item{tip_condlikes_of_data_on_each_state}{A numeric
  matrix with rows representing tips, and columns
  representing states/geographic ranges.  The cells give
  the likelihood of the observation data under the
  assumption that the tip has that state; typically this
  means that the known geographic range gets a '1' and all
  other states get a 0.}

  \item{phy}{A phylogeny object.  The function converts it
  to pruningwise order.}

  \item{Qmat}{A Q transition matrix representing the
  along-branch model for the evolution of geographic range,
  using parameters \emph{d} (dispersal/range expansion),
  \emph{e} (extinction/range contraction/local
  extirpation), and perhaps others (e.g. distance).  This
  matrix can be input in either dense or sparse (COO)
  format, as specified by \code{input_is_COO}.}

  \item{spPmat}{Default is \code{NULL}; users should
  usually use \code{spPmat_inputs}.  \code{spPmat} is A
  numeric matrix representing the probability of each
  ancestor range-->(Left range, Right range) transition at
  cladogenesis events.  There are different ways to
  represent this matrix.  In the simplest representation,
  this is just a rectangular matrix with numstates rows
  (representing the ancestral states) and numstates^2
  columns (representing all possible descendant pairs).
  Use of this type of matrix is specified by
  \code{cppSpMethod=1}. It is calculated from a textual
  speciation matrix (typically \code{spmat} in the code)
  via \code{\link{symbolic_to_relprob_matrix_sp}}. However,
  this matrix gets huge and slow for large numbers of
  states/ranges.  \code{cppSpMethod=2} and
  \code{cppSpMethod=3} implement successively more
  efficient and faster representation and processing of
  this matrix in COO-like formats.  See
  \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOprobs}}
  for the \code{cppSpMethod=2} method, and
  \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}}
  for the \code{cppSpMethod=3} method (the fastest).}

  \item{min_branchlength}{Nodes with branches below this
  branchlength will not be treated as cladogenesis events;
  instead, they will be treated as if an OTU had been
  sampled from an anagenetic lineage, i.e. as if you had a
  direct ancestor.  This is useful for putting fossils into
  the biogeography analysis, when you have fossil species
  that range through time. (Note: the proper way to obtain
  such trees, given that most phylogenetic methods force
  all OTUs to be tips rather than direct ancestors, is
  another question subject to active research.  However,
  one method might be to just set a branch-length cutoff,
  and treat any branches sufficiently small as direct
  ancestors.)}

  \item{return_what}{What should be returned to the user?
  Options are "loglike" (the log-likelihood of the data
  under the tree, model, and model parameters), "nodelikes"
  (the scaled conditional likelihoods at the nodes),
  "rootprobs" (the relative probability of the geographic
  ranges/states at the root), or "all" (all of the above in
  a list).  Typically the user will only want to return
  "loglike" while doing ML optimization, but then return
  "all" once the ML parameter values have been found.}

  \item{probs_of_states_at_root}{The prior probability of
  the states/geographic ranges at the root.  The default,
  \code{NULL}, effectively means an equal probability for
  each state (this is also what \code{LAGRANGE} assumes;
  and running with NULL will reproduce exactly the
  \code{LAGRANGE} parameter inferences and
  log-likelihood).}

  \item{rootedge}{Should the root edge be included in the
  calculation (i.e., calculate to the bottom of the root),
  if a root edge is present?  Default \code{FALSE}.}

  \item{sparse}{Should sparse matrix exponentiation be
  performed?  This should be faster for very large matrices
  (> 100-200 states), however, the calculations appear to
  be less accurate.  The function will transform a dense
  matrix to COO format (see
  \code{\link[rexpokit]{mat2coo}}) if necessary according
  to the \code{input_is_COO} parameter.}

  \item{printlevel}{If >= 1, various amounts of
  intermediate output will be printed to screen.  Note:
  Intermediate outputs from C++ and FORTRAN functions have
  been commented out, to meet CRAN guidelines.}

  \item{use_cpp}{Should the C++ routines from
  \code{\link[cladoRcpp]{cladoRcpp}} be used to speed up
  calculations?  Default \code{TRUE}.}

  \item{input_is_COO}{Is the input Q matrix a sparse,
  COO-formatted matrix (\code{TRUE}) or a standard dense
  matrix (\code{FALSE}). Default \code{FALSE}.}

  \item{spPmat_inputs}{A list of parameters so that
  \code{spPmat} (the speciation transition probability
  matrix) can be calculated on-the-fly, according to the
  method in \code{cppSpMethod}.  See example.}

  \item{cppSpMethod}{Three C++ methods from cladoRcpp for
  calculating and using the cladogenesis probability
  matrix.  1 is slowest but easiest to understand; 3 is
  fastest. If \code{spPmat_inputs} is given, the program
  will generate the appropriate spPmat on-the-fly, and the
  user does not have to input the full \code{spPmat}
  manually.}

  \item{cluster_already_open}{If the user wants to
  distribute the matrix exponentiation calculations from
  all the branches across a number of processors/nodes on a
  cluster, specify the cluster here.  E.g.
  \code{cluster_already_open =
  makeCluster(rep("localhost",num_cores_to_use), type =
  "SOCK")}.  Note: this will work on most platforms,
  including Macs running R from command line, but will NOT
  work on Macs running the R GUI \code{R.app}, because
  parallel processing functions like \code{MakeCluster}
  from e.g. \code{library(parallel)} for some reason crash
  R.app.  The program runs a check for R.app and will just
  run on 1 node if found.}

  \item{calc_ancprobs}{Should ancestral state estimation be
  performed (adds an uppass at the end).}

  \item{null_range_allowed}{Does the state space include
  the null range?#'}

  \item{fixnode}{If the state at a particular node is going
  to be fixed (e.g. for ML marginal ancestral states), give
  the node number.}

  \item{fixlikes}{The state likelihoods to be used at the
  fixed node.  I.e. 1 for the fixed state, and 0 for the
  others.}

  \item{stratified}{Default FALSE. If TRUE, you are running
  a stratified analysis, in which case uppass probs should
  be calculated elsewhere.}

  \item{states_allowed_TF}{Default NULL. If user gives a
  vector of TRUE and FALSE values, these states will be set
  to 0 likelihood throughout the calculations.}
}
\value{
  Return whatever is specified by \code{return_what}.
}
\description{
  This is the workhorse function of
  \code{\link[BioGeoBEARS]{BioGeoBEARS}}.  It calculates
  the likelihood of the tip data (the geographic ranges
  observed at the tips) given a phylogenetic tree, a Q
  transition matrix specifying the model of range evolution
  along branches, and a speciation probability matrix
  specifying the probability of the various possible
  ancestor-->(Left descendant, Right descendant) range
  evolution events at phylogenetic nodes/speciation events.
}
\details{
  This likelihood calculation will be repeated many
  hundreds or thousands of times in any ML (maximum
  likelihood) or Bayesian estimation procedure.  Thus, if
  the calculation of the log-likelihood of the data under
  one set of parameter values is too slow, inference takes
  days or becomes impossible.  However, by using fast
  matrix exponentiation (package
  \code{\link[rexpokit]{rexpokit}}) and fast C++ routines
  for calculating the probabilities of range inheritance
  scenarios at cladogenesis (package
  \code{\link[cladoRcpp]{cladoRcpp}}), major speed gains
  can be achieved. Most of the complexity in the input
  parameters and the code serves these more rapid
  alternatives.

  However, note that due to the explosion of the geographic
  range state space with more geographic areas (see
  \code{\link[cladoRcpp]{numstates_from_numareas}}), any
  computational method that explicitly calculates the
  likelihood of all states will eventually become unusable
  between 8-20 areas, depending on details.  An alternative
  method, which is fast for large numbers of areas, is
  BayArea, by Landis, Matzke, Moore, and Huelsenbeck; see
  \cite{Landis et al. (2013)}. However, BayArea does not
  currently implement cladogenesis models; it only has
  continuous-time model for evolutionary change along
  branches.  In effect, this means that the cladogenesis
  model is sympatric speciation with complete range copying
  with probability 1.
}
\note{
  Go BEARS!

  (COO = Coordinate list format for a matrix, see
  \url{http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29}
}
\examples{
testval=1
}
\author{
  Nicholas Matzke \email{matzke@berkeley.edu}
}
\references{
  Landis_Matzke_etal_2013_BayArea

  Matzke_2012_IBS

  ReeSmith2008
}
\seealso{
  \code{\link{calc_loglike_sp}},
  \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp}},
  \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOprobs}},
  \code{\link[cladoRcpp]{rcpp_calc_anclikes_sp_COOweights_faster}},
  \code{\link[rexpokit]{mat2coo}},
  \code{\link{rcpp_calc_anclikes_sp_COOweights_faster}}
}

