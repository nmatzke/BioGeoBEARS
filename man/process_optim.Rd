\name{process_optim}
\alias{process_optim}
\title{Extract \code{optim} results to a row}
\usage{
  process_optim(optim_results, max_num_params = NULL)
}
\arguments{
  \item{optim_results}{A results object from
  \code{\link[stats]{optim}}}

  \item{max_num_params}{Specify the number of parameters,
  if known. If \code{NULL}, the code will try to guess.}
}
\value{
  \code{tmprow3} A row holding the
  \code{\link[stats]{optim}} results, which can then be
  added to a table with \code{\link[base]{rbind}}.
}
\description{
  After running an ML (maximum likelihood) search with
  \code{\link[stats]{optim}}, \code{\link[stats]{optim}}
  returns a list with a variety of objects.  It is often
  handy to have the parameter values, log-likelihood, etc.,
  extracted to a table for comparison with other
  optimization runs.  \code{process_optim} does this.
}
\note{
  Go BEARS!
}
\examples{
testval=1
# Any optim() for a biogeography scenario would take too long to run for R CMD check.
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

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
}
\seealso{
  \code{\link[stats]{optim}}
}

