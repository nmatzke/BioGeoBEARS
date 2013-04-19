\name{remove_null_rowcols_from_mat}
\alias{remove_null_rowcols_from_mat}
\title{Remove rows or columns representing a null geographic range from a matrix}
\usage{
  remove_null_rowcols_from_mat(tmpmat, null_sym = "()")
}
\arguments{
  \item{tmpmat}{The matrix to check for null ranges.
  Function will only work if rows and columns have names,
  and one of the names matches \code{null_sym}.}

  \item{null_sym}{The character(s) denoting a null range.}
}
\value{
  \code{tmpmat3} The revised matrix.
}
\description{
  This function removes rows or columns representing a null
  geographic range from a matrix.
}
\details{
  LAGRANGE (\cite{Ree et al. (2008)}) and other models
  often assume that a null geographic range (the lineage
  inhabits no areas, i.e. is extinct) is a possible state.
  However, this is never a possible ancestral state (since
  an extinct lineage will never have descendants) so
  sometimes we must remove it.
}
\note{
  Go BEARS!
}
\examples{
testval=1
states_list = list("_", c("A"), c("B"), c("C"), c("A","B"), c("B","C"), c("A","C"), c("A","B","C"))

states_list = areas_list_to_states_list_new(areas=c("A","B","C"), include_null_range=TRUE, split_ABC=TRUE)
states_list

dedf = make_relprob_matrix_de(states_list=states_list, split_ABC=FALSE, split="", remove_simultaneous_events=TRUE, add_multiple_Ds=TRUE, dispersal_multiplier_matrix=make_dispersal_multiplier_matrix(states_list=states_list))

spmat_noNulls = remove_null_rowcols_from_mat(tmpmat=dedf, null_sym="()")
spmat_noNulls

spmat_noNulls = remove_null_rowcols_from_mat(tmpmat=dedf, null_sym="_")
spmat_noNulls
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
  \url{https://code.google.com/p/lagrange/}

  Ree RH and Smith SA (2008). "Maximum likelihood inference
  of geographic range evolution by dispersal, local
  extinction, and cladogenesis." _Systematic Biology_,
  *57*(1), pp. 4-14. <URL:
  http://dx.doi.org/10.1080/10635150701883881>, <URL:
  http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&dopt=Citation&list_uids=18253896>.

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
  \code{\link{areas_list_to_states_list_new}},
  \code{\link[cladoRcpp]{areas_list_to_states_list_old}},
  \code{\link{make_relprob_matrix_de}}
}

