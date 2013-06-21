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
states_list = list("_", c("A"), c("B"), c("C"), c("A","B"),
c("B","C"), c("A","C"), c("A","B","C"))

states_list = areas_list_to_states_list_new(areas=c("A","B","C"),
include_null_range=TRUE, split_ABC=TRUE)
states_list

dedf = make_relprob_matrix_de(states_list=states_list,
split_ABC=FALSE, split="", remove_simultaneous_events=TRUE,
add_multiple_Ds=TRUE,
dispersal_multiplier_matrix=make_dispersal_multiplier_matrix(states_list=states_list))

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

  ReeSmith2008

  Matzke_2012_IBS
}
\seealso{
  \code{\link{areas_list_to_states_list_new}},
  \code{\link[cladoRcpp]{areas_list_to_states_list_old}},
  \code{\link{make_relprob_matrix_de}}
}

