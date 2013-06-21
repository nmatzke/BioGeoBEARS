\name{nullsym_to_NA}
\alias{nullsym_to_NA}
\title{Convert a specified null range code to NA}
\usage{
  nullsym_to_NA(mat, nullsym = "-")
}
\arguments{
  \item{mat}{A matrix.}

  \item{nullsym}{A character specifying the null symbol.}
}
\value{
  \code{mat} The revised matrix
}
\description{
  Takes a matrix \code{mat}, converts any instances of the
  \code{nullsym} symbol to \code{NA}.
}
\note{
  Go BEARS!
}
\examples{
testval=1
mat = matrix(c("-",1,1,1,"-",1,1,1,"-"), nrow=3)
mat
mat2 = nullsym_to_NA(mat, nullsym="-")
mat2
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link{remove_null_rowcols_from_mat}}
}

