\name{expokit_dgpadm_Qmat2_prebyte}
\alias{expokit_dgpadm_Qmat2_prebyte}
\title{A version of expokit_dgpadm_Qmat to byte-compile}
\usage{
  expokit_dgpadm_Qmat2_prebyte(times, Qmat,
    transpose_needed = TRUE)
}
\arguments{
  \item{times}{one or more time values to exponentiate by}

  \item{Qmat}{an input Q transition matrix}

  \item{transpose_needed}{If TRUE (default), matrix will be
  transposed (apparently EXPOKIT needs the input matrix to
  be transposed compared to normal)}
}
\value{
  \code{tmpoutmat} the output matrix.
}
\description{
  Byte-compiling is supposed to speed up functions; this is
  an attempt to do this on the
  \code{\link[rexpokit]{rexpokit}} function
  \code{\link[rexpokit]{expokit_dgpadm_Qmat}}.  It is also
  possible to byte-compile everything during package
  installation (via \code{ByteCompile: true} in the
  DESCRIPTION file), which is implemented in
  \code{BioGeoBEARS}, so this may be redundant.
}
\details{
  \code{\link{expokit_dgpadm_Qmat2_prebyte}} gets
  byte-compiled into \code{\link{expokit_dgpadm_Qmat2}}.

  See \url{http://dirk.eddelbuettel.com/blog/2011/04/12/}
  for discussion of the \code{compiler}
  (\code{\link[compiler]{cmpfun}}) package.
}
\note{
  Go BEARS!
}
\examples{
testval=1
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link[rexpokit]{expokit_dgpadm_Qmat}},
  \code{\link{expokit_dgpadm_Qmat2}},
  \code{\link[compiler]{compile}},
  \code{\link[compiler]{cmpfun}}
}

