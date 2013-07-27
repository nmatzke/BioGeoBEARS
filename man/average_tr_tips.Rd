\name{average_tr_tips}
\alias{average_tr_tips}
\title{Average the heights of (non-fossil) tips to make ultrametric-ish.}
\usage{
  average_tr_tips(tr, fossils_older_than = 0.6)
}
\arguments{
  \item{tr}{An ape phylo object}

  \item{fossils_older_than}{Tips that are older than
  \code{fossils_older_than} will be excluded from the tips
  that are going to be averaged. This is not currently set
  to 0, because Newick files can have slight precision
  issues etc. that mean not all tips quite come to zero
  (which is why you need \code{\link{average_tr_tips}} in
  the first place!). Obviously you should be cautious about
  the value of , depending on the absolute timescale of
  your tree. Make sure you do not inappropriately average
  in fossils!!}
}
\value{
  \code{edge_times_bp} A 2-column matrix with the age (from
  the present) of the top and bottom of each edge.
}
\description{
  When you have a digitized tree, or other slightly uneven
  source tree, average the tips to get them all to line up
  at 0 my before present.  This makes an ultrametric tree
  if and only if there are no fossil tips in the tree.
}
\details{
  If the user includes fossils accidentally, this function
  can easily lead to pathological results (negative branch
  lengths etc.), so use with care!!
}
\note{
  Go BEARS!
}
\examples{
test=1
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link{prt}}, \code{\link{chainsaw2}},
  \code{\link{extend_tips_to_ultrametricize}}
}

