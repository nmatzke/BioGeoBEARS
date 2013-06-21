\name{AkaikeWeights_on_summary_table}
\alias{AkaikeWeights_on_summary_table}
\title{Calculate Akaike Weights, and add to table}
\usage{
  AkaikeWeights_on_summary_table(restable,
    colname_to_use = "AIC", add_to_table = TRUE)
}
\arguments{
  \item{restable}{A \code{\link[base]{data.frame}} with at
  least a column named as in \code{add_to_table}.}

  \item{colname_to_use}{The name of the column containing
  AIC values.}

  \item{add_to_table}{If TRUE, add to the main table and
  return the main table. If FALSE, return just the Akaike
  Weights results.}
}
\value{
  \code{restable}, the modified table, or \code{wt_vBest},
  the Akaike Weights results.
}
\description{
  This calculates Akaike Weights (relative probabilities on
  models explaining the same data) for the models in a
  column in a table.
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
  \url{http://www.brianomeara.info/tutorials/aic}

  Matzke_2012_IBS

  Burnham_Anderson_2002
}
\seealso{
  \code{\link{calc_AIC_column}},
  \code{\link{calc_AICc_column}}
}

