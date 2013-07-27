\name{prt}
\alias{prt}
\title{Print tree in table format}
\usage{
  prt(t, printflag = TRUE, relabel_nodes = FALSE,
    time_bp_digits = 7, add_root_edge = TRUE,
    get_tipnames = FALSE, fossils_older_than = 0.6)
}
\arguments{
  \item{t}{A \code{\link[ape]{phylo}} tree object.}

  \item{printflag}{Should the table be printed to screen?
  Default TRUE.}

  \item{relabel_nodes}{Manually renumber the internal
  nodes, if desired. Default FALSE.}

  \item{time_bp_digits}{The number of digits to print in
  the time_bp (time before present) column. Default=7.}

  \item{add_root_edge}{Should a root edge be added?
  Default \code{TRUE}.}

  \item{get_tipnames}{Should the list of tipnames
  descending from each node be printed as a string in
  another column? This is slow-ish, but useful for matching
  up nodes between differing trees. Default \code{FALSE}.}

  \item{fossils_older_than}{Tips that are older than
  \code{fossils_older_than} will be marked as \code{TRUE}
  in a column called \code{fossil}. This is not currently
  set to 0, because Newick files can have slight precision
  issues etc. that mean not all tips quite come to zero.
  You can attempt to fix this with
  \code{\link{average_tr_tips}} (but make sure you do not
  inappropriately average in fossils!!).}
}
\value{
  \code{dtf} A \code{\link[base]{data.frame}} holding the
  table. (Similar to the printout of a
  \code{\link[phylobase]{phylo4}} object.)
}
\description{
  Learning and using APE's tree structure can be difficult
  and confusing because much of the information is
  implicit.  This function prints the entire tree to a
  table, and makes much of the implicit information
  explicit.  It is not particularly fast, but it is useful.
}
\details{
  See \url{http://ape.mpl.ird.fr/ape_development.html} for
  the official documentation of R tree objects.
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
  \url{http://ape.mpl.ird.fr/ape_development.html}

  Matzke_2012_IBS
}
\seealso{
  \code{\link[ape]{phylo}}, \code{\link{average_tr_tips}}
}

