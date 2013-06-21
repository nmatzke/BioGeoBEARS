\name{plot_phylo3_nodecoords}
\alias{plot_phylo3_nodecoords}
\title{Modified plot.phylo to return the node coordinates}
\usage{
  plot_phylo3_nodecoords(x, type = "phylogram",
    use.edge.length = TRUE, node.pos = NULL,
    show.tip.label = TRUE, show.node.label = FALSE,
    edge.color = "black", edge.width = 1, edge.lty = 1,
    font = 3, cex = par("cex"), adj = NULL, srt = 0,
    no.margin = FALSE, root.edge = FALSE, label.offset = 0,
    underscore = FALSE, x.lim = NULL, y.lim = NULL,
    direction = "rightwards", lab4ut = "horizontal",
    tip.color = "black", plot = FALSE, rotate.tree = 0,
    ...)
}
\arguments{
  \item{x}{An ape phylo object}

  \item{type}{Default "phylogram"}

  \item{use.edge.length}{Default TRUE}

  \item{node.pos}{Default NULL}

  \item{show.tip.label}{Default TRUE}

  \item{show.node.label}{Default FALSE}

  \item{edge.color}{Default "black"}

  \item{edge.width}{Default 1}

  \item{edge.lty}{Default 1}

  \item{font}{Default 3}

  \item{cex}{Default par("cex")}

  \item{adj}{Default NULL}

  \item{srt}{Default 0}

  \item{no.margin}{Default FALSE}

  \item{root.edge}{Default FALSE}

  \item{label.offset}{Default 0}

  \item{underscore}{Default FALSE}

  \item{x.lim}{Default NULL}

  \item{y.lim}{Default NULL}

  \item{direction}{Default "rightwards"}

  \item{lab4ut}{Default "horizontal"}

  \item{tip.color}{Default "black"}

  \item{plot}{Default FALSE}

  \item{rotate.tree}{Default 0}

  \item{...}{Additional arguments to standard functions}
}
\value{
  A typical plot object, but with node coordinates added
}
\description{
  Just the standard plot.phylo (for phylo3 APE objects) In
  addition to the standard outputs, this function returns
  xx and yy, the plotting coordinates for those nodes.
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
  \code{\link[stats]{optim}}
}

