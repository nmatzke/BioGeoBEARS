\name{parse_lagrange_python_output}
\alias{parse_lagrange_python_output}
\title{Parse the output file from python \code{LAGRANGE}}
\usage{
  parse_lagrange_python_output(outfn = "output.results.txt",
    outputfiles = FALSE, results_dir = getwd(),
    new_splits_fn = TRUE, new_states_fn = FALSE,
    filecount = 0, append = FALSE)
}
\arguments{
  \item{outfn}{The python \code{LAGRANGE} output text
  file.}

  \item{outputfiles}{Should parsed output be written to
  files? Default FALSE.}

  \item{results_dir}{The directory \code{outfn} is in.}

  \item{new_splits_fn}{Should a text file containing a
  table of the splits and their probabilities be output?
  Default \code{TRUE}.}

  \item{new_states_fn}{Should a text file containing a
  table of the states and their probabilities be output?
  Default \code{FALSE}, as I don't believe python
  \code{LAGRANGE} will output the states at the nodes (C++
  \code{LAGRANGE} will, however).}

  \item{filecount}{The starting number for the filecount
  (relevant if one is processing many files).}

  \item{append}{Should results be appended to preexisting
  file? (default \code{FALSE})}
}
\value{
  sumstats A \code{\link[base]{data.frame}} containing the
  summary statistics (LnL, d and e rates, etc.)  The splits
  filename is output to screen.
}
\description{
  Parse the output of a python \code{LAGRANGE}.
}
\details{
  Python LAGRANGE is run from a UNIX/Terminal command-line
  with a command such as "\code{python
  lagrangefilename.py}".  You will need to have the
  "lagrange" python directory in your working directory.

  The input file can be obtained from
  \url{http://www.reelab.net/lagrange/configurator/index}
  (\cite{Ree (2009)}).

  Python comes installed on many machines, or can be
  downloaded from the Enthought Python Distribution
  (\url{https://www.enthought.com/products/epd/}).
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
  \url{https://code.google.com/p/lagrange/}
  \url{https://www.enthought.com/products/epd/}
  \url{http://www.reelab.net/lagrange/configurator/index}

  Ree2009configurator

  Matzke_2012_IBS

  ReeSmith2008
}
\seealso{
  \code{\link{get_lagrange_nodenums}},
  \code{\link{LGpy_splits_fn_to_table}}
}

