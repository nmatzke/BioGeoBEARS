\name{parse_lagrange_output_old}
\alias{parse_lagrange_output_old}
\title{Parse the output file from python \code{LAGRANGE} -- older version}
\usage{
  parse_lagrange_output_old(outfn, results_dir = getwd(),
    new_splits_fn = TRUE, new_states_fn = TRUE,
    filecount = 0)
}
\arguments{
  \item{outfn}{The C++ \code{LAGRANGE} output text file.}

  \item{results_dir}{The directory \code{outfn} is in.}

  \item{new_splits_fn}{Should a text file containing a
  table of the splits and their probabilities be output?
  Default \code{TRUE}.}

  \item{new_states_fn}{Should a text file containing a
  table of the splits and their probabilities be output?
  Default \code{TRUE}, unlike python \code{LAGRANGE}, C++
  \code{LAGRANGE} \emph{will} output the states at the
  nodes.}

  \item{filecount}{The starting number for the filecount
  (relevant if one is processing many files).}
}
\value{
  sumstats A \code{\link[base]{data.frame}} containing the
  summary statistics (LnL, d and e rates, etc.)  The splits
  filename is output to screen.
}
\description{
  Parse the output of a C++ \code{LAGRANGE} run.
}
\details{
  This function parses the output of \code{LAGRANGE},
  obtained by a command such as the following, run at a
  UNIX/Mac Terminal command line.  This is an older version
  useful for automating processing of many files.

  \code{cd
  /Users/nick/Desktop/__projects/_2011-07-15_Hannah_spider_fossils/_data/lagrange_for_nick}

  \code{./lagrange_cpp palp_no_Lacun_v1_2nd387.lg >
  lagrange_results_v1_2nd387.txt}

  C++ LAGRANGE can be obtained at
  \url{https://code.google.com/p/lagrange/}
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

  Matzke_2012_IBS

  ReeSmith2008
}
\seealso{
  \code{\link{get_lagrange_nodenums}},
  \code{\link{LGpy_splits_fn_to_table}}
}

