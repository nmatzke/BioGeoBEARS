\name{pdfit}
\alias{pdfit}
\title{Print a table to LaTeX format}
\usage{
  pdfit(table_vals, file_prefix = "tmptable",
    size = "\\\\tiny", gettex = FALSE, caption = NULL)
}
\arguments{
  \item{table_vals}{A table, hopefully produced by
  \code{\link{conditional_format_table}}.}

  \item{file_prefix}{The prefix for the output PDF and the
  intermediate files.}

  \item{size}{Font size, overriding
  \code{getOption("xtable.size")}. Default is "tiny" (with
  backslashes). You can also try "small".  Input
  \code{NULL} (without quotes or backslashes) for medium.
  (\code{NULL} is the options default.)}

  \item{gettex}{If TRUE, the \code{tex} code for the table
  is returned.}

  \item{caption}{A caption, if desired.}
}
\value{
  \code{texfile} The filename of the \code{tex} file.
}
\description{
  This function prints a table to PDF via
  \code{\link[xtable]{xtable}} and the LaTeX
  \code{pdflatex} function.  It will only work if you have
  command-line LaTeX installed.
}
\details{
  This function was inspired by
  \url{http://tex.stackexchange.com/questions/15013/generate-a-pdf-containing-r-output-inside-latex-table}.
}
\note{
  Go BEARS!
}
\examples{
test=1

# Setup data
\dontrun{
data = c(2.768443, 1.869964, 5.303702, 4.733483,  2.123816,
18.551051, 5.483625,  3.590745,  18.772389)
result = matrix(data, nrow=3, byrow=TRUE)
result = as.data.frame(result)
names(result) = c("CV", "LCB", "UCB")
rownames(result) = c("within", "between", "total")
result
pdfit(table_vals=result)#' }
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link{pdftable}}
}

