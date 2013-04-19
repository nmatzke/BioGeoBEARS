\name{getareas_from_tipranges_object}
\alias{getareas_from_tipranges_object}
\title{Get the names of the areas in a tipranges object}
\usage{
  getareas_from_tipranges_object(tipranges)
}
\arguments{
  \item{tipranges}{An object of class \code{tipranges}.}
}
\value{
  \code{areanames}, a list of the names of the areas
}
\description{
  This function extracts the names of the areas in a
  \code{tipranges} object.  Just a shortcut for
  \code{names(tipranges@df)}.
}
\note{
  Go BEARS!
}
\examples{
testval=1
tipranges_object = define_tipranges_object()
tipranges_object

areanames = getareas_from_tipranges_object(tipranges_object)
areanames
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\seealso{
  \code{\link{define_tipranges_object}},
  \code{\link[cladoRcpp]{areas_list_to_states_list_old}},
  \code{\link{areas_list_to_states_list_new}},
  \code{\link{tipranges_to_tip_condlikes_of_data_on_each_state}}
}

