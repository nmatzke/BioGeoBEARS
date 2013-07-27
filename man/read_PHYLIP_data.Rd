\name{read_PHYLIP_data}
\alias{read_PHYLIP_data}
\title{Read a PHYLIP-format file}
\usage{
  read_PHYLIP_data(lgdata_fn = "lagrange_area_data_file.data",
    regionnames = NULL)
}
\arguments{
  \item{lgdata_fn}{The filename to read.}

  \item{regionnames}{A list of the names of the areas. Only
  used if the names are NOT specified in the file.}
}
\value{
  \code{tmpdf} A \code{\link[base]{data.frame}} containing
  the data.
}
\description{
  This assumes data are interleaved, and that names are
  separated from data by a tab character; there is no
  10-character limit on names.
}
\details{
  This function is a precursor to
  \code{\link{getranges_from_LagrangePHYLIP}}.
}
\note{
  Go BEARS!
}
\examples{
testval=1

# Get the example files directory
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
# tmp hard code: extdata_dir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/"
# Set the filename (Hawaiian Psychotria from Ree & Smith 2008)
fn = np(paste(extdata_dir, "/Psychotria_geog.data", sep=""))

# Read in the file
tmpdf = read_PHYLIP_data(lgdata_fn=fn, regionnames=NULL)
tmpdf

# Read in the file
tmpdf = read_PHYLIP_data(lgdata_fn=fn,
regionnames=c("Kauai", "Oahu", "Maui-Nui","Big Island"))
tmpdf	# Note that regionnames are only
# used if they are NOT specified in the file.
# But, you could put them on manually
names(tmpdf) = c("Kauai", "Oahu", "Maui-Nui","Big Island")
tmpdf

# This one has no area names
fn = np(paste(extdata_dir, "/Psychotria_geog_noAreaNames.data", sep=""))
tmpdf = read_PHYLIP_data(lgdata_fn=fn,
regionnames=c("Kauai", "Oahu", "Maui-Nui","Big Island"))
tmpdf	# Note that regionnames are only
# used if they are NOT specified in the file.
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}

  Matzke_2012_IBS
}
\seealso{
  \code{\link{getranges_from_LagrangePHYLIP}}
}

