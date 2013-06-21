\name{get_probvals}
\alias{get_probvals}
\title{Calculate probability of ordered discrete states using a maxent distribution (equations 6.3-6.4 of Harte 2011)}
\usage{
  get_probvals(die_vals, meanval)
}
\arguments{
  \item{die_vals}{Values of the ordered discrete variable
  state (e.g., \code{seq(1,6)} for a six-sided die)}

  \item{meanval}{Mean value (the knowledge supplied to the
  MaxEnt function).}
}
\value{
  \code{Prob_nvals}, numeric values of the probability of
  each state from \code{die_vals}.
}
\description{
  This function is calculates the Maximum Entropy
  (\cite{Harte (2011)}) discrete probability distribution
  of a number of ordered states (e.g., faces of a 6-sided
  die) given the mean of many rolls.  Here, this is merely
  used so that a single parameter can control the
  probability distribution of small versus large descendant
  areas during cladogenesis.  This function could then used
  by \code{\link{relative_probabilities_of_subsets}} in
  BioGeoBEARS to weight different descendant range sizes
  (although, currently, the function
  \code{\link[FD]{maxent}} from the \code{\link[FD]{FD}}
  package is used).
}
\details{
  This calculation is based on Equations 6.3-6.4 of
  \cite{Harte (2011)}.

  See also: Maximum Entropy probability distribution for
  discrete variable with given mean (and discrete uniform
  flat prior)
  \url{http://en.wikipedia.org/wiki/Maximum_entropy_probability_distribution}
}
\note{
  Go BEARS!
}
\examples{
testval=1
# Examples
# Set up subplots
par(mfrow=c(3,2))

# Flat distribution (equal prob of any descendent size)
N = 6
# n = die_vals
die_vals = seq(1,N)
meanval = 3.5
probvals = get_probvals(die_vals, meanval)
probvals
barplot(height=probvals, width=1, names.arg=die_vals, ylim=c(0,1))
title(paste("Probabilities of each state, mean val=", meanval, sep=""))

# Descendents tend to have large ranges
N = 6
# n = die_vals
die_vals = seq(1,N)
meanval = 5.999
probvals = get_probvals(die_vals, meanval)
probvals
barplot(height=probvals, width=1, names.arg=die_vals, ylim=c(0,1))
title(paste("Probabilities of each state, mean val=", meanval, sep=""))

# Flat distribution (equal prob of any descendent size)
N = 6
# n = die_vals
die_vals = seq(1,N)
meanval = 5
probvals = get_probvals(die_vals, meanval)
probvals
barplot(height=probvals, width=1, names.arg=die_vals, ylim=c(0,1))
title(paste("Probabilities of each state, mean val=", meanval, sep=""))

# Flat distribution (equal prob of any descendent size)
N = 6
# n = die_vals
die_vals = seq(1,N)
meanval = 4
probvals = get_probvals(die_vals, meanval)
probvals
barplot(height=probvals, width=1, names.arg=die_vals, ylim=c(0,1))
title(paste("Probabilities of each state, mean val=", meanval, sep=""))

# Flat distribution (equal prob of any descendent size)
N = 6
# n = die_vals
die_vals = seq(1,N)
meanval = 2
probvals = get_probvals(die_vals, meanval)
probvals
barplot(height=probvals, width=1, names.arg=die_vals, ylim=c(0,1))
title(paste("Probabilities of each state, mean val=", meanval, sep=""))

# This produces the LAGRANGE default
# (all smaller descendents are of size 1)
N = 6
# n = die_vals
die_vals = seq(1,N)
meanval = 1.0001
probvals = get_probvals(die_vals, meanval)
probvals
barplot(height=probvals, width=1, names.arg=die_vals, ylim=c(0,1))
title(paste("LAGRANGE 'default', mean val=", meanval, sep=""))

# This is stopped by the error check
# (all smaller descendents are of size 1)
# N = 6
# # n = die_vals
# die_vals = seq(1,N)
# meanval = 0.5
# probvals = get_probvals(die_vals, meanval)
# probvals
# barplot(height=probvals, width=1, names.arg=die_vals, ylim=c(0,1))
# title(paste("Probabilities of each state, mean val=", meanval, sep=""))
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\references{
  \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
  \url{http://en.wikipedia.org/wiki/Maximum_entropy_probability_distribution}

  Harte2011

  Matzke_2012_IBS
}
\seealso{
  \code{\link{calcZ_part}}, \code{\link{calcP_n}},
  \code{\link[FD]{maxent}},
  \code{\link{symbolic_to_relprob_matrix_sp}},
  \code{\link{relative_probabilities_of_subsets}}
}

