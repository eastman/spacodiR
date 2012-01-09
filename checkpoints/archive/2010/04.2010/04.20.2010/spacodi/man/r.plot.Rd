\name{r.plot}

\alias{r.plot}

\title{simulating community phylogenetics samples}

\description{A utility for simulating species-by-plots matrices}

\usage{r.plot(species = 100, plots = 30, missing.prop = 0.6, mean.abund = 15, sim.tree = FALSE)}

\arguments{
  \item{species}{number of species sampled globally, across all plots}
  \item{plots}{number of community samples}
  \item{missing.prop}{a proportion of plots for which species are unsampled on average}
  \item{mean.abund}{a mean abundance (i.e., number of individuals) of species across all plots in which the species are present}
  \item{sim.tree}{a logical flag dictating whether a birth-death tree is to be simulated as well as the dataset}
}
\details{\code{r.plot} simulates a species-by-plots matrix (and a tree if \code{sim.tree=TRUE}). The dimensions of the dataset are determined by \code{species}
and \code{plots}. The proportion of 'empty' cells in the dataset (\code{missing.prop}) is a measure of the distinctiveness of community samples. Where species are 
present, \code{mean.abund} determines the average number of individuals expected to be sampled in a given plot: abundances are assumed to be log-normal distributed 
across plots.
}

\value{A list of two components

\item{out[[1]] }{a simulated dataset}
\item{out[[2]] }{a simulated tree (only if \code{sim.tree=TRUE})}
}

\author{Jonathan Eastman}

\seealso{\code{\link{bd.tree}} for the function used to simulate trees}


\examples{
# simulate a dataset where many species are shared among plots 
r.plot(species=25, plots=5, missing.prop=0.02, mean.abund=10, sim.tree=FALSE)

# simulate a dataset where many species are unshared among plots 
r.plot(species=25, plots=5, missing.prop=0.90, mean.abund=10, sim.tree=FALSE)

# compute phylogenetic diversity measures for a simulated dataset and tree
foo <- r.plot(species=25, plots=5, missing.prop=0.20, mean.abund=10, sim.tree=TRUE)
sp_plot <- foo[[1]]
phy <- foo[[2]]

## OUTPUT ##
spacodi.calc(sp_plot, phy)

}