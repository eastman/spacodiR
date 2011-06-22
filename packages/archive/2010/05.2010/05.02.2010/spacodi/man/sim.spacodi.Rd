\name{sim.spacodi}

\alias{sim.spacodi}

\title{simulating community phylogenetics samples}

\description{A utility for simulating species-by-plots matrices}

\usage{sim.spacodi(species = 100, plots = 30, missing.prop = 0.6, mean.abund = 15, sim.tree = FALSE, ...)}

\arguments{
  \item{species}{number of species sampled globally, across all plots}
  \item{plots}{number of community samples}
  \item{missing.prop}{a proportion of plots for which species are unsampled on average}
  \item{mean.abund}{a mean abundance (i.e., number of individuals) of species across all plots in which the species are present}
  \item{sim.tree}{a logical flag dictating whether a birth-death tree is to be simulated as well as the dataset}
  \item{\dots}{additional parameters to be directed to embedded functions (e.g., \code{\link{bd.tree}})}

}
\details{\code{sim.spacodi} simulates a community-phylogenetic dataset in \code{SPACoDi} format (and a tree if \code{sim.tree=TRUE}). The dimensions of the dataset are determined by \code{species}
and \code{plots}. The proportion of 'empty' cells in the dataset (\code{missing.prop}) is one measure of the distinctiveness of community samples. Where species are 
present, \code{mean.abund} determines the average number of individuals expected to be sampled in a given plot: abundances are assumed to be log-normal distributed 
amongst plots.
}

\value{A list of two components:
\itemize{
\item \code{sp_plot}{: a simulated dataset}
\item \code{sp_tree}{: a simulated tree (only if \code{sim.tree=TRUE})}
}
}

\author{Jonathan Eastman}

\seealso{\code{\link{bd.tree}} for the function used to simulate trees}


\examples{
# simulate a dataset where many species are shared among plots 
sim.spacodi(species=25, plots=5, missing.prop=0.02, mean.abund=10, sim.tree=FALSE)

# simulate a dataset where many species are unshared among plots 
sim.spacodi(species=25, plots=5, missing.prop=0.90, mean.abund=10, sim.tree=FALSE)

# compute phylogenetic diversity measures for a simulated dataset and tree
foo <- sim.spacodi(species=25, plots=5, missing.prop=0.20, mean.abund=10, sim.tree=TRUE)
sp.plot <- foo[[1]]
phy <- foo[[2]]

## OUTPUT ##
spacodi.calc(sp.plot=sp.plot, phy=phy)

}