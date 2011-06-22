\name{match.spacodi.data}

\alias{match.spacodi.data}

\title{pruning a tree to match datasets used in SPACoDi}

\description{
 \code{match.spacodi.data} is a utility for automating the process of pruning datasets for use in \pkg{spacodi}. Note that nearly all \pkg{spacodi} functions require that
 trees have present (at least) all species that are contained within the \code{spacodi.object}.}

\usage{match.spacodi.data(sp.plot, phy=NULL, sp.traits=NULL, prune=TRUE, verbose=FALSE)}

\arguments{
  \item{sp.plot}{a community-phylogenetic dataset in \code{SPACoDi} format (see \code{\link{as.spacodi}})}
  \item{phy}{a fully resolved phylogenetic tree, most sensibly ultrametric}
  \item{sp.traits}{a species-by-trait(s) dataframe}
  \item{prune}{a logical switch indicating whether plots with fewer than two species sampled are pruned}
  \item{verbose}{a logical switch determining whether dropped plots are written to the console}
}

\value{
A list of pruned dataset(s): 
\itemize{
    \item \code{sp.plot}
	
    \item \code{sp.tree}
	
	\item \code{sp.traits}
	}
	Note: if \code{prune=TRUE}, the list of plots pruned will be printed to the console.
}

\author{Jonathan Eastman}

\examples{
# load a species-by-plots matrix, along with a tree
data(sp.example)
attributes(sp.example)
attach(sp.example)
sp.plot
sp.tree
sp.traits

# prune out undersampled plots
sp.plot[,2]=0
match.spacodi.data(sp.plot) -> sp.plot.new
sp.plot.new


# match datasets where sp.traits is smaller than the remainder
match.spacodi.data(sp.plot=sp.plot, phy=sp.tree, sp.traits=sp.traits[1:6,])
							
}


