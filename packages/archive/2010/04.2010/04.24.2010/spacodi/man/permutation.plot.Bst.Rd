\name{permutation.plot.Bst}
\alias{permutation.plot.Bst}
\title{plotting observed and expected Bst across branching times of a phylogeny}

\usage{permutation.plot.Bst(phy, sp_plot, n.rep = 10, method = "1s", all.points = TRUE, envelope = TRUE, add.id = FALSE, cex = list(2, 0.2, 1), outfile = NULL, ...)}

\arguments{
  \item{phy}{a fully resolved phylogenetic tree, most sensibly ultrametric}
  \item{sp_plot}{dataset of species abundances among community plots}
  \item{n.rep}{number of randomization of \code{sp_plot} to perform}
  \item{method}{one of several permutation methods for the dataset; see \code{\link{Bst.permutation}} and the default method (\code{\link{resamp.1s}})}
  \item{all.points}{a logical indicator determining whether all expected \code{Bst}s are to be plotted}
  \item{envelope}{a logical indicator specifying whether a smooth-spline envelope is drawn around expected \code{Bst}s}
  \item{add.id}{a logical indicator dictating whether \code{node.labels} are placed near observed \code{Bst}s; see \code{\link[ape]{nodelabels}}}
  \item{cex}{a list of three elements: size of points for observed \code{Bst}s; for expected \code{Bst}s; and for size of text for \code{node.labels}}
  \item{outfile}{an optional .pdf file to which to write output}
  \item{\dots}{additional plotting parameters to be optionally supplied}
}

\details{
  This function generates a plot of observed and expected phylogenetic turnover (\code{Bst}) evaluated for 
  each subtree of the supplied phylogeny. If \code{all.points=TRUE}, empirical \code{Bst}s as well as those 
  expected from permuted community datasets are plotted together. Several methods are available for data 
  permutation (see \code{\link{Bst.permutation}} for further details). If \code{envelope=TRUE}, an envelope for 
  a credible range of \code{Bst}s is additionally plotted.  
}

\value{a plot of \code{Bst} values for a phylogeny, with the option to direct output to a .pdf file}

\author{Jonathan Eastman}


\seealso{\code{\link{spacodi.calc}} for further details on \code{Bst}; \code{\link[stats]{smooth.spline}} for generating the envelope}

\examples{
# load a species-by-plots matrix, along with a tree
data(spacodi.example)
attributes(spacodi.example)
attach(spacodi.example)
sp_plot
sp_tree

# generate a plot of observed and expected Bst
permutation.plot.Bst(phy=sp_tree, sp_plot=sp_plot, n.rep=50, add.id=TRUE)
}