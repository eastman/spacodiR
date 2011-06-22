\name{K.by.nodes}

\alias{K.by.nodes}

\title{conducting randomization tests of phylogenetic trait heritability through time}

\usage{K.by.nodes(sp.traits, phy, obs.only = FALSE, return.all = TRUE, n.rep = 10, rand.test = TRUE, r.rep = 10000)}

\arguments{
  \item{phy}{a fully resolved phylogenetic tree, most likely ultrametric}
  \item{sp.traits}{a dataframe of trait(s), with species as row.names}
  \item{obs.only}{a logical indicator determining whether only estimated \code{K}s for the empirical data are returned}
  \item{return.all}{a logical indicator determing whether results for all nodes are returned for the empirical data (regardless of whether \code{K} is non-null)}
  \item{n.rep}{number of simulation of trait evolution are to be performed}
  \item{rand.test}{a logical indicator for whether a randomization test of observed and expected \code{K}s is performed}
  \item{r.rep}{number of randomization iterations to perform if \code{rand.test=TRUE}}
}

\details{
  \code{K.by.nodes} is a function to conduct randomization tests on estimated and observed \code{K}s through time, for all
  subtrees of a phylogeny (if the number of tips subtended by a particular node exceeds two species).   
  
  Blomberg's \code{K} represents a measure of trait similarity among tips of a phylogeny.  For traits evolving under Brownian motion, Blomberg's \code{K} 
  is expected to take a value near 1.  If common ancestry is less predictive of species' traits, \code{K} is near 0.  Values of Blomberg's \code{K} near 0 may result from a 
  scenario where species are strongly 'drawn' to the same trait optimum (e.g., under an Ornstein-Uhlenbeck process), or where species' trait-values are not 
  phylogenetically heritable (that is, traits do not evolve in a treeward fashion, where related species would otherwise be most likely to exhibit similar trait values).
}

\value{  
  \code{K.by.nodes} returns a list of three elements:
  \itemize{
    \item \code{observed.K}{: estimated \code{K} for all nodes of the tree for which \code{K} is calculable}
	\item \code{expected.K}{: expected \code{K} for all nodes in the form of an array, with nodes as row names}
	\item \code{randomization.test}{: results from randomization tests of \code{K} for each node. Note: if randomization test is used, 
	some randomization iterations may return NA for particular nodes (see details above). The tally of valid expected 
	\code{K}s from permutation is recorded as a column in this element (i.e., \code{randomization.test$valid.comparisons}).  
	Be wary of significant results that are based off few valid comparisons.}
  }
}

\references{BLOMBERG SP, T GARLAND, and AR IVES. 2003. Testing for phylogenetic signal in comparative data: Behavioral traits are more labile. Evolution 57:717-745.}

\author{Jonathan Eastman and Timothy Paine}

\seealso{
   \code{\link[geiger]{fitContinuous}} for fitting Pagel's lambda to 
   species' traits as an alternative to Blomberg's K; \code{\link[geiger]{deltaTree}} for evolving characters under 
   various models; see \code{\link[picante]{Kcalc}} for the wrapped function to compute \code{K}; see \code{\link{spacodi.tree.plot}} and \code{\link{K.simulation.plot}}
   for plotting output 
}

\examples{
# load a species-by-plots matrix, along with a tree
data(sp.example)
attach(sp.example)
sp.traits
sp.tree

# define a single trait
trt=as.data.frame(sp.traits[,1])
row.names(trt)=row.names(sp.traits)
trt

# generate a plot of observed and expected K
K.simulation.plot(K.by.nodes(sp.traits=trt, phy=sp.tree, n.rep=10)->sp.permut, cex=list(2,0.25,0.75), bg=list("white",1), col=list(1,1),add.id=TRUE, envelope=FALSE, lwd=0.4)
}
