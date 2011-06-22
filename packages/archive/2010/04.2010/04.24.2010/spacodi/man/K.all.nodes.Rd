\name{K.all.nodes}

\alias{K.all.nodes}

\title{a wrapper for computing Blomberg's K for all subtrees of a phylogeny}

\description{\code{K.all.nodes} calculates Blomberg's K through time, given a phylogeny and trait data.  This function is most sensibly used for continuous data.}

\usage{K.all.nodes(phy, traits, return.all = TRUE)}

\arguments{
 \item{phy}{a fully resolved phylogenetic tree, most likely ultrametric}
 \item{traits}{a dataframe of trait(s), with species as row.names}
 \item{return.all}{a logical identifier -- whether to return estimates for all nodes or for only those for which K is non-null}
}

\details{This function generates estimates of Blomberg's \code{K} for all subtrees of a phylogeny (if the number of tips subtended by a particular node exceeds two species.)

Blomberg's \code{K} represents a measure of trait similarity among tips of a phylogeny.  For traits evolving under Brownian motion, Blomberg's \code{K} 
is expected to take a value near 1.  If common ancestry is less predictive of species' traits, \code{K} is near 0.  Values of Blomberg's \code{K} near 0 may result from a 
scenario where species are strongly 'drawn' to the same trait optimum (e.g., under an Ornstein-Uhlenbeck process), or where species' trait-values are not 
phylogenetically heritable (that is, traits do not evolve in a treeward fashion, where related species would be most likely to exhibit similar trait values).
}

\value{
 \code{K.all.nodes} returns an array of phylogenetic heritabilities (i.e., phylogenetic signal) for continuous traits, across a phylogeny.
 The first dimension of the array is determined by is the number of internal nodes of the supplied phylogeny. The second dimension 
  of the output array are the measures returned to the user: \code{Blomberg.K} -- computed for each node (if nodes subtend more than two 
  sampled species; \code{tips} -- the number of contained species at each node; \code{node.time} -- the branching time of the node in the 
  supplied phylogeny (in same units as those supplied); and \code{node.ID} -- the identifier for each node (as interpretable by \code{\link[ape]{nodelabels}}). 
}

\references{BLOMBERG SP, T GARLAND, and AR IVES. 2003. Testing for phylogenetic signal in comparative data: Behavioral traits are more labile. Evolution 57:717-745.}

\author{Timothy Paine and Jonathan Eastman}

\seealso{
   \code{\link[geiger]{fitContinuous}} for fitting Pagel's lambda to 
   species' traits, as an alternative to Blomberg's K; see 
   \code{\link[picante]{Kcalc}} for the wrapped function; see 
   \code{\link{K.exp.nodes}} for generating expected values of \code{K} 
   at nodes
}

\examples{

# define a random tree
phy	<- prune.last.split(bd.tree(taxa.stop=61->tt))
phy$tip.label=paste("sp",seq(1:(tt-1)),sep="")

trt.BM <- as.data.frame(sim.char(phy,as.matrix(0.5))) 
names(trt.BM) = "BM.trait"

# simulated trait evolving under Brownian motion							
trt.BM

# compute K for all subtrees, expecting most values to be near 1
K.all.nodes(phy=phy, traits=trt.BM)

# generate random trait
trt.NULL <- data.frame(runif(phy$tip.label))
row.names(trt.NULL) <- phy$tip.label
names(trt.NULL) = "null.trait"

# compute K for all subtrees, expecting most values to be near 0
K.all.nodes(phy=phy, traits=trt.NULL)

}