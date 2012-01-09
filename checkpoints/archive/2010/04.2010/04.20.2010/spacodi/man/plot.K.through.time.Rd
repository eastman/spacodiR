\name{plot.K.through.time}

\alias{plot.K.through.time}

\title{
	plotting Blomberg's K through time
}

\description{
	\code{plot.K.through.time} calculates and plots Blomberg's K through time, given a phylogeny and trait data.  This function is most sensibly used for continuous data.
}

\usage{
	plot.K.through.time(phy, traits, outfile = NULL)
}

\arguments{
  \item{phy}{a fully resolved phylogenetic tree, most likely ultrametric}
  \item{traits}{a dataframe of a trait, with species as row.names}
  \item{outfile}{an optional .pdf file to which to write output}
}

\details{This function generates estimates of Blomberg's \code{K} for all subtrees of a phylogeny (if the number of tips subtended by a particular node exceeds two species.)

Blomberg's \code{K} represents a measure of trait similarity among tips of a phylogeny.  For traits evolving under Brownian motion, Blomberg's \code{K} 
is expected to take a value near 1.  If common ancestry is less predictive of species' traits, \code{K} is near 0.  Values of Blomberg's \code{K} near 0 may result from a 
scenario where species are strongly 'drawn' to the same trait optimum (e.g., under an Ornstein-Uhlenbeck process), or where species' trait-values are not 
phylogenetically heritable (that is, traits do not evolve in a treeward fashion, where related species would be most likely to exhibit similar trait values).
}

\value{
	A plot of Blomberg's \code{K} against splitting times of the phylogeny.  The leftmost value is that for the root node.  Values closer to \code{branching.time=0} 
	are those closer to the present.
}

\references{
BLOMBERG SP, T GARLAND, and AR IVES. 2003. Testing for phylogenetic signal in comparative data: Behavioral traits are more labile. Evolution 57:717-745.
}
\author{
Jonathan Eastman}

\seealso{
	\code{geiger:::fitContinuous()} for fitting trait-evolution models to datasets; \code{\link{K.all.nodes}} for the underlying function}

\examples{

foo	<- r.plot(species=200,plots=2,missing.prop=0.15,sim.tree=TRUE)	

# define a random tree
phy	<- foo[[2]]

trt.BM <- as.data.frame(sim.char(phy,as.matrix(0.5))) 
names(trt.BM) = "BM.trait"

# simulated trait evolving under Brownian motion							
trt.BM

# plot K through time, expecting most values to be near 1
plot.K.through.time(phy, trt.BM)
}