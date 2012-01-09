\name{plot.Bst.through.time}

\alias{plot.Bst.through.time}

\title{
	plotting phylogenetic turnover through time
}

\description{
	\code{plot.Bst.through.time} calculates and plots phylogenetic turnover through time, given a phylogeny and a species-by-plots matrix.
}

\usage{
	plot.Bst.through.time(phy, sp_plot, outfile = NULL)
}

\arguments{
  \item{phy}{a fully resolved phylogenetic tree, most likely ultrametric}
  \item{sp_plot}{a species-by-plots matrix or dataframe}
  \item{outfile}{an optional .pdf file to which to write output}
}

\details{
  Bst is a measure of local phyletic proximity excess between individuals of distinct species, 
  thereby expressing phylogenetic turnover. This measure considers the balance between 
  phylogenetic diversity within and among plots. 
    
  For significance testing of \code{Bst}, see \code{\link{Bst.permutation}}
}

\value{
	A plot of phylogenetic turnover against splitting times of the phylogeny.  The leftmost value is that for the root node.  Values closer to \code{branching.time=0} 
	are those closer to the present.
}

\references{
  HARDY OJ and B SENTERRE. 2007. Characterizing the 
  phylogenetic structure of communities by an additive partitioning of 
  phylogenetic diversity. Journal of Ecology 95:493-506.
  
  HARDY OJ. 2008. Testing the spatial phylogenetic 
  structure of local communities: statistical performances of 
  different null models and test statistics on a locally neutral 
  community. Journal of Ecology 96:914-926.
}

\author{
Jonathan Eastman}

\seealso{
	\code{\link{Bst.all.nodes}} and \code{\link{spacodi.calc}} for the underlying functions; 
	see see \code{\link{plot.Bst}} for color-coded plotting of \code{Bst} on a tree}

\examples{

foo	<- r.plot(species=200,plots=2,missing.prop=0.15,sim.tree=TRUE)	

# define a random species-by-plots matrix
sp_plot <- foo[[1]]

# define a random tree
phy	<- foo[[2]]

# plot Bst through time
plot.Bst.through.time(phy, sp_plot)
}