\name{plot.Bst}
\alias{plot.Bst}
\title{
plotting color-coded phylogenetic turnover for nodes of a phylogeny}
\description{
  Bst is a measure of local phyletic proximity excess between individuals of distinct species, 
  thereby expressing phylogenetic turnover. This measure considers the balance between 
  phylogenetic diversity within and among plots. 
  
  \code{Bst.all.nodes} computes phylogenetic community diversity for nodes of 
  a phylogeny that subtend greater than two species and where nodes subtend species 
  that are sampled in multiple plots. 
}

\usage{
plot.Bst(phy, sp_plot, outfile = NULL)
 }
 \arguments{
 \item{phy}{a fully resolved phylogenetic tree, most sensibly ultrametric}
 \item{sp_plot}{a species-by-plots matrix or dataframe}
 \item{outfile}{an optional .pdf file to which to write output
 }
}

\details{
  Bst is a measure of local phyletic proximity excess between individuals of distinct species, 
  thereby expressing phylogenetic turnover. This measure considers the balance between 
  phylogenetic diversity within and among plots.
  
  Based on observed measures of \code{Bst}, this function plots phylogenetic turnover for nodes of a phylogeny, with color-coded values}
\value{a plot of \code{Bst}-values for a phylogeny.
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
Timothy Paine and Jonathan Eastman}

\seealso{
\code{\link{spacodi.calc}} for community diversity measures; 
see \code{\link{Bst.all.nodes}} for computing \code{Bst} for subtrees; 
see \code{\link{plot.Bst.through.time}} for plotting \code{Bst} against splitting times of a tree}

\examples{
# generate a species-by-plots matrix, along with a random tree
foo	<- r.plot(species=60,plots=10,missing.prop=0.05,sim.tree=TRUE)	

# define the data matrix
sp_plot	<- foo[[1]]
sp_plot

# define the tree
phy	<- foo[[2]]

# plot Bst-estimates across the phylogeny							
plot.Bst(phy, sp_plot)

}
