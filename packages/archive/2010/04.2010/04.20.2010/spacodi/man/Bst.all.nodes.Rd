\name{Bst.all.nodes}

\alias{Bst.all.nodes}

\title{computing phylogenetic turnover for all nodes of a phylogeny}

\usage{Bst.all.nodes(phy, sp_plot, return.all = TRUE)}

\arguments{
  \item{phy}{a fully resolved phylogenetic tree, most sensibly ultrametric}
  \item{sp_plot}{a species-by-plots matrix or dataframe}
  \item{return.all}{a logical switch to determine whether Bst should be returned for all nodes or only those for which Bst is sensibly calculated}
}

\details{
  Bst is a measure of local phyletic proximity excess between individuals of distinct species, 
  thereby expressing phylogenetic turnover.  This measure considers the balance between 
  phylogenetic diversity within and among plots. 
  
  \code{Bst.all.nodes} computes phylogenetic community diversity for nodes of 
  a phylogeny that subtend greater than two species and where nodes subtend species 
  that are sampled in multiple plots. 
    
  For significance testing, see \code{\link{Bst.permutation}}
}
  
\value{  
  \code{Bst.all.nodes} returns a two-dimensional array of community diversity measures (Bst) 
  from permuted datasets.  The first dimension of the array is determined
  by is the number of internal nodes of the supplied phylogeny. The second dimension 
  of the output array are the measures returned to the user: \code{Bst} -- computed for each node (if nodes 
  subtend more than two sampled species and if subtended species were sampled in multiple 
  plots); \code{tips} -- the number of contained species at each node; \code{node.time} -- the branching time of the
  node in the supplied phylogeny (in same units as those supplied); and \code{node.ID} -- the identifier for 
  each node (as interpreted by \code{ape:::nodelabels()}). 
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

\author{Timothy Paine and Jonathan Eastman}

\seealso{\code{\link{spacodi.calc}} for interpretation of results; see \code{\link{Bst.permutation}} for computing \code{Bst} for permuted datasets.}

\examples{
# generate a random dataset and phylogeny 
foo <- r.plot(species=35,plots=10,missing.prop=0.02,sim.tree=TRUE)
sp_plot <- foo[[1]]
phy <- foo[[2]]

# compute Bst for all nodes
Bst.all.nodes(sp_plot=sp_plot, phy=phy, return.all=TRUE)

# compute Bst for all nodes for which Bst is calculable
Bst.all.nodes(sp_plot=sp_plot, phy=phy, return.all=FALSE)
}
