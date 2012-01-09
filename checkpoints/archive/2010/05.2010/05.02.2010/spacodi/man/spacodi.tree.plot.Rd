\name{spacodi.tree.plot}
\alias{spacodi.tree.plot}
\alias{K.plot}
\alias{Bst.plot}

\title{plotting phylogenetic heritibility or phylogenetic turnover on trees}

\usage{spacodi.tree.plot(spacodi.object, phy, sig.plot = TRUE, cut.off = 0.05, cols = list("white", "gray", "black"), main = TRUE, outfile = NULL, add.id = FALSE, ...)}

\arguments{
  \item{spacodi.object}{a set of \code{K} or \code{Bst} estimates, from \code{\link{K.by.nodes}} or \code{\link{Bst.by.nodes}}}
  \item{phy}{a phylogenetic tree of class \code{phylo}, most sensibly ultrametric; see \code{\link[ape]{read.tree}}}
  \item{sig.plot}{a logical indicator dictating whether nodes are colored by significance of observed and expected \code{K}s or \code{Bst}s}
  \item{cut.off}{a value, if \code{sig.plot=TRUE}, distinguishing observed \code{K}s from expected \code{K}s or \code{Bst}s}
  \item{cols}{a list of three elements if using \code{sig.plot}: the first color is for values not different than expected; the second and third for values greater and less than expected}
  \item{main}{a logical indicator determining whether a title is plotted}
  \item{outfile}{an optional .pdf file to which to write output}
  \item{add.id}{a logical indicator dictating whether \code{node.labels} are placed near nodes; see \code{\link[ape]{nodelabels}}}
  \item{\dots}{additional plotting parameters to be optionally supplied}
}

\details{
Based on observed measures of Blomberg's \code{K} for a supplied trait, this function will plot phylogenetic heritability 
of a trait for nodes of a phylogeny, with color-coded \code{K}-values at nodes. Alternatively, this function will compute and plot \code{Bst} on the tree, with color-coded values.
Bst is a measure of local phyletic proximity excess between individuals of distinct species, thereby expressing phylogenetic turnover. 
This measure considers the balance between phylogenetic diversity within and among plots. 

Note: this function requires a \code{spacodi.object} with data for all nodes.  If using, for instance, \code{Bst.by.nodes} to
generate these data, the option \code{obs.only=TRUE} will ensure that data are returned for all nodes (see example below).

\code{spacodi.tree.plot} will return an error if \code{sig.plot=TRUE}, a plot of departure of \code{Bst} or \code{K} from expectations, if the \code{spacodi.object}
does not include expected values and result from a \code{randomization.test} as elements in the list (see \code{\link{Bst.by.nodes}} and \code{\link{K.by.nodes}}).  
}

\value{a plot of \code{K}- or \code{Bst}-values for a phylogeny, with the option to direct output to a .pdf file}

\references{
  BLOMBERG SP, T GARLAND, and AR IVES. 2003. Testing for phylogenetic signal in comparative data: Behavioral traits are more labile. Evolution 57:717-745.

  HARDY OJ and B SENTERRE. 2007. Characterizing the 
  phylogenetic structure of communities by an additive partitioning of 
  phylogenetic diversity. Journal of Ecology 95:493-506.
  
  HARDY OJ. 2008. Testing the spatial phylogenetic 
  structure of local communities: statistical performances of 
  different null models and test statistics on a locally neutral 
  community. Journal of Ecology 96:914-926.
}

\author{Timothy Paine and Jonathan Eastman}

\seealso{\code{\link{K.by.nodes}} for computing Blomberg's \code{K} for subtrees; 
see \code{\link{Bst.by.nodes}} for estimating phylogenetic turnover on trees with community-level sampling; 
\code{\link{K.simulation.plot}}; \code{\link{Bst.permutation.plot}}; \code{\link{phy.dot.plot}} }


\examples{
## Blomberg's K ##
# define a random tree
phy=prune.last.split(bd.tree(taxa.stop=41))

# evolve a character across the tree
trt.BM <- as.data.frame(sim.char(phy=phy,as.matrix(0.5))) 
names(trt.BM) = "BM.trait"

# plot K with color
K=K.by.nodes(sp.traits=trt.BM, phy=phy, obs.only=TRUE)
spacodi.tree.plot(spacodi.object=K, phy=phy, sig.plot=FALSE)

# plot Bst 
spl=sim.spacodi(species=40, plots=4, sim.tree=FALSE)
phy$tip.label=row.names(spl)
B=Bst.by.nodes(sp.plot=spl, phy=phy, return.all=TRUE, method="1s")
Bst.permutation.plot(B, sig.plot=TRUE, add.id=FALSE, edge.width=0.1)

}
