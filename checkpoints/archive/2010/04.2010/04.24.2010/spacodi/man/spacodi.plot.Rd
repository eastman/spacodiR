\name{spacodi.plot}
\alias{spacodi.plot}
\alias{K.plot}
\alias{Bst.plot}

\title{plotting phylogenetic heritibility or phylogenetic turnover on trees}

\usage{spacodi.plot(phy = NULL, sp_plot = NULL, trait = NULL, color.plot = FALSE, time.plot = FALSE, outfile = NULL, ...)}

\arguments{
  \item{phy}{a fully resolved phylogenetic tree, most sensibly ultrametric; required for all plotting}
  \item{sp_plot}{dataset provided if plotting \code{Bst}}
  \item{trait}{dataset to be supplied if plotting Blomberg's \code{K}}
  \item{color.plot}{a logical identifier indicating whether to generate a color-coded plot of \code{K} or \code{Bst} on the tree}
  \item{time.plot}{a logical identifier indicating whether to generate a scatter plot of \code{K} or \code{Bst} against time}
  \item{outfile}{an optional .pdf file to which to write output}
  \item{\dots}{additional plotting parameters to be optionally supplied}
}

\details{
These plotting functions require a tree and either trait data (for Blomberg's \code{K}) or a species-by-plots matrix (for \code{Bst}).
One of two plots may be generated: a color-coded depiction of values (\code{K} or \code{Bst}) or a scatter plot of values against time.

Based on observed measures of Blomberg's \code{K} for a supplied trait, this function will plot phylogenetic heritability 
of a trait for nodes of a phylogeny, either with color-coded \code{K}-values at nodes or \code{K} against splitting times of the tree. 
The underlying function (\code{\link{K.all.nodes}}) generates estimates of Blomberg's \code{K} for all subtrees of a phylogeny (if the number of tips subtended by a particular node 
exceeds two species.)

If supplied a species-by-plots matrix, this function will compute and plot \code{Bst} either on the tree or against splitting times of nodes.
Bst is a measure of local phyletic proximity excess between individuals of distinct species, thereby expressing phylogenetic turnover. 
This measure considers the balance between phylogenetic diversity within and among plots. The underlying function, \code{\link{Bst.all.nodes}}, 
computes phylogenetic community diversity for nodes of a phylogeny that subtend greater than two species and where nodes subtend species 
that are sampled in multiple plots. 
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

\seealso{\code{\link{K.all.nodes}} for computing Blomberg's \code{K} for subtrees; see \code{\link[picante]{Kcalc}} for the wrapped function; 
see \code{\link{Bst.all.nodes}} and \code{\link{spacodi.calc}} for estimating phylogenetic turnover on trees with community-level sampling}


\examples{
## Blomberg's K ##
# define a random tree
phy=prune.last.split(bd.tree(taxa.stop=201))

# evolve a character across the tree
trt.BM <- as.data.frame(sim.char(phy=phy,as.matrix(0.5))) 
names(trt.BM) = "BM.trait"

# plot K through time, expecting most values to be near 1
dev.new()
spacodi.plot(phy=phy, trait=trt.BM, time.plot=TRUE, cex=0.5)

# plot K with color
dev.new()
spacodi.plot(phy=phy, trait=trt.BM, color.plot=TRUE)

## Bst ##
# generate a species-by-plots matrix, along with a random tree
foo	<- sim.spacodi(species=40,plots=10,missing.prop=0.01,sim.tree=TRUE)	

# define the data matrix
s.plot	<- foo[[1]]
s.plot

# define the tree
phy	<- foo[[2]]

# plot Bst-estimates across the phylogeny							
dev.new()
spacodi.plot(phy=phy, sp_plot=s.plot, color.plot=TRUE)

# plot Bst-estimates against time							
dev.new()
spacodi.plot(phy=phy, sp_plot=s.plot, time.plot=TRUE, cex=0.5)

}
