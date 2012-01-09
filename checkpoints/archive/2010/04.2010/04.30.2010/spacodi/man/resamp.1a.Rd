\name{resamp.1a}

\alias{resamp.1a}

\title{randomizing a species-by-plots matrix: '1a' of Hardy (2008)}

\usage{resamp.1a(obj, abund.class.ratio = 4)}

\arguments{
  \item{obj}{a community-phylogenetic dataset in \code{SPACoDi} format (see \code{\link{as.spacodi}})}
  \item{abund.class.ratio}{a ratio defining abundance classes}
}

\details{A resampling procedure for a species-by-plots matrix, where species are shuffled within abundance classes. 
Species are grouped into distinct abundance classes characterized by a fixed ratio: 
\code{abund.class.ratio} = maximal abundance / minimal abundance.  For instance, if \code{abund.class.ratio = 4},
the limits between abundance classes could be 1, 4, 16, ... .  Species are randomly permuted within each 
class, which maintains most of the abundance phylogenetic structure originally present in a dataset.  
}

\value{A shuffled dataset}

\references{
  HARDY OJ. 2008. Testing the spatial phylogenetic 
  structure of local communities: statistical performances of 
  different null models and test statistics on a locally neutral 
  community. Journal of Ecology 96:914-926.
}

\author{Timothy Paine and Jonathan Eastman}

\seealso{\code{\link{Bst.by.nodes}} for permutation tests of community diversity and additional randomization methods; 
see \code{\link{spacodi.calc}} for the main underlying function}

\examples{
sp.plot <- sim.spacodi(species=10,plots=4,missing.prop=0.15,sim.tree=FALSE)
sp.plot

# shuffle dataset
resamp.1a(obj=sp.plot, abund.class.ratio=3) 
}
