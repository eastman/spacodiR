\name{resamp.2x}

\alias{resamp.2x}

\title{randomizing a community phylogenetics matrix: '2x' of Hardy (2008)}

\usage{resamp.2x(obj, level=0.1)}

\arguments{
  \item{obj}{a community-phylogenetic dataset in \code{SPACoDi} format (see \code{\link{as.spacodi}})}
  \item{level}{a proportion specifying the extent of data shuffling}
}

\details{
  A resampling procedure for a species-by-plots matrix, based on Gotelli swapping. 
  Shuffles abundances within a pair of plots and for a pair of species.  The \code{level} defines 
  the degree of sampling, with larger values dictating a higher level of reshuffling.  For instance, 
  if \code{level = 0.4} and the dataset involves 5 species and 10 plots, a total of 20 (\code{0.4*5*10}) 
  Gotelli swaps are performed.
}

\value{A shuffled dataset}

\references{
  HARDY OJ. 2008. Testing the spatial phylogenetic 
  structure of local communities: statistical performances of 
  different null models and test statistics on a locally neutral 
  community. Journal of Ecology 96:914-926.
  
  GOTELLI NJ. 2000. Null model analysis of species co-occurrence patterns. Ecology 81:2606–2621.
}

\author{Jonathan Eastman}

\seealso{\code{\link{Bst.by.nodes}} for permutation tests of community diversity and additional randomization methods; 
see \code{\link{spacodi.calc}} for the main underlying function; \code{\link[picante]{randomizeMatrix}}}

\examples{
sp_plot <- sim.spacodi(species=10,plots=4,missing.prop=0.15,sim.tree=FALSE)
sp_plot

# shuffle dataset
resamp.2x(obj=sp_plot, level=0.2) 
}
