\name{resamp.2s}

\alias{resamp.2s}

\title{randomizing a community phylogenetics matrix: '2s' of Hardy (2008)}
\description{\code{resamp.2s} is used for resampling data within a community dataset}

\usage{resamp.2s(obj)}

\arguments{
  \item{obj}{a community dataset in \pkg{spacodiR} format (see \code{\link{as.spacodi}})}
}

\details{A resampling procedure for a species-by-plots matrix, where observed abundances are shuffled across species but within plots.}

\value{A shuffled dataset}

\references{
  HARDY OJ. 2008. Testing the spatial phylogenetic 
  structure of local communities: statistical performances of 
  different null models and test statistics on a locally neutral 
  community. Journal of Ecology 96:914-926.
}

\author{Jonathan Eastman}

\seealso{\code{\link{spacodi.by.nodes}} for permutation tests of community diversity and additional randomization methods; 
see \code{\link{spacodi.calc}} for the main underlying function; \code{\link[picante]{randomizeMatrix}}}

\examples{
data(sp.example)
attach(sp.example)
spl

# shuffle dataset
resamp.2s(obj=spl) 
}
