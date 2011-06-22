\name{spacodi.matrices}

\alias{spacodi.matrices}

\title{generating pairwise-distance matrices of phylogenetic diversity: Ist, Pst, Bst, and PIst}

\description{
  This function automates the process of pairwise computation of measures of phylogenetic (or trait)
  diversity between all plots of a dataset. Average measures of diversity, from each plot to every other, 
  are assembled into a distance matrix for use in downstream analyses (e.g., Mantel tests). 
}

\usage{spacodi.matrices(sp.plot, phy = NULL, sp.traits = NULL)}

\arguments{
  \item{sp.plot}{a community-phylogenetic dataset in \code{SPACoDi} format (see \code{\link{as.spacodi}})}
  \item{phy}{a fully resolved phylogenetic tree, most sensibly ultrametric}
  \item{sp.traits}{a species-by-trait(s) dataframe}
}

\details{
  Either (but not both) a \code{phy} or \code{sp.traits} object is to be supplied along with the community dataset, \code{sp.plot}.
  If a tree is given, pairwise values for \code{Ist}, \code{Pst}, \code{Bst}, \code{PIst} are estimated for every 
  pair of unique plots using \code{\link{spacodi.calc}}. If fewer than two species are sampled for
  particular plots, these plots are dynamically removed from the dataset (as phylogenetic diversities are 
  incalculable in such cases). Analogously, matrices of pairwise \code{Ist}, \code{Tst}, \code{T*st}, \code{TAUst} are generated 
  if a traits matrix is supplied. 
}

\value{A list of four distance matrices (the following measures are further explained in the \code{\link{spacodi.calc}} documentation):
 \itemize{
  \item{\code{Ist}}{: a distance matrix of average species turnover between plots}
  \item{\code{Pst}}{: a distance matrix of average phylogenetic+species turnover between plots}
  \item{\code{Bst}}{: a distance matrix of average phylogenetic turnover between plots}
  \item{\code{PIst}}{: a distance matrix of \code{Pst}s, for presence | absence of species in plots}
 }
 See \code{\link{as.spacodi}} for the trait-based analogues of these measures. 
}

\references{
  HARDY OJ and B SENTERRE. 2007. Characterizing the 
  phylogenetic structure of communities by an additive partitioning of 
  phylogenetic diversity. Journal of Ecology 95:493-506.
}

\author{Jonathan Eastman}

\seealso{\code{\link{spacodi.calc}} for interpretation of results; see \code{\link[vegan]{mantel}} for 
statistical comparison of dissimilarity matrices}

\examples{
# load a community data matrix, along with a tree
data(sp.example)
attributes(sp.example)
attach(sp.example)
sp.plot
sp.tree
sp.traits

# randomly generate an 'environmental' variable
env.fake <- data.frame(runif(ncol(sp.plot))) 
names(env.fake) = "fake.env"
row.names(env.fake)=names(sp.plot)
e.f=dist(env.fake)

# generate phylogenetic diversity matrices
s.out=spacodi.matrices(sp.plot=sp.plot, phy=sp.tree)
t.out=spacodi.matrices(sp.plot=sp.plot, sp.traits=sp.traits)

# conduct Mantel tests, comparing pairwise distances between phylogenetic and environmental diversity
for(ii in 1:length(s.out)) {
	print(paste("### a test of environmental variation and ", names(s.out)[ii], sep=""))
	print(mantel(e.f, s.out[[ii]]))
	}
	
for(ii in 2:length(s.out)) {
	print(paste("### testing trait variation and phylogenetic diversity: ", names(s.out)[ii], " and ", names(t.out)[ii], sep=""))
	print(mantel(t.out[[ii]], s.out[[ii]]))
	}
	
	
}