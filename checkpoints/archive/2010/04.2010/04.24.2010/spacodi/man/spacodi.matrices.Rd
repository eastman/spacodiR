\name{spacodi.matrices}

\alias{spacodi.matrices}
\alias{Bst.matrix.for.plots}

\title{generating pairwise-distance matrices of phylogenetic diversity: Ist, Pst, Bst, and PIst}

\description{
  This function automates the process of pairwise computation of measures of phylogenetic 
  diversity between all plots of a dataset. Average measures of diversity, from each plot to every other, 
  are assembled into a distance matrix for use in downstream analyses (e.g., Mantel tests). 
}

\usage{spacodi.matrices(sp_plot, phy)}

\arguments{
  \item{sp_plot}{a species-by-plots matrix or dataframe}
  \item{phy}{a fully resolved phylogenetic tree, most sensibly ultrametric}
}

\details{
  Pairwise values for \code{Ist}, \code{Pst}, \code{Bst}, \code{PIst} are estimated for every 
  pair of unique plots using \code{\link{spacodi.calc}}. If fewer than two species are sampled for
  particular plots, these plots are dynamically removed from the dataset (as phylogenetic diversities are 
  incalculable in such cases).
}

\value{A list of four distance matrices (measures are further explained in the \code{\link{spacodi.calc}} documentation):
 \itemize{
 \item{\code{Ist} }{a distance matrix of average species turnover between plots}
 \item{\code{Pst} }{a distance matrix of average phylogenetic+species turnover between plots}
 \item{\code{Bst} }{a distance matrix of average phylogenetic turnover between plots}
 \item{\code{PIst} }{a distance matrix of \code{Pst}s, for presence|absence of species in plots}
 }
}


\author{Jonathan Eastman}

\seealso{\code{\link{spacodi.calc}} for interpretation of results; see \code{\link[vegan]{mantel}} for 
statistical comparison of dissimilarity matrices}

\examples{
# load a community data matrix, along with a tree
data(spacodi.example)
attributes(spacodi.example)
attach(spacodi.example)
sp_plot
sp_tree

# randomly generate an 'environmental' variable
env.fake <- data.frame(runif(ncol(sp_plot))) 
names(env.fake) = "fake.env"
row.names(env.fake)=names(sp_plot)
e.f=dist(env.fake)

# generate phylogenetic diversity matrices
s.out=spacodi.matrices(sp_plot=sp_plot, phy=sp_tree)

# conduct Mantel tests, comparing pairwise distances between phylogenetic and environmental diversity
for(ii in 1:length(s.out)) {
	print(paste("### a test of environmental variation and ", names(s.out)[ii], sep=""))
	print(mantel(e.f, s.out[[ii]]))
	}
}