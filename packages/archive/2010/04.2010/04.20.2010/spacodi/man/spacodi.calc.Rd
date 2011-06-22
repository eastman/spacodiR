\name{spacodi.calc}
\alias{spacodi.calc}

\title{measuring spatial and phylogenetic diversity of communities}

\description{To describe community structure, this function will compute 
various statistics describing phylogenetic (or trait) diversity of plots 
and phylogenetic (or trait) distances between plots in a pairwise fashion}

\usage{spacodi.calc(sp_plot, phy = NULL, sptraits = NA, which.plots = NA)}

\arguments{
  \item{sp_plot}{a species-by-plots matrix or dataframe}
  \item{phy}{a fully resolved phylogenetic tree, most sensibly ultrametric}
  \item{sptraits}{species-by-trait(s) matrix or dataframe}
  \item{which.plots}{a subset of plots to be included}
}

\details{
  \code{spacodi.calc} returns statistics of pairwise differentiation of plots
  (subcommunities), based either on trait or phylogenetic distances.    
  
  The following statistics may be computed:
  \itemize{
  
    \item \code{Ist} a measure of local species identity excess between 
	individuals, expressing species turnover. It is a form of spatial 
	partition of Gini-Simpson diversity (equivalent to \code{Fst} in population 
	genetics). \code{Ist} considers only the species-by-plot matrix. 
	
    \item \code{Pst}: a measure of local phyletic proximity excess between 
	individuals, expressing species + phylogenetic turnover. It is a form 
	of spatial partition of Rao's quadratic entropy (equivalent to \code{Nst} in 
	population genetics). \code{Tst} is the analogue for trait data.
	
	\item \code{Bst}: a measure of local phyletic proximity excess between 
	individuals of distinct species, expressing phylogenetic turnover 
	(independently of species turnover). \code{T*st} is the analogue for 
	trait data.
	
	\item \code{PIst}: Pst analogue for presence/absence data, expressing
	phylogenetic turnover (independently of species turnover). \code{TAUst}
	is the analogue for trait data.
  }
  \code{spatial phylogenetic clustering}: species within plots are more 
  related on average than species from distinct plots where \code{ Pst > Ist},
  \code{ Bst > 0}, or  \code{ PIst > 0}. 
  
  \code{spatial phylogenetic overdispersion}: species within plots are less
  related on average than species from distinct plots where \code{ Pst < Ist},
  \code{ Bst < 0}, or  \code{ PIst < 0}.
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

\author{Olivier Hardy and Timothy Paine}

\examples{
# generate a species-by-plots matrix, along with a random tree
foo	<- r.plot(species=60,plots=10,missing.prop=0.15,sim.tree=TRUE)	

# define the data matrix
sp_plot	<- foo[[1]]
sp_plot

# define the tree
phy	<- foo[[2]]

trt.BM <- as.data.frame(sim.char(phy,as.matrix(10))) 
names(trt.BM) = "BM.trait"
# simulated trait evolving under Brownian motion							
trt.BM

# community diversity Hardy and Senterre (2007) statistics: tree-based
spacodi.calc(sp_plot = sp_plot, phy = phy)

# community diversity Hardy and Senterre (2007) statistics: trait-based
spacodi.calc(sp_plot = sp_plot, sptraits = trt.BM)

# Ist: using abundance data only				
spacodi.calc(sp_plot = sp_plot)								
}
