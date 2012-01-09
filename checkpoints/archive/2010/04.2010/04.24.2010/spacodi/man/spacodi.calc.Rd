\name{spacodi.calc}
\alias{spacodi.calc}
\alias{Ist}
\alias{Pst}
\alias{PIst}
\alias{Bst}
\alias{Tst}
\alias{TAUst}
\alias{T*st}


\title{measuring spatial and phylogenetic diversity of communities }

\description{To describe community structure, this function will compute 
various statistics describing phylogenetic (or trait) diversity of plots 
and phylogenetic (or trait) distances between plots in a pairwise fashion}

\usage{spacodi.calc(sp_plot, phy = NULL, traits = NA, all.together = TRUE, prune.tree = FALSE, prune.plot = TRUE)}

\arguments{
  \item{sp_plot}{a species-by-plots matrix or dataframe}
  \item{phy}{a fully resolved phylogenetic tree, most sensibly ultrametric}
  \item{traits}{species-by-trait(s) array or dataframe}
  \item{all.together}{a logical identifier indicating whether to treat all traits together or separately}
  \item{prune.tree}{a logical identifier indicating whether to prune tree of species unsampled in \code{sp_plot}}
  \item{prune.plot}{a logical identifier indicating whether to dynamically prune plot data if fewer than two species were sampled}
}

\details{
  An \code{sp_plot} of numerical abundance, relative abundance, or presence|absence 
  can be handled by \code{spacodi.calc}, which returns statistics of pairwise 
  differentiation of plots (subcommunities), based either on trait or phylogenetic distances.
  
  If a phylogeny is supplied with species that are not present in the community dataset (i.e., \code{sp_plot}), 
  the user has the option to dynamically prune the tree to match the taxomomic set in \code{sp_plot}; use 
  \code{prune.tree=TRUE} for this option. If \code{prune.tree=FALSE}, the function will abort when unmatched 
  datasets are detected. 
  
  If fewer than two species are sampled in any particular plot, \code{spacodi.calc} will fail.  To allow the 
  dataset to be pruned of offending plots, choose \code{prune.plot=TRUE} (the default setting).
    
  The following statistics may be computed:
  \itemize{
  
    \item \code{Ist} a measure of local species identity excess between 
	individuals, expressing species turnover. It is a form of spatial 
	partition of Gini-Simpson diversity (equivalent to \code{Fst} in population 
	genetics). \code{Ist} considers only the species-by-plots matrix. 
	
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
  PHYLOGENETIC DIVERSITY
  \itemize{
    \item \code{spatial phylogenetic clustering}: species within plots are more 
	related on average than species from distinct plots where \code{ Pst > Ist},
	\code{ Bst > 0}, or  \code{ PIst > 0}. 
  
    \item \code{spatial phylogenetic overdispersion}: species within plots are less 
	related on average than species from distinct plots where \code{ Pst < Ist}, 
	\code{ Bst < 0}, or  \code{ PIst < 0}.
  }
  
  TRAIT DIVERSITY
  \itemize{
    \item If \code{all.together=TRUE}, all traits will be used to generate distance a distance 
	matrix for sampled species.  Where \code{all.together=FALSE} is used, a resulting 
	dataframe will be output, with measures of \code{Ist}, \code{Tst}, \code{T*st}, and 
	\code{TAUst} for each trait separately. 
  }
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

\author{Olivier Hardy, Timothy Paine, and Jonathan Eastman}

\examples{
# load a species-by-plots matrix, along with a tree
data(spacodi.example)
attributes(spacodi.example)
attach(spacodi.example)
sp_plot
sp_tree

# evolve a trait
trt.BM <- as.data.frame(sim.char(sp_tree,as.matrix(10))) 
names(trt.BM) = "BM.trait"
# simulated trait evolving under Brownian motion							
trt.BM

# community diversity Hardy and Senterre (2007) statistics: tree-based
spacodi.calc(sp_plot = sp_plot, phy = sp_tree)

# community diversity Hardy and Senterre (2007) statistics: trait-based
spacodi.calc(sp_plot = sp_plot, traits = trt.BM, all.together=TRUE)

# community diversity for a pair of traits
spacodi.calc(sp_plot = sp_plot, traits = sp_traits, all.together=TRUE)

# community diversity for a pair of traits, each singly
spacodi.calc(sp_plot = sp_plot, traits = sp_traits, all.together=FALSE)

# Ist: using abundance data only				
spacodi.calc(sp_plot = sp_plot)	

# calculations with missing taxa between tree and sp_plot
# excluding the last five species in sp_plot
spacodi.calc(sp_plot = sp_plot[1:15,], phy = sp_tree, prune.tree=TRUE) # carried out with warning

							
}
