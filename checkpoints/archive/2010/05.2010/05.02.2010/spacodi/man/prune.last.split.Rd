\name{prune.last.split}
\alias{prune.last.split}
\title{
pruning the most recent split from a phylogeny }
\description{
\code{prune.last.split} is intended for use with \code{\link{bd.tree}} or \code{\link[geiger]{birthdeath.tree}}. The function finds the last split 
in a tree and randomly selects one of either descendants to prune from the tree.  This last split is usually of zero-length, if simulated with 
the \code{taxa.stop} stopping criterion for tree simulation in either of the above two functions. }
\usage{
prune.last.split(phy)
}
\arguments{
  \item{phy}{a phylogenetic tree of class \code{phylo}; see \code{\link[ape]{read.tree}}}
}

\value{A pruned phylogeny}

\author{Jonathan Eastman}
\seealso{
	\code{\link[geiger]{birthdeath.tree}} or \code{\link{bd.tree}} for tree simulation
}
\examples{
# generate a tree with 101 species as tips
foo <- bd.tree(taxa.stop=101)

# prune the last split from a simulated tree
phy <- prune.last.split(phy=foo)
phy
plot(phy)

}