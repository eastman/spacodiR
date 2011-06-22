\name{node.time.extractor}

\alias{node.time.extractor}

\title{finding splitting times for nodes}

\description{\code{node.time.extractor} is a function to collect all nodes that are within a particular time fraction of the total temporal extent of a phylogeny}

\usage{node.time.extractor(phy, time.fraction = 0.5, return.times = FALSE)}

\arguments{
  \item{phy}{a fully resolved phylogenetic tree, most sensibly ultrametric}
  \item{time.fraction}{the temporal fraction of the tree for which to find contained nodes}
  \item{return.times}{a logical identifier specifying whether times are to be returned in addition to node names}
}

\details{Given a temporal fraction of a phylogeny (e.g., \code{time.fraction=0.5}, implying half the tree's duration), this function will return all contained nodes.}

\value{
A list of two components
  \item{output[[1]] }{a vector of node names (as interpreted by \code{ape:::nodelabels()}}
  \item{output[[2]] }{a corresponding vector of splitting times for the above nodes}
}

\author{Jonathan Eastman}

\seealso{\code{ape:::nodelabels()} for the function that names nodes}

\examples{
# generate a random tree
foo <- bd.tree(b=0.02, d=0, taxa.stop=25)

# get nodes and times
node.time.extractor(foo, time.fraction=0.35, return.times=TRUE)

# visualize the output
plot(foo)
nodelabels()

}