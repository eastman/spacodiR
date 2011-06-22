\name{node.time.extractor}

\alias{node.time.extractor}

\title{finding splitting times for nodes}

\description{\code{node.time.extractor} is a function to collect all nodes (and times) that are within a particular time fraction of the total temporal extent of a phylogeny}

\usage{node.time.extractor(phy, start.time = 0, stop.time = 1, return.times = TRUE, proportion = TRUE)}

\arguments{
  \item{phy}{a fully resolved phylogenetic tree, most sensibly ultrametric}
  \item{start.time}{the temporal slice of the tree defining the lower bound for extracted nodes}
  \item{stop.time}{the temporal slice of the tree defining the upper bound for extracted nodes}
  \item{return.times}{a logical identifier specifying whether times are to be returned in addition to node names}
  \item{proportion}{a logical identifier indicating whether \code{start.time} and \code{stop.time} are to be regarded as proportions of tree height}
}

\details{Given a temporal fraction of a phylogeny, this function will return all contained nodes (and 
branching times if desired).}

\value{
A list of two components
  \item{nodes}{a vector of node names (as interpreted by \code{\link[ape]{nodelabels}})}
  \item{times}{a corresponding vector of splitting times for the above nodes}
}

\author{Jonathan Eastman}

\seealso{\code{\link[ape]{nodelabels}} for the function used to label nodes}

\examples{
# generate a random tree
foo=bd.tree(b=0.02, d=0, taxa.stop=25)

# get nodes and times
node.time.extractor(foo, start.time=0.35, stop.time=0.80, return.times=TRUE, proportion=TRUE)

# visualize the correspondence
plot(foo)
nodelabels()

## an example of comparing observed and expected Bst at nodes ##
sim=sim.spacodi(missing=0.01,sim.tree=TRUE)
s=sim$sp_plot
p=sim$sp_tree

# get nodes for a time-slice
n.t=node.time.extractor(phy=p,start.time=0.25,stop.time=0.75)

# get Bsts for nodes
o=Bst.all.nodes(sp_plot=s, phy=p)
obs=o[!is.na(match(o[,4], n.t$nodes)),1]

# compile expected Bst-values from data reshuffling
Bst.permutation(sp_plot=s, phy=p, n.rep=5, method="1s")->exps
ee.b=cbind(c(exps[,"Bst",]),c(exps[,"node.ID",]))
exp=ee.b[!is.na(match(ee.b[,2], n.t$nodes)),1]

# randomization test of observed and expected Bst
randomization.test(obs=obs, exp=exp, iter=100, two.tailed=TRUE)
}