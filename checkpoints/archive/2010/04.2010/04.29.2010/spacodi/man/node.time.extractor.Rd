\name{node.time.extractor}

\alias{node.time.extractor}

\title{finding splitting times for nodes}

\description{\code{node.time.extractor} is a function to collect all nodes (and times) that are within a particular time fraction of the total temporal extent of a phylogeny}

\usage{node.time.extractor(phy, start.time = 0, stop.time = 1, return.times = TRUE, proportion = TRUE, ape.ID = TRUE)}

\arguments{
  \item{phy}{a fully resolved phylogenetic tree, most sensibly ultrametric}
  \item{start.time}{the temporal slice of the tree defining the lower bound for extracted nodes}
  \item{stop.time}{the temporal slice of the tree defining the upper bound for extracted nodes}
  \item{return.times}{a logical identifier specifying whether times are to be returned in addition to node names}
  \item{proportion}{a logical identifier indicating whether \code{start.time} and \code{stop.time} are to be regarded as proportions of tree height}
  \item{ape.ID}{a logical flag indicating whether node labels are returned as numbers (from \code{\link[ape]{nodelabels}}) or as preexisting node labels (\code{ape.ID=FALSE})}
}

\details{Given a temporal fraction of a phylogeny, this function will return all contained nodes (and 
branching times if desired). If used with \code{\link{Bst.by.nodes}}, one could conduct significance testing of \code{Bst} through 
particular time slices of the phylogeny. Note: when using \code{proportion=TRUE}, '0' is regarded as the present, '1' as the root of the tree (e.g., \code{start.time=1.0} and 
\code{stop.time=0.8} gathers nodes for the basalmost or most rootward 20 percent of nodes).}

\value{
A list of two components:
 \itemize{
  \item \code{nodes}{a vector of node names (as interpreted by \code{\link[ape]{nodelabels}})}
  \item \code{times}{a corresponding vector of splitting times for the above nodes}
 }
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
s=sim$sp.plot
p=sim$sp.tree

# get nodes for a time-slice
n.t=node.time.extractor(phy=p,start.time=0.25,stop.time=0.75)

# get Bsts for nodes
o=Bst.by.nodes(sp.plot=s, phy=p, obs.only=TRUE)[[1]]
obs=o[!is.na(match(o[,"node.ID"], n.t$nodes)),1]

# compile expected Bst-values from data reshuffling
Bst.by.nodes(sp.plot=s, phy=p, n.rep=5, method="1s")[[2]]->exps
e=exps[!is.na(match(row.names(exps), n.t$nodes)),]
exp=unname(unlist(e))

# randomization test of observed and expected Bst
randomization.test.sp(obs=obs, exp=exp, iter=100, two.tailed=TRUE)
}