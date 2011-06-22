\name{randomization.test}
\alias{randomization.test}
\title{testing significant by randomization}
\description{\code{randomization.test} is a function for conducting significance-testing by randomization.}
\usage{randomization.test(obs = obs, exp = exp, iter = 10000, verbose = FALSE, two.tailed = TRUE)}
\arguments{
  \item{obs}{an array of observed data}
  \item{exp}{an array of expected values}
  \item{iter}{number of randomization iterations to perform}
  \item{verbose}{a logical identifier controlling verbosity of output}
  \item{two.tailed}{a logical identifier indicating whether test should be two-tailed}
}
\details{Given a set of observed and expected values, \code{iter} comparisons between the vectors are performed.  
Expecting the sign of these comparisons to be random, if a large proportion of comparisons are in a particular 
direction, support is lent for a true difference in means.  For instance, if 90 percent of comparisons suggest 
that observed data are larger than expected, an approximate randomization \code{p}-value is interpreted to be 0.20.}

\value{an approximate \code{p}-value}
\author{Jonathan Eastman}

\examples{

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