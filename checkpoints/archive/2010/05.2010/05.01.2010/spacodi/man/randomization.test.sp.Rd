\name{randomization.test.sp}
\alias{randomization.test.sp}
\title{testing significant by randomization}
\description{\code{randomization.test.sp} is a function for conducting significance-testing by randomization.}
\usage{randomization.test.sp(obs = obs, exp = exp, iter = 10000, return.all = FALSE, two.tailed = TRUE)}

\arguments{
  \item{obs}{an vector of observed data}
  \item{exp}{an vector of expected values}
  \item{iter}{number of randomization iterations to perform}
  \item{return.all}{a logical identifier controlling whether all randomization comparisons are output into a vector}
  \item{two.tailed}{a logical identifier indicating whether randomization test should be two-tailed}
}

\details{Given a set of observed and expected values, \code{iter} comparisons between the vectors are performed.  
Expecting the sign of these comparisons to be random, if a large proportion of comparisons are in a particular 
direction, support is lent for a true difference in means.  For instance, if 90 percent of comparisons suggest 
that observed data are larger than expected, an approximate randomization \code{p}-value is interpreted to be 0.20. If 
a single value is supplied for the observed data, this method will be akin to finding the quantile (\code{q})in the expected values
where the observed datum appears. Under \code{two.tailed=TRUE}, the resulting \code{p.value} from this function will be twice (\code{q}), 
for an observed datum in the left tail of the expected distribution, or twice (\code{1-q}) for the opposite tail. }

\value{ A list of at most two primary elements:
    \itemize{
    \item \code{differences}{: the vector of all comparisons between \code{obs} and \code{exp} if \code{return.all=TRUE}}
	\item \code{p.value}{: result(s) from the randomization test; if \code{two.tailed=FALSE}, \code{p}-values from both tails are returned}
	}
	Note: be wary of significant results that are based off few comparisons. If the number of iterations exceeds the length of the
	empirical data, a warning will be issued.
}

\author{Jonathan Eastman}

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