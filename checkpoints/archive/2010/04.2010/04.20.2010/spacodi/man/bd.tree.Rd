\name{bd.tree}
\alias{bd.tree}
\title{
simulating trees under a birth-death process}
\description{
\code{bd.tree} is a modification of \code{geiger:::birthdeath.tree()}}
\usage{
bd.tree(b = 1, d = 0, time.stop = 0, taxa.stop = 15, return.all.extinct = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b}{
speciation rate}
  \item{d}{extinction rate}
  \item{time.stop}{a stopping criterion -- simulation stops as a function of time elapsed in the tree}
  \item{taxa.stop}{a stopping criterion -- simulation stops as a function of size of simulated tree}
  \item{return.all.extinct}{a logical flag -- whether extinct lineages are included in assessing the \code{taxa.stop} criterion 
}
}

\seealso{
	\code{geiger:::birthdeath.tree} for the original function; 
	see \code{\link{r.plot}} for simulating community phylogenetic datasets (i.e., trees and species-by-plots matrices)
}

