\name{plot.K}
\alias{plot.K}
\title{
plotting color-coded Blomberg's K for nodes of a phylogeny}
\description{This function generates estimates of Blomberg's \code{K} for all subtrees of a phylogeny (if the number of tips subtended by a particular node exceeds two species.)

Blomberg's \code{K} represents a measure of trait similarity among tips of a phylogeny.  For traits evolving under Brownian motion, Blomberg's \code{K} 
is expected to take a value near 1.  If common ancestry is less predictive of species' traits, \code{K} is near 0.  Values of Blomberg's \code{K} near 0 may result from a 
scenario where species are strongly 'drawn' to the same trait optimum (e.g., under an Ornstein-Uhlenbeck process), or where species' trait-values are not 
phylogenetically heritable (that is, traits do not evolve in a treeward fashion, where related species would be most likely to exhibit similar trait values). 

}
\usage{
plot.K(phy, traits, outfile = NULL)
}
\arguments{
 \item{phy}{a fully resolved phylogenetic tree, most sensibly ultrametric}
 \item{traits}{a dataframe of a trait (as the column), with species as row.names}
 \item{outfile}{an optional .pdf file to which to write output
}
}
\details{
Based on observed measures of Blomberg's \code{K}, this function plots phylogenetic heritability of the trait for nodes of a phylogeny, with color-coded values}

\value{a plot of \code{K}-values for a phylogeny
}
\references{
BLOMBERG SP, T GARLAND, and AR IVES. 2003. Testing for phylogenetic signal in comparative data: Behavioral traits are more labile. Evolution 57:717-745.
}

\author{
Timothy Paine and Jonathan Eastman}

\seealso{
\code{\link{K.all.nodes}} for computing Blomberg's \code{K} for subtrees; see \code{picante:::Kcalc()} for the wrapped function}
