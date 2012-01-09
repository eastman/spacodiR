\name{r.sp.plot}

\alias{r.sp.plot}

\title{randomly generating species-by-plots matrices by simulation from empirical data}

\usage{r.sp.plot(obj, n.rep = 10)}

\arguments{
  \item{obj}{a community-phylogenetic dataset in \code{SPACoDi} format (see \code{\link{as.spacodi}})}
  \item{n.rep}{a number of simulation iterations to perform}
}

\details{A simulation procedure for generating species-by-plots matrices, based upon the original dataset.  
Simulated datasets are equal in dimension to the original dataset. The average proportion of unsampled species 
among all plots is maintained in simulation. Simulated abundances are informed by an assumption of log-normal distributed 
observed abundances.}

\value{a simulated dataset}

\author{Jonathan Eastman}

\seealso{
  \code{\link{Bst.by.nodes}} for permutation tests of community diversity; 
  see \code{\link{sim.spacodi}} for simulating community phylogenetic data without an empirical dataset
}

\examples{
sp.plot <- sim.spacodi(species=10,plots=4,missing.prop=0.15,sim.tree=FALSE)
sp.plot

# shuffle dataset
r.sp.plot(obj=sp.plot, n.rep=3) 
}
