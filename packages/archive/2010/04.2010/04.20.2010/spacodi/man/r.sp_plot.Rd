\name{r.sp_plot}

\alias{r.sp_plot}

\title{randomly generating species-by-plots matrices by simulation}

\usage{r.sp_plot(obj, n.rep = 10)}

\arguments{
  \item{obj}{a species-by-plots matrix}
  \item{n.rep}{a number of simulation iterations to perform}
}

\details{A simulation procedure for generating species-by-plots matrices, based upon the original dataset.  
Simulated datasets are equal in dimension to the original dataset. The average extent of unsampled species 
among all plots is maintained. Simulated abundances are informed by an assumption of log-normal distributed 
actual abundances.}

\value{a simulated dataset}

\author{Jonathan Eastman}

\seealso{
  \code{\link{Bst.permutation}} for permutation tests of community diversity; 
  see \code{\link{r.plot}} for simulating community phylogenetic data without an empirical dataset
}

\examples{
sp_plot <- r.plot(species=10,plots=4,missing.prop=0.15,sim.tree=FALSE)
sp_plot

# shuffle dataset
r.sp_plot(obj=sp_plot, n.rep=3) 
}
