\name{write.spacodi.data}

\alias{write.spacodi.data}

\title{saving community phylogenetics datasets to a file}

\usage{write.spacodi.data(data, outfile)}
\arguments{
  \item{data}{a community-phylogenetics dataset formatted for the R-package \pkg{spacodi} (see \code{\link{as.spacodi}})}
  \item{outfile}{a formatted file for the Windows executable \code{SPACoDi}}
}

\details{This utility writes a species-by-plots matrix into a format readable 
by the external program \code{SPACoDi}, a Windows executable (available at 
  \url{http://ebe.ulb.ac.be/ebe/Software.html}).  
}

\author{Jonathan Eastman}

\seealso{\code{\link{as.spacodi}} for converting from \code{phylocom} to \pkg{spacodi} formats;
see \code{\link{spacodi.calc}} for community diversity measures of empirical data
}

\examples{
# generate a community-phylogenetics dataset
ss.plot	<- sim.spacodi(species=60,plots=10,missing.prop=0.15,sim.tree=FALSE)	

# save the dataset to working directory
write.spacodi.data(ss.plot, outfile="spacodi.formatted.txt")

}