\name{as.spacodi}

\alias{as.spacodi}

\title{converting between data formats for community phylogenetics}

\usage{as.spacodi(data, outfile = NULL)}

\arguments{
  \item{data}{a species-by-plots matrix}
  \item{outfile}{an optional text file to which to write output}
}

\details{This utility converts a species-by-plots matrix into a format readable by the external program \code{SPACoDi}, a Windows executable. 
The user has the option to save an output file, defined by \code{outfile}. 
}

\value{an array, formatted for use in \code{SPACoDi} for Windows}

\author{Jonathan Eastman}

\seealso{\code{\link{as.phylocom}} for converting between \code{phylocom} and \code{SPACoDi} formats;
see \code{\link{spacodi.calc}} for community diversity measures of empirical data}

\examples{
# generate a species-by-plots matrix
foo <- r.plot(species=10,plots=6,missing.prop=0.15,sim.tree=FALSE)	

# convert to phylocom format
as.phylocom(foo) -> p.foo
p.foo

# convert back to spacodi format and write to file
as.spacodi(p.foo, outfile="spacodi.formatted.txt") -> s.foo
s.foo		# re-converted SPACoDi dataset
foo		# comparison with the original dataset

}