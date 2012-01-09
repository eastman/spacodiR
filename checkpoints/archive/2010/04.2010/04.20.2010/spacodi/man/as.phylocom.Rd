\name{as.phylocom}

\alias{as.phylocom}

\title{converting between data formats for community phylogenetics}

\usage{as.phylocom(data, outfile = NULL)}
\arguments{
  \item{data}{a species-by-plots matrix}
  \item{outfile}{an optional text file to which to write output}
}

\details{This utility converts a species-by-plots matrix into \code{triplet} format, which is readable by the external program \code{phylocom}. 
The user has the option to save an output file, defined by \code{outfile}. 
}

\value{an array, formatted for use in \code{phylocom}}

\references{WEBB CO, DD ACKERLY and SW KEMBEL. 2008. Phylocom: software for the analysis of phylogenetic community structure and trait evolution. Bioinformatics 24:2098-2100.}

\author{Jonathan Eastman}

\seealso{\code{\link{as.spacodi}} for converting between \code{phylocom} and \code{SPACoDi} formats.}

\examples{
# generate a species-by-plots matrix
foo <- r.plot(species=10,plots=6,missing.prop=0.15,sim.tree=FALSE)	

# convert to phylocom format
as.phylocom(foo) -> p.foo
p.foo

# convert back to spacodi format
as.spacodi(p.foo) -> s.foo
s.foo

# save the converted dataset to working directory
as.phylocom(s.foo, outfile="phylocom.in.triplet.txt")


}