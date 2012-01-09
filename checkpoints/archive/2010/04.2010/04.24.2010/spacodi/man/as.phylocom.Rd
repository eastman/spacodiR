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

\seealso{\code{\link{as.spacodi}} for converting between \code{phylocom} 
and \code{SPACoDi} formats; see \code{\link[picante]{picante-package}} 
for an R-port of \code{phylocom}.}

\examples{
# call example data
data(spacodi.example)
attach(spacodi.example)

# convert to phylocom format
as.phylocom(sp_plot) -> p.sp_plot
p.sp_plot

# convert back to spacodi format
as.spacodi(p.sp_plot) -> s.sp_plot
s.sp_plot

# save the converted dataset to working directory
as.phylocom(s.sp_plot, outfile="phylocom.in.triplet.txt")


}