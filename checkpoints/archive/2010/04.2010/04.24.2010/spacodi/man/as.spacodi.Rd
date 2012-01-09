\name{as.spacodi}

\alias{as.spacodi}

\title{converting between data formats for community phylogenetics}

\usage{as.spacodi(data, outfile = NULL)}

\arguments{
  \item{data}{a species-by-plots matrix}
  \item{outfile}{an optional text file to which to write output}
}

\details{This utility converts a species-by-plots matrix into a format 
readable by the external program \code{SPACoDi}, a Windows executable 
(available at \url{http://ebe.ulb.ac.be/ebe/Software.html} ). 
The user has the option to save an output file, defined by \code{outfile}. 
}

\value{an array, formatted for use in \code{SPACoDi} for Windows}

\author{Jonathan Eastman}

\seealso{\code{\link{as.phylocom}} for converting between \code{phylocom} and \code{SPACoDi} formats;
see \code{\link{spacodi.calc}} for community diversity measures of empirical data}

\examples{

# call example data
data(spacodi.example)
attach(spacodi.example)

# save the converted dataset to working directory
as.phylocom(sp_plot, outfile="phylocom.in.triplet.txt") -> p.sp_plot

# convert back to spacodi format and write to file
as.spacodi(p.sp_plot, outfile="spacodi.formatted.txt") -> s.sp_plot
s.sp_plot		# re-converted SPACoDi dataset
sp_plot			# comparison with the original dataset

}