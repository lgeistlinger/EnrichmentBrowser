\name{isAvailable}
\alias{isAvailable}
\title{
   Is a required package available? 
}
\description{
   Convenience function for checking and installing required packages. 
}
\usage{
    isAvailable( pkg, type=c("annotation", "software", "data") )
}
\arguments{
    \item{pkg}{
        Character vector of length 1.
        A valid name of an existing R package.
    }
    \item{type}{
        Character vector of length 1.
        What type of package is this? 
        Choose one out of 'annotation', 'software', or 'data' package.
    }
}
\value{
    None. See details.
}
\details{
    Checks whether a required package is available in the library.
    If yes, the package is loaded via \code{\link{require}}.
    If not, the package is optionally installed via 
    \code{\link[BiocInstaller:biocLite]{biocLite}}
    and, subsequently, loaded via \code{\link{require}}.
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@sph.cuny.edu>
}
\seealso{  
    \code{require}, \code{biocLite}
}
\examples{
    \donttest{
        isAvailable("EnrichmentBrowser", type="software")
    }     
}
