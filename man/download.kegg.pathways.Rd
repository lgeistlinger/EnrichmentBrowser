\name{download.kegg.pathways}
\alias{download.kegg.pathways}
\title{
    Download of KEGG pathways for a particular organism
}
\description{
    The function downloads all metabolic and non-metabolic pathways in
    KEGG XML format for a specified organism.
}
\usage{
    download.kegg.pathways( org, cache=TRUE, out.dir = NULL, zip = FALSE )
}
\arguments{
    \item{org}{
        Organism in KEGG three letter code, 
        e.g. \sQuote{hsa} for \sQuote{homo sapiens}.
    }
    \item{cache}{
        Logical.
        Should a locally cached version used if available?
        Defaults to \code{TRUE}.
    }
    \item{out.dir}{
        Output directory. 
        If not null, pathways are written to files in the specified directory.  
    }
    \item{zip}{
        Logical.
        In case pathways are written to file (\sQuote{out.dir} is not null): 
        should output files be zipped? 
    }
}
\value{
    if(is.null(out.dir)): a list of KEGGPathway objects
    else: none, as pathways are written to file
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@sph.cuny.edu>
}
\seealso{
    \code{\link{keggList}}, \code{\link{keggGet}}, 
    \code{\link[KEGGgraph:KEGGPathway-class]{KEGGPathway-class}}, 
    \code{\link[KEGGgraph:parseKGML]{parseKGML}}
}
\examples{
    \donttest{
        pwys <- download.kegg.pathways("hsa")
    }
}
