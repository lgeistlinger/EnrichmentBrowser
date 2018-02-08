\name{parseGMT}
\alias{parse.genesets.from.GMT}
\alias{writeGMT}
\title{
    Reading from and writing to GMT file format for gene sets  
}
\description{
    Functionaliy for parsing and writing a list of gene sets 
    from/to a flat text file in GMT format.
}
\usage{
    writeGMT( gs, gmt.file )

    parse.genesets.from.GMT( gmt.file )
}
\arguments{
 \item{gs}{
        A list of gene sets (character vectors of gene IDs).
    }   
    \item{gmt.file}{
        Gene set file in GMT format. See details.
    }
}
\value{
    For \code{parse.genesets.from.GMT}: a list of gene sets (vectors of gene IDs).
}
\details{
    The GMT (Gene Matrix Transposed) file format
    is a tab delimited file format that describes gene sets. 
    In the GMT format, each row represents a gene set.
    Each gene set is described by a name, a description, 
    and the genes in the gene set. See references.
}
\references{
    GMT file format
    \url{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@sph.cuny.edu>
}
\examples{
    # parsing gene sets from GMT
    gmt.file <- system.file("extdata/hsa_kegg_gs.gmt", package="EnrichmentBrowser")
    gs <- parse.genesets.from.GMT(gmt.file)     
    
    # writing gene sets to file
    writeGMT(gs, gmt.file)
}
