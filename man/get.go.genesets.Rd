\name{get.go.genesets}
\alias{get.go.genesets}
\title{
    Definition of gene sets according to the Gene Ontology (GO)
}
\description{
    This function retrieves GO gene sets for an organism under investigation
    either via download from BioMart or based on BioC annotation packages.
}
\usage{
    get.go.genesets( org, onto = c("BP", "MF", "CC"), mode = c("GO.db","biomart") )
}
\arguments{
    \item{org}{
        An organism in (KEGG) three letter code, 
        e.g. \sQuote{hsa} for \sQuote{Homo sapiens}.
    }
    \item{onto}{
        Character. Specifies one of the three GO ontologies:
        'BP' (biological process), 'MF' (molecular function), 
        'CC' (cellular component). Defaults to 'BP'.
    }
    \item{mode}{
        Character. Determines in which way the gene sets are retrieved.
        This can be either 'GO.db' or 'biomart'.
        The 'GO.db' mode creates the gene sets based on BioC annotation
        packages - which is fast, but represents not necessarily the most
        up-to-date mapping. In addition, this option is only available for
        the currently supported model organisms in BioC.
        The 'biomart' mode downloads the mapping from BioMart - which can be 
        time consuming, but allows to select from a larger range of organisms
        and contains the latest mappings.
        Defaults to 'GO.db'.
    }
}
\value{
    A list of gene sets (vectors of gene IDs).
}
\references{
    \url{http://geneontology.org/}
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@bio.ifi.lmu.de>
}
\seealso{  
    \code{\link{annFUN}} for general GO2gene mapping used in the 'GO.db' mode,
    \code{\link{getBM}} for general queries to BioMart.
    \code{\link{get.kegg.genesets}} for defining gene sets according to KEGG, 
    \code{\link{parse.genesets.from.GMT}} to parse user-def. gene sets from file.
}
\examples{
    # Typical usage for gene set enrichment analysis:
    # Biological process terms based on BioC annotation (for human)
    gs <- get.go.genesets("hsa")
    
    # eq.:  
    # gs <- get.go.genesets(org="hsa", onto="BP", mode="GO.db")

    # Alternatively:
    # downloading from BioMart 
    # this may take a few minutes ...
    \donttest{
    gs <- get.go.genesets(org="hsa", mode="biomart")
    }
}
