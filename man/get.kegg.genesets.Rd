\name{get.kegg.genesets}
\alias{get.kegg.genesets}
\alias{parse.genesets.from.GMT}
\title{
    Definition of gene sets according to KEGG pathways for a specified organism
}
\description{
    To perform a gene set enrichment analysis on KEGG pathways,
    it is necessary to build up the gene set database in a format
    that the GSEA method can read. 
    Parsing a list of gene sets from a flat text file in GMT format.
    This function performs the necessary steps, 
    including the retrieval of the participating gene IDs
    for each pathway and the conversion to GMT format.
}
\usage{
    get.kegg.genesets( pwys, gmt.file = NULL )

    parse.genesets.from.GMT( gmt.file )
}
\arguments{
    \item{pwys}{
        Either a list of \code{\linkS4class{KEGGPathway}} objects or an 
        absolute file path of a zip compressed archive of pathway xml files 
        in KGML format. 
        Alternatively, an organism in KEGG three letter code, 
        e.g. \sQuote{hsa} for \sQuote{Homo sapiens}.
    }
    \item{gmt.file}{
        Gene set file in GMT format. See details.
    }
}
\value{
    A list of gene sets (vectors of gene IDs).
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

    KEGG Organism code
    \url{http://www.genome.jp/kegg/catalog/org_list.html}
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@bio.ifi.lmu.de>
}
\seealso{  
    \code{\link{keggList}}, \code{\link{keggLink}},
    \code{\link[KEGGgraph:KEGGPathway-class]{KEGGPathway-class}}, 
    \code{\link[KEGGgraph:parseKGML]{parseKGML}}
}
\examples{
    # WAYS TO DEFINE GENE SETS ACCORDING TO HUMAN KEGG PATHWAYS

    # (1) from scratch: via organism ID 
    \donttest{
    gs <- get.kegg.genesets("hsa")
    }

    # (2) extract from pathways
    # download human pathways via: 
    # pwys <- download.kegg.pathways("hsa")
    pwys <- system.file("extdata/hsa_kegg_pwys.zip", package="EnrichmentBrowser")
    gs <- get.kegg.genesets(pwys)

    # (3) parsing gene sets from GMT
    gmt.file <- system.file("extdata/hsa_kegg_gs.gmt", package="EnrichmentBrowser")
    gs <- parse.genesets.from.GMT(gmt.file)     
}
