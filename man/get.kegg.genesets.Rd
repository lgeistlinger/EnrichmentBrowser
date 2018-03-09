\name{get.kegg.genesets}
\alias{get.kegg.genesets}
\title{
    Definition of gene sets according to KEGG pathways for a specified organism
}
\description{
    This function retrieves KEGG gene sets for an organism under investigation
    either via download from the KEGG REST server.
}
\usage{
    get.kegg.genesets( pwys, cache = TRUE, gmt.file = NULL )
}
\arguments{
    \item{pwys}{
        Either a list of \code{\linkS4class{KEGGPathway}} objects or an 
        absolute file path of a zip compressed archive of pathway xml files 
        in KGML format. 
        Alternatively, an organism in KEGG three letter code, 
        e.g. \sQuote{hsa} for \sQuote{Homo sapiens}.
    }
    \item{cache}{
        Logical.
        Should a locally cached version used if available?
        Defaults to \code{TRUE}.
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
    Ludwig Geistlinger <Ludwig.Geistlinger@sph.cuny.edu>
}
\seealso{  
    \code{\link{keggList}}, \code{\link{keggLink}},
    \code{\linkS4class{KEGGPathway}}, 
    \code{\link{parseKGML}}
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
}
