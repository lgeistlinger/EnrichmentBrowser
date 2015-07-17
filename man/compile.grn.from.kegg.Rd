\name{compile.grn.from.kegg}
\alias{compile.grn.from.kegg}
\title{
    Compilation of a gene regulatory network from KEGG pathways
}
\description{
    To perform network-based enrichment analysis a gene regulatory network (GRN)
    is required. There are well-studied processes and organisms for which 
    comprehensive and well-annotated regulatory networks are available, e.g. the
    RegulonDB for E. coli and Yeastract for S. cerevisiae. 
    However, in many cases such a network is missing. 
    A first simple workaround is to compile a network from regulations in the 
    KEGG database.
}
\usage{
    compile.grn.from.kegg( pwys, out.file = NULL )
}
\arguments{
    \item{pwys}{
        Either a list of \code{\linkS4class{KEGGPathway}} objects or an absolute
        file path of a zip compressed archive of pathway xml files in KGML format. 
        Alternatively, you can specify an organism in KEGG three letter code, 
        e.g. \sQuote{hsa} for \sQuote{Homo sapiens}, and the pathways will be
        downloaded automatically.
    }
    \item{out.file}{
        Optional output file the gene regulatory network will be written to.
    }
}
\value{
    if(is.null(out.file)): the gene regulatory network;
    else: none, as the gene regulatory network is written to file
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@bio.ifi.lmu.de>
}
\seealso{
    \code{\link[KEGGgraph:KEGGPathway-class]{KEGGPathway-class}}, 
    \code{\link[KEGGgraph:parseKGML]{parseKGML}}, 
    \code{\link{download.kegg.pathways}}
}
\examples{
    # (1) download human pathways
    #   pwys <- download.kegg.pathways("hsa")
    # (2) compile gene regulatory network
    #   grn <- compile.grn.from.kegg(pwys)
    
    pwys <- system.file("extdata/hsa_kegg_pwys.zip", package="EnrichmentBrowser")
    hsa.grn <- compile.grn.from.kegg(pwys)
}
