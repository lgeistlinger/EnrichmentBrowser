\name{ggea.graph}
\alias{ggea.graph}
\alias{ggea.graph.legend}
\title{
    GGEA graphs of consistency between regulation and expression
}
\description{
    Gene graph enrichment analysis (GGEA) is a network-based enrichment 
    analysis method implemented in the EnrichmentBrowser package.
    The idea of GGEA is to evaluate the consistency of known regulatory 
    interactions with the observed gene expression data.
    A GGEA graph for a gene set of interest displays the consistency 
    of each interaction in the network that involves a gene set member.
    Nodes (genes) are colored according to expression (up-/down-regulated) and
    edges (interactions) are colored according to consistency, i.e. how well 
    the interaction type (activation/inhibition) is reflected in the correlation
    of the expression of both interaction partners.     
}
\usage{
    ggea.graph(gs, grn, eset, org=NA, alpha=0.05, beta=1, max.edges=50, cons.thresh=0.7)
    ggea.graph.legend()
}
\arguments{
    \item{gs}{
        Gene set under investigation.
        This should be a character vector of KEGG gene IDs.
    }
    \item{grn}{
        Gene regulatory network.
        Character matrix with exactly *THREE* cols; 
        1st col = IDs of regulating genes; 
        2nd col = corresponding regulated genes;
        3rd col = regulation effect; Use '+' and '-' for activation/inhibition.
    }
    \item{eset}{
        Expression set. 
        An object of class 'ExpressionSet' containing the gene expression set. 
        See 'read.eset' and 'probe.2.gene.eset' for required annotations 
        in the pData and fData slot.
    }
    \item{org}{
        Organism under investigation in KEGG three letter code, 
        e.g. \sQuote{hsa} for \sQuote{homo sapiens}.
        Used to map gene IDs to gene symbols as displayed in 
        KEGG pathway maps.
        If not annotated, gene IDs are plotted as default.
    }
    \item{alpha}{
        Statistical significance level. Defaults to 0.05.
    }
    \item{beta}{
        Log2 fold change significance level. Defaults to 1 (2-fold).
    }  
    \item{max.edges}{
        Maximum number of edges that should be displayed. 
        Defaults to 50.
    }  
    \item{cons.thresh}{
        Consistency threshold.
        Graphical parameter that correspondingly increases line width 
        of edges with a consistency above the chosen threshold (defaults to 0.7).
    }
}
\value{
    None, plots to a graphics device.
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@bio.ifi.lmu.de>
}
\seealso{
    \code{\link{nbea}} to perform network-based enrichment analysis.
    \code{\link{ea.browse}} for exploration of resulting gene sets. 
}
\examples{
    # (1) reading the expression data from file
    exprs.file <- system.file("extdata/ALL_exprs.tab", package="EnrichmentBrowser")
    pdat.file <- system.file("extdata/ALL_pData.tab", package="EnrichmentBrowser")
    fdat.file <- system.file("extdata/ALL_fData.tab", package="EnrichmentBrowser")
    probe.eset <- read.eset(exprs.file, pdat.file, fdat.file)

    # (2) summarizing probe expression on gene level
    gene.eset <- probe.2.gene.eset(probe.eset) 

    # (3a) getting all human KEGG gene sets
    # hsa.gs <- get.kegg.genesets("hsa")
    gs.file <- system.file("extdata/hsa_kegg_gs.gmt", package="EnrichmentBrowser")
    hsa.gs <- parse.genesets.from.GMT(gs.file)

    # (3b) compiling gene regulatory network from KEGG pathways
    # hsa.grn <- compile.grn.from.kegg("hsa")
    pwys <- system.file("extdata/hsa_kegg_pwys.zip", package="EnrichmentBrowser")
    hsa.grn <- compile.grn.from.kegg(pwys)

    # plot consistency graph
    ggea.graph(gs=hsa.gs[["hsa05416_Viral_myocarditis"]], 
                                                grn=hsa.grn, eset=gene.eset)

    # get legend
    ggea.graph.legend()
}
