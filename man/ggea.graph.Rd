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
    ggea.graph( gs, grn, eset, 
        alpha = 0.05, beta = 1, max.edges = 50, cons.thresh = 0.7 )

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
        An object of class \code{\linkS4class{ExpressionSet}} containing the 
        gene expression set. See \code{\link{read.eset}} and 
        \code{\link{probe.2.gene.eset}} for required annotations 
        in the pData and fData slot.
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
    # (1) expression set: 
    # simulated expression values of 100 genes
    # in two sample groups of 6 samples each
    eset <- make.example.data(what="eset")
    eset <- de.ana(eset)

    # (2) gene sets:
    # draw 10 gene sets with 15-25 genes
    gs <- make.example.data(what="gs", gnames=featureNames(eset))

    # (3) compiling artificial regulatory network 
    grn <- make.example.data(what="grn", nodes=featureNames(eset))

    # (4) plot consistency graph
    ggea.graph(gs=gs[[1]], grn=grn, eset=eset)

    # (5) get legend
    ggea.graph.legend()
}
