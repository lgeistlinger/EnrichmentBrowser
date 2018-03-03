\name{ggeaGraph}
\alias{ggeaGraph}
\alias{ggeaGraphLegend}
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
    ggeaGraph( gs, grn, se, 
        alpha = 0.05, beta = 1, max.edges = 50, cons.thresh = 0.7 )

    ggeaGraphLegend()
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
    \item{se}{
        Expression data given as an object of class 
            \code{\linkS4class{SummarizedExperiment}}. 
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
    Ludwig Geistlinger <Ludwig.Geistlinger@sph.cuny.edu>
}
\seealso{
    \code{\link{nbea}} to perform network-based enrichment analysis.
    \code{\link{ea.browse}} for exploration of resulting gene sets. 
}
\examples{
    # (1) expression data: 
    # simulated expression values of 100 genes
    # in two sample groups of 6 samples each
    se <- makeExampleData(what="SE")
    se <- deAna(se)

    # (2) gene sets:
    # draw 10 gene sets with 15-25 genes
    gs <- makeExampleData(what="gs", gnames=names(se))

    # (3) compiling artificial regulatory network 
    grn <- makeExampleData(what="grn", nodes=names(se))

    # (4) plot consistency graph
    ggeaGraph(gs=gs[[1]], grn=grn, se=se)

    # (5) get legend
    ggeaGraphLegend()
}
