\name{ea.browse}
\alias{ea.browse}
\alias{gs.ranking}
\title{
    Exploration of enrichment analysis results
}
\description{
    Functions to extract a flat gene set ranking from an 
    enrichment analysis result object and to detailedly explored it. 
}
\usage{
    ea.browse(res, nr.show=-1, set.view=TRUE, graph.view=NULL)

    gs.ranking(res, signif.only=TRUE)
}
\arguments{
    \item{res}{
        Enrichment analysis result object as returned by the functions 
        \sQuote{sbea} and \sQuote{nbea}.
    } 
    \item{nr.show}{
        Number of gene sets to show. 
        As default all statistical significant gene sets are displayed.
    }
    \item{set.view}{
        Logical.
        Should a set-based summary (includes a report, plots, and browsing in 
        KEGG) be created for the result? Defaults to TRUE. 
    }
    \item{graph.view}{
        Optional. 
        Should a graph-based summary (reports and visualizes consistency of 
        regulations) be created for the result?
        If specified, it needs to be a gene regulatory network, i.e. either 
        an absolute file path to a tabular file or a character matrix with 
        exactly *THREE* cols; 1st col = IDs of regulating genes; 
        2nd col = corresponding regulated genes;
        3rd col = regulation effect; Use '+' and '-' for activation/inhibition.
    }
    \item{signif.only}{
        Logical.
        Display only those gene sets in the ranking which satisfy the 
        significance level? Defaults to TRUE.
    }
}
\value{
    gs.ranking: data.frame with gene sets ranked by the corresponding p-value;

    ea.browse: none, opens the browser to explore results.
} 
\seealso{
    \code{\link{sbea}}, \code{\link{nbea}}, \code{\link{comb.ea.results}}}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@bio.ifi.lmu.de>
}
\examples{
    # (1) reading the expression data from file
    exprs.file <- system.file("extdata/ALL_exprs.tab", package="EnrichmentBrowser")
    pdat.file <- system.file("extdata/ALL_pData.tab", package="EnrichmentBrowser")
    fdat.file <- system.file("extdata/ALL_fData.tab", package="EnrichmentBrowser")
    probe.eset <- read.eset(exprs.file, pdat.file, fdat.file)

    # (2) summarizing probe expression on gene level
    gene.eset <- probe.2.gene.eset(probe.eset) 

    # (3) getting all human KEGG gene sets
    # hsa.gs <- get.kegg.genesets("hsa")
    gs.file <- system.file("extdata/hsa_kegg_gs.gmt", package="EnrichmentBrowser")
    hsa.gs <- parse.genesets.from.GMT(gs.file)

    # (4) performing the enrichment analysis
    ea.res <- sbea(method="ora", eset=gene.eset, gs=hsa.gs, perm=0)

    # (5) result visualization and exploration
    gs.ranking(ea.res)
    ea.browse(ea.res)
}

