\name{ea.browse}
\alias{ea.browse}
\alias{gs.ranking}
\title{
    Exploration of enrichment analysis results
}
\description{
    Functions to extract a flat gene set ranking from an 
    enrichment analysis result object and to detailedly explore it. 
}
\usage{
    ea.browse( res, nr.show = -1, graph.view = NULL, html.only = FALSE )

    gs.ranking( res, signif.only = TRUE )
}
\arguments{
    \item{res}{
        Enrichment analysis result list (as returned by the functions 
        \code{\link{sbea}} and \code{\link{nbea}}).
    } 
    \item{nr.show}{
        Number of gene sets to show. 
        As default all statistically significant gene sets are displayed.
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
    \item{html.only}{
        Logical.
        Should the html file only be written (without opening the browser to 
        view the result page)? Defaults to FALSE.
    }
    \item{signif.only}{
        Logical.
        Display only those gene sets in the ranking, which satisfy the 
        significance level? Defaults to TRUE.
    }
}
\value{
    gs.ranking: DataFrame with gene sets ranked by the corresponding p-value;

    ea.browse: none, opens the browser to explore results.
} 
\seealso{
    \code{\link{sbea}}, \code{\link{nbea}}, \code{\link{comb.ea.results}}}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@sph.cuny.edu>
}
\examples{
    
    # real data
    # (1) reading the expression data from file
    exprs.file <- system.file("extdata/exprs.tab", package="EnrichmentBrowser")
    cdat.file <- system.file("extdata/colData.tab", package="EnrichmentBrowser")
    rdat.file <- system.file("extdata/rowData.tab", package="EnrichmentBrowser")
    probeSE <- readSE(exprs.file, cdat.file, rdat.file)
    geneSE <- probe2gene(probeSE) 
    geneSE <- de.ana(geneSE)
    metadata(geneSE)$annotation <- "hsa"

    # artificial enrichment analysis results
    gs <- make.example.data(what="gs", gnames=names(geneSE))
    ea.res <- make.example.data(what="ea.res", method="ora", eset=geneSE, gs=gs)

    # (5) result visualization and exploration
    gs.ranking(ea.res)
    ea.browse(ea.res)
}

