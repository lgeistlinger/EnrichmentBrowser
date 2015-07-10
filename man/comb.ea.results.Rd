\name{comb.ea.results}
\alias{comb.ea.results}
\title{
    Combining enrichment analysis results
}
\description{
    Different enrichment analysis methods usually result in 
    different gene set rankings for the same dataset.
    This function allows to combine results from the different 
    set-based and network-based enrichment analysis methods.
    This includes the computation of average gene set ranks 
    across methods and of combined gene set p-values. 
}
\usage{
    comb.ea.results(res.list)
}
\arguments{
    \item{res.list}{
        A list of enrichment analysis result lists as returned 
        by the functions \sQuote{sbea} and \sQuote{nbea}.
    } 
}
\value{
    An enrichment analysis result list 
    that can be detailedly explored by calling 'ea.browse' and from 
    which a flat gene set ranking can be extracted by calling 'gs.ranking'.
}
\seealso{
    \code{\link{sbea}}, \code{\link{nbea}}, \code{\link{ea.browse}}}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@bio.ifi.lmu.de>
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

    # (3) make artificial enrichment analysis results:
    # 2 ea methods with 5 significantly enriched gene sets each 
    ora.res <- make.example.data(what="ea.res", method="ora", eset=eset, gs=gs) 
    gsea.res <- make.example.data(what="ea.res", method="gsea", eset=eset, gs=gs)
    
    # (4) combining the results
    res.list <- list(ora.res, gsea.res)
    comb.res <- comb.ea.results(res.list)

    # (5) result visualization and exploration
    gs.ranking(comb.res)
}

