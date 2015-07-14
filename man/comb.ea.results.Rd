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
    comb.ea.results(res.list, 
        rank.fun=c("rel.ranks", "abs.ranks"), 
        comb.fun=c("mean", "median", "min", "max", "sum"))
}
\arguments{
    \item{res.list}{
        A list of enrichment analysis result lists as returned 
        by the functions \sQuote{sbea} and \sQuote{nbea}.
    }
    \item{rank.fun}{
        Ranking function.
        Used to rank gene sets according to the result table of individual
        enrichment methods (as returned from the \code{\link{gs.ranking}} 
        function). This is typically done according to gene set p-values,
        but can also take into account gene set scores/statistics, especially
        in case of gene sets with equal p-value.
        Can be either one of the predefined functions (rel.ranks, abs.ranks)
        or a user-defined function.
        See details.
        Default to 'rel.ranks', i.e. relative ranks are computed in 2 steps:
        \enumerate{
            \item Ranks are assigned according to distinct gene set p-value
                *categories*, i.e. gene sets with equal p-value obtain the *same*
                rank. Thus the gene sets with lowest p-value obtain rank 1, 
                and so on.
            \item As opposed to absolute ranks (rank.fun='abs.ranks'), which
                are returned from step 1, relative ranks are then computed by 
                dividing the absolute rank by number of distinct p-value 
                categories and multiplying with 100 (= percentile rank).
        }
    }
    \item{comb.fun}{
        Rank combination function.
        Used to combine gene set ranks across methods.
        Can be either one of the predefined functions (mean, median,
        max, min, sum) or a user-defined function.
        See details.
        Defaults to 'mean', i.e. the average rank across methods is
        computed.  
    } 
}
\value{
    An enrichment analysis result list 
    that can be detailedly explored by calling 'ea.browse' and from 
    which a flat gene set ranking can be extracted by calling 'gs.ranking'.
}
\seealso{
    \code{\link{sbea}}, \code{\link{nbea}}, \code{\link{ea.browse}}
}
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

    # user-defined ranking and combination functions
    # (a) dummy ranking, give 1:nrow(res.tbl)
    dummy.rank <- function(res.tbl) seq_len(nrow(res.tbl))

    # (b) weighted average for combining ranks
    wavg <- function(r) mean(c(1,2) * r)

    comb.res <- comb.ea.results(res.list, rank.fun=dummy.rank, comb.fun=wavg)
}

