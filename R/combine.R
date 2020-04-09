###########################################################
#
# author: Ludwig Geistlinger
# date: 28 Feb 2011
#
# Combine result files of enrichment analyses
#
# Update:
#   27 May 2015: rank-based only combination
#
###########################################################

#' Combining enrichment analysis results
#' 
#' Different enrichment analysis methods usually result in different gene set
#' rankings for the same dataset.  This function allows to combine results from
#' the different set-based and network-based enrichment analysis methods.  This
#' includes the computation of average gene set ranks across methods.
#' 
#' 
#' @aliases comb.ea.results
#' @param res.list A list of enrichment analysis result lists (as returned by
#' the functions \code{\link{sbea}} and \code{\link{nbea}}).
#' @param rank.col Rank column.  Column name of the enrichment analysis result
#' table that should be used to rank the gene sets.  Defaults to the gene set
#' p-value column, i.e. gene sets are ranked according to gene set
#' significance.
#' @param decreasing Logical.  Should smaller (decreasing=FALSE, default) or
#' larger (decreasing=TRUE) values in rank.col be ranked better?  In case of
#' gene set p-values the smaller the better, in case of gene set scores the
#' larger the better.
#' @param rank.fun Ranking function.  Used to rank gene sets according to the
#' result table of individual enrichment methods (as returned from the
#' \code{\link{gsRanking}} function). This is typically done according to gene
#' set p-values, but can also take into account gene set scores/statistics,
#' especially in case of gene sets with equal p-value.  Can be either one of
#' the predefined functions ('comp.ranks', 'rel.ranks', 'abs.ranks') or a
#' user-defined function.  Defaults to 'comp.ranks', i.e. competitive
#' (percentile) ranks are computed by calculating for each gene set the
#' percentage of gene sets with a p-value as small or smaller.  Alternatively,
#' 'rel.ranks', i.e. relative ranks are computed in 2 steps: \enumerate{ \item
#' Ranks are assigned according to distinct gene set p-value *categories*, i.e.
#' gene sets with equal p-value obtain the *same* rank. Thus, the gene sets
#' with lowest p-value obtain rank 1, and so on.  \item As opposed to absolute
#' ranks (rank.fun = 'abs.ranks'), which are returned from step 1, relative
#' ranks are then computed by dividing the absolute rank by number of distinct
#' p-value categories and multiplying with 100 (= percentile rank).  }
#' @param comb.fun Rank combination function.  Used to combine gene set ranks
#' across methods.  Can be either one of the predefined functions (mean,
#' median, max, min, sum) or a user-defined function.  Defaults to 'sum', i.e.
#' the rank sum across methods is computed.
#' @return An enrichment analysis result list that can be detailedly explored
#' by calling \code{\link{eaBrowse}} and from which a flat gene set ranking
#' can be extracted by calling \code{\link{gsRanking}}.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{\link{sbea}}, \code{\link{nbea}}, \code{\link{eaBrowse}}
#' @examples
#' 
#' 
#'     # (1) expression data: 
#'     # simulated expression values of 100 genes
#'     # in two sample groups of 6 samples each
#'     se <- makeExampleData(what="SE")
#'     se <- deAna(se)
#' 
#'     # (2) gene sets:
#'     # draw 10 gene sets with 15-25 genes
#'     gs <- makeExampleData(what="gs", gnames=names(se))
#' 
#'     # (3) make artificial enrichment analysis results:
#'     # 2 ea methods with 5 significantly enriched gene sets each 
#'     ora.res <- makeExampleData(what="ea.res", method="ora", se=se, gs=gs) 
#'     gsea.res <- makeExampleData(what="ea.res", method="gsea", se=se, gs=gs)
#'     
#'     # (4) combining the results
#'     res.list <- list(ora.res, gsea.res)
#'     comb.res <- combResults(res.list)
#' 
#'     # (5) result visualization and exploration
#'     gsRanking(comb.res)
#' 
#'     # user-defined ranking and combination functions
#'     # (a) dummy ranking, give 1:nrow(res.tbl)
#'     dummy.rank <- function(res.tbl) seq_len(nrow(res.tbl))
#' 
#'     # (b) weighted average for combining ranks
#'     wavg <- function(r) mean(c(1,2) * r)
#' 
#'     comb.res <- combResults(res.list, rank.fun=dummy.rank, comb.fun=wavg)
#' 
#' @export combResults
combResults <- function(res.list, 
    rank.col=configEBrowser("PVAL.COL"),
    decreasing=FALSE,
    rank.fun=c("comp.ranks", "rel.ranks", "abs.ranks"), 
    comb.fun=c("mean", "median", "min", "max", "sum"))
{
    GS.COL <- configEBrowser("GS.COL")
    PVAL.COL <- configEBrowser("PVAL.COL")

    # requires min. 2 results
    stopifnot(length(res.list) > 1)   

    # compute the intersect of gene sets in all result tables
    nr.res <- length(res.list)  
    gs <- res.list[[1]]$res.tbl[,GS.COL]
    for(i in seq_len(nr.res-1)) 
        gs <- intersect(gs, res.list[[i+1]]$res.tbl[,GS.COL])
    for(i in seq_len(nr.res))
    {
        curr <- res.list[[i]]$res.tbl 
        res.list[[i]]$res.tbl <- curr[curr[,GS.COL] %in% gs,]
    }
    nr.gs <- length(gs)

    # determine average ranks
    rankm <- matrix(0, nrow=nr.gs, ncol=nr.res)
    pvalm <- matrix(0.0, nrow=nr.gs, ncol=nr.res)
    
    if(length(rank.col) != nr.res) rank.col <- rep(rank.col[1], nr.res)
    if(length(decreasing) != nr.res) decreasing <- rep(decreasing[1], nr.res)

    for(i in seq_len(nr.res))
    {
        res.tbl <- res.list[[i]]$res.tbl
        ranks <- .getRanks(res.tbl, 
            rank.fun=rank.fun, rank.col=rank.col[i], decreasing=decreasing[i])
        ord <- match(gs, res.tbl[,GS.COL])
        pvalm[,i] <- res.tbl[ord, PVAL.COL]
        rankm[,i] <- ranks[ord] 
    }

    # compute average ranks
    av.ranks <- .combRanks(rankm, comb.fun)
    if(is.character(comb.fun)) avr.cname <- toupper(match.arg(comb.fun))
    else avr.cname <- "COMB"
    avr.cname <- paste(avr.cname, "RANK", sep=".")

    # construct combined table
    res.tbl <- data.frame(rankm, av.ranks, pvalm)
    methods <- vapply(res.list, function(x) toupper(x$method), character(1))
    colnames(res.tbl) <- c(paste0(methods, ".RANK"), avr.cname, 
                paste0(methods, ".PVAL"))
    rownames(res.tbl) <- gs

    # order and format results
    res.tbl <- res.tbl[order(av.ranks),]
    res.tbl[,seq_len(nr.res+1)] <- round(res.tbl[,seq_len(nr.res+1)], digits=1)
    res.tbl[,(nr.res+2):ncol(res.tbl)] <- 
        signif(res.tbl[,(nr.res+2):ncol(res.tbl)], digits=3)
    res.tbl <- DataFrame(rownames(res.tbl), res.tbl)
    colnames(res.tbl)[1] <- GS.COL
    rownames(res.tbl) <- NULL

    res <- res.list[[1]]
    res$res.tbl <- res.tbl
    res$method <- "comb"
    res$nr.sigs <- min(nr.gs, configEBrowser("NR.SHOW"))

    return(res)
}

# Stouffer's combination of pvalues 
.stoufferComb <- function(pvals, two.sided=TRUE)
{
    ## remove NA
    pvals <- pvals[!is.na(pvals)]
    ## correct for two sided tests
    if(two.sided) pvals <- pvals/2
    ## z transformation
    z <- sum(-qnorm(pvals)) / sqrt(length(pvals))
    ## backtransformation
    p <- 1 - pnorm(abs(z))
    ## correct for two sided tests
    if(two.sided) p <- p * 2
    return(p)
}

# Fisher's combination of pvalues
.fisherComb <- function(pvals) {
    Xsq <- -2*sum(log(pvals))
    p <- 1 - pchisq(Xsq, df = 2*length(pvals))
    return(p)
}

# input: named ranking measure vector
# e.g. pvalues/scores of gene sets
# res <- c(0.01, 0.02, ..)
# names(res) <- c("gs1", "gs2", ..)
.getRanks <- function(res, 
    rank.fun=c("comp.ranks", "rel.ranks", "abs.ranks"), 
    rank.col=configEBrowser("PVAL.COL"),
    decreasing=FALSE)
{
    if(is.function(rank.fun)) ranks <- rank.fun(res)
    else{
        rank.fun <- match.arg(rank.fun)
        GS.COL <- configEBrowser("GS.COL")
        rcol <- res[, rank.col]
        if(decreasing) rcol <- -rcol
        names(rcol) <- res[, GS.COL] 

        if(rank.fun == "comp.ranks")
            ranks <- vapply(rcol, function(p) mean(rcol <= p) * 100, numeric(1))
        else
        {
            ucats <- unique(rcol)
            ranks <- match(rcol, ucats)
            if(rank.fun == "rel.ranks") ranks <- ranks / length(ucats) * 100
        }
        return(ranks)
    }
    return(ranks)
}

.combRanks <- function(rankm, 
    comb.fun=c("mean", "median", "min", "max", "sum"))
{   
    if(is.function(comb.fun)) ranks <- apply(rankm, 1, comb.fun)
    else
    {   
        comb.fun <- match.arg(comb.fun)
        ranks <- apply(rankm, 1, 
            function(x) do.call(comb.fun, list(x, na.rm=TRUE)))
    }
    return(ranks)
}


