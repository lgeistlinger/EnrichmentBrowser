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

# Stouffer's combination of pvalues 
stouffer.comb <- function(pvals, two.sided=TRUE)
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
fisher.comb <- function(pvals) {
    Xsq <- -2*sum(log(pvals))
    p <- 1 - pchisq(Xsq, df = 2*length(pvals))
    return(p)
}

###
# combine results of different enrichment analyises
# and compute average measures (rank & p-value)
###
# TODO: (1) weighted average ranks (to prioritize methods)
#       (2) what to do with equal average ranks (stouffer?)
comb.ea.results <- function(res.list)
{
    GS.COL <- config.ebrowser("GS.COL")
    GSP.COL <- config.ebrowser("GSP.COL")

    # requires min. 2 results
    if(length(res.list) < 2) stop("Length of \'res.list\' should be > 1")   

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
    
    for(i in seq_len(nr.res))
    {
        res.tbl <- res.list[[i]]$res.tbl
        ranking <- sapply(gs, function(s) which(res.tbl[,GS.COL] == s))
        rankm[,i] <- ranking 
        pvalm[,i] <- as.numeric(res.tbl[rankm[,i], GSP.COL])
    }

    # compute average ranks
    av.ranks <- rowMeans(rankm, na.rm=TRUE)

    # construct combined table
    res.tbl <- data.frame(rankm, av.ranks, pvalm)
    methods <- sapply(res.list, function(x) toupper(x$method))
    colnames(res.tbl) <- c(paste0(methods, ".RANK"), "AVG.RANK", 
                paste0(methods, ".PVAL"))
    rownames(res.tbl) <- gs

    # order and format results
    res.tbl <- res.tbl[order(av.ranks),]
    res.tbl[,seq_len(nr.res+1)] <- floor(res.tbl[,seq_len(nr.res+1)])
    res.tbl[,(nr.res+2):ncol(res.tbl)] <- 
        signif(res.tbl[,(nr.res+2):ncol(res.tbl)], digits=3)
    res.tbl <- DataFrame(rownames(res.tbl), res.tbl)
    colnames(res.tbl)[1] <- GS.COL
    rownames(res.tbl) <- NULL

    res <- res.list[[1]]
    res$res.tbl <- res.tbl
    res$method <- "comb"
    res$nr.sigs <- min(nr.gs, config.ebrowser("NR.SHOW"))

    return(res)
}
