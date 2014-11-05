###########################################################
#
# author: Ludwig Geistlinger
# date: 28 Feb 2011
#
# Combine result files of enrichment analyses
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
comb.ea.results <- function(    res.list, 
                                pcomb.meth=c("fisher","stouffer"), 
                                out.file=NULL)
{
    # requires min. 2 results
    pcomb.meth <- pcomb.meth[1]
    if(length(res.list) < 2) stop("Length of \'res.list\' should be > 1")   
    if(!(pcomb.meth %in% c("stouffer", "fisher"))) 
        stop("\'pcomb.meth\' must be either \'stouffer\' or \'fisher\'")

    # compute the intersect of gene sets in all result tables
    nr.res <- length(res.list)  
    gs <- res.list[[1]]$res.tbl[,"GENE.SET"]
    for(i in seq_len(nr.res-1)) 
        gs <- intersect(gs, res.list[[i+1]]$res.tbl[,"GENE.SET"])
    for(i in seq_len(nr.res))
    {
        curr <- res.list[[i]]$res.tbl 
        res.list[[i]]$res.tbl <- curr[curr[,"GENE.SET"] %in% gs,]
    }
    nr.gs <- length(gs)

    # determine average ranks & merged pvals
    rankm <- matrix(0, nrow=nr.gs, ncol=nr.res)
    pvalm <- matrix(0.0, nrow=nr.gs, ncol=nr.res)

    for(i in seq_len(nr.res))
    {
        res.tbl <- res.list[[i]]$res.tbl
        # compute average ranks
        ranking <- sapply(gs, function(s) which(res.tbl[,"GENE.SET"] == s))
        rankm[,i] <- ranking 

        # compute merged pvals
        pvalm[,i] <- as.numeric(res.tbl[rankm[,i], "P.VALUE"])
    }

    av.ranks <- rowMeans(rankm, na.rm=TRUE)
    if(pcomb.meth[1] == "fisher") merged.ps <- apply(pvalm, 1, fisher.comb)
    else merged.ps <- apply(pvalm, 1, stouffer.comb)

    # construct combined table
    res.tbl <- cbind(rankm, av.ranks, pvalm, merged.ps)
    methods <- sapply(res.list, function(x) toupper(x$method))
    colnames(res.tbl) <- c(paste(methods, ".RANK", sep=""), "AVG.RANK", 
                paste(methods, ".PVAL", sep=""), "COMB.PVAL")
    rownames(res.tbl) <- gs

    # order and format res.tblults
    res.tbl <- res.tbl[order(merged.ps),]
    res.tbl[,seq_len(nr.res+1)] <- floor(res.tbl[,seq_len(nr.res+1)])
    res.tbl[,(nr.res+2):ncol(res.tbl)] <- 
        signif(res.tbl[,(nr.res+2):ncol(res.tbl)], digits=3)
    res.tbl <- cbind(rownames(res.tbl), res.tbl)
    colnames(res.tbl)[1] <- "GENE.SET"
    rownames(res.tbl) <- NULL

    res <- res.list[[1]]
    res$res.tbl <- res.tbl
    res$method <- "comb"
    res$nr.sigs <- min(nr.gs, NROW.TOP.TABLE)

    # write an output table if desired
    if(!is.null(out.file)) 
        write.table(res.tbl, file=out.file, 
            quote=FALSE, row.names=FALSE, sep="\t")
    return(res)
}
