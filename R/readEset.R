# TODO: Add comment
# 
# Author: ludwig geistlinger
# Date: May 27th, 2010
#
# reads in data from plain text files and stores it in an expression
# set and performs simultaneously a de primary analysis
###############################################################################

read.eset <- function(
    exprs.file, 
    pdat.file, 
    fdat.file, 
    value.type=c("log2count", "log2ratio"),
    out.file=NULL, 
    heatm.file=NULL, 
    distr.file=NULL, 
    volc.file=NULL)
{
    value.type <- value.type[1]
    if(!(value.type %in% c("log2count", "log2ratio")))
        stop("\'value.type\' must be either \'log2count\' or \'log2ratio\'")

    # read features
    fDat <- scan(fdat.file, what="character", quiet=TRUE)
    nr.features <- length(fDat) / 2
    fDat <- matrix(fDat, nrow=nr.features, ncol=2, byrow=TRUE)

    # read samples
    pDat <- scan(pdat.file, what="character", quiet=TRUE)
    nr.samples <- length(pDat) / 2
    pDat <- matrix(pDat, nrow=nr.samples, ncol=2, byrow=TRUE)
    pDat[,2] <- gsub("/", "-", pDat[,2]) 

    
    # read expression values
    expr <- matrix(scan(exprs.file, quiet=TRUE), 
        nrow=nr.features, ncol=nr.samples, byrow=TRUE)
    rownames(expr) <- fDat[,1]
    colnames(expr) <- pDat[,1]

    # exclude NAs
    na.indr <- which(apply(expr, 1, function(x) any(is.na(x))))
    for(i in seq_along(na.indr))
    {
        cexpr <- expr[na.indr[i],]
        na.indc <- is.na(cexpr)
        cexpr[na.indc] <- mean(cexpr[!na.indc])
        expr[na.indr[i],] <- cexpr
    }
    
    # create the eset
    eset <- new("ExpressionSet", exprs=expr)
    pData(eset)[, SMPL.COL] <- pDat[,1]

    # is phenotype binary
    pData(eset)[, GRP.COL] <- pDat[,2]

    fData(eset)[, PRB.COL] <- fDat[,1]
    fData(eset)[, GN.COL] <- fDat[,2]

    eset <- de.ana(eset, value.type=value.type)

    # output the eset
    if(!is.null(out.file)) save(eset, file=out.file)
    
    # plot de
    # (a) heatmap
    if(!is.null(heatm.file)) exprs.heatmap(eset, heatm.file)
    
    # (b) pval distribution
    if(!is.null(distr.file)) pdistr(eset, distr.file)   
    
    # (c) volcano plot
    if(!is.null(volc.file)) volcano(eset, volc.file)
    
    return(eset)
}

# perform de analysis
de.ana <- function(eset, value.type=c("log2count", "log2ratio"))
{
    value.type <- value.type[1]
    groups <- sort(unique(pData(eset)[,GRP.COL]))
    if(!all(groups == GRPS)) 
        stop(paste("Group classification in pData is not binary:\n",
            "Expected (0, 1) but found (", 
                paste(groups, collapse=", "), ")", sep=""))

    if(value.type == "log2count")
    {
        fc.tt <- get.fold.change.and.t.test(eset, GRP.COL, groups)
        fData(eset)[, FC.COL] <- fc.tt@fc
        fData(eset)[, RAWP.COL] <- fc.tt@tt
    }
    # log2ratio
    else{
        group <- pData(eset)[, GRP.COL] == "1"
        fData(eset)[, FC.COL] <- apply(exprs(eset), 1, 
            function(x) mean(x[group]) - mean(x[!group]))
        fData(eset)[, RAWP.COL] <- rowttests(
            exprs(eset), factor(pData(eset)[, GRP.COL]))$p.value
    }

    fData(eset)[, ADJP.COL] <- p.adjust(
        fData(eset)[, RAWP.COL], method=ADJ.METH)
    return(eset)
}

# P-Value Distribution
pdistr <- function(eset, out.file=NULL, use.adjp=TRUE)
{
    if(use.adjp) p <- fData(eset)[, ADJP.COL]
    else p <- fData(eset)[, RAWP.COL]
    if(!is.null(out.file)) png(filename=out.file) 
    truehist(p, nbins=100, prob=TRUE,
            main="P-Value Distribution", xlab="P-Value", ylab="Frequency")
    if(!is.null(out.file)) dev.off()
}

# Volcano Plot (fold change vs. p-value)
volcano <- function(eset, out.file=NULL)
{
    if(!is.null(out.file)) png(filename=out.file)
    plot(x=fData(eset)[, FC.COL], 
        y=-log(fData(eset)[,ADJP.COL], base=10), col="red", 
            main="Volcano Plot", xlab="log2(foldChange)", ylab="-log10(p)")
    if(!is.null(out.file)) dev.off()
}

# Heatmap
exprs.heatmap <- function(eset, out.file=NULL)
{
    HMcols <- rev(brewer.pal(10,"RdBu"))
    cols <- brewer.pal(10, "BrBG")

    expr <- exprs(eset)
    grp <- pData(eset)[, GRP.COL]
    ord <- order(grp)
    expr <- expr[,ord]
    grp <- grp[ord] 
    grpColors <- ifelse(grp==1, "brown", "grey")  
    if(!is.null(out.file)) png(filename=out.file)
    heatmap(exprs(eset), scale="row", col = HMcols, Colv=NA,
        ColSideColors=grpColors, keep.dendro=TRUE,
        #distfun=onecor,
        labRow=FALSE, xlab="Samples", ylab="Features")
    if(!is.null(out.file)) dev.off()
}

onecor <- function(x) as.dist(1-cor(t(x)))


