# 
# 
# Author: ludwig geistlinger
# Date: May 27th, 2010
#
# reads in data from plain text files and stores it in an expression
# set and performs simultaneously a de primary analysis
#
# UPDATE: Oct 21th, 2014
#
# de analysis is now based on limma rather than simpleaffy
# supported data types now include RNA-seq
#
###############################################################################

read.eset <- function(
    exprs.file, 
    pdat.file, 
    fdat.file, 
    de=TRUE,
    heatm.file=NULL, 
    distr.file=NULL, 
    volc.file=NULL)
{
    # read features
    anno <- ifelse(file.exists(fdat.file), NA, fdat.file)
    if(is.na(anno))
    {
        ncol.fdat <- length(scan(fdat.file, what="character", nlines=1, quiet=TRUE))
        fDat <- scan(fdat.file, what="character", quiet=TRUE)
        nr.features <- length(fDat) / ncol.fdat
        fDat <- matrix(fDat, nrow=nr.features, ncol=ncol.fdat, byrow=TRUE)
    }
    else 
    {
        fDat <- anno.p2g(fdat.file)
        ncol.fdat <- ncol(fDat)
        nr.features <- nrow(fDat)
    }

    # read samples
    ncol.pdat <- length(scan(pdat.file, what="character", nlines=1, quiet=TRUE))
    pDat <- scan(pdat.file, what="character", quiet=TRUE)
    nr.samples <- length(pDat) / ncol.pdat
    pDat <- matrix(pDat, nrow=nr.samples, ncol=ncol.pdat, byrow=TRUE)
    
    # read expression values
    expr <- matrix(scan(exprs.file, quiet=TRUE), 
        nrow=nr.features, ncol=nr.samples, byrow=TRUE)
    rownames(expr) <- fDat[,1]
    colnames(expr) <- pDat[,1]

    # replace NAs by mean expression values
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
    if(!is.na(anno)) annotation(eset) <- anno
    for(i in seq_len(ncol.pdat)) pData(eset)[, i] <- pDat[,i]
    for(i in seq_len(ncol.fdat)) fData(eset)[, i] <- fDat[,i]
    
    colnames(pData(eset))[1:2] <- c(SMPL.COL, GRP.COL)
    if(ncol.pdat > 2) colnames(pData(eset))[3] <- BLK.COL

    colnames(fData(eset))[1:2] <- c(PRB.COL, GN.COL)
    
    # (a) heatmap
    if(!is.null(heatm.file)) exprs.heatmap(eset, heatm.file)
   
    if(de)
    { 
        eset <- de.ana(eset)
        # (b) pval distribution
        if(!is.null(distr.file)) pdistr(eset, distr.file)   
    
        # (c) volcano plot
        if(!is.null(volc.file)) volcano(eset, volc.file)
    }
    
    return(eset)
}

# perform de analysis
de.ana <- function(eset)
{
    groups <- sort(unique(pData(eset)[,GRP.COL]), decreasing=TRUE)
    if(!all(groups == GRPS)) 
        stop(paste0("Group classification in pData is not binary:\n",
            "Expected (0, 1) but found (", paste(groups, collapse=", "), ")"))
    
    block <- 0
    paired <- BLK.COL %in% colnames(pData(eset))
    if(paired) block <- factor(pData(eset)[,BLK.COL])
    
    g <- factor(ifelse(pData(eset)[,GRP.COL] == "1", "d", "c"))
    if(paired) design <- model.matrix(~0 + g + block)
    else design <- model.matrix(~0 + g)
    colnames(design)[1:2] <- levels(g)
    fit <- lmFit(exprs(eset), design)
    cont.matrix <- makeContrasts(contrasts="d-c", levels=design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2)
    aT1 <- topTable(fit2, coef=1, 
        number=nrow(eset), sort.by="none", adjust.method="none")
    fData(eset)[, FC.COL] <- aT1[, "logFC"]
    fData(eset)[, RAWP.COL] <- aT1[, "P.Value"]
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



