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
#
# UPDATE: Oct 21th, 2014
#
###############################################################################

read.eset <- function(exprs.file, pdat.file, fdat.file, 
    data.type=c(NA, "ma", "rseq"), NA.method=c("mean", "rm", "keep"))
{
    data.type <- match.arg(data.type)
    NA.method <- match.arg(NA.method)

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

    # deal with NAs
    expr <- na.treat(expr, NA.method) 
    if(NA.method=="rm") fDat <- fDat[fDat[,1] %in% rownames(expr),]
   
    # create the eset
    eset <- new("ExpressionSet", exprs=expr)
    if(!is.na(anno)) annotation(eset) <- anno
    pData(eset) <- data.frame(pDat, stringsAsFactors=FALSE, row.names=sampleNames(eset))
    fData(eset) <- data.frame(fDat, stringsAsFactors=FALSE, row.names=featureNames(eset))
    
    colnames(pData(eset))[1:2] <- sapply(c("SMPL.COL", "GRP.COL"), config.ebrowser)
    if(ncol.pdat > 2) colnames(pData(eset))[3] <- config.ebrowser("BLK.COL")
    if(ncol.fdat == 1) colnames(fData(eset))[1] <- config.ebrowser("EZ.COL") 
    else colnames(fData(eset))[1:2] <- sapply(c("PRB.COL", "EZ.COL"), config.ebrowser)
   
    # ma or rseq?
    if(!(data.type %in% c("ma", "rseq"))) 
        data.type <- auto.detect.data.type(expr)
    experimentData(eset)@other$dataType <- data.type
 
    return(eset)
}

na.treat <- function(expr, NA.method=c("mean", "rm", "keep"))
{
    NA.method <- match.arg(NA.method)
    # replace NAs by mean expression values
    if(!(NA.method %in% c("mean", "rm"))) return(expr)

    na.indr <- which(apply(expr, 1, function(x) any(is.na(x))))
    if(length(na.indr) == 0) return(expr)
    if(NA.method == "rm") return(expr[-na.indr,])

    data.type <- auto.detect.data.type(expr)

    for(i in seq_along(na.indr))
    {
        cexpr <- expr[na.indr[i],]
        na.indc <- is.na(cexpr)
        m <- mean(cexpr[!na.indc])
        if(data.type == "rseq") m <- round(m)
        cexpr[na.indc] <- m         
        expr[na.indr[i],] <- cexpr
    }    
    return(expr)
}

auto.detect.data.type <- function(expr) 
    ifelse(all(is.wholenumber(expr), na.rm=TRUE), "rseq", "ma")
    
is.wholenumber <- function(x, tol=.Machine$double.eps^0.5) abs(x-round(x)) < tol
