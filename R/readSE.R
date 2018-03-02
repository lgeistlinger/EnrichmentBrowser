# 
# 
# Author: ludwig geistlinger
# Date: May 27th, 2010
#
# reads in data from plain text files and stores it in an expression set
#
###############################################################################

read.eset <- function(exprs.file, pdat.file, fdat.file, 
    data.type=c(NA, "ma", "rseq"), NA.method=c("mean", "rm", "keep"))
{
    .Deprecated("readSE")
    readSE(assay.file=exprs.file, cdat.file=pdat.file,
        rdat.file=fdat.file, data.type=data.type, NA.method=NA.method)
}

readSE <- function(assay.file, cdat.file, rdat.file, 
    data.type=c(NA, "ma", "rseq"), NA.method=c("mean", "rm", "keep"))
{
    data.type <- match.arg(data.type)
    NA.method <- match.arg(NA.method)

    # read features
    anno <- ifelse(file.exists(rdat.file), NA, rdat.file)
    if(is.na(anno))
    {
        ncol.rdat <- length(scan(rdat.file, what="character", nlines=1, quiet=TRUE))
        rdat <- scan(rdat.file, what="character", quiet=TRUE)
        nr.features <- length(rdat) / ncol.rdat
        rdat <- matrix(rdat, nrow=nr.features, ncol=ncol.rdat, byrow=TRUE)
        rdat <- as.data.frame(rdat, stringsAsFactors=FALSE)
    }
    else 
    {
        rdat <- .annoP2G(rdat.file)
        ncol.rdat <- ncol(rdat)
        nr.features <- nrow(rdat)
    }

    # read samples
    ncol.cdat <- length(scan(cdat.file, what="character", nlines=1, quiet=TRUE))
    cdat <- scan(cdat.file, what="character", quiet=TRUE)
    nr.samples <- length(cdat) / ncol.cdat
    cdat <- matrix(cdat, nrow=nr.samples, ncol=ncol.cdat, byrow=TRUE)
    cdat <- as.data.frame(cdat, stringsAsFactors=FALSE)    
    cdat[,2] <- as.integer(cdat[,2])
    
    # read expression values
    expr <- matrix(scan(assay.file, quiet=TRUE), 
        nrow=nr.features, ncol=nr.samples, byrow=TRUE)
    rownames(expr) <- rdat[,1]
    colnames(expr) <- cdat[,1]

    # deal with NAs
    expr <- .naTreat(expr, NA.method) 
    if(NA.method=="rm") rdat <- rdat[rdat[,1] %in% rownames(expr),]
   
    # create the eset
    eset <- SummarizedExperiment(assays=list(exprs=expr))
    if(!is.na(anno)) metadata(eset)$annotation <- anno
    colData(eset) <- DataFrame(cdat, row.names=colnames(eset))
    rowData(eset) <- DataFrame(rdat, row.names=rownames(eset))
    
    colnames(colData(eset))[1:2] <- sapply(c("SMPL.COL", "GRP.COL"), config.ebrowser)
    if(ncol.cdat > 2) colnames(colData(eset))[3] <- config.ebrowser("BLK.COL")
    if(ncol.rdat == 1) colnames(rowData(eset))[1] <- config.ebrowser("EZ.COL") 
    else colnames(rowData(eset))[1:2] <- sapply(c("PRB.COL", "EZ.COL"), config.ebrowser)
   
    # ma or rseq?
    if(!(data.type %in% c("ma", "rseq"))) 
        data.type <- .detectDataType(expr)
    metadata(eset)$dataType <- data.type
 
    return(eset)
}

.naTreat <- function(expr, NA.method=c("mean", "rm", "keep"))
{
    NA.method <- match.arg(NA.method)
    # replace NAs by mean expression values
    if(!(NA.method %in% c("mean", "rm"))) return(expr)

    na.indr <- which(apply(expr, 1, function(x) any(is.na(x))))
    if(length(na.indr) == 0) return(expr)
    if(NA.method == "rm") return(expr[-na.indr,])

    data.type <- .detectDataType(expr)

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

.detectDataType <- function(expr) 
    ifelse(all(.isWholenumber(expr), na.rm=TRUE), "rseq", "ma")
    
.isWholenumber <- function(x, tol=.Machine$double.eps^0.5) abs(x-round(x)) < tol
