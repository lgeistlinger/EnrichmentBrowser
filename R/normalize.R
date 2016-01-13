############################################################
# 
# author: Ludwig Geistlinger
# date: 2015-03-10 17:26:50
# 
# descr: normalization of microarray and RNAseq data
# 
############################################################

normalize <- function(eset, norm.method="quantile", within=FALSE, data.type=c(NA, "ma", "rseq"))
{
    # dealing with an eset?
    if(class(eset) != "ExpressionSet") 
    {
        if(class(eset) == "matrix") eset <- new("ExpressionSet", exprs=eset)
        else stop("\'eset\' must be either a matrix or an ExpressionSet")
    }

    # ma or rseq data?
    
    if("dataType" %in% names(experimentData(eset)@other))
            data.type <- experimentData(eset)@other$dataType
    else
    {
        data.type <- match.arg(data.type)
        if(is.na(data.type))
            data.type <- auto.detect.data.type(exprs(eset))
        experimentData(eset)@other$dataType <- data.type  
    }

    # rseq normalization with EDASeq
    if(data.type == "rseq")
    {
        # remove genes with low read count
	is.too.low <- rowSums(exprs(eset)) < ncol(eset)
	nr.too.low <- sum(is.too.low)
        if(nr.too.low > 0) message(paste("Removing",
        	nr.too.low, "genes with low read count ..."))
	eset <- eset[!is.too.low,]
 
        if(norm.method == "quantile") norm.method <- "full"
        
        eset <- EDASeq::newSeqExpressionSet(
                    counts=exprs(eset), 
                    phenoData=pData(eset), 
                    featureData=fData(eset),
                    annotation=annotation(eset), 
                    experimentData=experimentData(eset),
                    protocolData=protocolData(eset))
        
        # also within lane? -> gc content normalization
        if(within)
        {
            gc.col <- grep("^[gG][cC]$", colnames(fData(eset)), value=TRUE)
            if(length(gc.col) == 0)
            {
                org <- annotation(eset)
                if(!length(org)) stop(paste("Please provide organism under", 
                    "investigation in the annotation slot. See man page for details."))
                MODEL.ORGS <- c("cel", "dme", "hsa", "mmu", "rno",  "sce")
                mode <- ifelse(org %in% MODEL.ORGS, "org.db", "biomart") 
                lgc <- EDASeq::getGeneLengthAndGCContent(
                    id=featureNames(eset), org=org, mode=mode)
                fData(eset) <- cbind(fData(eset), lgc)
                gc.col <- "gc"
            }
            message("Normalizing for GC content ...")
            na.gc <- is.na(fData(eset)[,"gc"])
            nr.na.gc <- sum(na.gc)
            if(nr.na.gc > 0) message(paste("Removing", 
                nr.na.gc, "genes due to missing GC content ..."))
            eset <- eset[!na.gc,]
            eset <- EDASeq::withinLaneNormalization(eset, "gc", which=norm.method)
        }
        eset <- EDASeq::betweenLaneNormalization(eset, which=norm.method)
    } 
    # ma normalization with limma
    else exprs(eset) <- limma::normalizeBetweenArrays(exprs(eset), method=norm.method)
    return(eset)
}


