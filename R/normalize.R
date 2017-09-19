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
    if(is(eset, "ExpressionSet")) eset <- as(eset, "RangedSummarizedExperiment")
    
    if(!is(eset, "SummarizedExperiment")) 
    {
        if(is.matrix(eset)) 
            eset <- new("SummarizedExperiment", assays=list(exprs=eset))
        else stop(paste("\'eset\' must be either",
            "a matrix, a SummarizedExperiment, or an ExpressionSet"))
    }

    # ma or rseq data?
    if("dataType" %in% names(metadata(eset)))
            data.type <- metadata(eset)$dataType
    else
    {
        data.type <- match.arg(data.type)
        if(is.na(data.type)) data.type <- .detectDataType(assay(eset))
        metadata(eset)$dataType <- data.type  
    }

    # rseq normalization with EDASeq
    if(data.type == "rseq")
    {
        # remove genes with low read count
	    is.too.low <- rowSums(assay(eset)) < ncol(eset)
	    nr.too.low <- sum(is.too.low)
        if(nr.too.low > 0) message(paste("Removing",
        	nr.too.low, "genes with low read count ..."))
	    eset <- eset[!is.too.low,]
 
        if(norm.method == "quantile") norm.method <- "full"
        
        # also within lane? -> gc content normalization
        if(within)
        {
            gc.col <- grep("^[gG][cC]$", colnames(rowData(eset)), value=TRUE)
            if(length(gc.col) == 0)
            {
                org <- metadata(eset)$annotation
                if(!length(org)) stop(paste("Please provide organism under", 
                    "investigation in the annotation slot. See man page for details."))
                MODEL.ORGS <- c("cel", "dme", "hsa", "mmu", "rno",  "sce")
                mode <- ifelse(org %in% MODEL.ORGS, "org.db", "biomart") 
                lgc <- EDASeq::getGeneLengthAndGCContent(
                        id=rownames(eset), org=org, mode=mode)
                rowData(eset)$gc <- lgc
                gc.col <- "gc"
            }
            message("Normalizing for GC content ...")
            na.gc <- is.na(rowData(eset)[,gc.col])
            nr.na.gc <- sum(na.gc)
            if(nr.na.gc > 0) message(paste("Removing", 
                nr.na.gc, "genes due to missing GC content ..."))
            eset <- eset[!na.gc,]
            assay(eset) <- EDASeq::withinLaneNormalization(
                assay(eset), rowData(eset)[,gc.col], which=norm.method)
        }
        assay(eset) <- EDASeq::betweenLaneNormalization(assay(eset), which=norm.method)
    } 
    # ma normalization with limma
    else assay(eset) <- limma::normalizeBetweenArrays(assay(eset), method=norm.method)
    
    return(eset)
}


