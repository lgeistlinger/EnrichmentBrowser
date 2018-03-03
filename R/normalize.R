############################################################
# 
# author: Ludwig Geistlinger
# date: 2015-03-10 17:26:50
# 
# descr: normalization of microarray and RNAseq data
# 
############################################################

normalize <- function(se, norm.method="quantile", within=FALSE, data.type=c(NA, "ma", "rseq"))
{
    # dealing with an se?
    if(is(se, "ExpressionSet")) se <- as(se, "RangedSummarizedExperiment")
    
    if(!is(se, "SummarizedExperiment")) 
    {
        if(is.matrix(se)) 
            se <- new("SummarizedExperiment", assays=list(exprs=se))
        else stop(paste("\'se\' must be either",
            "a matrix, a SummarizedExperiment, or an ExpressionSet"))
    }

    # ma or rseq data?
    if("dataType" %in% names(metadata(se)))
            data.type <- metadata(se)$dataType
    else
    {
        data.type <- match.arg(data.type)
        if(is.na(data.type)) data.type <- .detectDataType(assay(se))
        metadata(se)$dataType <- data.type  
    }

    # rseq normalization with EDASeq
    if(data.type == "rseq")
    {
        # remove genes with low read count
	    is.too.low <- rowSums(assay(se)) < ncol(se)
	    nr.too.low <- sum(is.too.low)
        if(nr.too.low > 0) message(paste("Removing",
        	nr.too.low, "genes with low read count ..."))
	    se <- se[!is.too.low,]
 
        if(norm.method == "quantile") norm.method <- "full"
        
        # also within lane? -> gc content normalization
        if(within)
        {
            gc.col <- grep("^[gG][cC]$", colnames(rowData(se)), value=TRUE)
            if(length(gc.col) == 0)
            {
                org <- metadata(se)$annotation
                if(!length(org)) stop(paste("Please provide organism under", 
                    "investigation in the annotation slot. See man page for details."))
                MODEL.ORGS <- c("cel", "dme", "hsa", "mmu", "rno",  "sce")
                mode <- ifelse(org %in% MODEL.ORGS, "org.db", "biomart") 
                lgc <- EDASeq::getGeneLengthAndGCContent(
                        id=rownames(se), org=org, mode=mode)
                rowData(se)$gc <- lgc
                gc.col <- "gc"
            }
            message("Normalizing for GC content ...")
            na.gc <- is.na(rowData(se)[,gc.col])
            nr.na.gc <- sum(na.gc)
            if(nr.na.gc > 0) message(paste("Removing", 
                nr.na.gc, "genes due to missing GC content ..."))
            se <- se[!na.gc,]
            assay(se) <- EDASeq::withinLaneNormalization(
                assay(se), rowData(se)[,gc.col], which=norm.method)
        }
        assay(se) <- EDASeq::betweenLaneNormalization(assay(se), which=norm.method)
    } 
    # ma normalization with limma
    else assay(se) <- limma::normalizeBetweenArrays(assay(se), method=norm.method)
    
    return(se)
}


