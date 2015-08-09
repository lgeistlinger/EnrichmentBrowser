# perform de analysis
de.ana <- function(expr, grp=NULL, blk=NULL, 
    de.method=c("limma", "edgeR", "DESeq"), padj.method="BH", stat.only=FALSE)
{
    is.eset <- is(expr, "eSet")
    if(is.eset) 
    { 
        GRP.COL <- config.ebrowser("GRP.COL")    
        BLK.COL <- config.ebrowser("BLK.COL")

        eset <- expr
        if(is(expr, "SeqExpressionSet")) expr <- EDASeq::counts(eset)
        else expr <- exprs(eset)
        
        # check for group annotation
        if(!(GRP.COL %in% colnames(pData(eset))))
            stop(paste("Expression data \'expr\' is an ExpressionSet",
                "but contains no group assignment in the pData slot"))
        grp <- pData(eset)[,GRP.COL]

        # check for block annotation
        if(BLK.COL %in% colnames(pData(eset))) blk <- pData(eset)[,BLK.COL] 
    }
    if(class(expr) != "matrix")
        stop("Expression data in \'expr\' must be either a matrix or a (Seq)ExpressionSet")

    if(is.null(grp)) stop("Group assignment 'grp' must be specified") 
    groups <- sort(unique(grp, decreasing=TRUE))
    if(!all(groups == c(0, 1))) 
        stop(paste0("Group classification is not binary:\n",
            "Expected (0, 1) but found (", paste(groups, collapse=", "), ")"))
    
    group <- factor(grp)
    paired <- !is.null(blk)
    f <- "~"
    if(paired) 
    {
        block <- factor(blk)
        f <- paste0(f, "block + ") 
    }
    f <- formula(paste0(f, "group"))
 
    # EDGER
    de.method <- match.arg(de.method)
    if(de.method == "edgeR")
    {
        design <- model.matrix(f)
        y <- edgeR::DGEList(counts=expr,group=group)
        y <- edgeR::calcNormFactors(y)
        y <- edgeR::estimateGLMCommonDisp(y,design)
        y <- edgeR::estimateGLMTrendedDisp(y,design)
        y <- edgeR::estimateGLMTagwiseDisp(y,design)
        fit <- edgeR::glmFit(y,design)
        lrt <- edgeR::glmLRT(fit)
        if(stat.only) return(lrt$table[,"LR"])
        de.tbl <- lrt$table[, c("logFC", "PValue", "LR")] 
        colnames(de.tbl)[3] <- "edgeR.STAT"
    }
    # DESEQ
    else if(de.method == "DESeq")
    {
        colData <- data.frame(group=group)
        if(paired) colData$block <- block
        dds <- DESeq2::DESeqDataSetFromMatrix(
            countData=expr, colData=colData, design=f)
        dds <- suppressMessages(DESeq2::DESeq(dds))
        res <- DESeq2::results(dds, pAdjustMethod="none")
        if(stat.only) return(res[,"stat"])
        de.tbl <- data.frame(res[,c("log2FoldChange","pvalue","stat")])
        colnames(de.tbl)[3] <- "DESeq.STAT"
    }
    # LIMMA  
    else if(de.method == "limma")
    {
        design <- model.matrix(f)
        data.type <- auto.detect.data.type(expr)
        if(data.type == "rseq") 
        {
            dge <- edgeR::DGEList(counts=expr)
            dge <- edgeR::calcNormFactors(dge)
            expr <- limma::voom(dge, design)
        }
        fit <- limma::lmFit(expr, design)
        fit <- limma::eBayes(fit)
        aT1 <- limma::topTable(fit, number=nrow(expr), coef="group1", 
            sort.by="none", adjust.method="none")
        if(stat.only) return(aT1[,"t"])
        de.tbl <- aT1[, c("logFC", "P.Value", "t")]
        colnames(de.tbl)[3] <- "limma.STAT"
    }
    else stop(paste(de.method, "is not supported. See man page for supported de.method."))

    de.tbl[,2] <- p.adjust(de.tbl[,2], method=padj.method)
    colnames(de.tbl)[1:2] <- sapply(c("FC.COL", "ADJP.COL"), config.ebrowser)

    if(is.eset)
    {
        fData(eset)[,colnames(de.tbl)] <- de.tbl
        return(eset)
    }
    return(de.tbl)
}


