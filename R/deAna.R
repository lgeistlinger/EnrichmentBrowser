de.ana <- function(expr, grp=NULL, blk=NULL, 
                    de.method=c("limma", "edgeR", "DESeq"), 
                    padj.method="BH", stat.only=FALSE, min.cpm=2)
{
    .Deprecated("deAna")
    deAna(expr=expr, grp=grp, blk=blk, de.method=de.method,
            padj.method=padj.method, stat.only=stat.only, min.cpm=min.cpm)
}

deAna <- function(expr, grp=NULL, blk=NULL, 
                    de.method=c("limma", "edgeR", "DESeq"), 
                    padj.method="BH", stat.only=FALSE, min.cpm=2)
{
    ### TEMPORARY: will be replaced by as(eSet,SummarizedExperiment)
    if(is(expr, "ExpressionSet")) expr <- as(expr, "RangedSummarizedExperiment")
    ###

    isSE <- is(expr, "SummarizedExperiment")
    if(isSE) 
    { 
        GRP.COL <- config.ebrowser("GRP.COL")    
        BLK.COL <- config.ebrowser("BLK.COL")

        se <- expr
        expr <- assay(se)
        
        # check for group annotation
        if(!(GRP.COL %in% colnames(colData(se))))
            stop("Group assignment must be specified")
		grp <- colData(se)[,GRP.COL]

        # check for block annotation
        if(BLK.COL %in% colnames(colData(se))) blk <- colData(se)[,BLK.COL] 
    }
    if(!is.matrix(expr))
        stop(paste("Expression data in \'expr\' must be either", 
            "a matrix, a SummarizedExperiment, or an ExpressionSet"))

    if(is.null(grp)) stop("Group assignment 'grp' must be specified") 
    groups <- sort(unique(grp))
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
    
	de.method <- match.arg(de.method)
	data.type <- .detectDataType(expr)
	if(data.type == "rseq")
	{
		# filter low-expressed genes
		rs <- rowSums(edgeR::cpm(expr) > min.cpm)
		keep <-  rs >= ncol(expr) / 2
		nr.low <- sum(!keep)
		if(nr.low)
		{ 
			message(paste("Excluding", nr.low, 
                "genes not satisfying min.cpm threshold")) 
			expr <- expr[keep,]	
			if(isSE) se <- se[keep,]
		}	
	}
	else
	{
		# check for appropriate choice of DE method
		if(de.method %in% c("edgeR", "DESeq"))
			stop(paste(de.method, "only applicable to integer read counts"))
	} 


    # EDGER
    if(de.method == "edgeR")
    {
        # TODO: wait for edgeR_3.18.1 to remove this
        isAvailable("edgeR", type="software")

        y <- edgeR::DGEList(counts=expr,group=group)
        y <- edgeR::calcNormFactors(y)
        design <- model.matrix(f)
        if(length(group) == 2)
        { 
            message("Calling edgeR without replicates")
            message("Using default BCV (square-root-dispersion) of 0.4")
            fit <- edgeR::glmQLFit(y, design, dispersion=0.4)
        }
        else
        { 
            y <- edgeR::estimateDisp(y, design, robust=TRUE)
            fit <- edgeR::glmQLFit(y, design, robust=TRUE)
        }
        qlf <- edgeR::glmQLFTest(fit)
        if(stat.only) return(qlf$table[,"F"])
        de.tbl <- qlf$table[, c("logFC", "PValue", "F")] 
        colnames(de.tbl)[3] <- "edgeR.STAT"
    }
    # DESEQ
    else if(de.method == "DESeq")
    {
        colData <- data.frame(group=group)
        if(paired) colData$block <- block
        suppressMessages({
            dds <- DESeq2::DESeqDataSetFromMatrix(
                countData=expr, colData=colData, design=f)
            dds <- DESeq2::DESeq(dds)
        })
        res <- DESeq2::results(dds, pAdjustMethod="none")
        if(stat.only) return(res[,"stat"])
        de.tbl <- data.frame(res[,c("log2FoldChange","pvalue","stat")])
        colnames(de.tbl)[3] <- "DESeq.STAT"
    }
    # LIMMA  
    else if(de.method == "limma")
    {
        design <- model.matrix(f)
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
    else stop(paste(de.method, "is not supported.", 
                    "See man page for supported de.method."))

    de.tbl[,2] <- p.adjust(de.tbl[,2], method=padj.method)
    colnames(de.tbl)[1:2] <- sapply(c("FC.COL", "ADJP.COL"), config.ebrowser)

    if(isSE)
    {
        i <- grep("STAT$", colnames(rowData(se)))
        if(length(i)) rowData(se) <- rowData(se)[,-i]
        rowData(se)[,colnames(de.tbl)] <- DataFrame(de.tbl)
        return(se)
    }
    return(de.tbl)
}


