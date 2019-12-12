#' Differential expression analysis between two sample groups
#' 
#' The function carries out a differential expression analysis between two
#' sample groups. Resulting fold changes and derived p-values are returned.
#' Raw p-values are corrected for multiple testing.
#' 
#' 
#' @aliases de.ana
#' @param expr Expression data.  A numeric matrix. Rows correspond to genes,
#' columns to samples.  Alternatively, this can also be an object of class
#' \code{\linkS4class{SummarizedExperiment}}.
#' @param grp *BINARY* group assignment for the samples.  Use '0' and '1' for
#' unaffected (controls) and affected (cases) samples, respectively.  If NULL,
#' this is assumed to be defined via a column named 'GROUP' in the
#' \code{\link{colData}} slot if 'expr' is a
#' \code{\linkS4class{SummarizedExperiment}}.
#' @param blk Optional. For paired samples or sample blocks.  This can also be
#' defined via a column named 'BLOCK' in the \code{\link{colData}} slot if
#' 'expr' is a \code{\linkS4class{SummarizedExperiment}}.
#' @param de.method Differential expression method.  Use 'limma' for microarray
#' and RNA-seq data.  Alternatively, differential expression for RNA-seq data
#' can be also calculated using edgeR ('edgeR') or DESeq2 ('DESeq2').  Defaults
#' to 'limma'.
#' @param padj.method Method for adjusting p-values to multiple testing.  For
#' available methods see the man page of the stats function
#' \code{\link{p.adjust}}.  Defaults to 'BH'.
#' @param stat.only Logical. Should only the test statistic be returned?  This
#' is mainly for internal use, in order to carry out permutation tests on the
#' DE statistic for each gene.  Defaults to FALSE.
#' @param filter.by.expr Logical. For RNA-seq data: include only genes with
#' sufficiently large counts in the DE analysis? If TRUE, excludes genes not 
#' satisfying a minimum number of read counts across samples using the 
#' \code{\link{filterByExpr}} function from the edgeR package.
#' Defaults to TRUE.
#' @return A DE-table with measures of differential expression for each
#' gene/row, i.e. a two-column matrix with log2 fold changes in the 1st column
#' and derived p-values in the 2nd column.  If 'expr' is a
#' \code{\linkS4class{SummarizedExperiment}}, the DE-table will be
#' automatically appended to the \code{\link{rowData}} slot.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{\link{readSE}} for reading expression data from file,
#' \code{\link{normalize}} for normalization of expression data,
#' \code{\link{voom}} for preprocessing of RNA-seq data, \code{\link{p.adjust}}
#' for multiple testing correction, \code{\link{eBayes}} for DE analysis with
#' limma, \code{\link{glmFit}} for DE analysis with edgeR, and
#' \code{DESeq} for DE analysis with DESeq2.
#' @examples
#' 
#'     # (1) microarray data: intensity measurements
#'     maSE <- makeExampleData(what="SE", type="ma")
#'     maSE <- deAna(maSE)
#'     rowData(maSE)
#'     
#'     # (2) RNA-seq data: read counts
#'     rseqSE <- makeExampleData(what="SE", type="rseq")
#'     rseqSE <- deAna(rseqSE, de.method="DESeq2")
#'     rowData(rseqSE)
#' 
#' @export deAna
deAna <- function(expr, grp=NULL, blk=NULL, 
                    de.method=c("limma", "edgeR", "DESeq2"), 
                    padj.method="BH", stat.only=FALSE, filter.by.expr=TRUE)
{
    if(is(expr, "ExpressionSet")) expr <- as(expr, "SummarizedExperiment")

    isSE <- is(expr, "SummarizedExperiment")
    if(isSE) 
    { 
        se <- expr
        info <- .extractInfoFromSE(se)
        expr <- info$expr
        grp <- info$grp
        blk <- info$blk
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

    if(data.type != "rseq" && de.method %in% c("edgeR", "DESeq2"))
        stop(paste(de.method, "only applicable to integer read counts"))

    # filter low-expressed genes
    if(data.type == "rseq" && !stat.only && filter.by.expr)
    {
        expr <- .filterRSeq(expr)
        if(isSE) se <- se[rownames(expr),]
	}

    # DE analysis
    if(de.method == "edgeR") de.tbl <- .edger(expr, group, f, stat.only)
    else if(de.method == "DESeq2") 
        de.tbl <- .deseq(expr, group, paired, block, f, stat.only)
    else if(de.method == "limma") de.tbl <- .limma(expr, f, data.type, stat.only)

    if(stat.only) return(de.tbl)

    # format result table
    colnames(de.tbl)[1:2] <- sapply(c("FC.COL", "PVAL.COL"), configEBrowser)
    de.tbl <- de.tbl[,c(1,3,2)]

    # multiple testing correction
    if(padj.method != "none")
    { 
        adjp <- p.adjust(de.tbl[,3], method=padj.method)
        de.tbl <- cbind(de.tbl, adjp)
        colnames(de.tbl)[4] <- configEBrowser("ADJP.COL")
    }

    if(isSE)
    {
        i <- grep("STAT$", colnames(rowData(se)))
        if(length(i)) rowData(se) <- rowData(se)[,-i]
        rowData(se)[,colnames(de.tbl)] <- DataFrame(de.tbl)
        return(se)
    }
    return(de.tbl)
}

#' @export
#' @keywords internal
de.ana <- function(expr, grp=NULL, blk=NULL, 
                    de.method=c("limma", "edgeR", "DESeq"), 
                    padj.method="BH", stat.only=FALSE)
{
    .Deprecated("deAna")
    deAna(expr=expr, grp=grp, blk=blk, de.method=de.method,
            padj.method=padj.method, stat.only=stat.only)
}

.limma <- function(expr, f, data.type, stat.only)
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
    return(de.tbl)
}

.edger <- function(expr, group, f, stat.only)
{
    dge <- edgeR::DGEList(counts=expr, group=group)
    dge <- edgeR::calcNormFactors(dge)
    design <- model.matrix(f)
    if(length(group) == 2)
    { 
        message("Calling edgeR without replicates")
        message("Using default BCV (square-root-dispersion) of 0.4")
        fit <- edgeR::glmQLFit(dge, design, dispersion=0.4)
    }
    else
    { 
        dge <- edgeR::estimateDisp(dge, design, robust=TRUE)
        fit <- edgeR::glmQLFit(dge, design, robust=TRUE)
    }
    qlf <- edgeR::glmQLFTest(fit)
    if(stat.only) return(qlf$table[,"F"])
    de.tbl <- qlf$table[, c("logFC", "PValue", "F")] 
    colnames(de.tbl)[3] <- "edgeR.STAT"
    return(de.tbl)
}

.deseq <- function(expr, group, paired, block, f, stat.only)
{
    DESeqDataSetFromMatrix <- DESeq <- results <- NULL
    isAvailable("DESeq2", type="software")
    
    colData <- data.frame(group=group)
    if(paired) colData$block <- block
    suppressMessages({
        dds <- DESeqDataSetFromMatrix(
            countData=expr, colData=colData, design=f)
        dds <- DESeq(dds)
    })
    res <- results(dds, pAdjustMethod="none")
    if(stat.only) return(res[,"stat"])
    de.tbl <- data.frame(res[,c("log2FoldChange","pvalue","stat")])
    colnames(de.tbl)[3] <- "DESeq2.STAT"
    return(de.tbl)
}

.filterRSeq <- function(expr)
{
    keep <- edgeR::filterByExpr(expr)
    nr.low <- sum(!keep)
    if(nr.low)
    { 
        message(paste("Excluding", nr.low, 
            "genes not satisfying min.cpm threshold")) 
        expr <- expr[keep,]	
    }
    return(expr)
}

.extractInfoFromSE <- function(se)
{
    GRP.COL <- configEBrowser("GRP.COL")    
    BLK.COL <- configEBrowser("BLK.COL")

    expr <- assay(se)
    
    # check for group annotation
    if(!(GRP.COL %in% colnames(colData(se))))
        stop("Group assignment must be specified")
    grp <- colData(se)[,GRP.COL]

    # check for block annotation
    blk <- NULL
    if(BLK.COL %in% colnames(colData(se))) blk <- colData(se)[,BLK.COL] 

    res <- list(expr=expr, grp=grp, blk=blk)
    return(res)
}
