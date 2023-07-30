############################################################
# 
# author: Ludwig Geistlinger
# date: 2020-06-22 16:38:50
# 
# descr: functionality for importing datasets from limma,
#        edgeR, and DESeq2
# 
############################################################

#' Import results from differential expression (DE) analysis 
#' 
#' This function imports fully processed expression datasets and results of 
#' differential expression (DE) analysis from limma, edgeR, and DESeq2. 
#' The imported data is converted to a \code{\linkS4class{SummarizedExperiment}},
#' annotating experimental design and genewise differential expression, which
#' then allows straightforward application of enrichment analysis methods. 
#' 
#' The expression data object (argument \code{obj}) is expected to be fully
#' processed (including normalization and dispersion estimation) and to have
#' the experimental design annotated.
#' The experimental design is expected to describe *a comparison of two groups*
#' with an optional blocking variable for paired samples / sample batches (i.e.
#' \code{design = ~ group} or \code{design = ~ batch + group}.)
#'
#' The differential expression result (argument \code{res}) is expected to have
#' the same number of rows as the expression data object,
#' and also that the order of the rows is the same / consistent, i.e. that there
#' is a 1:1 correspondence between the rownames of \code{obj} and the rownames 
#' of \code{res}. Note that the expression dataset is automatically restricted
#' to the genes for which DE results are available. However, for an appropriate
#' estimation of the size of the universe for competitive gene set tests, it is
#' recommended to provide DE results for all genes in the expression data object
#' whenever possible (see examples).
#'
#' @param obj Expression data object. Supported options include 
#' \code{\linkS4class{EList}} (voom/limma),
#' \code{\linkS4class{DGEList}} (edgeR), and \code{DESeqDataSet}
#' (DESeq2). See details.
#' @param res Differential expression results. Expected to match the provided
#' expression data object type, i.e. should be an object of class \itemize{ 
#' \item \code{\link{data.frame}} if \code{obj} is provided as an
#' \code{\linkS4class{EList}} (voom/limma),
#' \item \code{\linkS4class{TopTags}} if \code{obj} is provided as a
#' \code{\linkS4class{DGEList}} (edgeR), and
#' \item \code{DESeqResults} if \code{obj} is provided as a
#' \code{DESeqDataSet} (DESeq2). See details.}
#' @param from Character. Differential expression method from which to import 
#' results from. Defaults to \code{"auto"}, which automatically determines
#' the import type based on the expression data object provided. 
#' Can also be explicitly chosen as either \code{'limma'}, \code{'edgeR'} or
#' \code{'DESeq2'}.
#' @param anno Character. The organism under study in KEGG three letter
#' code (e.g. \sQuote{hsa} for \sQuote{Homo sapiens}). For microarray probe-level
#' data: the ID of a recognized microarray platform. See references. 
#' @return An object of class \code{\linkS4class{SummarizedExperiment}} that 
#' stores \itemize{
#' \item the expression matrix in the \code{assay} slot,
#' \item information about the samples, including the experimental design, in the
#' \code{colData} slot, and
#' \item information about the genes, including measures of differential expression,
#' in the \code{rowData} slot.}
#'
#' Mandatory annotations:
#' \itemize{ \item colData column storing binary group assignment (named
#' "GROUP") \item rowData column storing (log2) fold changes of differential
#' expression between sample groups (named "FC") \item rowData column storing
#' adjusted (corrected for multiple testing) p-values of differential
#' expression between sample groups (named "ADJ.PVAL"). } Additional optional
#' annotations: \itemize{ \item colData column defining paired samples or
#' sample blocks (named "BLOCK") \item metadata slot named "annotation" giving
#' the organism under investigation in KEGG three letter code (e.g. "hsa" for
#' Homo sapiens) \item metadata slot named "dataType" indicating the expression
#' data type ("ma" for microarray, "rseq" for RNA-seq). }
#' @author Ludwig Geistlinger
#' @seealso \code{\link{readSE}} for reading expression data from file,
#' \code{\link{normalize}} for normalization of expression data,
#' \code{\link{voom}} for preprocessing of RNA-seq data, \code{\link{p.adjust}}
#' for multiple testing correction, \code{\link{eBayes}} for DE analysis with
#' limma, \code{\link{glmQLFit}} for DE analysis with edgeR, and
#' \code{DESeq} for DE analysis with DESeq2.
#' @references KEGG Organism code
#' \url{http://www.genome.jp/kegg/catalog/org_list.html}
#' 
#' Recognized microarray platforms
#' \url{http://www.bioconductor.org/packages/release/BiocViews.html#___ChipName}
#' @examples
#' 
#'   # Setup
#'   ## i) Generate count data
#'   nsamples <- 4
#'   ngenes <- 1000
#'   dispers <- 1 / rchisq(ngenes, df = 10)
#'   rdesign <- model.matrix(~factor(rep(c(1, 2), each = 2)))
#'     
#'   counts <- rnbinom(ngenes * nsamples, mu = 20, size = 1 / dispers)
#'   counts <- matrix(counts, nrow = ngenes, ncol = nsamples)
#'   
#'   ## ii) Generate intensity data
#'   sd <- 0.3 * sqrt(4 / rchisq(100, df = 4))
#'   intens <- matrix(rnorm(100 * 6, sd = sd), nrow = 100, ncol = 6)
#'   rownames(intens) <- paste0("Gene", 1:100)
#'   intens[1:2, 4:6] <- intens[1:2, 4:6] + 2
#'   mdesign <- cbind(Grp1 = 1, Grp2vs1 = rep(c(0,1), each = 3))
#'
#'   # (1) import from edgeR (RNA-seq count data)
#'   # (1a) create the expression data object
#'   library(edgeR)
#'   d <- DGEList(counts)
#'   d <- calcNormFactors(d)
#'   d <- estimateDisp(d, rdesign)
#'   
#'   # (1b) obtain differential expression results 
#'   fit <- glmQLFit(d, rdesign)
#'   qlf <- glmQLFTest(fit)
#'   res <- topTags(qlf, n = nrow(d), sort.by = "none")
#'
#'   # (1c) import
#'   se <- import(d, res)
#'
#'   # (2) import from DESeq2 (RNA-seq count data)
#'   # (2a) create the expression data object
#'   library(DESeq2)
#'   condition <- factor(rep(c("A", "B"), each = 2))
#'   dds <- DESeqDataSetFromMatrix(counts, 
#'                                 colData = DataFrame(condition = condition),
#'                                 design = ~ condition)
#'
#'   # (2b) obtain differential expression results 
#'   dds <- DESeq(dds)
#'   res <- results(dds)
#'
#'   # (2c) import
#'   se <- import(dds, res)
#'
#'   # (3) import from voom/limma (RNA-seq count data)
#'   # (3a) create the expression data object
#'   library(limma)
#'   keep <- filterByExpr(counts, rdesign)
#'   el <- voom(counts[keep,], rdesign)
#'   
#'   # (3b) obtain differential expression results 
#'   fit <- lmFit(el, rdesign)
#'   fit <- eBayes(fit, robust = TRUE) 
#'   res <- topTable(fit, coef = 2, number = nrow(counts), sort.by = "none")
#'
#'   # (3c) import
#'   se <- import(el, res)
#'   
#'   # (4) import from limma-trend (RNA-seq count data)
#'   # (4a) create the expression data object
#'   logCPM <- edgeR::cpm(counts[keep,], log = TRUE, prior.count = 3)
#'   el <- new("EList", list(E = logCPM, design = rdesign))
#'   
#'   # (4b) obtain differential expression results 
#'   fit <- lmFit(el, rdesign)
#'   fit <- eBayes(fit, trend = TRUE) 
#'   res <- topTable(fit, coef = 2, number = nrow(el), sort.by = "none")
#'
#'   # (4c) import
#'   se <- import(el, res)
#'
#'   # (5) import from limma (microarray intensity data)
#'   # (5a) create the expression data object
#'   el <- new("EList", list(E = intens, design = mdesign))
#'   
#'   # (5b) obtain differential expression results 
#'   fit <- lmFit(el, mdesign)
#'   fit <- eBayes(fit, robust = TRUE) 
#'   res <- topTable(fit, coef = 2, number = nrow(el), sort.by = "none")
#'
#'   # (5c) import
#'   se <- import(el, res)
#' 
#' @export
import <- function(obj, res,
                    from = c("auto", "limma", "edgeR", "DESeq2"), 
                    anno = NA)
{
    from <- match.arg(from)
    if(from == "auto") from <- .detectImportType(obj)

    if(nrow(obj) != nrow(res) || !all(rownames(obj) == rownames(res)))
    {
        n <- intersect(rownames(obj), rownames(res))
        len <- length(n)
        if(!len) stop("\'obj\' and \'res\' have no gene IDs in common")
        if(len < nrow(obj))
        { 
            message("Restricting \'obj\' and \'res\' to ", 
                    len, " common gene IDs")
            message("It is recommended to provide DE results",
                    " for all genes in \'obj\'")
        }
        obj <- obj[n,]
        res <- res[n,]
    }
    
    if(from == "limma") se <- .importFromLimma(obj, res)
    else if(from == "edgeR") se <- .importFromEdgeR(obj, res) 
    else se <- .importFromDESeq2(obj, res)

    if(!is.na(anno)) metadata(se)$annotation <- anno
 
    return(se)
}

.detectImportType <- function(obj)
{
    type <- NA

    # limma
    if(is(obj, "EList")) type <- "limma"
    # voom / edgeR
    else if(is(obj, "DGEList")) type <- "edgeR"
    # DESeq2
    else if(is(obj, "DESeqDataSet")) type <- "DESeq2"

    if(is.na(type)) stop("Invalid import object. Supported options include",
                         " EList (voom/limma), DGEList (edgeR),",
                         " and DESeqDataSet (DESeq2).")
    return(type)
}



# voom:
# keep <- filterByExpr(counts, design)
# el <- voom(counts[keep,], design)

# limma:
# el <- normalizeBetweenArrays(expr)
# fit <- lmFit(el, design=c(-1,1,-1,1))
# fit <- eBayes(fit)
# res <- topTable(fit, coef = "grouptrt", number = nrow(expr), sort.by = "none")
.importFromLimma <- function(obj, res)
{
    stopifnot(is(obj, "EList") || is(obj, "DGEList"))
    stopifnot(is.data.frame(res))
    rnames <- c("logFC", "t", "P.Value", "adj.P.Val")
    stopifnot(all(rnames %in% colnames(res)))

    se <- SummarizedExperiment(assays = list(exprs = obj$E))
    if("weights" %in% names(obj)) 
        assay(se, "weights", withDimnames = FALSE) <- obj$weights

    # (1) rowData
    rrnames <- c("FC.COL", "PVAL.COL", "ADJP.COL")
    rrnames <- vapply(rrnames, configEBrowser, character(1), USE.NAMES = FALSE)
    ind <- match(rnames, colnames(res))
    colnames(res)[ind] <- c(rrnames[1], "limma.STAT", rrnames[2:3])
    if("genes" %in% names(obj)) res <- cbind(res, obj$genes)
    rowData(se) <- res

    # (2) colData
    if("targets" %in% names(obj)) colData(se) <- DataFrame(obj$targets)
    vars <- names(attr(obj$design, "contrasts"))
    design <- obj$design[, -1, drop = FALSE]
    gcol <- configEBrowser("GRP.COL")
    grp.col <- design[,ncol(design)]
    ngrps <- length(unique(grp.col))
    if(length(vars) > 2 || ngrps > 2) 
        stop("Supported experimental designs include binary",
             " group comparisons with an optional blocking",
             " variable for paired samples / sample batches")

    # (2i) group
    se[[gcol]] <- unname(grp.col)
    
    # (2ii) block
    if(length(vars) > 1)
    {
        blk.cols <- design[,seq_len(ncol(design) - 1)]
        colnames(blk.cols) <- sub(vars[1], "", colnames(blk.cols))
        for(i in seq_len(ncol(blk.cols))) 
            blk.cols[,i] <- ifelse(blk.cols[,i] == 1, i, 0)
        blk.col <- rowSums(blk.cols)
        blk.col[blk.col == 0] <- ncol(blk.cols) + 1
        blk.col <- paste0("block", blk.col) 
        bcol <- configEBrowser("BLK.COL")
        se[[bcol]] <- blk.col 
    }

    # (3) metadata 
    metadata(se)$dataType <- "ma"
    return(se)
}



# obj:
# group <- dds$dex
# y <- DGEList(counts=assay(dds),group=group)
# keep <- filterByExpr(y)
# y <- y[keep,,keep.lib.sizes=FALSE]
# y <- calcNormFactors(y)
# block <- dds$cell
# design <- model.matrix(~ block + group)
# y <- estimateDisp(y, design)

# res:
# fit <- glmQLFit(y,design)
# qlf <- glmQLFTest(fit)
# res <- topTags(qlf, n = nrow(y), sort.by = "none")
.importFromEdgeR <- function(obj, res)
{
    stopifnot(is(obj, "DGEList") && is(res, "TopTags"))
    res <- res$table
    padj <- grep("^F[WED]{1,2}R$", colnames(res), value = TRUE)
    rnames <- c("logFC", "F", "PValue", padj)
    stopifnot(all(rnames %in% colnames(res)))

    se <- SummarizedExperiment(assays = list(counts = obj$counts),
                               colData = obj$samples)

    assay(se, "norm") <- limma::voom(obj)$E
    names(assays(se))[1] <- "raw" 

    # (1) rowData
    rrnames <- c("FC.COL", "PVAL.COL", "ADJP.COL")
    rrnames <- vapply(rrnames, configEBrowser, character(1), USE.NAMES = FALSE)
    ind <- match(rnames, colnames(res))
    colnames(res)[ind] <- c(rrnames[1], "edgeR.STAT", rrnames[2:3])
    rowData(se) <- res

    # (2) colData
    vars <- names(attr(obj$design, "contrasts"))
    design <- obj$design[, -1, drop = FALSE]
    gcol <- configEBrowser("GRP.COL")
    grp.col <- design[,ncol(design)]
    ngrps <- length(unique(grp.col))
    if(length(vars) > 2 || ngrps > 2) 
        stop("Supported experimental designs include binary",
             " group comparisons with an optional blocking",
             " variable for paired samples / sample batches")

    # (2i) group
    se[[gcol]] <- unname(grp.col)
    
    # (2ii) block
    if(length(vars) > 1)
    {
        blk.cols <- design[,seq_len(ncol(design) - 1)]
        colnames(blk.cols) <- sub(vars[1], "", colnames(blk.cols))
        for(i in seq_len(ncol(blk.cols))) 
            blk.cols[,i] <- ifelse(blk.cols[,i] == 1, i, 0)
        blk.col <- rowSums(blk.cols)
        blk.col[blk.col == 0] <- ncol(blk.cols) + 1
        blk.col <- paste0("block", blk.col) 
        bcol <- configEBrowser("BLK.COL")
        se[[bcol]] <- blk.col 
    }

    # (3) metadata 
    metadata(se)$dataType <- "rseq"
    metadata(se)$common.dispersion <- obj$common.dispersion
    return(se)
}

# obj <- DESeq(dds)
# res <- results(obj)
.importFromDESeq2 <- function(obj, res)
{
    stopifnot(is(obj, "DESeqDataSet") && is(res, "DESeqResults"))
    rnames <- c("log2FoldChange", "stat", "pvalue", "padj")
    stopifnot(all(rnames %in% colnames(res)))

    varianceStabilizingTransformation <- NULL
    isAvailable("DESeq2", type = "software")
    se <- SummarizedExperiment(assays = assays(obj)[1],
                               colData = colData(obj),
                               rowData = rowData(obj),
                               metadata = metadata(obj))

    dts <- varianceStabilizingTransformation(obj, blind = FALSE)
    assay(se, "norm", withDimnames = FALSE) <- assay(dts)
    names(assays(se))[1] <- "raw" 

    # (1) rowData
    if(length(unlist(rowRanges(obj)))) rowRanges(se) <- rowRanges(obj)
    snames <- setdiff(colnames(rowData(se)), colnames(res))
    if(length(snames)) rowData(se) <- rowData(se)[,snames]
     
    rrnames <- c("FC.COL", "PVAL.COL", "ADJP.COL")
    rrnames <- vapply(rrnames, configEBrowser, character(1), USE.NAMES = FALSE)
    ind <- match(rnames, colnames(res))
    colnames(res)[ind] <- c(rrnames[1], "DESeq2.STAT", rrnames[2:3])
    rowData(se) <- cbind(res, rowData(se))
    
    # metadata 
    design <- obj@design
    metadata(se)$design <- design 
    metadata(se)$dataType <- "rseq"

    # (2) colData
    design <- as.character(design)[2]
    design <- unlist(strsplit(design, " \\+ "))
    if(length(design) > 2) stop("Supported experimental designs include binary",
                                " group comparisons with an optional blocking",
                                " variable for paired samples / sample batches")

    # (2i) group
    grp.col <- design[length(design)]
    ref <- levels(obj[[grp.col]])[1]
    gcol <- configEBrowser("GRP.COL")
    se[[gcol]] <- ifelse(obj[[grp.col]] == ref, 0, 1)
    
    # (2ii) block
    if(length(design) > 1)
    {
        blk.col <- design[1]
        bcol <- configEBrowser("BLK.COL")
        se[[bcol]] <- se[[blk.col]] 
    }
    
    return(se)
} 

