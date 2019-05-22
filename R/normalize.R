############################################################
# 
# author: Ludwig Geistlinger
# date: 2015-03-10 17:26:50
# 
# descr: normalization of microarray and RNAseq data
# 
############################################################



#' Normalization of microarray and RNA-seq expression data
#' 
#' This function wraps commonly used functionality from limma for microarray
#' normalization and from EDASeq for RNA-seq normalization.
#' 
#' Normalization of high-throughput expression data is essential to make
#' results within and between experiments comparable.  Microarray (intensity
#' measurements) and RNA-seq (read counts) data exhibit typically distinct
#' features that need to be normalized for.  For specific needs that deviate
#' from standard normalizations, the user should always refer to more
#' specific functions/packages.  See also the limma's user guide
#' \url{http://www.bioconductor.org/packages/limma} for definition and
#' normalization of the different expression data types.
#' 
#' Microarray data is expected to be single-channel.  For two-color arrays, it
#' is expected here that normalization within arrays has been already carried
#' out, e.g. using \code{\link{normalizeWithinArrays}} from limma.
#' 
#' RNA-seq data is expected to be raw read counts.  Please note that
#' normalization for downstream DE analysis, e.g. with edgeR and DESeq2, is not
#' ultimately necessary (and in some cases even discouraged) as many of these
#' tools implement specific normalization approaches.  See the vignette of
#' EDASeq, edgeR, and DESeq2 for details.
#' 
#' @param se An object of class \code{\linkS4class{SummarizedExperiment}}.
#' @param norm.method Determines how the expression data should be normalized.
#' For available microarray normalization methods see the man page of the limma
#' function \code{\link{normalizeBetweenArrays}}.  For available RNA-seq
#' normalization methods see the man page of the EDASeq function
#' \code{\link{betweenLaneNormalization}}.  Defaults to 'quantile', i.e.
#' normalization is carried out so that quantiles between arrays/lanes/samples
#' are equal. See details.
#' @param within Logical.  Is only taken into account if data.type='rseq'.
#' Determine whether GC content normalization should be carried out (as
#' implemented in the EDASeq function \code{\link{withinLaneNormalization}}).
#' Defaults to FALSE. See details.
#' @param data.type Expression data type.  Use 'ma' for microarray and 'rseq'
#' for RNA-seq data.  If NA, data.type is automatically guessed.  If the
#' expression values in 'se' are decimal numbers they are assumed to be
#' microarray intensities.  Whole numbers are assumed to be RNA-seq read
#' counts.  Defaults to NA.
#' @return An object of class \code{\linkS4class{SummarizedExperiment}}.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{\link{readSE}} for reading expression data from file;
#' 
#' \code{\link{normalizeWithinArrays}} and \code{\link{normalizeBetweenArrays}}
#' for normalization of microarray data;
#' 
#' \code{\link{withinLaneNormalization}} and
#' \code{\link{betweenLaneNormalization}} for normalization of RNA-seq data.
#' @examples
#' 
#'     #
#'     # (1) simulating expression data: 100 genes, 12 samples
#'     #
#'     
#'     # (a) microarray data: intensity measurements
#'     maSE <- makeExampleData(what="SE", type="ma")
#'     
#'     # (b) RNA-seq data: read counts
#'     rseqSE <- makeExampleData(what="SE", type="rseq")
#' 
#'     #
#'     # (2) Normalization
#'     #
#'     
#'     # (a) microarray ... 
#'     normSE <- normalize(maSE) 
#' 
#'     # (b) RNA-seq ... 
#'     normSE <- normalize(rseqSE) 
#' 
#'     # ... normalize also for GC content
#'     gc.content <- rnorm(100, 0.5, sd=0.1)
#'     rowData(rseqSE)$gc <- gc.content 
#' 
#'     normSE <- normalize(rseqSE, within=TRUE)
#' 
#' @export normalize
normalize <- function(se, norm.method="quantile", within=FALSE, data.type=c(NA, "ma", "rseq"))
{
    # dealing with an se?
    if(is(se, "ExpressionSet")) se <- as(se, "SummarizedExperiment")
    
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

#' Variance-stabilizing transformation for RNA-seq expression data
#' 
#' This function implements a variance-stabilizing transformation (VST) for 
#' RNA-seq read count data. It accounts for differences in sequencing depth
#' between samples and over-dispersion of read count data. Permutation-based
#' enrichment methods can then be applied as for microarray data. 
#' 
#' The VST uses the cpm function implemented in the edgeR package to compute 
#' moderated log2 read counts. Using edgeR's estimate of the common dispersion 
#' phi, the prior.count parameter of the cpm function is chosen as 0.5 / phi as 
#' previously suggested (Harrison, 2015).
#' 
#' @param se An object of class \code{\linkS4class{SummarizedExperiment}}.
#' @return An object of class \code{\linkS4class{SummarizedExperiment}}.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{\link{cpm}} and \code{\link{estimateDisp}}
#' @references
#' Harrison (2015) Anscombe's 1948 variance stabilizing transformation for
#' the negative binomial distribution is well suited to RNA-seq expression
#' data. doi:10.7490/f1000research.1110757.1
#'
#' Anscombe (1948) The transformation of Poisson, binomial and
#' negative-binomial data. Biometrika 35(3-4):246-54.
#'
#' Law et al. (2014) voom: precision weights unlock linear model analysis tools 
#' for RNA-seq read counts. Genome Biol 15:29.
#' 
#' @examples
#' 
#'     se <- makeExampleData(what="SE", type="rseq")
#'     vstSE <- vst(se) 
#' 
#' @export vst
vst <- function(se)
{
    expr <- assay(se)
    GRP.COL <- configEBrowser("GRP.COL")
    BLK.COL <- configEBrowser("BLK.COL")

    grp <- colData(se)[,GRP.COL]
    blk <- NULL
    if(BLK.COL %in% colnames(colData(se))) blk <- colData(se)[,BLK.COL]

    group <- factor(grp)
    paired <- !is.null(blk)
    f <- "~" 
    if(paired) 
    {   
        block <- factor(blk)
        f <- paste0(f, "block + ") 
    }   
    f <- formula(paste0(f, "group"))
    design <- model.matrix(f)

    dge <- edgeR::DGEList(counts=expr, group=group)
    dge <- .filterRSeq(dge)
    dge <- edgeR::calcNormFactors(dge)
    dge <- edgeR::estimateDisp(dge, design, robust=TRUE)
    
    pc <- 0.5 / dge$common.dispersion
    cpms <- edgeR::cpm(dge, log=TRUE, prior.count=pc)

    se <- se[rownames(dge),]
    assay(se) <- cpms
    metadata(se)$dataType <- "ma"
    return(se)
}
