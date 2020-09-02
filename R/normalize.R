############################################################
# 
# author: Ludwig Geistlinger
# date: 2015-03-10 17:26:50
# 
# descr: normalization of microarray and RNAseq data
# 
# EDIT 02 Mar 2020: VST for RNA-seq data  
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
#' is expected that normalization within arrays has been already carried
#' out, e.g. using \code{\link{normalizeWithinArrays}} from limma.
#' 
#' RNA-seq data is expected to be raw read counts.  Please note that
#' normalization for downstream DE analysis, e.g. with edgeR and DESeq2, is not
#' ultimately necessary (and in some cases even discouraged) as many of these
#' tools implement specific normalization approaches.  See the vignette of
#' EDASeq, edgeR, and DESeq2 for details.
#'
#' Using \code{norm.method = "vst"} invokes a variance-stabilizing 
#' transformation (VST) for RNA-seq read count data. This accounts for differences 
#' in sequencing depth between samples and over-dispersion of read count data. 
#' The VST uses the \code{\link{cpm}} function implemented in the edgeR package 
#' to compute moderated log2 read counts. Using edgeR's estimate of the common 
#' dispersion phi, the \code{prior.count} parameter of the \code{\link{cpm}} 
#' function is chosen as 0.5 / phi as previously suggested (Harrison, 2015).
#'
#' @param se An object of class \code{\linkS4class{SummarizedExperiment}}.
#' @param norm.method Determines how the expression data should be normalized.
#' For available microarray normalization methods see the man page of the limma
#' function \code{\link{normalizeBetweenArrays}}.  For available RNA-seq
#' normalization methods see the man page of the EDASeq function
#' \code{betweenLaneNormalization}.  Defaults to \code{'quantile'}, i.e.
#' normalization is carried out so that quantiles between arrays/lanes/samples
#' are equal. For RNA-seq data, this can also be \code{'vst'}, 
#' \code{'voom'}, or \code{'deseq2'} to invoke a variance-stabilizing transformation
#' that allows statistical modeling as for microarry data. See details.
#' @param data.type Expression data type.  Use \code{'ma'} for microarray and 
#' \code{'rseq'} for RNA-seq data.  If \code{NA}, the data type is automatically 
#' guessed: if the expression values in \code{se} are decimal (float) numbers, 
#' they are assumed to be microarray intensities;  whole (integer) numbers are 
#' assumed to be RNA-seq read counts.  Defaults to \code{NA}.
#' @param filter.by.expr Logical. For RNA-seq data: include only genes with
#' sufficiently large counts in the DE analysis? If TRUE, excludes genes not 
#' satisfying a minimum number of read counts across samples using the 
#' \code{\link{filterByExpr}} function from the edgeR package.
#' Defaults to TRUE.

#' @return An object of class \code{\linkS4class{SummarizedExperiment}}.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#'
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
#' @seealso \code{\link{readSE}} for reading expression data from file;
#' 
#' \code{\link{normalizeWithinArrays}} and \code{\link{normalizeBetweenArrays}}
#' for normalization of microarray data;
#' 
#' \code{withinLaneNormalization} and \code{betweenLaneNormalization} from the 
#' EDASeq package for normalization of RNA-seq data;
#'
#' \code{\link{cpm}}, \code{\link{estimateDisp}}, \code{\link{voom}}, and
#' \code{varianceStabilizingTransformation} from the DESeq2 package.
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
#'     maSE <- normalize(maSE) 
#'     assay(maSE, "raw")[1:5,1:5] 
#'     assay(maSE, "norm")[1:5,1:5] 
#' 
#'     # (b) RNA-seq ... 
#'     normSE <- normalize(rseqSE, norm.method = "vst") 
#'     assay(maSE, "raw")[1:5,1:5] 
#'     assay(maSE, "norm")[1:5,1:5] 
#'
#' @export normalize
normalize <- function(se, 
                      norm.method = "quantile", 
                      data.type = c(NA, "ma", "rseq"),
                      filter.by.expr = TRUE)
{
    # dealing with an SE?
    if(is.matrix(se)) se <- new("SummarizedExperiment", assays = list(raw = se))
    if(is(se, "ExpressionSet")) se <- as(se, "SummarizedExperiment")
    if(!is(se, "SummarizedExperiment")) 
        stop(paste("\'se\' must be either",
                   "a matrix, a SummarizedExperiment, or an ExpressionSet"))

    # ma or rseq data?
    if("dataType" %in% names(metadata(se)))
            data.type <- metadata(se)$dataType
    else
    {
        data.type <- match.arg(data.type)
        if(is.na(data.type)) data.type <- .detectDataType(assay(se))
        metadata(se)$dataType <- data.type  
    }

    # rseq normalization
    if(data.type == "rseq")
    {
        if(filter.by.expr)
        {
            # remove genes with low read count
            GRP.COL <- configEBrowser("GRP.COL")
            if(GRP.COL %in% colnames(colData(se))) grp <- se[[GRP.COL]]
            else grp <- sample(c(0,1), ncol(se), replace = TRUE, prob = c(0.5, 0.5))
        
            keep <- .filterRSeq(assay(se), group = grp, index.only = TRUE)
            se <- se[keep,]
        }
        if(norm.method %in% c("voom", "vst", "deseq2")) se <- .vst(se, norm.method)
        else
        {
            betweenLaneNormalization <- NULL
            isAvailable("EDASeq", type = "software")
            if(norm.method == "quantile") norm.method <- "upper"  
            assays(se)[[2]] <- betweenLaneNormalization(assay(se), 
                                                            which = norm.method)
        }
    } 
    # ma normalization with limma
    else assays(se)[[2]] <- limma::normalizeBetweenArrays(assay(se), 
                                                            method = norm.method)
    names(assays(se)) <- c("raw", "norm")
    return(se)
}

.vst <- function(se, method = c("vst", "voom", "deseq2"))
{
    method <- match.arg(method)
    design <- .getDesign(se)

    if(method == "deseq2")
    { 
        varianceStabilizingTransformation <- DESeqDataSet <- NULL
        isAvailable("DESeq2", type = "software")
        dds <- DESeqDataSet(se, design = design)
        dds <- varianceStabilizingTransformation(dds, blind = FALSE)
        cpms <- assay(dds) 
    }

    else
    {
        dge <- edgeR::DGEList(counts = assay(se), group = design[,"group1"])
        dge <- edgeR::calcNormFactors(dge)

        # anscombe
        if(method == "vst")
        {
            dge <- edgeR::estimateDisp(dge, design, robust = TRUE)
            pc <- 0.5 / dge$common.dispersion
            cpms <- edgeR::cpm(dge, log = TRUE, prior.count = pc)
        }
        # voom 
        else cpms <- limma::voom(dge, design)$E
        se <- se[rownames(dge),]
    }
    assays(se)[[2]] <- cpms
    return(se)
}

.getDesign <- function(se)
{
    GRP.COL <- configEBrowser("GRP.COL")
    BLK.COL <- configEBrowser("BLK.COL")

    grp <- se[[GRP.COL]]
    blk <- NULL
    if(BLK.COL %in% colnames(colData(se))) blk <- se[[BLK.COL]]

    group <- factor(grp)
    paired <- !is.null(blk)
    f <- "~" 
    if(paired) 
    {   
        block <- factor(blk)
        f <- paste0(f, "block + ") 
    }   
    f <- formula(paste0(f, "group"))
    model.matrix(f)
}
