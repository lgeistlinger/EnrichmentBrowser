# 
# 
# Author: ludwig geistlinger
# Date: May 27th, 2010
#
# reads in data from plain text files and stores it in an expression set
#
###############################################################################

#' Reading gene expression data from file
#' 
#' The function reads in plain expression data from file with minimum
#' annotation requirements for the colData and rowData slots.
#' 
#' 
#' @aliases read.eset
#' @param assay.file Expression matrix.  A tab separated text file containing
#' expression values.  Columns = samples/subjects; rows =
#' features/probes/genes; NO headers, row or column names.  See details.
#' @param cdat.file Column (phenotype) data.  A tab separated text file
#' containing annotation information for the samples in either *two or three*
#' columns.  NO headers, row or column names.  The number of rows/samples in
#' this file should match the number of columns/samples of the expression
#' matrix.  The 1st column is reserved for the sample IDs; The 2nd column is
#' reserved for a *BINARY* group assignment.  Use '0' and '1' for unaffected
#' (controls) and affected (cases) sample class, respectively.  For paired
#' samples or sample blocks a third column is expected that defines the blocks.
#' @param rdat.file Row (feature) data.  A tab separated text file containing
#' annotation information for the features.  In case of probe level data:
#' exactly *TWO* columns; 1st col = probe/feature IDs; 2nd col = corresponding
#' gene ID for each feature ID in 1st col; In case of gene level data: The list
#' of gene IDs newline-separated (i.e. just one column).  It is recommended to
#' use *ENTREZ* gene IDs (to benefit from downstream visualization and
#' exploration functionality of the enrichment analysis).  NO headers, row or
#' column names.  The number of rows (features/probes/genes) in this file
#' should match the number of rows/features of the expression matrix.
#' Alternatively, this can also be the ID of a recognized platform such as
#' 'hgu95av2' (Affymetrix Human Genome U95 chip) or 'ecoli2' (Affymetrix E.
#' coli Genome 2.0 Array).
#' @param data.type Expression data type.  Use 'ma' for microarray and 'rseq'
#' for RNA-seq data.  If NA, data.type is automatically guessed.  If the
#' expression values in the expression matrix are decimal numbers, they are
#' assumed to be microarray intensities.  Whole numbers are assumed to be
#' RNA-seq read counts.  Defaults to NA.
#' @param NA.method Determines how to deal with NA's (missing values).  This
#' can be one out of: \itemize{ \item mean: replace NA's by the row means for a
#' feature over all samples.  \item rm: rows (features) that contain NA's are
#' removed.  \item keep: do nothing. Missing values are kept (which, however,
#' can then cause several issues in the downstream analysis) } Defaults to
#' 'mean'.
#' @return An object of class \code{\linkS4class{SummarizedExperiment}}.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{\linkS4class{SummarizedExperiment}}
#' @examples
#' 
#'     # reading the expression data from file
#'     assay.file <- system.file("extdata/exprs.tab", package="EnrichmentBrowser")
#'     cdat.file <- system.file("extdata/colData.tab", package="EnrichmentBrowser")
#'     rdat.file <- system.file("extdata/rowData.tab", package="EnrichmentBrowser")
#'     se <- readSE(assay.file, cdat.file, rdat.file)
#' 
#' @export readSE
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
   
    # create the se
    se <- SummarizedExperiment(assays=list(exprs=expr))
    if(!is.na(anno)) metadata(se)$annotation <- anno
    colData(se) <- DataFrame(cdat, row.names=colnames(se))
    rowData(se) <- DataFrame(rdat, row.names=rownames(se))
    
    colnames(colData(se))[1:2] <- sapply(c("SMPL.COL", "GRP.COL"), configEBrowser)
    if(ncol.cdat > 2) colnames(colData(se))[3] <- configEBrowser("BLK.COL")
    if(ncol.rdat == 1) colnames(rowData(se))[1] <- configEBrowser("EZ.COL") 
    else colnames(rowData(se))[1:2] <- sapply(c("PRB.COL", "EZ.COL"), configEBrowser)
   
    # ma or rseq?
    if(!(data.type %in% c("ma", "rseq"))) data.type <- .detectDataType(expr)
    metadata(se)$dataType <- data.type
 
    return(se)
}

#' @export
#' @keywords internal
read.eset <- function(exprs.file, pdat.file, fdat.file, 
    data.type=c(NA, "ma", "rseq"), NA.method=c("mean", "rm", "keep"))
{
    .Deprecated("readSE")
    readSE(assay.file=exprs.file, cdat.file=pdat.file,
        rdat.file=fdat.file, data.type=data.type, NA.method=NA.method)
}

.naTreat <- function(expr, NA.method=c("mean", "rm", "keep"))
{
    NA.method <- match.arg(NA.method)
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
