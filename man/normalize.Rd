\name{normalize}
\alias{normalize}
\title{ 
    Normalization of microarray and RNA-seq expression data 
}
\description{
    This function wraps commonly used functionality from limma for microarray
    normalization and from EDASeq for RNA-seq normalization.
}
\usage{
    normalize( eset, 
        norm.method = "quantile", within = FALSE, data.type = c(NA, "ma", "rseq") )
}
\arguments{
    \item{eset}{
        Expression set. 
        An object of \code{\link[Biobase:ExpressionSet-class]{ExpressionSet-class}}.
        See the man page of \code{\link{read.eset}} for prerequisites 
        for the expression data.
    }
    \item{norm.method}{
        Determines how the expression data should be normalized.
        For available microarray normalization methods see the man page of the 
        limma function \code{\link{normalizeBetweenArrays}}.
        For available RNA-seq normalization methods see the man page of the 
        EDASeq function \code{\link{betweenLaneNormalization}}.
        Defaults to 'quantile', i.e. normalization is carried out so that 
        quantiles between arrays/lanes/samples are equal. See details.
    }
    \item{within}{
        Logical.
        Is only taken into account if data.type='rseq'.
        Determine whether GC content normalization should be carried out (as 
        implemented in the EDASeq function \code{\link{withinLaneNormalization}}).
        Defaults to FALSE. See details.
    }
    \item{data.type}{
        Expression data type.
        Use 'ma' for microarray and 'rseq' for RNA-seq data.
        If NA, data.type is automatically guessed. 
        If the expression values in 'eset' are decimal numbers 
        they are assumed to be microarray intensities.
        Whole numbers are assumed to be RNA-seq read counts. 
        Defaults to NA.
    }
}
\value{ 
    An object of \code{\link[Biobase:ExpressionSet-class]{ExpressionSet-class}}.
    For RNA-seq data, an object of 
    \code{\link[EDASeq:SeqExpressionSet-class]{SeqExpressionSet-class}}
    to conform with downstream DE analysis.
}
\details{
    Normalization of high-throughput expression data is essential
    to make results within and between experiments comparable.
    Microarray (intensity measurements) and RNA-seq (read counts) data
    exhibit typically distinct features that need to be normalized for.
    For specific needs that deviate from these standard normalizations, the user
    should always refer to more specific functions/packages.        

    Microarray data is expected to be single-channel.
    For two-color arrays, it is expected here that normalization within arrays
    has been already carried out, e.g. using \code{\link{normalizeWithinArrays}}
    from limma.

    RNA-seq data is expected to be raw read counts.
    Please note that normalization for downstream DE analysis, e.g. with edgeR 
    and DESeq, is not ultimately necessary (and in some cases even discouraged)
    as many of these tools implement specific normalization approaches.
    See the vignette of EDASeq, edgeR, and DESeq for details.
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@bio.ifi.lmu.de>
}
\seealso{
    \code{\link{read.eset}} describes prerequisites for the expression data;

    \code{\link{normalizeWithinArrays}} and
    \code{\link{normalizeBetweenArrays}} for normalization of microarray data;

    \code{\link{withinLaneNormalization}} and
    \code{\link{betweenLaneNormalization}} for normalization of RNA-seq data.
}
\examples{
    #
    # (1) simulating expression data: 100 genes, 12 samples
    #
    
    # (a) microarray data: intensity measurements
    ma.eset <- make.example.data(what="eset", type="ma")
    
    # (b) RNA-seq data: read counts
    rseq.eset <- make.example.data(what="eset", type="rseq")

    #
    # (2) Normalization
    #
    
    # (a) microarray ... 
    norm.eset <- normalize(ma.eset) 

    # (b) RNA-seq ... 
    norm.eset <- normalize(rseq.eset) 

    # ... normalize also for GC content
    gc.content <- rnorm(100, 0.5, sd=0.1)
    fData(rseq.eset)$gc <- gc.content 

    norm.eset <- normalize(rseq.eset, within=TRUE)
}

