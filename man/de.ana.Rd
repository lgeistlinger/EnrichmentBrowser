\name{de.ana}
\alias{de.ana}
\title{ 
    Differential expression analysis between two sample groups 
}
\description{
    The function carries out a differential expression analysis between two
    sample groups. Resulting fold changes and derived p-values are returned.
    Raw p-values are corrected for multiple testing.
}
\usage{
    de.ana( expr, grp = NULL, blk = NULL, 
        de.method = c("limma", "edgeR", "DESeq"), padj.method = "BH", stat.only=FALSE, min.cpm=2 )
}
\arguments{
    \item{expr}{
        Expression data. 
        A numeric matrix. Rows correspond to genes, columns to samples.
        Alternatively, this can also be an object of class 
        \code{\linkS4class{ExpressionSet}} (in case of microarray data) 
        or an object of class \code{\linkS4class{SeqExpressionSet}} 
        (in case of RNA-seq data).
        See the man page of \code{\link{read.eset}} for prerequisites 
        for the expression data.
    }
    \item{grp}{
        *BINARY* group assignment for the samples. 
        Use '0' and '1' for unaffected (controls) and affected (cases) 
        samples, respectively.
        If NULL, this is assumed to be defined via a column named 'GROUP' 
        in the pData slot if 'expr' is a (Seq)ExpressionSet.
    }
    \item{blk}{
        Optional. For paired samples or sample blocks.
        This can also be defined via a column named 'BLOCK' 
        in the pData slot if 'expr' is a (Seq)ExpressionSet.
    }
    \item{de.method}{
        Differential expression method.
        Use 'limma' for microarray and RNA-seq data. 
        Alternatively, differential expression for RNA-seq data can be also
        calculated using edgeR ('edgeR') or DESeq2 ('DESeq').
        Defaults to 'limma'.
    }
    \item{padj.method}{
        Method for adjusting p-values to multiple testing.
        For available methods see the man of page the of the 
        stats function \code{\link{p.adjust}}.
        Defaults to 'BH'.
    }
    \item{stat.only}{
        Logical. Should only the test statistic be returned?
        This is mainly for internal use, in order to carry out
        permutation tests on the DE statistic for each gene.
        Defaults to FALSE.
    }
    \item{min.cpm}{
		In case of RNA-seq data: should genes not satisfying
		a minimum counts-per-million (cpm) threshold be excluded
		from the analysis? This is typically recommended. See the
		edgeR vignette for details. The default filter is to exlcude
		genes with cpm < 2 in more than half of the samples.  
	}
}
\value{ 
    A DE-table with measures of differential expression for each gene/row, 
    i.e. a two-column matrix with log2 fold changes in the 1st column 
    and derived p-values in the 2nd column.
    If 'expr' is a (Seq)ExpressionSet, the DE-table will be automatically
    appended to the fData slot. 
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@bio.ifi.lmu.de>
}
\seealso{
    \code{\link{read.eset}} describes prerequisites for the expression data,
    \code{\link{normalize}} for normalization of expression data,
    \code{\link{voom}} for preprocessing of RNA-seq data,
    \code{\link{p.adjust}} for multiple testing correction,
    \code{\link{eBayes}} for DE analysis with limma,
    \code{\link{glmFit}} for DE analysis with edgeR, and
    \code{\link{DESeq}} for DE analysis with DESeq. 
}
\examples{
    # (1) microarray data: intensity measurements
    ma.eset <- make.example.data(what="eset", type="ma")
    ma.eset <- de.ana(ma.eset)
    head(fData(ma.eset))
    
    # (2) RNA-seq data: read counts
    rseq.eset <- make.example.data(what="eset", type="rseq")
    rseq.eset <- de.ana(rseq.eset, de.method="DESeq")
    head(fData(rseq.eset))
}

