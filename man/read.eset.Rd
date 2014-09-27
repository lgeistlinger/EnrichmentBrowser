\name{read.eset}
\alias{read.eset}
\title{ 
    Reading gene expression data from file into an expression set 
}
\description{
    The function reads in plain expression data with minimum annotation 
    requirements for the pData and fData slots. A simultaneous primary analysis
    of differential expression of single genes is performed and annotated to 
    the expression set.
    If desired, the distribution of the t-test derived p-values and the 
    relation between them and the fold change, can be plotted in a histogram
    and a volcano plot, respectively. 
}
\usage{
    read.eset(exprs.file, pdat.file, fdat.file, 
        value.type=c("log2count", "log2ratio"), out.file=NULL, 
        heatm.file=NULL, distr.file = NULL, volc.file = NULL)
}
\arguments{
    \item{exprs.file}{
        Expression matrix.
        A tab separated text file containing the numerical expression values.
        Columns = samples/subjects; rows = features/probes/genes; 
        NO headers, row or column names. 
    }
    \item{pdat.file}{
        Phenotype data.
        A tab separated text file containing annotation information 
        for the samples. Exactly *TWO* cols; 1st col = sample names; 
        2nd column = *BINARY* group assignment. Use '0' and '1' for unaffected 
        and affected sample class, respectively. 
        NO headers, row or column names.
        The number of rows/samples in this file should match the number of 
        columns/samples of the expression matrix. 
    }
    \item{fdat.file}{
        Feature data.
        A tab separated text file containing annotation information 
        for the features.
        Exactly *TWO* columns; 1st col = feature IDs; 
        2nd col = corresponding KEGG gene ID for each feature ID in 1st col; 
        NO headers, row or column names.    
        The number of rows/features in this file should match the number of 
        rows/features of the expression matrix. 
    }
    \item{value.type}{
        Type of values in the expression matrix: 
        either 'log2count' or 'log2ratio' (defaults to 'log2count').   
    }
    \item{out.file}{
        Optional. If specified, the processed expression set is saved to file 
        (using \code{\link{save}}).
    }
    \item{heatm.file}{
        Optional. If specified, a heatmap is plotted in PNG format to file.
    }
    \item{distr.file}{
        Optional. If specified, the p-value distribution is plotted 
        in PNG format to file.
    }
    \item{volc.file}{
        Optional. If specified, the volcano plot is plotted in 
        PNG format to file.
    }
}
\value{
    An \code{\link[Biobase:ExpressionSet-class]{ExpressionSet-class}} 
    with measures of differential expression annotated in the fData slot.
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@bio.ifi.lmu.de>
}

\seealso{
    \code{\link[Biobase:ExpressionSet-class]{ExpressionSet-class}}, 
    \code{\link{get.fold.change.and.t.test}}, 
    \code{\link{p.adjust}}
}
\examples{
    # reading the expression data from file
    exprs.file <- system.file("extdata/ALL_exprs.tab", package="EnrichmentBrowser")
    pdat.file <- system.file("extdata/ALL_pData.tab", package="EnrichmentBrowser")
    fdat.file <- system.file("extdata/ALL_fData.tab", package="EnrichmentBrowser")
    eset <- read.eset(exprs.file, pdat.file, fdat.file)
    head(fData(eset))
}

