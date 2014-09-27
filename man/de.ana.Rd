\name{de.ana}
\alias{de.ana}
\title{ 
    Differential expression analysis between two sample groups 
}
\description{
    The function carries out a differential expression analysis between two
    sample groups. Resulting fold changes and t-test derived p-values for each
    feature are appended to the fData slot.
    Raw p-values are corrected for multiple testing using the method from 
    Benjamini and Hochberg.
}
\usage{
    de.ana(eset, value.type=c("log2count", "log2ratio"))
}
\arguments{
    \item{eset}{
        An object of class 'ExpressionSet'    
    }
    \item{value.type}{
        Type of values in the expression matrix: 
        either 'log2count' or 'log2ratio' (defaults to 'log2count').   
    }
}
\value{
    An object of \code{\link[Biobase:ExpressionSet-class]{ExpressionSet-class}} 
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
    eset <- de.ana(eset)
    head(fData(eset))
}

