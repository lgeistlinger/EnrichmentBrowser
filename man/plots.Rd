\name{plots}
\alias{pdistr}
\alias{volcano}
\alias{exprs.heatmap}
\title{
    Visualization of gene expression
}
\description{
    Visualization of differential gene expression via 
    heatmap, p-value histogram and volcano plot (fold change vs. p-value).
}
\usage{
    pdistr( p )
    volcano( fc, p )
    exprs.heatmap( expr, grp )
}
\arguments{
    \item{p}{
        Numeric vector of p-values for each gene.
    }
    \item{fc}{
        Numeric vector of fold changes (typically on log2 scale).
    }
    \item{expr}{
        Expression matrix. Rows correspond to genes, columns to samples.    
    }
    \item{grp}{
        *BINARY* group assignment for the samples. 
        Use '0' and '1' for unaffected (controls) and affected (cases) 
        samples, respectively.
    }
}
\value{
    None, plots to a graphics device.
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@bio.ifi.lmu.de>
}
\seealso{
    \code{\link{de.ana}} for differential expression analysis,
    \code{\link{heatmap}} and \code{\link{truehist}} for generic plotting. 
}
\examples{
    # (1) simulating expression data: 100 genes, 12 samples
    eset <- make.example.data(what="eset") 
    
    # plot heatmap
    exprs.heatmap(expr=exprs(eset), grp=as.factor(pData(eset)$GROUP))

    # (2) DE analysis
    eset <- de.ana(eset)
    pdistr(fData(eset)$ADJ.PVAL)
    volcano(fc=fData(eset)$FC, p=fData(eset)$ADJ.PVAL)
}
