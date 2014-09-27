\name{plots}
\alias{pdistr}
\alias{volcano}
\alias{exprs.heatmap}
\title{
    Visualization of gene expression
}
\description{
    Visualization of differential gene expression via 
    heatmaps, p-value histograms and volcano plots (fold change vs. p-value).
}
\usage{
    pdistr(eset, out.file=NULL, use.adjp=TRUE)
    volcano(eset, out.file=NULL)
    exprs.heatmap(eset, out.file=NULL)
}
\arguments{
    \item{eset}{
        Expression set. An object of class 'ExpressionSet' containing the 
        gene expression set. See 'read.eset' and 'probe.2.gene.eset' for 
        required annotations in the pData and fData slot.
    }
    \item{use.adjp}{
        Logical. Use differential expression p-values that are already corrected
        for multiple testing? Defaults to TRUE.
    }
    \item{out.file}{
        Optional output file to which graphics are plotted to.
    }
}
\value{
    None, plots to a graphics device.
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@bio.ifi.lmu.de>
}
\seealso{
    \code{\link{heatmap}}, \code{\link{truehist}} 
}
\examples{
    # (1) reading the expression data from file
    exprs.file <- system.file("extdata/ALL_exprs.tab", package="EnrichmentBrowser")
    pdat.file <- system.file("extdata/ALL_pData.tab", package="EnrichmentBrowser")
    fdat.file <- system.file("extdata/ALL_fData.tab", package="EnrichmentBrowser")
    eset <- read.eset(exprs.file, pdat.file, fdat.file)

    # (2) plot
    exprs.heatmap(eset)
    pdistr(eset)
    volcano(eset)
}
