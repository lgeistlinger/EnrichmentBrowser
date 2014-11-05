\name{probe.2.gene.eset}
\alias{probe.2.gene.eset}
\title{
    Transformation of probe level expression to gene level expression 
}
\description{
    Reads expression data at probe level and summarizes gene expression 
    behavior by averaging over all probes that are annotated to a particular 
    gene. A subsequent differential expression analysis at the gene level 
    is automatically performed.
}
\usage{
    probe.2.gene.eset(probe.eset, use.mean=TRUE,  
        heatm.file=NULL, distr.file=NULL, volc.file=NULL)
}
\arguments{
    \item{probe.eset}{
        Probe expression set of class 'ExpressionSet'.         
        The fData slot of the expression set must contain a 'GENE' column 
        that lists for each probe the corresponding KEGG gene ID.
    }
    \item{use.mean}{
        Logical.
        Determining, in case of multiple probes for one gene,
        whether a mean value is computed or the probe
        that discriminate the most between the two sample
        group is kept.
        Defaults to TRUE.
    }
    \item{heatm.file}{
        Optional. If specified, a heatmap is plotted in PNG format to file.
    }
    \item{distr.file}{
        Optional. If specified, the p-value distribution is plotted in 
        PNG format to file.
    }
    \item{volc.file}{
        Optional. If specified, the volcano plot is plotted in PNG format 
        to file.
    }
}
\value{
    An \code{\link[Biobase:ExpressionSet-class]{ExpressionSet-class}} 
    on gene level with measures of differential expression annotated.
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@bio.ifi.lmu.de>
}
\seealso{
    \code{\link[Biobase:ExpressionSet-class]{ExpressionSet-class}}, 
    \code{\link{eBayes}}, 
    \code{\link{p.adjust}}
}
\examples{
    # (1) reading the expression data from file
    exprs.file <- system.file("extdata/ALL_exprs.tab", package="EnrichmentBrowser")
    pdat.file <- system.file("extdata/ALL_pData.tab", package="EnrichmentBrowser")
    fdat.file <- system.file("extdata/ALL_fData.tab", package="EnrichmentBrowser")
    probe.eset <- read.eset(exprs.file, pdat.file, fdat.file)
    gene.eset <- probe.2.gene.eset(probe.eset) 
    head(fData(gene.eset))
}
