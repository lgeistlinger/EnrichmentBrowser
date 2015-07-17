\name{probe.2.gene.eset}
\alias{probe.2.gene.eset}
\title{
    Transformation of probe level expression to gene level expression 
}
\description{
    Reads expression data at probe level and summarizes gene expression 
    behavior by averaging over all probes that are annotated to a particular 
    gene. 
}
\usage{
    probe.2.gene.eset( probe.eset, use.mean = TRUE )
}
\arguments{
    \item{probe.eset}{
        Probe expression set of class \code{\linkS4class{ExpressionSet}}.
        The fData slot of the expression set must contain a 'GENE' column 
        that lists for each probe the corresponding KEGG gene ID.
    }
    \item{use.mean}{
        Logical.
        Determining, in case of multiple probes for one gene,
        whether a mean value is computed (use.mean=TRUE), or the probe
        that discriminate the most between the two sample
        group is kept (use.mean=FALSE).
        Defaults to TRUE.
    }
}
\value{
    An \code{\linkS4class{ExpressionSet}} on gene level.
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@bio.ifi.lmu.de>
}
\seealso{
    \code{\link[Biobase:ExpressionSet-class]{ExpressionSet-class}}, 
    \code{\link{read.eset}} for reading expression data from file,
    \code{\link{de.ana}} for differential expression analysis. 
}
\examples{
    # (1) reading the expression data from file
    exprs.file <- system.file("extdata/exprs.tab", package="EnrichmentBrowser")
    pdat.file <- system.file("extdata/pData.tab", package="EnrichmentBrowser")
    fdat.file <- system.file("extdata/fData.tab", package="EnrichmentBrowser")
    probe.eset <- read.eset(exprs.file, pdat.file, fdat.file)
    gene.eset <- probe.2.gene.eset(probe.eset) 
}
