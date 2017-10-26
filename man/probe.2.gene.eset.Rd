\name{probe.2.gene.eset}
\alias{probe.2.gene.eset}
\title{
    Transformation of probe level expression to gene level expression 
}
\description{
    Transforms expression data on probe level to gene level expression 
    by summarizing all probes that are annotated to a particular gene. 
}
\usage{
    probe.2.gene.eset( probe.eset, use.mean = TRUE )
}
\arguments{
    \item{probe.eset}{
        Probe expression data of class \code{\linkS4class{SummarizedExperiment}}.
        The \code{rowData} slot of the expression dataset must contain 
        a 'GENE' column that lists for each probe the corresponding KEGG gene ID.
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
    A \code{\linkS4class{SummarizedExperiment}} on gene level.
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@sph.cuny.edu>
}
\seealso{
    \code{\link{read.eset}} for reading expression data from file,
    \code{\link{de.ana}} for differential expression analysis. 
}
\examples{
    # (1) reading the expression data from file
    exprs.file <- system.file("extdata/exprs.tab", package="EnrichmentBrowser")
    pdat.file <- system.file("extdata/colData.tab", package="EnrichmentBrowser")
    fdat.file <- system.file("extdata/rowData.tab", package="EnrichmentBrowser")
    probe.eset <- read.eset(exprs.file, pdat.file, fdat.file)
    gene.eset <- probe.2.gene.eset(probe.eset) 
}
