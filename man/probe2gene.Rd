\name{probe2gene}
\alias{probe2gene}
\alias{probe.2.gene.eset}
\title{
    Transformation of probe level expression to gene level expression 
}
\description{
    Transforms expression data on probe level to gene level expression 
    by summarizing all probes that are annotated to a particular gene. 
}
\usage{
    probe2gene( probeSE, use.mean = TRUE )
}
\arguments{
    \item{probeSE}{
        Probe expression data. 
        An object of class \code{\linkS4class{SummarizedExperiment}}.
        Make sure that the \code{\link{metadata}} contains an element named
        \code{annotation} that provides the corresponding ID of a 
        recognized platform such as \code{hgu95av2} (Affymetrix Human Genome 
        U95 chip).
        This requires that a corresponding \code{.db} package exists (see 
        \url{http://www.bioconductor.org/packages/release/BiocViews.html#___ChipName}
        for available chips/packages) and that you have it installed.
        Alternatively, The mapping from probe to gene can also be defined in the
        \code{\link{rowData}} slot via two columns named (i) \code{PROBEID} for the 
        platform-specific probe ID, and (ii) \code{ENTREZID} for the corresponding
        NCBI Entrez Gene ID.
    }
    \item{use.mean}{
        Logical.
        Determining, in case of multiple probes for one gene,
        whether a mean value is computed (\code{use.mean=TRUE}), or the probe
        that discriminate the most between the two sample
        group is kept (\code{use.mean=FALSE}).
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
    \code{\link{readSE}} for reading expression data from file,
    \code{\link{deAna}} for differential expression analysis. 
}
\examples{
    # (1) reading the expression data from file
    exprs.file <- system.file("extdata/exprs.tab", package="EnrichmentBrowser")
    pdat.file <- system.file("extdata/colData.tab", package="EnrichmentBrowser")
    fdat.file <- system.file("extdata/rowData.tab", package="EnrichmentBrowser")
    probeSE <- readSE(exprs.file, pdat.file, fdat.file)
    geneSE <- probe2gene(probeSE) 
}
