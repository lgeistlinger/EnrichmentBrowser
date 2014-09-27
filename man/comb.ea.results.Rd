\name{comb.ea.results}
\alias{comb.ea.results}
\title{
    Combining enrichment analysis results
}
\description{
    Different enrichment analysis methods usually result in 
    different gene set rankings for the same dataset.
    This function allows to combine results from the different 
    set-based and network-based enrichment analysis methods.
    This includes the computation of average gene set ranks 
    across methods and of combined gene set p-values. 
}
\usage{
    comb.ea.results(res.list, pcomb.meth=c("fisher","stouffer"), out.file=NULL)
}
\arguments{
    \item{res.list}{
        A list of enrichment analysis result objects as returned 
        by the functions \sQuote{sbea} and \sQuote{nbea}.
    } 
    \item{pcomb.meth}{
        P-value combination method.
        See details. Defaults to \sQuote{fisher}.
    }
    \item{out.file}{
        Optional output file the combined ranking is written to.
    }
}
\value{
    if(is.null(out.file)): an enrichment analysis result object 
    that can be detailedly explored by calling 'ea.browse' and from 
    which a flat gene set ranking can be extracted by calling 'gs.ranking'.
    If 'out.file' is given, the ranking is written to the specified file. 
}
\details{
    Fisher's method (pcomb.meth="fisher") combines the statistical significance 
    from several independent tests.
    Stouffer's method (pcomb.meth="stouffer") allows the combination from 
    dependent tests. See Kim et al., 2013.
}
\references{
    Kim et al. (2013)
    Stouffer's test in a large scale simultaneous hypothesis testing
    PLoS One, 8(5), e63290.
}
\seealso{
    \code{\link{sbea}}, \code{\link{nbea}}, \code{\link{ea.browse}}}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@bio.ifi.lmu.de>
}
\examples{
    # (1) reading the expression data from file
    exprs.file <- system.file("extdata/ALL_exprs.tab", package="EnrichmentBrowser")
    pdat.file <- system.file("extdata/ALL_pData.tab", package="EnrichmentBrowser")
    fdat.file <- system.file("extdata/ALL_fData.tab", package="EnrichmentBrowser")
    probe.eset <- read.eset(exprs.file, pdat.file, fdat.file)

    # (2) summarizing probe expression on gene level
    gene.eset <- probe.2.gene.eset(probe.eset) 

    # (3a) getting all human KEGG gene sets
    # hsa.gs <- get.kegg.genesets("hsa")
    gs.file <- system.file("extdata/hsa_kegg_gs.gmt", package="EnrichmentBrowser")
    hsa.gs <- parse.genesets.from.GMT(gs.file)

    # (3b) compiling gene regulatory network from KEGG pathways
    # hsa.grn <- compile.grn.from.kegg("hsa")
    pwys <- system.file("extdata/hsa_kegg_pwys.zip", package="EnrichmentBrowser")
    hsa.grn <- compile.grn.from.kegg(pwys)

    # (4) performing the set-based and network-based enrichment analysis
    # Note: reduced permutations for demonstration
    #       recommended default is 1000 permutations
    # sbea.res <- sbea(method="gsea", eset=gene.eset, gs=hsa.gs)
    # nbea.res <- nbea(method="ggea", eset=gene.eset, gs=hsa.gs, grn=hsa.grn)
    sbea.res <- sbea(method="ora", eset=gene.eset, gs=hsa.gs, perm=0)
    nbea.res <- nbea(method="ggea", 
                        eset=gene.eset, gs=hsa.gs, grn=hsa.grn, perm=100)
    
    # (5) combining the results
    res.list <- list(sbea.res, nbea.res)
    comb.res <- comb.ea.results(res.list)

    # (6) result visualization and exploration
    gs.ranking(comb.res)
    \donttest{
    ea.browse(comb.res)
    }
}

