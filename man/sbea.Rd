\name{sbea}
\alias{sbea}
\alias{ora}
\alias{gsea}
\alias{sbea.methods}
\title{
    Set-based enrichment analysis (SBEA)
}
\description{
    This is the main function for the enrichment analysis of gene sets.
    It implements and uses existing implementations of several frequently 
    used state-of-art methods and allows a flexible inspection of 
    resulting gene set rankings. 
}
\usage{
    sbea(method=c("ora", "safe", "gsea", "samgs"), eset, gs, 
            alpha=0.05, perm=1000, out.file=NULL, browse=FALSE)

    sbea.methods()
}
\arguments{
    \item{method}{
        Set-based enrichment analysis method.
        Currently, the following set-based enrichment analysis methods 
        are supported:
        \sQuote{ora}, \sQuote{safe}, \sQuote{gsea}, and \sQuote{samgs}.
        See Details. For basic ora also set 'perm=0'. Default is \sQuote{ora}.
        This can also be the name of a tailored function implementing set-based
        enrichment. See Details. 
    } 
    \item{eset}{
        Expression set. Either an object of class 'ExpressionSet' or an absolute
        file path to an RData file containing the gene expression set. 
        See 'read.eset' and 'probe.2.gene.eset' for required annotations 
        in the pData and fData slot.
    }
    \item{gs}{
        Gene sets.
        Either a list of gene sets (vectors of KEGG gene IDs) or a text file 
        in GMT format storing all gene sets under investigation.
    }
    \item{alpha}{
        Statistical significance level. Defaults to 0.05.
    }
    \item{perm}{
        Number of permutations of the expression matrix to estimate the null 
        distribution. Defaults to 1000. For basic ora set 'perm=0'.
    }  
    \item{out.file}{
        Optional output file the gene set ranking will be written to.
    }
    \item{browse}{
        Logical. Should results be displayed in the browser for interactive 
        exploration? Defaults to FALSE.
    }
}
\value{
    sbea.methods: a character vector of currently supported methods;

    sbea: if(is.null(out.file)): an enrichment analysis result object 
    that can be detailedly explored by calling 'ea.browse' and from 
    which a flat gene set ranking can be extracted by calling 'gs.ranking'.
    If 'out.file' is given, the ranking is written to the specified file. 
}
\details{
    'ora': overrepresentation analysis, simple and frequently used test based 
    on the hypergeometric distribution (see Goeman and Buhlmann, 2007, for a 
    critical review). 'safe': significance analysis of function and expression,
    generalization of ORA, includes other test statistics, e.g. Wilcoxon's rank
    sum, and allows to estimate the significance of gene sets by sample 
    permutation; implemented in the safe package (Barry et al., 2005).
    'gsea': gene set enrichment analysis, frequently used and widely accepted, 
    uses a Kolmogorov-Smirnov statistic to test whether the ranks of the 
    p-values of genes in a gene set resemble a uniform distribution 
    (Subramanian et al., 2005). 'samgs': significance analysis of microarrays 
    on gene sets, extends the SAM method for single genes to gene set analysis
    (Dinu et al., 2007).

    It is also possible to use additional set-based enrichment methods.
    This requires to implement a function that takes 'eset', 'gs', 'alpha', and
    'perm' as arguments and returns a numeric vector 'ps' storing the resulting
    p-value for each gene set in 'gs'. This vector must be named accordingly 
    (i.e. names(ps) == names(gs)). See examples. 
}
\references{
    Goeman and Buhlmann (2007) 
    Analyzing gene expression data in terms of gene sets: methodological issues.
    Bioinformatics, 23, 980-7.

    Barry et al. (2005)
    Significance Analysis of Function and Expression.
    Bioinformatics, 21:1943-9.

    Subramanian et al. (2005) 
    Gene Set Enrichment Analysis: a knowledge-based approach for interpreting 
    genome-wide expression profiles.
    Proc Natl Acad Sci USA, 102:15545-50.

    Dinu et al. (2007)
    Improving gene set analysis of microarray data by SAM-GS.
    BMC Bioinformatics, 8:242
}
\seealso{
    Input:    
    \code{\link{read.eset}}, \code{\link{probe.2.gene.eset}}
      \code{\link{get.kegg.genesets}} to retrieve gene sets from KEGG.

    Output:
    \code{\link{gs.ranking}} to retrieve the ranked list of gene sets.  
    \code{\link{ea.browse}} for exploration of resulting gene sets. 

    Other:
    \code{\link{nbea}} to perform network-based enrichment analysis.
    \code{\link{comb.ea.results}} to combine results from different methods.
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@bio.ifi.lmu.de>
}
\examples{
    # currently supported methods
    sbea.methods()

    # (1) reading the expression data from file
    exprs.file <- system.file("extdata/ALL_exprs.tab", package="EnrichmentBrowser")
    pdat.file <- system.file("extdata/ALL_pData.tab", package="EnrichmentBrowser")
    fdat.file <- system.file("extdata/ALL_fData.tab", package="EnrichmentBrowser")
    probe.eset <- read.eset(exprs.file, pdat.file, fdat.file)

    # (2) summarizing probe expression on gene level
    gene.eset <- probe.2.gene.eset(probe.eset) 

    # (3) getting all human KEGG gene sets
    # hsa.gs <- get.kegg.genesets("hsa")
    gs.file <- system.file("extdata/hsa_kegg_gs.gmt", package="EnrichmentBrowser")
    hsa.gs <- parse.genesets.from.GMT(gs.file)

    # (4) performing the enrichment analysis
    ea.res <- sbea(method="ora", eset=gene.eset, gs=hsa.gs, perm=0)

    # (5) result visualization and exploration
    gs.ranking(ea.res)
    \donttest{
    ea.browse(ea.res)
    }

    # using your own tailored function as enrichment method
    dummy.sbea <- function(eset, gs, alpha, perm)
    {
        sig.ps <- sample(seq(0,0.05, length=1000),5)
        insig.ps <- sample(seq(0.1,1, length=1000), length(gs)-5)
        ps <- sample(c(sig.ps, insig.ps), length(gs))
        names(ps) <- names(gs)
        return(ps)
    }

    sbea.res2 <- sbea(method="dummy.sbea", eset=gene.eset, gs=hsa.gs)
    gs.ranking(sbea.res2) 
}

