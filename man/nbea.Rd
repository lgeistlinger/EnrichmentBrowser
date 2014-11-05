\name{nbea}
\alias{nbea}
\alias{nbea.methods}
\alias{ggea}
\title{
    Network-based enrichment analysis (NBEA)
}
\description{
    This is the main function for network-based enrichment analysis.
    It implements and uses existing implementations of several frequently 
    used state-of-art methods and allows a flexible inspection of 
    resulting gene set rankings. 
}
\usage{
    nbea(method=c("ggea", "nea", "spia"), eset, gs, grn, 
                    alpha=0.05, perm=1000, out.file=NULL, browse=FALSE, ...)

    nbea.methods()
}
\arguments{
    \item{method}{
        Network-based enrichment analysis method.
        Currently, the following network-based enrichment analysis methods 
        are supported:
        \sQuote{ggea}, \sQuote{nea}, \sQuote{spia}.
        See Details. Default is 'ggea'.
        This can also be the name of a tailored function implementing 
        network-based enrichment. See Details. 
    } 
    \item{eset}{
        Expression set. Either an object of class 'ExpressionSet' or an 
        absolute file path to an RData file containing the gene expression set. 
        See 'read.eset' and 'probe.2.gene.eset'
        for required annotations in the pData and fData slot.
    }
    \item{gs}{
        Gene sets.
        Either a list of gene sets (vectors of KEGG gene IDs) or 
        a text file in GMT format storing all gene sets under investigation.
    }
    \item{grn}{
        Gene regulatory network.
        Either an absolute file path to a tabular file or a character matrix 
        with exactly *THREE* cols; 1st col = IDs of regulating genes; 
        2nd col = corresponding regulated genes;
        3rd col = regulation effect; Use '+' and '-' for activation/inhibition.
    }
    \item{alpha}{
        Statistical significance level. Defaults to 0.05.
    }
    \item{perm}{
        Number of permutations of the expression matrix to estimate the 
        null distribution. Defaults to 1000. If using method=\sQuote{ggea},
        it is possible to set perm < 1 to use a fast approximation of
        gene set significance to avoid permutation testing. See Details.
    }  
    \item{out.file}{
        Optional output file the gene set ranking will be written to.
    }
    \item{browse}{
        Logical. Should results be displayed in the browser for interactive 
        exploration? Defaults to FALSE.
    }
    \item{...}{
        Additional arguments passed to individual nbea methods.
        This includes currently for GGEA:
        \itemize{
        \item beta: Log2 fold change significance level. Defaults to 1 (2-fold).
        \item cons.thresh: consistency threshold. Defaults to -1.
        \item gs.edges: Decides which edges of the grn are considered for a
                gene set under investigation. Should be one out of c('&', '|'), 
                denoting logical AND and OR. respectively. Accordingly, this 
                either includes edges for which regulator AND / OR target gene
                are members of the investigated gene set.            
        }
    }
}
\value{
    nbea.methods: a character vector of currently supported methods;

    nbea: if(is.null(out.file)): an enrichment analysis result object 
    that can be detailedly explored by calling 'ea.browse' and from 
    which a flat gene set ranking can be extracted by calling 'gs.ranking'.
    If 'out.file' is given, the ranking is written to the specified file. 
}
\details{
    'ggea': gene graph enrichment analysis, scores gene sets according to 
    consistency within the given gene regulatory network, i.e. checks 
    activating regulations for positive correlation and repressing regulations 
    for negative correlation of regulator and target gene expression 
    (Geistlinger et al., 2011). When using 'ggea' it is possible to estimate
    the statistical significance of the consistency score of each gene set in 
    two different ways: (1) based on sample permutation as described in the 
    original publication (Geistlinger et al., 2011) or (2) using an 
    approximation based on Bioconductor's npGSEA package that is much faster.  
    
    'nea': network enrichment analysis, implemented in Bioconductor's neaGUI 
    package. 
    
    'spia': signaling pathway impact analysis, implemented in Bioconductor's 
    SPIA package.

    It is also possible to use additional network-based enrichment methods.
    This requires to implement a function that takes 'eset', 'gs', 'grn', 
    'alpha', and 'perm' as arguments and returns a numeric matrix 'res.tbl' 
    with a mandatory column named 'P.VALUE' storing the resulting p-value for 
    each gene set in 'gs'. The rows of this matrix must be named 
    accordingly (i.e. rownames(res.tbl) == names(gs)). See examples. 

}
\references{
    Geistlinger et al. (2011) 
    From sets to graphs: towards a realistic enrichment analysis of 
    transcriptomic systems. 
    Bioinformatics, 27(13), i366--73.
}
\seealso{
    Input:    
    \code{\link{read.eset}}, \code{\link{probe.2.gene.eset}}
    \code{\link{get.kegg.genesets}} to retrieve gene set definitions from KEGG.
    \code{\link{compile.grn.from.kegg}} to construct a GRN from KEGG pathways.
  
    Output:
    \code{\link{gs.ranking}} to rank the list of gene sets.  
    \code{\link{ea.browse}} for exploration of resulting gene sets. 

    Other:
    \code{\link{sbea}} to perform set-based enrichment analysis.
    \code{\link{comb.ea.results}} to combine results from different methods.
    \code{\link{spia}} for more information on signaling pathway impact analysis.
    \code{\link{nea}} for more information on network enrichment analysis.    
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@bio.ifi.lmu.de>
}
\examples{
    # currently supported methods
    nbea.methods()

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

    # (4) performing the enrichment analysis
    # Note: reduced permutations for demonstration
    #       recommended default is 1000 permutations
    # ea.res <- nbea(method="ggea", eset=gene.eset, gs=hsa.gs, grn=hsa.grn)
    ea.res <- nbea(method="ggea", 
                    eset=gene.eset, gs=hsa.gs, grn=hsa.grn, perm=100)

    # (5) result visualization and exploration
    gs.ranking(ea.res)
    \donttest{
    ea.browse(ea.res, graph.view=hsa.grn)
    }

    # using your own tailored function as enrichment method
    dummy.nbea <- function(eset, gs, grn, alpha, perm)
    {
        sig.ps <- sample(seq(0,0.05, length=1000),5)
        insig.ps <- sample(seq(0.1,1, length=1000), length(gs)-5)
        ps <- sample(c(sig.ps, insig.ps), length(gs))
        score <- sample(1:100, length(gs), replace=TRUE)
        res.tbl <- cbind(score, ps)
        colnames(res.tbl) <- c("SCORE", "P.VALUE")
        rownames(res.tbl) <- names(gs)
        return(res.tbl[order(ps),])
    }

    nbea.res2 <- nbea(method="dummy.nbea", 
        eset=gene.eset, gs=hsa.gs, grn=hsa.grn)
    gs.ranking(nbea.res2) 
}

