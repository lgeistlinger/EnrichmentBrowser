\name{nbea}
\alias{nbea}
\alias{nbea.methods}
\alias{ggea}
\title{
    Network-based enrichment analysis (NBEA)
}
\description{
    This is the main function for network-based enrichment analysis.
    It implements and wraps existing implementations of several frequently 
    used state-of-art methods and allows a flexible inspection of 
    resulting gene set rankings. 
}
\usage{
    nbea( method = EnrichmentBrowser::nbea.methods(), eset, gs, grn, 
            prune.grn=TRUE, alpha = 0.05, perm = 1000, padj.method = "none", 
            out.file = NULL, browse = FALSE, ... )

    nbea.methods()
}
\arguments{
    \item{method}{
        Network-based enrichment analysis method.
        Currently, the following network-based enrichment analysis methods 
        are supported:
        \sQuote{ggea}, \sQuote{spia}, \sQuote{pathnet}, \sQuote{degraph},
		\sQuote{topologygsa}, \sQuote{ganpa}, \sQuote{cepa}, \sQuote{netgsa},
		and \sQuote{nea}. 
        Default is 'ggea'.
        This can also be the name of a user-defined function implementing 
        network-based enrichment. See Details. 
    } 
    \item{eset}{
        Expression set. 
        An object of class \code{\linkS4class{ExpressionSet}}.
        See \code{\link{read.eset}} and \code{\link{probe.2.gene.eset}}
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
    \item{prune.grn}{
        Logical. 
        Should the GRN be pruned? This removes duplicated, self, and
        reversed edges. Defaults to TRUE.
    }
    \item{alpha}{
        Statistical significance level. Defaults to 0.05.
    }
    \item{perm}{
        Number of permutations of the expression matrix to estimate the 
        null distribution. Defaults to 1000. If using method=\sQuote{ggea},
        it is possible to set 'perm=0' to use a fast approximation of
        gene set significance to avoid permutation testing. See Details.
    } 
    \item{padj.method}{
        Method for adjusting nominal gene set p-values to multiple testing.
        For available methods see the man of page the of the 
        stats function \code{\link{p.adjust}}.
        Defaults to 'none', i.e. leaves the nominal gene set p-values unadjusted.
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
        This includes currently:
        \itemize{
            \item beta: Log2 fold change significance level. Defaults to 1 (2-fold).
        }
        For SPIA and NEA:
        \itemize{
            \item sig.stat: decides which statistic is used for determining 
                significant DE genes. 
                Options are:
                    \itemize{
                        \item 'p' (Default): genes with p-value below alpha.
                        \item 'fc': genes with abs(log2(fold change)) above beta
                        \item '&': p & fc (logical AND)
                        \item '|': p | fc (logical OR)
                    } 
        }
        For GGEA:
        \itemize{
        \item cons.thresh: edge consistency threshold between -1 and 1.
            Defaults to 0.2, i.e. only edges of the GRN with consistency >= 0.2 
            are included in the analysis. Evaluation on real datasets has shown
            that this works best to distinguish relevant gene sets. 
            Use consistency of -1 to include all edges.
        \item gs.edges: decides which edges of the grn are considered for a
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
    that can be detailedly explored by calling \code{\link{ea.browse}} and from 
    which a flat gene set ranking can be extracted by calling \code{\link{gs.ranking}}.
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
    approximation in the spirit of Bioconductor's npGSEA package that is much faster.  
    
    'spia': signaling pathway impact analysis, combines ORA with the probability that 
	expression changes are propagated across the pathway topology;    
	implemented in Bioconductor's SPIA package.

    'pathnet': pathway analysis using network information, applies ORA on combined 
	evidence for the observed signal for gene nodes and the signal implied by connected 
	neighbors in the network; implemented in Bioconductor's PathNet package.

	'degraph': differential expression testing for gene graphs,
    multivariate testing of differences in mean incorporating underlying
    graph structure; implemented in Bioconductor's DEGraph package

	'topologygsa': topology-based gene set analysis,
    uses Gaussian graphical models to incorporate the dependence structure
	among genes as implied by pathway topology;
    implemented in CRAN's topologyGSA package.

	'ganpa': gene association network-based pathway analysis,
    incorporates network-derived gene weights in the enrichment analysis;
    implemented in CRAN's GANPA package.

	'cepa': centrality-based pathway enrichment,
    incorporates network centralities as node weights mapped from
    differentially expressed genes in pathways;
    implemented in CRAN's CePa package.

	'netgsa': network-based gene set analysis,
    incorporates external information about interactions among genes
    as well as novel interactions learned from data;
    implemented in CRAN's NetGSA package.

	'nea': network enrichment analysis, applies ORA on interactions instead of genes; 
	implemented in Bioconductor's neaGUI package. 

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
    \code{\link{read.eset}}, \code{\link{probe.2.gene.eset}}, 
    \code{\link{get.kegg.genesets}} to retrieve gene set definitions from KEGG.
    \code{\link{compile.grn.from.kegg}} to construct a GRN from KEGG pathways.
  
    Output:
    \code{\link{gs.ranking}} to rank the list of gene sets.  
    \code{\link{ea.browse}} for exploration of resulting gene sets. 

    Other:
    \code{\link{sbea}} to perform set-based enrichment analysis.
    \code{\link{comb.ea.results}} to combine results from different methods.
    the SPIA package for more information on signaling pathway impact analysis.
    the neaGUI package for more information on network enrichment analysis.    
    the PathNet package for more information on pathway analysis using network information.
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@bio.ifi.lmu.de>
}
\examples{
    # currently supported methods
    nbea.methods()

    # (1) expression data: 
    # simulated expression values of 100 genes
    # in two sample groups of 6 samples each
    eset <- make.example.data(what="eset")
    eset <- de.ana(eset)

    # (2) gene sets:
    # draw 10 gene sets with 15-25 genes
    gs <- make.example.data(what="gs", gnames=featureNames(eset))

    # (3) make 2 artificially enriched sets:
    sig.genes <- featureNames(eset)[fData(eset)$ADJ.PVAL < 0.1]
    gs[[1]] <- sample(sig.genes, length(gs[[1]])) 
    gs[[2]] <- sample(sig.genes, length(gs[[2]]))   
   
    # (4) gene regulatory network 
    grn <- make.example.data(what="grn", nodes=featureNames(eset))
    
    # (5) performing the enrichment analysis
    ea.res <- nbea(method="ggea", eset=eset, gs=gs, grn=grn)

    # (6) result visualization and exploration
    gs.ranking(ea.res, signif.only=FALSE)

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

    ea.res2 <- nbea(method=dummy.nbea, eset=eset, gs=gs, grn=grn)
    gs.ranking(ea.res2) 
}

