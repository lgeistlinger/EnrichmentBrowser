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
    It implements and wraps existing implementations of several frequently 
    used state-of-art methods and allows a flexible inspection of 
    resulting gene set rankings. 
}
\usage{
    sbea( method = EnrichmentBrowser::sbea.methods(), eset, gs, alpha = 0.05, 
        perm = 1000, padj.method = "none", out.file = NULL, browse = FALSE, ... )

    sbea.methods()
}
\arguments{
    \item{method}{
        Set-based enrichment analysis method.
        Currently, the following set-based enrichment analysis methods 
        are supported:
        \sQuote{ora}, \sQuote{safe}, \sQuote{gsea}, \sQuote{padog},     
		\sQuote{roast}, \sQuote{camera}, \sQuote{gsa}, \sQuote{gsva},
		\sQuote{globaltest}, \sQuote{samgs}, \sQuote{ebm}, and \sQuote{mgsa}.
        For basic ora also set 'perm=0'. Default is \sQuote{ora}.
        This can also be the name of a user-defined function implementing 
        set-based enrichment. See Details. 
    } 
    \item{eset}{
        Expression set. An object of class \code{\linkS4class{ExpressionSet}}. 
        See \code{\link{read.eset}} and \code{\link{probe.2.gene.eset}} for 
        required annotations in the pData and fData slots.
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
        Using method="gsea" and 'perm=0' invokes the permutation
        approximation from the npGSEA package.
    }  
    \item{padj.method}{
        Method for adjusting nominal gene set p-values to multiple testing.
        For available methods see the man page of the 
        stats function \code{\link{p.adjust}}.
        Defaults to'none', i.e. leaves the nominal gene set p-values unadjusted.
    }
    \item{out.file}{
        Optional output file the gene set ranking will be written to.
    }
    \item{browse}{
        Logical. Should results be displayed in the browser for interactive 
        exploration? Defaults to FALSE.
    }
    \item{...}{
        Additional arguments passed to individual sbea methods.
        This includes currently for ORA and MGSA:
        \itemize{
            \item beta: Log2 fold change significance level. Defaults to 1 (2-fold).
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
    }
}
\value{
    sbea.methods: a character vector of currently supported methods;

    sbea: if(is.null(out.file)): an enrichment analysis result object 
    that can be detailedly explored by calling \code{\link{ea.browse}} and from 
    which a flat gene set ranking can be extracted by calling \code{\link{gs.ranking}}.
    If 'out.file' is given, the ranking is written to the specified file. 
}
\details{
    'ora': overrepresentation analysis, simple and frequently used test based 
    on the hypergeometric distribution (see Goeman and Buhlmann, 2007, for a 
    critical review). 

    'safe': significance analysis of function and expression,
    generalization of ORA, includes other test statistics, e.g. Wilcoxon's rank
    sum, and allows to estimate the significance of gene sets by sample 
    permutation; implemented in the safe package (Barry et al., 2005).
    
    'gsea': gene set enrichment analysis, frequently used and widely accepted, 
    uses a Kolmogorov-Smirnov statistic to test whether the ranks of the 
    p-values of genes in a gene set resemble a uniform distribution 
    (Subramanian et al., 2005). 

	'padog': pathway analysis with down-weighting of overlapping genes, 
    incorporates gene weights to favor genes appearing in few pathways versus genes 
    that appear in many pathways; implemented in the PADOG package.

	'roast': rotation gene set test, 
    uses rotation instead of permutation for assessment of gene set significance; 
    implemented in the limma and edgeR packages for microarray and RNA-seq data, respectively.

	'camera': correlation adjusted mean rank gene set test,
    accounts for inter-gene correlations as implemented in the limma 
    and edgeR packages for microarray and RNA-seq data, respectively.

	'gsa': gene set analysis,
    differs from GSEA by using the maxmean statistic, i.e. the mean of the positive or  
    negative part of gene scores in the gene set; implemented in the GSA package.

	'gsva': gene set variation analysis,
    transforms the data from a gene by sample matrix to a gene set by sample matrix, 
    thereby allowing the evaluation of gene set enrichment for each sample;
    implemented in the GSVA package.

	'globaltest': global testing of groups of genes,
    general test of groups of genes for association with a response variable;
    implemented in the globaltest package.

    'samgs': significance analysis of microarrays 
    on gene sets, extends the SAM method for single genes to gene set analysis
    (Dinu et al., 2007).

    'ebm': empirical Brown's method, combines $p$-values of genes in a 
    gene set using Brown's method to combine $p$-values from dependent tests; 
    implemented in the EmpiricalBrownsMethod package. 

	'mgsa': model-based gene set analysis,
    Bayesian modeling approach taking set overlap into account by working on all sets 
    simultaneously, thereby reducing the number of redundant sets; 
	implemented in the mgsa package.

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

    # (4) performing the enrichment analysis
    ea.res <- sbea(method="ora", eset=eset, gs=gs, perm=0)

    # (5) result visualization and exploration
    gs.ranking(ea.res)

    # using your own tailored function as enrichment method
    dummy.sbea <- function(eset, gs, alpha, perm)
    {
        sig.ps <- sample(seq(0, 0.05, length=1000), 5)
        nsig.ps <- sample(seq(0.1, 1, length=1000), length(gs)-5)
        ps <- sample(c(sig.ps, nsig.ps), length(gs))
        names(ps) <- names(gs)
        return(ps)
    }

    ea.res2 <- sbea(method=dummy.sbea, eset=eset, gs=gs)
    gs.ranking(ea.res2) 
}

