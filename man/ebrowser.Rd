\name{ebrowser}
\alias{ebrowser}
\title{
    Seamless navigation through enrichment analysis results 
}
\description{
    This is the all-in-one wrapper function to perform 
    the standard enrichment analysis pipeline implemented 
    in the EnrichmentBrowser package.
    
    Given flat gene expression data, the 
    data is read in and subsequently
    subjected to chosen enrichment analysis methods.

    The results from different methods can be combined and
    investigated in detail in the default browser. 
}
\usage{
    ebrowser( meth, exprs, pdat, fdat, org, data.type = c(NA, "ma", "rseq"), 
            norm.method = "quantile", de.method = "limma",
            gs, grn = NULL, perm = 1000, alpha = 0.05, beta = 1,
            comb = FALSE, browse = TRUE, nr.show = -1 )
}
\arguments{
    \item{meth}{
        Enrichment analysis method.
        Currently, the following enrichment analysis methods are supported:
        \sQuote{ora}, \sQuote{safe}, \sQuote{gsea}, \sQuote{samgs}, 
        \sQuote{ggea}, \sQuote{spia}, \sQuote{nea}, and \sQuote{pathnet}.
        See \code{\link{sbea}} and \code{\link{nbea}} for details.
    }
    \item{exprs}{
      Expression matrix.
        A tab separated text file containing *normalized* expression values 
        on a *log* scale.
        Columns = samples/subjects; rows = features/probes/genes; 
        NO headers, row or column names.
        Supported data types are log2 counts (microarray single-channel),
        log2 ratios (microarray two-color), and log2-counts per million 
        (RNA-seq logCPMs).
        See limma's user guide for definition and normalization of the 
        different data types.       
        Alternatively, this can be an object of 
        \code{\link[Biobase:ExpressionSet-class]{ExpressionSet-class}},
        assuming the expression matrix in the 'exprs' slot.
    }
    \item{pdat}{
        Phenotype data.
        A tab separated text file containing annotation information 
        for the samples in either *two or three* columns. 
        NO headers, row or column names.
        The number of rows/samples in this file should match the number of 
        columns/samples of the expression matrix.
        The 1st column is reserved for the sample IDs; 
        The 2nd column is reserved for a *BINARY* group assignment. 
        Use '0' and '1' for unaffected (controls) and affected (cases) sample 
        class, respectively. 
        For paired samples or sample blocks a third column is expected that
        defines the blocks. 
        If 'exprs' is an object of 
        \code{\link[Biobase:ExpressionSet-class]{ExpressionSet-class}},
        the 'pdat' argument can be left unspecified, which then expects
        group and optional block assignments in respectively named columns
        'GROUP' (mandatory) and 'BLOCK' (optional) in the 'pData' slot of the 
        ExpressionSet.
    }
    \item{fdat}{
        Feature data.
        A tab separated text file containing annotation information 
        for the features.
        Exactly *TWO* columns; 1st col = feature IDs; 
        2nd col = corresponding KEGG gene ID
        for each feature ID in 1st col; NO headers, row or column names.    
        The number of rows/features in this file should match the number of 
        rows/features of the expression matrix. 
        If 'exprs' is an object of 
        \code{\link[Biobase:ExpressionSet-class]{ExpressionSet-class}},
        the 'fdat' argument can be left unspecified, which then expects
        feature and gene IDs in respectively named columns
        'PROBE' and 'GENE' in the 'fData' slot of the ExpressionSet.
    }
    \item{org}{
        Organism under investigation in KEGG three letter code, 
        e.g. \sQuote{hsa} for \sQuote{Homo sapiens}.
        See also \code{\link{kegg.species.code}} to convert your 
        organism of choice to KEGG three letter code.
    }
    \item{data.type}{
        Expression data type.
        Use 'ma' for microarray and 'rseq' for RNA-seq data.
        If NA, data.type is automatically guessed. 
        If the expression values in 'eset' are decimal numbers 
        they are assumed to be microarray intensities.
        Whole numbers are assumed to be RNA-seq read counts. 
        Defaults to NA.
    }
    \item{norm.method}{
        Determines whether and how the expression data should be normalized.
        For available microarray normalization methods see the man page of the 
        limma function \code{\link{normalizeBetweenArrays}}.
        For available RNA-seq normalization methods see the man page of the 
        EDASeq function \code{\link{betweenLaneNormalization}}.
        Defaults to 'quantile', i.e. normalization is carried out so that 
        quantiles between arrays/lanes/samples are equal. 
        Use 'none' to indicate that the data is already normalized and should
        not be normalized by ebrowser.
        See the man page of \code{\link{normalize}} for details.
    }
    \item{de.method}{
        Determines which method is used for per-gene differential expression 
        analysis. See the man page of \code{\link{de.ana}} for details.
        Defaults to 'limma', i.e. differential expression is calculated based 
        on the typical limma \code{\link{lmFit}} procedure.    
    }
    \item{gs}{
        Gene sets.
        Either a list of gene sets (vectors of KEGG gene IDs) 
        or a text file in GMT format storing all gene sets under investigation.
    }
    \item{grn}{
        Gene regulatory network.
        Either an absolute file path to a tabular file or a character matrix 
        with exactly *THREE* cols; 1st col = IDs of regulating genes; 
        2nd col = corresponding regulated genes;
        3rd col = regulation effect; Use '+' and '-' for activation/inhibition.
    }
    \item{perm}{
        Number of permutations of the expression matrix to estimate the null 
        distribution. Defaults to 1000.
    } 
    \item{alpha}{
        Statistical significance level. Defaults to 0.05.
    }
    \item{beta}{
        Log2 fold change significance level. Defaults to 1 (2-fold).
    } 
    \item{comb}{
        Logical. Should results be combined if more then one enrichment method 
        is selected? Defaults to FALSE.
    }
    \item{browse}{
        Logical. Should results be displayed in the browser for interactive 
        exploration? Defaults to TRUE.
    }
    \item{nr.show}{
        Number of gene sets to show. 
        As default all statistical significant gene sets are displayed.
    }
}   
\value{
    None, opens the browser to explore results.
}
\seealso{
    \code{\link{read.eset}} to read expression data from file; 
    \code{\link{probe.2.gene.eset}} to transform probe to gene level expression;
    \code{\link{kegg.species.code}} maps species name to KEGG code.
    \code{\link{get.kegg.genesets}} to retrieve gene set definitions from KEGG;
    \code{\link{compile.grn.from.kegg}} to construct a GRN from KEGG pathways;
    \code{\link{sbea}} to perform set-based enrichment analysis;
    \code{\link{nbea}} to perform network-based enrichment analysis;
    \code{\link{comb.ea.results}} to combine results from different methods;
    \code{\link{ea.browse}} for exploration of resulting gene sets    
}
\references{
Limma User's guide:
\url{http://www.bioconductor.org/packages/limma}
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@bio.ifi.lmu.de>
}
\examples{
    # expression data from file
    exprs.file <- system.file("extdata/exprs.tab", package="EnrichmentBrowser")
    pdat.file <- system.file("extdata/pData.tab", package="EnrichmentBrowser")
    fdat.file <- system.file("extdata/fData.tab", package="EnrichmentBrowser")
    
    # getting all human KEGG gene sets
    # hsa.gs <- get.kegg.genesets("hsa")
    gs.file <- system.file("extdata/hsa_kegg_gs.gmt", package="EnrichmentBrowser")
    hsa.gs <- parse.genesets.from.GMT(gs.file)

    # set-based enrichment analysis
    ebrowser(   meth="ora", 
            exprs=exprs.file, pdat=pdat.file, fdat=fdat.file, 
            gs=hsa.gs, org="hsa", nr.show=3)

    # compile a gene regulatory network from KEGG pathways
    # hsa.grn <- compile.grn.from.kegg("hsa")
    pwys <- system.file("extdata/hsa_kegg_pwys.zip", package="EnrichmentBrowser")
    hsa.grn <- compile.grn.from.kegg(pwys)
   
    # network-based enrichment analysis
    ebrowser(   meth="ggea", 
            exprs=exprs.file, pdat=pdat.file, fdat=fdat.file, 
            gs=hsa.gs, grn=hsa.grn, org="hsa", nr.show=3 )

    # combining results
    ebrowser(   meth=c("ora", "ggea"), comb=TRUE,
            exprs=exprs.file, pdat=pdat.file, fdat=fdat.file, 
            gs=hsa.gs, grn=hsa.grn, org="hsa", nr.show=3 )
}

