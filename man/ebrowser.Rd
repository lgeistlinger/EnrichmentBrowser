\name{ebrowser}
\alias{ebrowser}
\title{
    Seamless navigation through enrichment analysis results 
}
\description{
    This is the all-in-one wrapper function to perform 
    the standard enrichment analysis pipeline implemented 
    in the EnrichmentBrowser package.
    
    Given flat gene expression data on probe level, the 
    data is read in, transformed to gene level and subsequently
    subjected to chosen enrichment analysis methods.

    The results from different methods can be combined and
    investigated in detail in the default browser. 
}
\usage{
    ebrowser( meth, exprs, pdat, fdat, value.type="log2count",
            gs, grn=NULL, perm=1000, alpha=0.05, beta=1,
            comb=FALSE, browse=TRUE, nr.show=-1, out.dir=NULL )
}
\arguments{
    \item{meth}{
        Enrichment analysis method.
        Currently, the following enrichment analysis methods are supported:
        \sQuote{ora}, \sQuote{safe}, \sQuote{gsea}, \sQuote{samgs}, 
        \sQuote{ggea}, and \sQuote{nea}.
        See \sQuote{sbea} and \sQuote{nbea} for details.
    }
    \item{exprs}{
        Expression matrix.
        A tab separated text file containing the numerical expression values.
        Colums = samples/subjects; rows = features/probes/genes; NO headers, 
        row or column names. 
    }
    \item{pdat}{
        Phenotype data.
        A tab separated text file containing annotation information 
        for the samples. 
        Exactly *TWO* cols; 1st col = sample names; 
        2nd column = *BINARY* group assignment. 
        Use '0' and '1' for unaffected and affected sample class, respectively. 
        NO headers, row or column names.
        The nr of rows/samples in this file should match the number of 
        columns/samples of the expression matrix. 
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
    }
    \item{value.type}{
        Type of values in the expression matrix: 
        either 'log2count' or 'log2ratio' (defaults to 'log2count').   
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
        Log2 fold change signifcance level. Defaults to 1 (2-fold).
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
    \item{out.dir}{
        Optional output directory all results will be written to.
    }
}   
\value{
    None, opens the browser to explore results.
}
\seealso{
    \code{\link{read.eset}} to read expression data from file; 
    \code{\link{probe.2.gene.eset}} to transform probe to gene level expression;
    \code{\link{get.kegg.genesets}} to retrieve gene set definitions from KEGG;
    \code{\link{compile.grn.from.kegg}} to construct a GRN from KEGG pathways;
    \code{\link{sbea}} to perform set-based enrichment analysis;
    \code{\link{nbea}} to perform network-based enrichment analysis;
    \code{\link{comb.ea.results}} to combine results from different methods;
    \code{\link{ea.browse}} for exploration of resulting gene sets    
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@bio.ifi.lmu.de>
}
\examples{
    # expression data from file
    exprs.file <- system.file("extdata/ALL_exprs.tab", package="EnrichmentBrowser")
    pdat.file <- system.file("extdata/ALL_pData.tab", package="EnrichmentBrowser")
    fdat.file <- system.file("extdata/ALL_fData.tab", package="EnrichmentBrowser")
    
    # getting all human KEGG gene sets
    # hsa.gs <- get.kegg.genesets("hsa")
    gs.file <- system.file("extdata/hsa_kegg_gs.gmt", package="EnrichmentBrowser")
    hsa.gs <- parse.genesets.from.GMT(gs.file)

    \donttest{
    # set-based enrichment analysis
    ebrowser(   meth="ora", 
            exprs=exprs.file, 
            pdat=pdat.file, 
            fdat=fdat.file, 
            gs=hsa.gs)
    }

    # compile a gene regulatory network from KEGG pathways
    # hsa.grn <- compile.grn.from.kegg("hsa")
    pwys <- system.file("extdata/hsa_kegg_pwys.zip", package="EnrichmentBrowser")
    hsa.grn <- compile.grn.from.kegg(pwys)
   
    \donttest{
    # network-based enrichment analysis
    ebrowser(   meth="ggea", 
            exprs=exprs.file, 
            pdat=pdat.file, 
            fdat=fdat.file, 
            gs=hsa.gs,
            grn=hsa.grn)
    }

    # combining results
    ebrowser(   meth=c("ora", "ggea"), 
            perm=100,
            exprs=exprs.file, 
            pdat=pdat.file, 
            fdat=fdat.file, 
            gs=hsa.gs,
            grn=hsa.grn)
}

