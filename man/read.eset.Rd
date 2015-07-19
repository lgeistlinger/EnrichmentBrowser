\name{read.eset}
\alias{read.eset}
\title{ 
    Reading gene expression data from file into an expression set 
}
\description{
    The function reads in plain expression data from file with minimum 
    annotation requirements for the pData and fData slots. 
}
\usage{
    read.eset( exprs.file, pdat.file, fdat.file, 
        data.type = c(NA, "ma", "rseq"), NA.method = c("mean", "rm", "keep") )
}
\arguments{
    \item{exprs.file}{
        Expression matrix.
        A tab separated text file containing expression values. 
        Columns = samples/subjects; rows = features/probes/genes; 
        NO headers, row or column names.
        See details.
    }
    \item{pdat.file}{
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
    }
    \item{fdat.file}{
        Feature data.
        A tab separated text file containing annotation information 
        for the features.
        In case of probe level data: exactly *TWO* columns; 
        1st col = probe/feature IDs; 
        2nd col = corresponding gene ID for each feature ID in 1st col; 
        In case of gene level data: The list of gene IDs newline-separated 
        (i.e. just one column).
        It is recommended to use *ENTREZ* gene IDs (to benefit from downstream 
        visualization and exploration functionality of the enrichment analysis).
        NO headers, row or column names.    
        The number of rows (features/probes/genes) in this file should match 
        the number of rows/features of the expression matrix. 
        Alternatively, this can also be the ID of a recognized platform
        such as 'hgu95av2' (Affymetrix Human Genome U95 chip) or 'ecoli2'
        (Affymetrix E. coli Genome 2.0 Array). 
        See details.
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
    \item{NA.method}{
        Determines how to deal with NA's (missing values).
        This can be one out of:
        \itemize{
        \item mean: replace NA's by the row means for a feature over all samples.
        \item rm: rows (features) that contain NA's are removed.
        \item keep: do nothing. Missing values are kept (which, however, can 
            then cause several issues in the downstream analysis)          
        }
        Defaults to 'mean'.
    }
}
\value{
    An object of \code{\linkS4class{ExpressionSet}}.
}
\details{
See the limma's user guide 
\url{http://www.bioconductor.org/packages/limma}
for definition and normalization of the different expression data types.        

In case of microarry data the feature IDs typically correspond to probe IDs.
Thus, the fdat.file should define a mapping from probe ID (1st column) to
corresponding KEGG gene ID (2nd column).
The mapping can be defined automatically by providing the ID of a 
recognized platform such as 'hgu95av2' (Affymetrix Human Genome U95 chip).
This requires that a corresponding '.db' package exists (see 
\url{http://www.bioconductor.org/packages/release/BiocViews.html#___ChipName}
for all available chips/packages) and that you have it installed.
*However, this option should be used with care*. Existing mappings might
be outdated and sometimes the KEGG gene ID does not correspond to the Entrez ID 
(e.g. for E. coli and S. cerevisae). In these cases probe identifiers are mapped
 twice (probe ID -> Entrez ID -> KEGG ID), which almost always results in loss
of information. Thus, mapping quality should always be checked and in case
properly defined with a 2-column fdat.file.
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@bio.ifi.lmu.de>
}
\seealso{
    \code{\link[Biobase:ExpressionSet-class]{ExpressionSet-class}}
}
\examples{
    # reading the expression data from file
    exprs.file <- system.file("extdata/exprs.tab", package="EnrichmentBrowser")
    pdat.file <- system.file("extdata/pData.tab", package="EnrichmentBrowser")
    fdat.file <- system.file("extdata/fData.tab", package="EnrichmentBrowser")
    eset <- read.eset(exprs.file, pdat.file, fdat.file)
}

