\name{read.eset}
\alias{read.eset}
\title{ 
    Reading gene expression data from file into an expression set 
}
\description{
    The function reads in plain expression data with minimum annotation 
    requirements for the pData and fData slots. A differential 
    expression analysis can be performed and annotated to the expression set.
    If desired, overview plots such as a heatmap, p-value distribution and 
    volcano plot (fold change vs. p-value) are created.
}
\usage{
    read.eset(exprs.file, pdat.file, fdat.file, de=TRUE,
                heatm.file=NULL, distr.file = NULL, volc.file = NULL)
}
\arguments{
    \item{exprs.file}{
        Expression matrix.
        A tab separated text file containing *normalized* expression values 
        on a *log* scale.
        Columns = samples/subjects; rows = features/probes/genes; 
        NO headers, row or column names.
        Supported data types are log2 counts (microarray single-channel),
        log2 ratios (microarray two-color), and log2-counts per million 
        (RNA-seq logCPMs).
        See details.
    }
    \item{pdat.file}{
        Phenotype data.
        A tab separated text file containing annotation information 
        for the samples in either *two or three* columns. 
        NO headers, row or column names.
        The number of rows/samples in this file should match the number of 
        columns/samples of the expression matrix.
        The 1st colum is reserved for the sample IDs; 
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
        Exactly *TWO* columns; 1st col = feature IDs; 
        2nd col = corresponding KEGG gene ID for each feature ID in 1st col; 
        NO headers, row or column names.    
        The number of rows/features in this file should match the number of 
        rows/features of the expression matrix. 
        Alternatively, this can also be the ID of a recognized platform
        such as 'hgu95av2' (Affymetrix Human Genome U95 chip) or 'ecoli2'
        (Affymetrix E. coli Genome 2.0 Array). 
        See details.
    }
    \item{de}{
        Logical.
        Should a simultanous differential expression (de) analysis be carried 
        out for each gene. Defaults to TRUE.
    }
    \item{heatm.file}{
        Optional. If specified, a heatmap is plotted in PNG format to file.
    }
    \item{distr.file}{
        Optional. If specified, the p-value distribution is plotted 
        in PNG format to file.
    }
    \item{volc.file}{
        Optional. If specified, the volcano is plotted in 
        PNG format to file.
    }
}
\value{
    An \code{\link[Biobase:ExpressionSet-class]{ExpressionSet-class}} 
    with measures of differential expression annotated in the fData slot.
}
\details{
See the limma's user guide 
http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
for definition and normalization of the different expression data types.        

In case of microarry data the feature IDs typically correspond to probe IDs.
Thus, the fdat.file should define a mapping from probe ID (1st column) to
corresponding KEGG gene ID (2nd column).
The mapping can be defined automatically by providing the ID of a 
recognized platform such as 'hgu95av2' (Affymetrix Human Genome U95 chip).
This requires that a corresponding '.db' package exists (see 
http://www.bioconductor.org/packages/release/BiocViews.html#___ChipName
for all available chips/packages) and that you have it installed.
*However, this option should be used with care*. Existing mappings might
be outdated and sometimes the KEGG gene ID does not correspond to the Entrez ID 
(e.g. for E. coli and S. cerevisae). In these cases probe identifiers are mapped
 twice (probe ID -> Entrez ID -> KEGG ID) which almost always results in loss
of information. Thus, mapping quality should always be checked and in case
properly defined with a 2 column fdat.file.
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@bio.ifi.lmu.de>
}
\seealso{
    \code{\link[Biobase:ExpressionSet-class]{ExpressionSet-class}}, 
    \code{\link{eBayes}},
    \code{\link{voom}}, 
    \code{\link{p.adjust}}
}
\examples{
    # reading the expression data from file
    exprs.file <- system.file("extdata/ALL_exprs.tab", package="EnrichmentBrowser")
    pdat.file <- system.file("extdata/ALL_pData.tab", package="EnrichmentBrowser")
    fdat.file <- system.file("extdata/ALL_fData.tab", package="EnrichmentBrowser")
    eset <- read.eset(exprs.file, pdat.file, fdat.file)
    head(fData(eset))
}

