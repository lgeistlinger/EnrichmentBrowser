######################################################################
#
# author: Ludwig Geistlinger
# date: 25 Feb 2011
#
#
# THE ENRICHMENT BROWSER
#
# update: 05 June 2014, eBrowser2.0 all-in-one wrapper function
#
######################################################################


# CONSTANTS: 
# TODO resolve:
# options: 
# 1) par-like style --> current solution
#
# 2) methods/functions for getting/setting, e.g. group(obj), adjp(obj), ... 
#
# 3) convert them all to function parameters 

.ebrowser_config_cache <- new.env(parent=emptyenv())



#' Configuring the EnrichmentBrowser
#' 
#' Function to get and set configuration parameters determining the default
#' behavior of the EnrichmentBrowser
#' 
#' Important colData, rowData, and result column names: \itemize{ \item
#' SMPL.COL: colData column storing the sample IDs (default: "SAMPLE") \item
#' GRP.COL: colData column storing binary group assignment (default: "GROUP")
#' \item BLK.COL: colData column defining paired samples or sample blocks
#' (default: "BLOCK")
#' 
#' \item PRB.COL: rowData column storing probe/feature IDs ("PROBEID",
#' read-only) \item EZ.COL: rowData column storing gene ENTREZ IDs ("ENTREZID",
#' read-only) \item SYM.COL: rowData column storing gene symbols ("SYMBOL",
#' read-only) \item GN.COL: rowData column storing gene names ("GENENAME",
#' read-only)
#' 
#' \item FC.COL: rowData column storing (log2) fold changes of differential
#' expression between sample groups (default: "FC") \item ADJP.COL: rowData
#' column storing adjusted (corrected for multiple testing) p-values of
#' differential expression between sample groups (default: "ADJ.PVAL")
#' 
#' \item GS.COL: result table column storing gene set IDs (default: "GENE.SET")
#' \item PVAL.COL: result table column storing gene set significance (default:
#' "PVAL") \item PMID.COL: gene table column storing PUBMED IDs ("PUBMED",
#' read-only) }
#' 
#' Important URLs (all read-only): \itemize{ \item NCBI.URL:
#' \url{http://www.ncbi.nlm.nih.gov/} \item PUBMED.URL:
#' \url{http://www.ncbi.nlm.nih.gov/pubmed/} \item GENE.URL:
#' \url{http://www.ncbi.nlm.nih.gov/gene/} \item KEGG.URL:
#' \url{http://www.genome.jp/dbget-bin/} \item KEGG.GENE.URL:
#' \url{http://www.genome.jp/dbget-bin/www_bget?} \item KEGG.SHOW.URL:
#' \url{http://www.genome.jp/dbget-bin/show_pathway?} \item GO.SHOW.URL:
#' \url{http://amigo.geneontology.org/amigo/term/} }
#' 
#' Default output directory: \itemize{ \item EBROWSER.HOME:
#' \code{rappdirs::user_data_dir("EnrichmentBrowser")} \item OUTDIR.DEFAULT:
#' \code{file.path(EBROWSER.HOME, "results")} }
#' 
#' Gene set size: \itemize{ \item GS.MIN.SIZE: minimum number of genes per gene
#' set (default: 5) \item GS.MAX.SIZE: maximum number of genes per gene set
#' (default: 500) }
#' 
#' Result appearance: \itemize{ \item RESULT.TITLE: (default: "Table of
#' Results") \item NR.SHOW: maximum number of entries to show (default: 20) }
#'
#' @aliases config.ebrowser 
#' @param key Configuration parameter.
#' @param value Value to overwrite the current value of key.
#' @return If is.null(value) this returns the value of the selected
#' configuration parameter.  Otherwise, it updates the selected parameter with
#' the given value.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @examples
#' 
#'     # getting config information
#'     configEBrowser("GS.MIN.SIZE") 
#' 
#'     # setting config information
#'     # WARNING: this is for advanced users only!
#'     # inappropriate settings will impair EnrichmentBrowser's functionality
#'     configEBrowser(key="GS.MIN.SIZE", value=3)  
#' 
#' @export configEBrowser
configEBrowser <- function(key, value=NULL) 
{
    .key_readonly <- c(
        "PRB.COL", "EZ.COL", "GN.COL", "SYM.COL", "PMID.COL", 
        "NCBI.URL", "PUBMED.URL", "GENE.URL", "KEGG.URL", "KEGG.GENE.URL",
        "KEGG.SHOW.URL", "GO.SHOW.URL", "SBEA.PKGS", "NBEA.PKGS")
 
    if(is.null(value)) .ebrowser_config_cache[[key]]
    else if(!(key %in% .key_readonly)) .ebrowser_config_cache[[key]] <- value
}

#' @export
#' @keywords internal
config.ebrowser <- function(key, value=NULL)     
{
    .Deprecated("configEBrowser")
    configEBrowser(key, value)
}

##
# check the output directory
##
.checkOutDir <- function(out.dir)
{
    if(!file.exists(dirname(out.dir)))
        stop(paste0("Not a valid output directory path \'",out.dir,"\'"))
    if(!file.exists(out.dir)){
        message(paste("Output directory", out.dir, 
            "does not exist.\nThus, the directory is going to be created."))
        dir.create(out.dir)
    }
}


##
##
# eBrowser functionality
##


#' Seamless navigation through enrichment analysis results
#' 
#' This is the all-in-one wrapper function to perform the standard enrichment
#' analysis pipeline implemented in the EnrichmentBrowser package.
#' 
#' Given flat gene expression data, the data is read in and subsequently
#' subjected to chosen enrichment analysis methods.
#' 
#' The results from different methods can be combined and investigated in
#' detail in the default browser.
#' 
#' *On data type and normalization:*
#'  
#' Normalization of high-throughput expression data is essential to make
#' results within and between experiments comparable.  Microarray (intensity
#' measurements) and RNA-seq (read counts) data exhibit typically distinct
#' features that need to be normalized for.  This function wraps commonly used 
#' functionality from limma for microarray normalization and from EDASeq for 
#' RNA-seq normalization. For specific needs that deviate
#' from standard normalizations, the user should always refer to more
#' specific functions/packages.  See also the limma's user guide
#' \url{http://www.bioconductor.org/packages/limma} for definition and
#' normalization of the different expression data types.
#' 
#' Microarray data is expected to be single-channel.  For two-color arrays, it
#' is expected here that normalization within arrays has been already carried
#' out, e.g. using \code{\link{normalizeWithinArrays}} from limma.
#' 
#' RNA-seq data is expected to be raw read counts.  Please note that
#' normalization for downstream DE analysis, e.g. with edgeR and DESeq2, is not
#' ultimately necessary (and in some cases even discouraged) as many of these
#' tools implement specific normalization approaches.  See the vignette of
#' EDASeq, edgeR, and DESeq2 for details.
#' 
#' 
#' @param meth Enrichment analysis method(s).  See \code{\link{sbeaMethods}} and
#' \code{\link{nbeaMethods}} for currently supported enrichment analysis
#' methods.  See also \code{\link{sbea}} and \code{\link{nbea}} for details.
#' @param exprs Expression matrix.  A tab separated text file containing
#' the expression values (microarray: intensity measurements, RNA-seq: read counts).  
#' Columns = samples/subjects; rows = features/probes/genes; NO headers, row or 
#' column names. Alternatively, this can be a
#' \code{\linkS4class{SummarizedExperiment}}, assuming the expression matrix in
#' the \code{\link{assays}} slot. See details.
#' @param cdat Column (phenotype) data.  A tab separated text file containing annotation
#' information for the samples in either *two or three* columns.  NO headers,
#' row or column names.  The number of rows/samples in this file should match
#' the number of columns/samples of the expression matrix.  The 1st column is
#' reserved for the sample IDs; The 2nd column is reserved for a *BINARY* group
#' assignment.  Use '0' and '1' for unaffected (controls) and affected (cases)
#' sample class, respectively.  For paired samples or sample blocks a third
#' column is expected that defines the blocks.  If 'exprs' is a
#' \code{\linkS4class{SummarizedExperiment}}, the 'cdat' argument can be left
#' unspecified, which then expects group and optional block assignments in
#' respectively named columns 'GROUP' (mandatory) and 'BLOCK' (optional) in the
#' \code{\link{colData}} slot.
#' @param rdat Row (feature) data.  A tab separated text file containing annotation
#' information for the features. In case of probe level data:
#' exactly *TWO* columns; 1st col = probe/feature IDs; 2nd col = corresponding
#' gene ID for each feature ID in 1st col. In case of gene level data: the
#' gene IDs newline-separated (i.e. just *one* column).  It is recommended to
#' use *ENTREZ* gene IDs (to benefit from downstream visualization and
#' exploration functionality of the EnrichmentBrowser).  NO headers, row or
#' column names.  The number of rows (features/probes/genes) in this file
#' should match the number of rows/features of the expression matrix.
#' Alternatively, this can also be the ID of a recognized platform such as
#' 'hgu95av2' (Affymetrix Human Genome U95 chip) or 'ecoli2' (Affymetrix E.
#' coli Genome 2.0 Array).
#' If 'exprs' is a \code{\linkS4class{SummarizedExperiment}}, the 'rdat' argument
#' can be left unspecified, which then expects probe and corresponding Entrez
#' Gene IDs in respectively named columns 'PROBEID' and 'ENTREZID' in the
#' \code{\link{rowData}} slot.
#' @param org Organism under investigation in KEGG three letter code, e.g.
#' \sQuote{hsa} for \sQuote{Homo sapiens}.  See also
#' \code{\link{kegg.species.code}} to convert your organism of choice to KEGG
#' three letter code.
#' @param data.type Expression data type.  Use 'ma' for microarray and 'rseq'
#' for RNA-seq data.  If NA, data.type is automatically guessed.  If the
#' expression values in 'exprs' are decimal numbers they are assumed to be
#' microarray intensities.  Whole numbers are assumed to be RNA-seq read
#' counts.  Defaults to NA.
#' @param norm.method Determines whether and how the expression data should be
#' normalized.  For available microarray normalization methods see the man page
#' of the limma function \code{\link{normalizeBetweenArrays}}.  For available
#' RNA-seq normalization methods see the man page of the EDASeq function
#' \code{\link{betweenLaneNormalization}}.  Defaults to 'quantile', i.e.
#' normalization is carried out so that quantiles between arrays/lanes/samples
#' are equal.  Use 'none' to indicate that the data is already normalized and
#' should not be normalized by ebrowser.  See the man page of
#' \code{\link{normalize}} for details.
#' @param de.method Determines which method is used for per-gene differential
#' expression analysis. See the man page of \code{\link{deAna}} for details.
#' Defaults to 'limma', i.e. differential expression is calculated based on the
#' typical limma \code{\link{lmFit}} procedure.
#' @param gs Gene sets.  Either a list of gene sets (character vectors of gene
#' IDs) or a text file in GMT format storing all gene sets under investigation.
#' @param grn Gene regulatory network.  Either an absolute file path to a
#' tabular file or a character matrix with exactly *THREE* cols; 1st col = IDs
#' of regulating genes; 2nd col = corresponding regulated genes; 3rd col =
#' regulation effect; Use '+' and '-' for activation/inhibition.
#' @param perm Number of permutations of the sample group assignments. 
#' Defaults to 1000. Can also be an integer vector matching
#' the length of 'meth' to assign different numbers of permutations for
#' different methods.
#' @param alpha Statistical significance level. Defaults to 0.05.
#' @param beta Log2 fold change significance level. Defaults to 1 (2-fold).
#' @param comb Logical. Should results be combined if more then one enrichment
#' method is selected? Defaults to FALSE.
#' @param browse Logical. Should results be displayed in the browser for
#' interactive exploration? Defaults to TRUE.
#' @param nr.show Number of gene sets to show.  As default all statistical
#' significant gene sets are displayed.  Note that this only influences the
#' number of gene sets for which additional visualization will be provided
#' (typically only of interest for the top / signifcant gene sets).  Selected
#' enrichment methods and resulting flat gene set rankings still include the
#' complete number of gene sets under study.
#' @param out.dir Output directory. If \code{NULL}, defaults to a 
#' timestamp-generated subdirectory of \code{configEBrowser("OUTDIR.DEFAULT")}. 
#' @param report.name Character. Name of the HTML report. Defaults to \code{"index"}.
#' @param ... Additional arguments passed on to the individual building blocks.
#' @return None, writes an HTML report and, if selected, opens the browser to 
#' explore results.
#
#' If not instructed otherwise (via argument \code{out.dir}), 
#' the main HTML report and associated files are written to 
#' \code{configEBrowser("OUTDIR.DEFAULT")}. 
#' See \code{?configEBrowser} to change the location. 
#' If \code{browse=TRUE}, the HTML report will automatically be opened in 
#' the your default browser.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{\link{readSE}} to read expression data from file;
#' \code{\link{probe2gene}} to transform probe to gene level expression;
#' \code{\link{kegg.species.code}} maps species name to KEGG code.
#' \code{\link{getGenesets}} to retrieve gene set databases such as GO or KEGG;
#' \code{\link{compileGRN}} to construct a GRN from pathway databases;
#' \code{\link{sbea}} to perform set-based enrichment analysis;
#' \code{\link{nbea}} to perform network-based enrichment analysis;
#' \code{\link{combResults}} to combine results from different methods;
#' \code{\link{eaBrowse}} for exploration of resulting gene sets
#' @references Limma User's guide:
#' \url{http://www.bioconductor.org/packages/limma}
#' @examples
#' 
#'     # expression data from file
#'     exprs.file <- system.file("extdata/exprs.tab", package="EnrichmentBrowser")
#'     cdat.file <- system.file("extdata/colData.tab", package="EnrichmentBrowser")
#'     rdat.file <- system.file("extdata/rowData.tab", package="EnrichmentBrowser")
#'     
#'     # getting all human KEGG gene sets
#'     # hsa.gs <- getGenesets(org="hsa", db="kegg")
#'     gs.file <- system.file("extdata/hsa_kegg_gs.gmt", package="EnrichmentBrowser")
#'     hsa.gs <- getGenesets(gs.file)
#'
#'     # output destination 
#'     out.dir <- configEBrowser("OUTDIR.DEFAULT") 
#'
#'     # set-based enrichment analysis
#'     ebrowser( meth="ora", perm=0,
#'             exprs=exprs.file, cdat=cdat.file, rdat=rdat.file, 
#'             gs=hsa.gs, org="hsa", nr.show=3,
#'             out.dir=out.dir, report.name="oraReport")
#' 
#'     # compile a gene regulatory network from KEGG pathways
#'     hsa.grn <- compileGRN(org="hsa", db="kegg")
#'    
#'     # network-based enrichment analysis
#'     ebrowser(   meth="ggea", 
#'             exprs=exprs.file, cdat=cdat.file, rdat=rdat.file, 
#'             gs=hsa.gs, grn=hsa.grn, org="hsa", nr.show=3,
#'             out.dir=out.dir, report.name="ggeaReport")
#' 
#'     # combining results
#'     ebrowser( meth=c("ora", "ggea"), perm=0, comb=TRUE,
#'             exprs=exprs.file, cdat=cdat.file, rdat=rdat.file, 
#'             gs=hsa.gs, grn=hsa.grn, org="hsa", nr.show=3,
#'             out.dir=out.dir, report.name="combReport")
#' 
#' @export ebrowser
ebrowser <- function(
    meth, exprs, cdat, rdat, org, data.type=c(NA, "ma", "rseq"),
    norm.method="quantile", de.method="limma",
    gs, grn=NULL, perm=1000, alpha=0.05, beta=1, 
    comb=FALSE, browse=TRUE, nr.show=-1, out.dir=NULL, report.name="index", ...)
{
    GRP.COL <- configEBrowser("GRP.COL")
    FC.COL <- configEBrowser("FC.COL")
    ADJP.COL <- configEBrowser("ADJP.COL")
    EZ.COL <- configEBrowser("EZ.COL")    

    METHODS <- c(sbeaMethods(), nbeaMethods())
    is.method <- (meth %in% METHODS) | is.function(meth)
    if(!all(is.method)) stop("No such method: ", meth[!is.method])
    
    if(any(nbeaMethods() %in% meth)) 
        if(is.null(grn))
            stop("\'grn\' must be not null")
  
    if(is.null(out.dir))
    { 
        out.dir <- configEBrowser("OUTDIR.DEFAULT")
        stamp <- format(Sys.time(), "%a_%b%d_%Y_%H%M%S")
        out.dir <- file.path(out.dir, stamp)
    }
    else out.dir <- path.expand(out.dir)
    if(!file.exists(out.dir)) dir.create(out.dir, recursive=TRUE)    

    # execution
    # read expression data
    data.type <- match.arg(data.type)
    if(is.character(exprs))
    {
        message("Read expression data ...")
        se <- readSE( assay.file=exprs, 
            cdat.file=cdat, rdat.file=rdat, data.type=data.type )
    }
    else
    { 
        if(is(exprs, "ExpressionSet")) exprs <- as(exprs, "SummarizedExperiment")
        se <- exprs
        if(is.na(data.type))
            data.type <- .detectDataType(assay(se))
        metadata(se)$dataType <- data.type
    }
    
    # normalize?
    if(norm.method != 'none')
    {
        message("Normalize ...")
        se <- normalize(se, norm.method=norm.method)
    }
    
    # probe 2 gene if data.type = ma
    # ... and it's not already a gene level se
    if(metadata(se)$dataType == "ma")
    {
        has.pcol <- configEBrowser("PRB.COL") %in% colnames(rowData(se))
        anno <- metadata(se)$annotation
        has.anno <- ifelse(length(anno), nchar(anno) > 3, FALSE)
        if(has.pcol || has.anno)
        {
            message("Transform probe expression to gene expression ...")    
            geneSE <- probe2gene(se)
        }
        else geneSE <- se
    }
    else geneSE <- se

    message("DE analysis ...")    
    geneSE <- deAna(geneSE, de.method=de.method)
    if(missing(org)) org <- metadata(geneSE)$annotation 
    else metadata(geneSE)$annotation <- org
        
    nr.meth <- length(meth)
	if(length(perm) != nr.meth) perm <- rep(perm[1], nr.meth)
    if(comb) res.list <- vector("list", length=nr.meth)
    for(i in seq_len(nr.meth))
    {
        m <- meth[i]
        message(paste("Execute", toupper(m), "..."))
        out.file <- file.path(out.dir, paste0(m, ".txt"))

        if(m %in% nbeaMethods()) 
            res <- nbea( method=m, se=geneSE, gs=gs, 
                    grn=grn, alpha=alpha, beta=beta, perm=perm[i], ... )

        else res <- sbea( method=m, se=geneSE, 
                    gs=gs, alpha=alpha, perm=perm[i], ... )

        write.table(res$res.tbl, file=out.file, 
            quote=FALSE, row.names=FALSE, sep="\t")

        # produce html reports, if desired
        if(browse) eaBrowse(res, nr.show, 
                            graph.view=grn, html.only=TRUE, out.dir=out.dir)
        
        # link gene statistics
        sam.file <- file.path(out.dir, "samt.RData")    
        if(m == "samgs" && file.exists(sam.file))
        {
            samt <- round(get(load(sam.file)), digits=2)
            rowData(geneSE, use.names=TRUE)[names(s2n), "SAM.T"] <- unname(samt) 
        }       

        s2n.file <- file.path(out.dir, "gsea_s2n.RData")
        if(m == "gsea" && file.exists(s2n.file))
        {
            s2n <- round(get(load(s2n.file)), digits=2)
            rowData(geneSE, use.names=TRUE)[names(s2n),"GSEA.S2N"] <- unname(s2n) 
        }
        if(comb) res.list[[i]] <- res
    }

    # write genewise differential expression
    gt.file <- file.path(out.dir, "de.txt")
	message("Annotating genes ...")
    gt <- .getGeneAnno(names(geneSE), org)
    gt <- cbind(gt, rowData(geneSE, use.names=TRUE))
    gt <- .sortGeneTable(gt)
    ind <- gt[,configEBrowser("EZ.COL")]
    geneSE <- geneSE[ind, ]
    rowData(geneSE) <- gt

    message(paste("Genewise differential expression written to", gt.file))    
    write.table(gt, file=gt.file, row.names=FALSE, quote=FALSE, sep="\t")

    if(comb)
    {
        # combine results (average ranks, ztrans p-vals)
        message("Combine results ...")
        if(nr.meth == 1) 
            message(paste("Only one method given, \'comb\' ignored."))
        else
        {
            # combine results in a specified out file
            out.file <- file.path(out.dir, "comb.txt")
            res <- combResults(res.list=res.list)
            write.table(res$res.tbl, 
                file=out.file, quote=FALSE, row.names=FALSE, sep="\t")

            # produce html report, if desired
            if(browse) eaBrowse(res, nr.show, 
                                graph.view=grn, html.only=TRUE, out.dir=out.dir)
        }
    }
    
    message(paste("Your output files are in", out.dir))
    
    # create INDEX.html
    if(browse)
    { 
        # plot DE
        if(length(geneSE) > 1000)
        {
            message("Restricting global view to the 1000 most significant genes")
            geneSE <- geneSE[1:1000,]
        }
        vs <- .viewSet(geneSE, out.prefix=file.path(out.dir, "global"))
        
        message("Produce html report ...")
        .createIndex(meth, comb, out.dir, report.name)
    }
}

