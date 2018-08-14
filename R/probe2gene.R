# 
# Author: ludwig geistlinger
# Date: May 27th, 2010
#
# transforms a probe eset to a gene eset
#
# UPDATE Oct 21th, 2014
# automatic probe 2 gene mapping with BioC annotation packages
# different summarization methods for probes (mean or most signif.)
#
#
###############################################################################

#' Transformation of probe level expression to gene level expression
#' 
#' Transforms expression data on probe level to gene level expression by
#' summarizing all probes that are annotated to a particular gene.
#' 
#' @aliases probe.2.gene.eset
#' @param probeSE Probe expression data.  An object of class
#' \code{\linkS4class{SummarizedExperiment}}.  Make sure that the
#' \code{\link{metadata}} contains an element named \code{annotation} that
#' provides the corresponding ID of a recognized platform such as
#' \code{hgu95av2} (Affymetrix Human Genome U95 chip).  This requires that a
#' corresponding \code{.db} package exists (see
#' \url{http://www.bioconductor.org/packages/release/BiocViews.html#___ChipName}
#' for available chips/packages) and that you have it installed.
#' Alternatively, the mapping from probe to gene can also be defined in the
#' \code{\link{rowData}} slot via two columns named (i) \code{PROBEID} for the
#' platform-specific probe ID, and (ii) \code{ENTREZID} for the corresponding
#' NCBI Entrez Gene ID.
#' @param chip Character. The ID of a recognized microarray platform. 
#' Only required if not provided in the \code{\link{metadata}} of \code{probeSE}
#' via an element named \code{annotation}. 
#' @param from Character. ID type from which should be mapped. Corresponds to the
#' ID type of the names of argument \code{se}, with the default \code{PROBEID}
#' being appropriate if the mapping is based on Bioconductor annotation packages. 
#' Note that \code{from} is ignored if \code{to} is a \code{\link{rowData}} column 
#' of \code{probeSE}. 
#' @param to Character. Gene ID type to which should be mapped. Corresponds to 
#' the gene ID type the rownames of argument \code{probeSE} should be updated with.
#' Note that this can also be the name of a column in the \code{\link{rowData}} 
#' slot of \code{probeSE} to specify user-defined mappings in which conflicts 
#' have been manually resolved. Defaults to \code{ENTREZID}.
#' @param multi.to How to resolve 1:many mappings, i.e. multiple gene IDs for a
#' single probe ID? This is passed on to the \code{multiVals} argument of
#' \code{\link{mapIds}} and can thus take several pre-defined values, but also
#' the form of a user-defined function. However, note that this requires that a
#' single gene ID is returned for each probe ID. Default is \code{"first"}, 
#' which accordingly returns the first gene ID mapped onto the respective probe ID.
#' @param multi.from How to resolve many:1 mappings, i.e. multiple probe IDs 
#' mapping to the same gene ID? Pre-defined options include:
#' \itemize{ \item 'mean' (Default): updates the respective gene expression with
#'  the average over the expression of all probes mapping to that gene, 
#' \item 'first': returns the first probe ID for each gene ID with 
#' multiple probe IDs,
#' \item 'minp' selects the probe ID with minimum p-value (according to the
#' \code{\link{rowData}} column \code{PVAL} of \code{probeSE}),
#' \item 'maxfc' selects the probe ID with maximum absolute log2 fold change 
#' (according to the \code{\link{rowData}} column \code{FC} of \code{probeSE}).}
#' @return A \code{\linkS4class{SummarizedExperiment}} on gene level.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{\link{readSE}} for reading expression data from file,
#' \code{\link{deAna}} for differential expression analysis.
#' @examples
#' 
#'     # (1) reading the expression data from file
#'     exprs.file <- system.file("extdata/exprs.tab", package="EnrichmentBrowser")
#'     cdat.file <- system.file("extdata/colData.tab", package="EnrichmentBrowser")
#'     rdat.file <- system.file("extdata/rowData.tab", package="EnrichmentBrowser")
#'     probeSE <- readSE(exprs.file, cdat.file, rdat.file)
#'     geneSE <- probe2gene(probeSE) 
#' 
#' @export probe2gene
probe2gene <- function(probeSE, chip=NA, from="PROBEID", to="ENTREZID", 
    multi.to="first", multi.from="mean")
{
    if(multi.from == "mean") 
        geneSE <- .probe2geneAverage(probeSE, chip, from, to, multi.to)
    else geneSE <- idMap(probeSE, chip, from, to, multi.to, multi.from)
    return(geneSE)
}

#' @export
#' @keywords internal
probe.2.gene.eset <- function(probe.eset, use.mean=TRUE)
{
    .Deprecated("probe2gene")
    multi.from <- ifelse(use.mean, "mean", "minp")
    probe2gene(probe.eset, multi.from=multi.from)
}


.probe2geneAverage <- function(se, chip, from, to, multi.to)
{
    if(is(se, "ExpressionSet")) 
        se <- as(se, "SummarizedExperiment")  

    if(!(to %in% colnames(rowData(se))))
    {
        anno <- chip 
        if(is.na(anno)) anno <- metadata(se)$annotation
        if(!length(anno)) stop("Chip under investigation not annotated")
        x <- .idmap(names(se), anno, from, to, excl.na=FALSE, 
                        multi.to=multi.to, resolve.multiFrom=FALSE)
        rowData(se)[[to]] <- unname(x)
        se <- .annotateSE(se, anno)
    }
    
    # remove probes without gene annotation
    not.na <- !is.na(rowData(se)[,to])
    if(sum(not.na) < nrow(se)) se <- se[not.na,]

    probe.exprs <- assay(se)
    p2g <- as.vector(rowData(se)[, to])
    
    # determine unique genes
    genes <- unique(p2g)
    
    # compute gene expression
    gene.grid <- seq_along(genes)
    names(gene.grid) <- genes
    gene.int.map <- gene.grid[p2g]

    .summarizeP2G <-
        function(g)
        {
            curr.probes <- which(gene.int.map == g)
            curr.exprs <- probe.exprs[curr.probes,]
            if(is.matrix(curr.exprs))
                curr.exprs <- colMeans(curr.exprs, na.rm=TRUE)
            return(curr.exprs)
        }


    gene.exprs <- vapply(gene.grid, .summarizeP2G, numeric(ncol(se))) 
    gene.exprs <- t(gene.exprs)    

    colnames(gene.exprs) <- colnames(se)
    rownames(gene.exprs) <- genes

    # create new SE
    geneSE <- SummarizedExperiment(assays=list(exprs=gene.exprs), 
        colData=colData(se), metadata=metadata(se))  

    return(geneSE)
}

.annoP2G <- function(se) 
{
    PRB.COL <- configEBrowser("PRB.COL")
    EZ.COL <- configEBrowser("EZ.COL")

    is.se <- is(se, "SummarizedExperiment")
    if(is.se) anno <- metadata(se)$annotation
    else anno <- se
    anno.pkg <- paste0(anno, ".db")
    
    # check whether annotation package is installed
    isAvailable(anno.pkg)
    anno.pkg <- get(anno.pkg)  
 
    # determine mapping
    org <- AnnotationDbi::species(anno.pkg) 
    data(korg, package="pathview")
    suppressMessages(org <- pathview::kegg.species.code(org))
    org.start <- paste0("^", org, ":")
    anno.keys <- AnnotationDbi::keys(anno.pkg)
    # check whether there is a mapping to Entrez (= NCBI gene-ID)
    if(EZ.COL %in% AnnotationDbi::keytypes(anno.pkg))
    {
        p2g.map <- suppressMessages( AnnotationDbi::mapIds(anno.pkg, 
            keys=anno.keys, keytype=PRB.COL, column=EZ.COL) )
        probes <- names(p2g.map)

        # check whether EntrezID is used by KEGG
#        first.nna <- p2g.map[!is.na(p2g.map)][1]
#        first.keggid <- KEGGREST::keggConv("genes", 
#            paste("ncbi-geneid", first.nna, sep=":"))
#        kegg.uses.entrez <- sub(org.start, "", first.keggid) == first.nna

        # otherwise convert from EntrezID to KEGGID
#        if(!kegg.uses.entrez)
#        {
#            message("KEGG does not use EntrezIDs for your organism")
#            message("Downloading mapping Entrez -> KEGG & converting accordingly")
#            map.e2k <- KEGGREST::keggConv(org, "ncbi-geneid")
#            names(map.e2k) <- sub("^ncbi-geneid:", "", names(map.e2k))
#            map.e2k <- sub(org.start, "", map.e2k)
#            p2g.map <- map.e2k[p2g.map]
#            names(p2g.map) <- probes
#        }     
    }
    else
    {
        message(paste("Did not found mapping: probe -> EntrezID in", anno, ".db"))
        message("Try to find mapping: probe -> KEGGID")
        avail.maps <-  AnnotationDbi::keytypes(anno.pkg)
        kegg.ids <- KEGGREST::keggLink(paste0("path:", org, "00010"))
        kegg.ids <- grep(org.start, kegg.ids[,2], value=TRUE)[1:3]
        kegg.ids <- sub(org.start, "", kegg.ids)
        is.map <- sapply(avail.maps, 
            function(m)
            {
                x <- AnnotationDbi::mapIds(anno.pkg, 
                    keys=anno.keys, keytype=PRB.COL, column=m)
                return(all(kegg.ids %in% x))
           })
        if(!any(is.map)) stop("Found no suitable mapping")
        rel.map <- avail.maps[which(is.map)[1]]
        p2g.map <- AnnotationDbi::mapIds(anno.pkg, 
            keys=anno.keys, keytype=PRB.COL, column=rel.map)
    }
    if(is.se)
    {
        p2g.map <- p2g.map[names(se)]
        rowData(se)[,PRB.COL] <- names(p2g.map)
        rowData(se)[,EZ.COL] <- unname(p2g.map)
        metadata(se)$annotation <- org
        return(se)
    }
    fDat <- data.frame(names(p2g.map), unname(p2g.map), stringsAsFactors=FALSE)
    colnames(fDat) <- c("PRB.COL", "EZ.COL")
    return(fDat)
}

# converts from an annotation package ID to an organism ID
# e.g. hgu95av2.db -> hsa
.annoPkg2Org <- function(anno.pkg)
{
    org <- AnnotationDbi::species(anno.pkg) 
    data(korg, package="pathview")
    suppressMessages(org <- pathview::kegg.species.code(org))
    return(unname(org))
}   

