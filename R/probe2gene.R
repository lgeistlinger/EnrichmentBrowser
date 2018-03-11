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

# fast conversion of probe 2 gene expression


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
#' Alternatively, The mapping from probe to gene can also be defined in the
#' \code{\link{rowData}} slot via two columns named (i) \code{PROBEID} for the
#' platform-specific probe ID, and (ii) \code{ENTREZID} for the corresponding
#' NCBI Entrez Gene ID.
#' @param use.mean Logical.  Determining, in case of multiple probes for one
#' gene, whether a mean value is computed (\code{use.mean=TRUE}), or the probe
#' that discriminate the most between the two sample group is kept
#' (\code{use.mean=FALSE}).  Defaults to TRUE.
#' @return A \code{\linkS4class{SummarizedExperiment}} on gene level.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{\link{readSE}} for reading expression data from file,
#' \code{\link{deAna}} for differential expression analysis.
#' @examples
#' 
#'     # (1) reading the expression data from file
#'     exprs.file <- system.file("extdata/exprs.tab", package="EnrichmentBrowser")
#'     pdat.file <- system.file("extdata/colData.tab", package="EnrichmentBrowser")
#'     fdat.file <- system.file("extdata/rowData.tab", package="EnrichmentBrowser")
#'     probeSE <- readSE(exprs.file, pdat.file, fdat.file)
#'     geneSE <- probe2gene(probeSE) 
#' 
#' @export probe2gene
probe2gene <- function(probeSE, use.mean=TRUE)
{
    ### TEMPORARY: will be replaced by as(eSet,SummarizedExperiment)
    if(is(probeSE, "ExpressionSet")) 
        probeSE <- as(probeSE, "RangedSummarizedExperiment")
    ### 

    EZ.COL <- config.ebrowser("EZ.COL")
    se <- probeSE
    if(!(EZ.COL %in% colnames(rowData(se)))) se <- .annoP2G(se)
    
    # remove probes without gene annotation
    not.na <- !is.na(rowData(se)[, EZ.COL])
    if(sum(not.na) < nrow(se)) se <- se[not.na,]

    probe.exprs <- assay(se)
    p2g <- as.vector(rowData(se)[, EZ.COL])
    
    # determine unique genes
    genes <- unique(p2g)
    
    # compute gene expression
    gene.grid <- seq_along(genes)
    names(gene.grid) <- genes
    gene.int.map <- gene.grid[p2g]
    gene.exprs <- t(sapply(gene.grid, 
        function(g)
        {
            curr.probes <- which(gene.int.map == g)
            curr.exprs <- probe.exprs[curr.probes,]
            if(is.matrix(curr.exprs))
            { 
                if(use.mean) curr.exprs <- colMeans(curr.exprs, na.rm=TRUE)
                else
                { 
                    FC.COL <- config.ebrowser("FC.COL") 
                    if(!(FC.COL %in% colnames(rowData(se))))
                        stop(paste("use.mean=FALSE, but did not find differential", 
                            "expression in rowData. Run deAna first."))
                    curr.exprs <- curr.exprs[
                        which.max(abs(rowData(se)[curr.probes, FC.COL])),]
                }
            }
            return(curr.exprs)
        })) 
    colnames(gene.exprs) <- colnames(se)
    rownames(gene.exprs) <- genes

    # create new SE
    geneSE <- SummarizedExperiment(assays=list(exprs=gene.exprs), 
        colData=colData(se), metadata=metadata(se))  
    
    return(geneSE)
}

#' @export
#' @keywords internal
probe.2.gene.eset <- function(probe.eset, use.mean=TRUE)
{
    .Deprecated("probe2gene")
    probe2gene(probe.eset, use.mean=use.mean)
}

.annoP2G <- function(se) 
{
    PRB.COL <- config.ebrowser("PRB.COL")
    EZ.COL <- config.ebrowser("EZ.COL")

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
    # check whether there is a mapping to Entrez (= NCBI gene-ID)
    if(EZ.COL %in% keytypes(anno.pkg))
    {
        p2g.map <- suppressMessages(mapIds(anno.pkg, 
            keys=keys(anno.pkg), keytype=PRB.COL, column=EZ.COL))
        probes <- names(p2g.map)

        # check whether EntrezID is used by KEGG
        first.nna <- p2g.map[!is.na(p2g.map)][1]
        first.keggid <- KEGGREST::keggConv("genes", 
            paste("ncbi-geneid", first.nna, sep=":"))
        kegg.uses.entrez <- sub(org.start, "", first.keggid) == first.nna

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
        avail.maps <-  keytypes(anno.pkg)
        kegg.ids <- KEGGREST::keggLink(paste0("path:", org, "00010"))
        kegg.ids <- grep(org.start, kegg.ids[,2], value=TRUE)[1:3]
        kegg.ids <- sub(org.start, "", kegg.ids)
        is.map <- sapply(avail.maps, 
            function(m)
            {
                x <- mapIds(anno.pkg, 
                    keys=keys(anno.pkg), keytype=PRB.COL, column=m)
                return(all(kegg.ids %in% x))
           })
        if(!any(is.map)) stop("Found no suitable mapping")
        rel.map <- avail.maps[which(is.map)[1]]
        p2g.map <- mapIds(anno.pkg, 
            keys=keys(anno.pkg), keytype=PRB.COL, column=rel.map)
    }
    if(is.se)
    {
        p2g.map <- p2g.map[rownames(se)]
        rowData(se)[,PRB.COL] <- names(p2g.map)
        rowData(se)[,EZ.COL] <- p2g.map
        metadata(se)$annotation <- org
        return(se)
    }
    fDat <- cbind(names(p2g.map), p2g.map)
    return(fDat)
}

# converts from an annotation package ID to an organism ID
# e.g. hgu95av2.db -> hsa
.annoPkg2Org <- function(anno.pkg)
{
    info <- suppressMessages(capture.output(get(anno.pkg)))
    org.line <- grep(" ORGANISM: ", info, value=TRUE)
    org <- tolower(unlist(strsplit(org.line, ": "))[2])
    org.spl <- unlist(strsplit(org, " "))
    org.id <- paste0(substring(org.spl[1],1,1), substring(org.spl[2],1,2))
    return(org.id)
}

