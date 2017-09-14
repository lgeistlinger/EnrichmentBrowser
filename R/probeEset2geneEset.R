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


### FUNCTIONS

# fast conversion of probe 2 gene expression
probe.2.gene.eset <- function(probe.eset, use.mean=TRUE)
{
    ### TEMPORARY: will be replaced by as(eSet,SummarizedExperiment)
    if(is(probe.eset, "ExpressionSet")) 
        probe.eset <- as(probe.eset, "RangedSummarizedExperiment")
    ### 

    EZ.COL <- config.ebrowser("EZ.COL")
    eset <- probe.eset
    if(!(EZ.COL %in% colnames(rowData(eset)))) eset <- .annoP2G(eset)
    
    # remove probes without gene annotation
    not.na <- !is.na(rowData(eset)[, EZ.COL])
    if(sum(not.na) < nrow(eset)) eset <- eset[not.na,]

    probe.exprs <- assay(eset)
    probe2gene <- as.vector(rowData(eset)[, EZ.COL])
    
    # determine unique genes
    genes <- unique(probe2gene)
    
    # compute gene expression
    gene.grid <- seq_len(length(genes))
    names(gene.grid) <- genes
    gene.int.map <- gene.grid[probe2gene]
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
                    if(!(FC.COL %in% colnames(rowData(eset))))
                        stop(paste("use.mean=FALSE, but did not find differential", 
                            "expression in rowData. Run de.ana first."))
                    curr.exprs <- curr.exprs[
                        which.max(abs(rowData(eset)[curr.probes, FC.COL])),]
                }
            }
            return(curr.exprs)
        })) 
    colnames(gene.exprs) <- colnames(eset)
    rownames(gene.exprs) <- genes

    # create new eset
    gene.eset <- SummarizedExperiment(assays=list(exprs=gene.exprs), 
        colData=colData(eset), metadata=metadata(eset))  
    
    return(gene.eset)
}

.annoP2G <- function(eset) 
{
    PRB.COL <- config.ebrowser("PRB.COL")
    EZ.COL <- config.ebrowser("EZ.COL")

    is.eset <- is(eset, "SummarizedExperiment")
    if(is.eset) anno <- metadata(eset)$annotation
    else anno <- eset
    anno.pkg <- paste0(anno, ".db")
    
    # check whether annotation package is installed
    .isAvailable(anno.pkg)
    anno.pkg <- get(anno.pkg)  
 
    # determine mapping
    org <- AnnotationDbi::species(anno.pkg) 
    data(korg, package="pathview")
    org <- pathview::kegg.species.code(org)
    org.start <- paste0("^", org, ":")
    # check whether there is a mapping to Entrez (= NCBI gene-ID)
    if(EZ.COL %in% keytypes(anno.pkg))
    {
        p2g.map <- suppressMessages(mapIds(anno.pkg, 
            keys=keys(anno.pkg), keytype=PRB.COL, column=EZ.COL))
        probes <- names(p2g.map)

        # check whether EntrezID is used by KEGG
        first.nna <- p2g.map[!is.na(p2g.map)][1]
        first.keggid <- keggConv("genes", 
            paste("ncbi-geneid", first.nna, sep=":"))
        kegg.uses.entrez <- sub(org.start, "", first.keggid) == first.nna

        # otherwise convert from EntrezID to KEGGID
#        if(!kegg.uses.entrez)
#        {
#            message("KEGG does not use EntrezIDs for your organism")
#            message("Downloading mapping Entrez -> KEGG & converting accordingly")
#            map.e2k <- keggConv(org, "ncbi-geneid")
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
        kegg.ids <- keggLink(paste0("path:", org, "00010"))
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
    if(is.eset)
    {
        p2g.map <- p2g.map[rownames(eset)]
        rowData(eset)[,PRB.COL] <- names(p2g.map)
        rowData(eset)[,EZ.COL] <- p2g.map
        metadata(eset)$annotation <- org
        return(eset)
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

