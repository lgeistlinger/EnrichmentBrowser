############################################################
#
# author: Ludwig Geistlinger
# date: 19 April 2011
#
# Compilation of GRNs from databases and literature
#
############################################################

#' Compilation of a gene regulatory network from pathway databases
#' 
#' To perform network-based enrichment analysis a gene regulatory network (GRN)
#' is required. There are well-studied processes and organisms for which
#' comprehensive and well-annotated regulatory networks are available, e.g. the
#' RegulonDB for E. coli and Yeastract for S. cerevisiae.  However, in many
#' cases such a network is missing.  A first simple workaround is to compile a
#' network from regulations in pathway databases such as KEGG.
#' 
#' 
#' @aliases compile.grn.from.kegg
#' @param org An organism in KEGG three letter code, e.g. \sQuote{hsa} for
#' \sQuote{Homo sapiens}.  Alternatively, and mainly for backward
#' compatibility, this can also be either a list of
#' \code{\linkS4class{KEGGPathway}} objects or an absolute file path of a zip
#' compressed archive of pathway xml files in KGML format.
#' @param db Pathway database.  This should be one or more DBs out of 'kegg',
#' 'reactome', 'pathbank', and 'wikipathways'.  See \code{\link{pathwayDatabases}} for
#' available DBs of the respective organism.  Default is 'kegg'.
#' Note: when dealing with non-model organisms, GRN compilation is currently
#' only supported directly from KEGG (the argument \code{kegg.native} should
#' accordingly be set to \code{TRUE}).
#' @param act.inh Should gene regulatory interactions be classified as
#' activating (+) or inhibiting (-)?  If TRUE, this will drop interactions for
#' which such a classification cannot be made (e.g. binding events).
#' Otherwise, all interactions found in the pathway DB will be included.
#' Default is \code{TRUE}.
#' @param map2entrez Should gene identifiers be mapped to NCBI Entrez Gene IDs?
#' This only applies to Reactome and PathBank as they both use UNIPROT IDs.
#' This is typically recommended when using the GRN for network-based enrichment
#' analysis with the EnrichmentBrowser.  Default is \code{TRUE}.
#' @param keep.type Should the original interaction type descriptions be kept?
#' If TRUE, this will keep the long description of interaction types as found
#' in the original KGML and BioPax pathway files.  Default is \code{FALSE}.
#' @param kegg.native For KEGG: should the GRN be compiled from the native KGML
#' files or should graphite's pathway topology conversion be used?  See the
#' vignette of the graphite package for details.  This is mostly for backward
#' compatibility.  Default is \code{FALSE}. Note: when dealing with non-model 
#' organisms (not supported by graphite) this argument should be set to 
#' \code{TRUE}.
#' @return The GRN in plain matrix format.  Two columns named \code{FROM} (the
#' regulator) and \code{TO} (the regulated gene) are guaranteed.  Additional
#' columns, named \code{TYPE} and \code{LONG.TYPE}, are included if option
#' \code{act.inh} or \code{keep.type} is activated.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{\link{pathwayDatabases}}, \code{\link{pathways}},
#' \code{\linkS4class{KEGGPathway}}, \code{\link{parseKGML}},
#' \code{\link{downloadPathways}}
#' @examples
#' 
#'     kegg.grn <- compileGRN(org="hsa", db="kegg")
#' 
#' @export compileGRN
compileGRN <- function(org, db="kegg",
    act.inh=TRUE, map2entrez=TRUE, keep.type=FALSE, kegg.native=FALSE)
{
    kegg.native <- length(db) == 1 && db == "kegg" && kegg.native
    if(kegg.native) grn <- .compileGRNFromKEGG(org) 
    else grn <- .compileGRNFromGraphite(org, db=db, act.inh=act.inh, 
                    map2entrez=map2entrez, keep.type=keep.type)
    return(grn)
}    

.compileGRNFromGraphite <- function(org, 
    db="kegg", act.inh=TRUE, 
    map2entrez=TRUE, keep.type=FALSE)
{
    valid.dbs <- c("kegg", "reactome", "pathbank", "wikipathways")
    if(!all(db %in% valid.dbs)) 
        stop(paste("Valid values of \'db\':", paste(valid.dbs, collapse=", ")))

    #isAvailable("graphite", type="software")
    all.dbs <- as.matrix(graphite::pathwayDatabases())

    # valid org?
    org.ind <- grep(org, all.dbs[,"species"])
    if(!length(org.ind)) 
        stop(paste0("No graphite data for organism \"", org, "\"")) 

    # which DBs are available for this org?
    gorg <- all.dbs[org.ind[1], "species"]
    org.dbs <- all.dbs[org.ind, "database"]

    # any invalid user dbs?
    na.dbs <- setdiff(db, org.dbs)
    if(length(na.dbs)) for(n in na.dbs) warning(paste("Invalid DB:", n)) 
    dbs <- intersect(org.dbs, db)

    # construct
    edge.coll <- lapply(dbs, function(d) 
        .processDB(graphite::pathways(gorg, d), act.inh, map2entrez, keep.type))
    grn <- do.call(rbind, edge.coll)
    grn <- unique(grn)
    grn <- grn[do.call(order, as.data.frame(grn)),]
    return(grn)
}    

#pdbs <- c("kegg", "reactome", "pathbank", "panther", "wikipathways")
#inspectPDBs <- function(db)
#{
#    message(db)
#    pwys <- pathways("hsapiens", db)
#    db.edges <- lapply(pwys, edges)
#    db.edges <- do.call(rbind, db.edges)
#    db.edges <- as.matrix(db.edges)
#    db.edges <- unique(db.edges)
#    return(table(db.edges[,"type"]))
#}
#lapply(pdbs, inspectPDBs)


.processDB <- function(pwys, act.inh=TRUE, map2entrez=TRUE, keep.type=FALSE)
{    
    d <- tolower(pwys@name)
    gorg <- pwys@species
    db.edges <- lapply(pwys, edges)
    db.edges <- do.call(rbind, db.edges)
    db.edges <- as.matrix(db.edges)
    db.edges <- unique(db.edges)
    
    if(act.inh) db.edges <- .filterEdges(db.edges, db=d, keep.type=keep.type)
    else
    {
        rel.cols <- c("src", "dest")
        if(keep.type) rel.cols <- c(rel.cols, "type") 
        db.edges <- db.edges[,rel.cols]
        cnames <- c("FROM", "TO")
        if(keep.type) cnames[3] <- "TYPE"
        colnames(db.edges) <- cnames
    }    

    # map to ENTREZ IDs necessary?
    if(!(d %in% c("kegg", "wikipathways")) && map2entrez)
    {
        org <- tolower(substring(gorg, 1, 3))
        ids <- unique(as.vector(db.edges[,1:2]))
        x <- .idmap(ids, org, from="UNIPROT", to="ENTREZID")
        for(i in 1:2) db.edges[,i] <- x[db.edges[,i]]
        ind.na <- apply(db.edges, 1, function(e) any(is.na(e)))
        db.edges <- db.edges[!ind.na,]
        db.edges <- unique(db.edges)
    }    
    rownames(db.edges) <- NULL
    return(db.edges)
}

.filterEdges <- function(db.edges, db, keep.type=FALSE)
{
    if(db=="kegg")
    {
        rel.types <- c("activation", "expression", "inhibition", "repression")
        rel.types <- paste0("Process(", rel.types, ")")
    }
    # else: reactome
    else if(db == "reactome")
    {
        rel.types <- paste0( "Control(Out: ",
                                rep(c("ACTIVATION", "INHIBITION"), each=3),
                                " of ",
                                rep(c("BiochemicalReaction", "TemplateReaction",
                                        "ACTIVATION"), 2),
                                ")")
    }
    
    # humancyc: no information on act/inh
    # panther:
    # no information on act/inh 
    # ID mapping UNIPROT -> ENTREZID fails greatly
    
    db.edges <- db.edges[db.edges[,"type"] %in% rel.types,]
    n <- length(rel.types) / 2
    ind <- seq_len(n)
    ltype <- db.edges[,"type"]
    db.edges[,"type"] <- ifelse(ltype %in% rel.types[ind], "+", "-")
    db.edges <- db.edges[, c("src", "dest", "type")]
    colnames(db.edges) <- c("FROM", "TO", "TYPE")
    if(keep.type)
    {
        db.edges <- cbind(db.edges, ltype)
        colnames(db.edges)[ncol(db.edges)] <- "LONG.TYPE" 
    }
    db.edges <- unique(db.edges)
    return(db.edges)
}

.compileGRNFromKEGG <- function(pwys, out.file=NULL)
{
    kegg.rels <- unique(.getKEGGRels(pwys)[,1:3])
    kegg.rels[,"TYPE"] <- ifelse(kegg.rels[,"TYPE"] == "-->", "+", "-") 
    if(is.null(out.file)) return(kegg.rels)
    write.table(kegg.rels, file=out.file, 
                row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
    message(paste("GRN written to", out.file))
}

# rels %in% c("ECrel", "GErel", "PCrel", "PPrel")
.getKEGGRels <- function(  pwys, 
                           out.file=NULL, 
                           types=c("-->", "--|"), 
                           rels=c("GErel", "PPrel"))
{
    if(is.character(pwys))
    {
        if(nchar(pwys) %in% c(3,4)) pwys <- downloadPathways(pwys)
        else pwys <- .extractPwys(pwys)
    }

    org <- getPathwayInfo(pwys[[1]])@org
    no.out <- FALSE
    if(is.null(out.file))
    {
        no.out <- TRUE
        out.dir <- configEBrowser("OUTDIR.DEFAULT")
        if(!file.exists(out.dir)) dir.create(out.dir, recursive=TRUE)
        out.file <- file.path(out.dir, paste(org, "rels.txt", sep="_"))
    }

    if(file.exists(out.file)) file.remove(out.file)
    con <- file(out.file, open="at")
    GRN.HEADER.COLS <- c("FROM", "TO", "TYPE", "REL", "PWY")
    writeLines(paste(GRN.HEADER.COLS, collapse="\t"), con)

    for(p in pwys)
    {
        nr <- getPathwayInfo(p)@number
        for(e in edges(p))
        {  
            # effect type: -->, --|, ...
            # relation type: GErel, PPrel, ...
            rt <- getType(e)
            if(rt %in% rels)
            {
                et <- .getEdgeType(e)
                if(et %in% types)
                {
                    entries <- getEntryID(e)
                    ids1 <- .getNodeName(entries[1], nodes(p))
                    ids2 <- .getNodeName(entries[2], nodes(p))
                    for(i in ids1) 
                        for(j in ids2)
                            writeLines(
                                paste(c(i,j,et,rt,nr), collapse="\t"), 
                            con)
                        
                }
            }
        }
        flush(con)
    }

    close(con)
    cont <- as.matrix(read.delim(out.file))
    cont <- gsub(" ", "", cont)
    cont <- unique(cont)
    if(no.out)
    { 
        file.remove(out.file)
        return(cont)
    }
    write.table(cont, file=out.file, sep="\t", row.names=FALSE, quote=FALSE)
    message(paste(org, "KEGG relations written to", out.file))
}

.getEdgeType <- function(edge)
{
    if(length(getSubtype(edge)) == 0) return(NA)
    return(getSubtype(edge)$subtype@value)
}

.getNodeName <- function(entry, nodes)
{
    n <- nodes[[entry]]
    ids <- getName(n)
    if((length(ids) == 1) && (ids == "undefined"))
    {
        ids <- NULL 
        ncomp <- getComponent(n)
        ids <- unlist(lapply(ncomp, function(comp) getName(nodes[[comp]])))
    }
    ids <- sub("^[a-z]{3,4}:", "", ids)
    return(ids)
}
