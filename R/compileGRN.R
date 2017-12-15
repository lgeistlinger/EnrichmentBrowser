############################################################
#
# author: Ludwig Geistlinger
# date: 19 April 2011
#
# Compilation of GRNs from databases and literature
#
############################################################

compileGRN <- function(org, mode=c("kegg", "graphite"), act.inh=TRUE)
{
    mode <- match.arg(mode)    
    if(mode == "kegg") grn <- compile.grn.from.kegg() 
}    

compileGRNFromGraphite <- function(org, dbs=c("kegg", "reactome"), act.inh=TRUE)
{
    pathways <- pathwayDatabases <- NULL
    .isAvailable("graphite", type="software")

    all.dbs <- as.matrix(pathwayDatabases())

    # valid org?
    org.ind <- grep(org, all.dbs[,"species"])
    if(!length(org.ind)) 
        stop(paste0("No graphite data for organism \"", org, "\"")) 

    # which DBs are available for this org?
    gorg <- all.dbs[org.ind[1], "species"]
    org.dbs <- all.dbs[org.ind, "database"]

    db.csv <- paste(org.dbs[seq_len(length(org.dbs) - 1)], collapse=", ")
    db.csv <- paste(db.csv, org.dbs[length(org.dbs)], sep=" and ")
    message(paste("Found data from", db.csv, "for", gorg, "in graphite"))

    # any invalid user dbs?
    na.dbs <- setdiff(dbs, org.dbs)
    if(length(na.dbs)) for(n in na.dbs) warning(paste("Ignoring", n)) 

    dbs <- intersect(org.dbs, dbs)
    db.csv <- paste(dbs[seq_len(length(dbs) - 1)], collapse=", ")
    db.csv <- paste(db.csv, dbs[length(dbs)], sep=" and ")
    message(paste("Constructing GRN from", db.csv))

    # construct
    edge.coll <- lapply(dbs, function(d) .processDB(d, gorg, act.inh))
    grn <- do.call(rbind, edge.coll)
    grn <- unique(grn)
    return(grn)
}    

.processDB <- function(d, gorg, act.inh)
{    
    message(paste("Processing", d, "...")) 
    pwys <- pathways(gorg, d)
    db.edges <- lapply(pwys, edges)
    db.edges <- do.call(rbind, db.edges)
    db.edges <- as.matrix(db.edges)
    db.edges <- unique(db.edges)
    
    if(act.inh) db.edges <- .filterEdges(db.edges, db=d, act.inh=act.inh)
    else
    {
        db.edges <- db.edges[,c("src", "dest")]
        db.edges <- unique(db.edges)
    }    

    colnames(db.edges)[1:2] <- c("FROM", "TO")
    colnames(db.edges) <- toupper(colnames(db.edges))

    # map to ENTREZ IDs necessary?
    if(!(d %in% c("kegg", "biocarta")))
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

.filterEdges <- function(db.edges, db, act.inh)
{
    if(db=="kegg")
    {
        rel.types <- c("activation", "inhibition", "expression", "repression")
        rel.types <- paste0("Process(", rel.types, ")")
        db.edges <- db.edges[db.edges[,"type"] %in% rel.types,]
        db.edges[,"type"] <- ifelse(db.edges[,"type"] %in% rel.types[c(1,3)], "+", "-")
    }
    else if(db=="reactome")
    {
            
    }
    else if(db=="biocarta")
    {
            
    }    
    else if(db=="humancyc")
    {
        # no information on act/inh
    }    
    else if(db=="nci")
    {

    }
    # panther
    else
    {
        # no information on act/inh 
        # ID mapping UNIPROT -> ENTREZID fails greatly
    }    
    db.edges <- db.edges[, c("src", "dest", "type")]
    db.edges <- unique(db.edges)
    return(db.edges)
}

compile.grn.from.kegg <- function(pwys, out.file=NULL)
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
        if(nchar(pwys) == 3) pwys <- download.kegg.pathways(pwys)
        else pwys <- .extractPwys(pwys)
    }

    org <- getPathwayInfo(pwys[[1]])@org
    no.out <- FALSE
    if(is.null(out.file))
    {
        no.out <- TRUE
        out.dir <- config.ebrowser("OUTDIR.DEFAULT")
        if(!file.exists(out.dir)) dir.create(out.dir)
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
    ids <- sub("^[a-z]{3}:", "", ids)
    return(ids)
}
