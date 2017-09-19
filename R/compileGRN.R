############################################################
#
# author: Ludwig Geistlinger
# date: 19 April 2011
#
# Compilation of GRNs from databases and literature
#
############################################################


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

    sapply(pwys, 
        function(p)
        {
            nr <- getPathwayInfo(p)@number
            sapply(edges(p),
                function(e)
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
                            sapply(ids1, function(i) 
                                sapply(ids2, function(j)
                                    writeLines(
                                        paste(c(i,j,et,rt,nr), collapse="\t"), con)))
                        }
                    }
                })
            flush(con)
        })

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
        ids <- unlist(sapply(ncomp, function(comp) getName(nodes[[comp]])))
    }
    ids <- sub("^[a-z]{3}:", "", ids)
    return(ids)
}
