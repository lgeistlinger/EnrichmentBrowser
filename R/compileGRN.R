############################################################
#
# author: Ludwig Geistlinger
# date: 19 April 2011
#
# Compilation of GRNs from databases and literature
#
############################################################

GRN.HEADER.COLS <- c("FROM", "TO", "TYPE", "REL", "PWY")

compile.grn.from.kegg <- function(pwys, out.file=NULL)
{
    kegg.rels <- unique(get.kegg.rels.of.organism(pwys)[,1:3])
    kegg.rels[,"TYPE"] <- ifelse(kegg.rels[,"TYPE"] == "-->", "+", "-") 
    if(is.null(out.file)) return(kegg.rels)
    write.table(kegg.rels, file=out.file, 
                row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
    message(paste("GRN written to", out.file))
}

# rels %in% c("ECrel", "GErel", "PCrel", "PPrel")
get.kegg.rels.of.organism <- function(  pwys, 
                                        out.file=NULL, 
                                        types=c("-->", "--|"), 
                                        rels=c("GErel", "PPrel"))
{
    if(class(pwys) == "character")
    {
        if(nchar(pwys) == 3) pwys <- download.kegg.pathways(pwys)
        else pwys <- extract.pwys(pwys)
    }

    org <- getPathwayInfo(pwys[[1]])@org
    no.out <- FALSE
    if(is.null(out.file))
    {
        no.out <- TRUE
        out.file <- file.path(system.file(package="EnrichmentBrowser"), 
                                paste(org, "_rels.txt", sep=""))
    }

    if(file.exists(out.file)) file.remove(out.file)
    con <- file(out.file, open="at")
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
                    et <- get.edge.type(e)
                    rt <- getType(e)

                    if((et %in% types) && (rt %in% rels))
                    {
                        entries <- getEntryID(e)
                        ids1 <- get.node.name(entries[1], nodes(p))
                        ids2 <- get.node.name(entries[2], nodes(p))
                        sapply(ids1, function(i) 
                            sapply(ids2, function(j)
                                writeLines(
                                    paste(c(i,j,et,rt,nr), collapse="\t"), con)))
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

get.edge.type <- function(edge)
        ifelse(length(getSubtype(edge)) == 0,
                        NA,
                        getSubtype(edge)$subtype@value)

get.node.name <- function(entry, nodes)
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
