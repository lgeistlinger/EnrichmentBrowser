###############################################################################
# 
# Author: ludwig geistlinger
# Date: 16 June 2010
#
# script for extracting genesets from kgml files (KEGG XML)
# of pathways. The genesets are subsequently written in gmt format
# that is suitable for input to GSEA.
#
# Update, 05 May 2014:  easy & fasr  geneset getter for organism of choice 
#           based on KEGGREST functionality 
#
#
###############################################################################

get.kegg.genesets <- function(pwys, gmt.file=NULL)
{
    if(class(pwys) == "character")
    {
        if(nchar(pwys) == 3)
            return(download.kegg.genesets(org=pwys, gmt.file=gmt.file))
        else pwys <- extract.pwys(pwys)
    }
    return(extract.kegg.genesets(pwys, gmt.file=gmt.file))
}

# 
download.kegg.genesets <- function(org, gmt.file=NULL)
{
    pwys <- keggList("pathway", org)
    pwy2gene <- keggLink(org, "pathway")

    gs <- sapply(names(pwys), 
            function(pwy)
            { 
                genes <- pwy2gene[names(pwy2gene) == pwy]
                genes <- sub("^[a-z]{3}:", "", genes)
                genes <- sort(genes)
                names(genes) <- NULL
                return(genes)
            }, simplify=FALSE)

    names(gs) <- make.gs.names(names(pwys), pwys)
    
    if(!is.null(gmt.file)) write.gmt(gs, gmt.file=gmt.file)
    return(gs)
}

# only preferred over 'download.kegg.genesets' 
# when pathway kgmls have already been download
extract.kegg.genesets <- function(pwys, gmt.file=NULL)
{
    # read in & parse pathways
    if(class(pwys) == "character") pwys <- extract.pwys(pwys)
    
    # get pathway annotations
    nn <- sapply(pwys, getName)
    tt <- sapply(pwys, getTitle)
    
    # extract genesets
    gs <- sapply(pwys, 
        function(pwy)
        {
            genes <- get.genes.by.pathway(pwy)
            genes <- sub("^[a-z]{3}:", "", genes)
            genes <- sort(genes)
            return(genes)
        }, simplify=FALSE)

    names(gs) <- make.gs.names(nn, tt)
    
    if(!is.null(gmt.file)) write.gmt(gs, gmt.file=gmt.file)
    return(gs)
}

## extract pwys from zip archive and parse KGML files
extract.pwys <- function(pwy.zip)
{
    pwy.dir <- dirname(pwy.zip)
    unzip(pwy.zip, exdir=pwy.dir, junkpaths=TRUE)
    pwy.files <- list.files(pwy.dir, pattern="*.xml", full.names=TRUE)
    pwys <- sapply(pwy.files, parseKGML)
    ## clean up
    sapply(pwy.files, file.remove)
    return(pwys)
}

get.genes.by.pathway <- function(pwy)
{
    ts <- sapply(nodes(pwy), getType)
    genes <- unique(unlist(sapply(nodes(pwy)[ts == "gene"], getName)))
    return(genes)
}

# build first gmt column: the ID (format: <pwy.nr>_<pwy.title>)
make.gs.names <- function(ids, titles)
{
    ids <- sub("path:", "", ids)
    titles <- sapply(titles, function(title) unlist(strsplit(title, " - "))[1])
    titles <- str_trim(titles)
    titles <- gsub(" ", "_", titles)
    ids <- paste(ids, titles, sep="_")
}

# write genesets to file in GMT format
write.gmt <- function(gs, gmt.file)
{
    ## collapse geneset members to one tab separated string 
    gs.strings <- sapply(gs, function(x) paste(x,collapse="\t"))
    
    ## paste an not annotated second column (artifact of gmt format)
    ann <- paste(names(gs), rep(NA,length(gs)), sep="\t")
    
    ## paste all together
    all <- paste(ann, gs.strings, sep="\t")
    
    ## collapse all together to a single newline separated string
    all.str <- paste(all, collapse="\n")
    all.str <- paste(all, "\n", sep="")
    
    ## write the gs in gmt format
    cat(all.str, file=gmt.file, sep="")
}

## parse geneset database
parse.genesets.from.GMT <- function(gmt.file)
{
    content <- readLines(gmt.file, warn=FALSE)
    le <- length(content)
    genesets <- vector("list", length=le)
    gs.names <- vector("character", length=le)
    for(i in seq_len(le))
    {
        line <- content[i]
        spl <- unlist(strsplit(line, "\t"))
        gs.names[i] <- spl[1]
        genesets[[i]] <- spl[-c(1,2)]
    }
    names(genesets) <- gs.names
    return(genesets)
}


