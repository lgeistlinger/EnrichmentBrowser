###############################################################################
# 
# Author: ludwig geistlinger
# Date: 16 June 2010
#
# script for extracting genesets from kgml files (KEGG XML)
# of pathways. The genesets are subsequently written in gmt format
# that is suitable for input to GSEA.
#
# Update, 05 May 2014:  easy & fast geneset getter for organism of choice 
#           based on KEGGREST functionality 
#
# Update, 10 March 2015: including getter for GO genesets
###############################################################################

#
# (1) GO
#
get.go.genesets <- function(org, 
    onto=c("BP", "MF", "CC"), mode=c("GO.db","biomart"))
{
    onto <- match.arg(onto)
    mode <- match.arg(mode)

    if(mode=="GO.db")
    {
        gs <- topGO::annFUN.org(whichOnto=onto, mapping=.org2pkg(org))
        GO2descr <- AnnotationDbi::as.list(GO.db::GOTERM)
        gs <- gs[intersect(names(gs), names(GO2descr))]
        GO2descr <- GO2descr[names(gs)]
        GO2descr <- sapply(GO2descr, AnnotationDbi::Term)
        names(gs) <- paste(names(gs), gsub(" ", "_", GO2descr), sep="_")
    }
    else
    {
        # setting mart
        ensembl <- biomaRt::useMart("ensembl")
        ds <- biomaRt::listDatasets(ensembl)[,"dataset"]
        ds <- grep(paste0("^", org), ds, value=TRUE)
        ensembl <- biomaRt::useDataset(ds, mart=ensembl)

        message("Downloading mapping from BioMart ...")
        message("This may take a few minutes ...")
        
        GO2descr <- biomaRt::getBM(attributes=
            c("go_id", "name_1006", "namespace_1003"), mart=ensembl)
        GO2descr <- GO2descr[GO2descr$go_id != "", ]
        GO2descr <- GO2descr[order(GO2descr[,"go_id"]),]
        ontos <- sapply(GO2descr[,3], 
            function(x)
            {
                spl <- unlist(strsplit(x, "_"))
                onto <- paste(substring(spl,1,1), collapse="")         
                return(toupper(onto))
            }) 
        GO2descr <- GO2descr[ontos==onto,1:2]
        
        gene2GO <- biomaRt::getBM(attributes = c("entrezgene", "go_id"), mart=ensembl)
        gene2GO <- gene2GO[apply(gene2GO, 1 , function(r) all(r != "")), ]
        gene2GO <- gene2GO[order(gene2GO[,"go_id"]),]
        gene2GO <- gene2GO[gene2GO[,"go_id"] %in% GO2descr[,"go_id"],]
        gs <- sapply(GO2descr[,"go_id"], 
            function(g) gene2GO[gene2GO[,"go_id"] == g, "entrezgene"])
        gs <- sapply(gs, function(s) as.character(sort(s)))
        names(gs) <- paste(GO2descr[,1], gsub(" ", "_", GO2descr[,2]), sep="_")
    }
    return(gs)
}

#
# (2) KEGG
#
get.kegg.genesets <- function(pwys, gmt.file=NULL)
{
    if(class(pwys) == "character")
    {
        if(length(pwys) == 1 && file.exists(pwys)) pwys <- extract.pwys(pwys)
        else return(download.kegg.genesets(pwys, gmt.file=gmt.file))
    }
    return(extract.kegg.genesets(pwys, gmt.file=gmt.file))
}

download.kegg.genesets <- function(pwys, gmt.file=NULL)
{
    # download all gs of organism
    if(length(pwys) == 1 && length(grep("^[a-z]{3}$", pwys)))
    {   
        org <- pwys
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
    }
    # download selected ids
    else
    {
        gs <- sapply(pwys, 
            function(pwy)
            { 
                info <- keggLink(paste0("path:", pwy))
                genes <- grep(paste0("^", 
                    substring(pwy, 1, 3), ":"), info, value=TRUE)
                genes <- sub("^[a-z]{3}:", "", genes)
                genes <- sort(genes)
                names(genes) <- NULL
                return(genes)
            }, simplify=FALSE)
        titles <- sapply(pwys, 
            function(pwy) keggList(paste0("map", sub("^[a-z]{3}", "", pwy))))
        names(titles) <- pwys
        pwys <- titles
    }
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


#
# (3) UTILS
#

# build first gmt column: the ID (format: <pwy.nr>_<pwy.title>)
make.gs.names <- function(ids, titles)
{
    ids <- sub("path:", "", ids)
    titles <- sapply(titles, function(title) unlist(strsplit(title, " - "))[1])
    titles <- stringr::str_trim(titles)
    titles <- gsub(" ", "_", titles)
    ids <- paste(ids, titles, sep="_")
    return(ids)
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

# prepare gene sets as gene set collection from a variety of input formats
prep.gsets <- function(gsets)
{
    if(class(gsets) != "GeneSetCollection")
	{	
        if(class(gsets) != "list")
        {
		    # gmt file or char vector of set names/ids
		    if(class(gsets) != "character")
                stop(paste(
                    "Invalid input: \'gsets\' need to be either a filename",
                    "(GMT file), a character vector (gene set IDs), a named",
                    "list of character vectors, or a GeneSetCollection.", 
                    "See the man page for details."
                ))
            else
            {
			    # gmt file
			    if(file.exists(gsets[1])) 
                    gsets <- parse.genesets.from.GMT(gsets)
			    # char vector
			    else gsets <- auto.create.gsets(gsets)
		    }
        }
        gsets <- gs.list.2.gs.coll(gsets)
	}
    return(gsets)
}

# coerces a list of gene sets (char vecs) into a GeneSetCollection
gs.list.2.gs.coll <- function(gs.list)
{
	gs.type <- auto.detect.gs.type(names(gs.list)[1])
    ctype <- paste0(gs.type, "Collection")
	gs.list.new <- sapply(names(gs.list), 
        function(s){
                spl <- unlist(strsplit(s,"_")) 
                args <- list()
                if(gs.type != "Computed") args <- list(spl[1])
                descr <- ""
                if(length(spl) > 1) 
                    descr <- paste(spl[2:length(spl)],collapse=" ")
                org <- ifelse(gs.type == "KEGG", substring(spl[1], 1, 3), "")
                
                sname <- gsub("[<>:\\?\\|\"\\*\\/]", "", spl[1]) 

                gset <- GeneSet(setName=sname, 
                        geneIds=gs.list[[s]],
                        collectionType=do.call(ctype, args),
                        shortDescription=descr,
                        organism=org)
                return(gset)
        })
    gsc <- GeneSetCollection(gs.list.new) 
    return(gsc)
}

# create gene sets from a list of gene set ids such as
# 'hsa00010' (KEGG ID) or 'GO:0000002' (GO ID)
auto.create.gsets <- function(gs.ids)
{
	gs.type <- auto.detect.gs.type(gs.ids[1])
	if(gs.type == "GO") gs <- get.go.genesets(gs.ids)
	else if(gs.type == "KEGG") gs <- get.kegg.genesets(gs.ids)
	else stop(paste("Automatic gene set recognition 
        is currently only supported for GO and KEGG"))
	return(gs)
}

# currently supported: GO, KEGG, user.def
auto.detect.gs.type <- function(gs.id)
{
	if(substring(gs.id, 1, 3) == "GO:") return("GO")
	else if(grepl("^[a-z]{3}[0-9]{5}", gs.id)) return("KEGG")
	else return("Computed") 
}


