###############################################################################
# 
# Author: ludwig geistlinger
# Date: 16 June 2010
#
# Getting gene sets for enrichment analysis
#
# Update, 05 May 2014: KEGG gene set getter for organism of choice 
# Update, 10 March 2015: including getter for GO genesets
# Update, 24 April 2018: unified getGenesets fassade for gene sets from
#                        different DBs       
###############################################################################

#' @name getGenesets
#'
#' @title Definition of gene sets according to different sources
#' 
#' @description Functionality for retrieving gene sets for an organism under
#' investigation from databases such as GO and KEGG. Parsing and writing a list
#' of gene sets from/to a flat text file in GMT format is also supported.
#'
#' The GMT (Gene Matrix Transposed) file format is a tab delimited file format
#' that describes gene sets.  In the GMT format, each row represents a gene
#' set.  Each gene set is described by a name, a description, and the genes in
#' the gene set. See references.
#'
#' @aliases get.go.genesets get.kegg.genesets parse.genesets.from.GMT
#' @param org An organism in (KEGG) three letter code, e.g. \sQuote{hsa} for
#' \sQuote{Homo sapiens}. Alternatively, this can also be a text file storing 
#' gene sets in GMT format. See details.
#' @param db Database from which gene sets should be retrieved. Currently, 
#' either 'go' (default) or 'kegg'. 
#' @param cache Logical.  Should a locally cached version used if available?
#' Defaults to \code{TRUE}.
#' @param go.onto Character. Specifies one of the three GO ontologies: 'BP'
#' (biological process), 'MF' (molecular function), 'CC' (cellular component).
#' Defaults to 'BP'.
#' @param go.mode Character. Determines in which way the gene sets are retrieved.
#' This can be either 'GO.db' or 'biomart'.  The 'GO.db' mode creates the gene
#' sets based on BioC annotation packages - which is fast, but represents not
#' necessarily the most up-to-date mapping. In addition, this option is only
#' available for the currently supported model organisms in BioC.  The
#' 'biomart' mode downloads the mapping from BioMart - which can be time
#' consuming, but allows to select from a larger range of organisms and
#' contains the latest mappings.  Defaults to 'GO.db'.
#' @param gs A list of gene sets (character vectors of gene IDs).
#' @param gmt.file Gene set file in GMT format. See details.
#' @return For \code{getGenesets}: a list of gene sets (vectors of gene IDs).
#' For \code{writeGMT}: none, writes to file.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{\link{annFUN}} for general GO2gene mapping used in the
#' 'GO.db' mode, and the biomaRt package for general queries to BioMart. 
#' 
#' \code{\link{keggList}} and \code{\link{keggLink}} for accessing the KEGG REST
#' server.
#' @references GO: \url{http://geneontology.org/}
#' 
#' KEGG Organism code \url{http://www.genome.jp/kegg/catalog/org_list.html}
#'
#' GMT file format
#' \url{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}
#' @examples
#' 
#'     # (1) Typical usage for gene set enrichment analysis with GO:
#'     # Biological process terms based on BioC annotation (for human)
#'     go.gs <- getGenesets(org="hsa", db="go")
#'     
#'     # eq.:  
#'     # go.gs <- getGenesets(org="hsa", db="go", go.onto="BP", go.mode="GO.db")
#' 
#'     # Alternatively:
#'     # downloading from BioMart 
#'     # this may take a few minutes ...
#'     \donttest{
#'      go.gs <- getGenesets(org="hsa", db="go", mode="biomart")
#'     }
#'
#'     # (2) Defining gene sets according to KEGG  
#'     kegg.gs <- getGenesets(org="hsa", db="kegg")
#'
#'     # (3) parsing gene sets from GMT
#'     gmt.file <- system.file("extdata/hsa_kegg_gs.gmt", package="EnrichmentBrowser")
#'     gs <- getGenesets(gmt.file)     
#'     
#'     # (4) writing gene sets to file
#'     writeGMT(gs, gmt.file)
#' 
NULL

#' @export
#' @rdname getGenesets
getGenesets <- function(org, 
    db=c("go", "kegg"), 
    cache=TRUE,
    go.onto=c("BP", "MF", "CC"),
    go.mode=c("GO.db", "biomart"))
{
    stopifnot(length(org) == 1 && is.character(org))
    if(file.exists(org)) gs <- .parseGMT(org)
    else
    {
        if(!grepl("^[a-z]{3}$", org))
            stop(paste("\'org\' must be an organism in KEGG three letter code",
                        "or a file in GMT format"))

        db <- match.arg(db)
        if(db == "go") gs <- .getGO(org, go.onto, go.mode, cache)
        else gs <- .getKEGG(org, cache) 
    }
    return(gs)
}

#' @export
#' @keywords internal
get.go.genesets <- function(org, 
    onto=c("BP", "MF", "CC"), mode=c("GO.db","biomart"), cache=TRUE)
{
    .Deprecated("getGenesets")
    .getGO(org, onto, mode, cache)
}

.getGO <- function(org, 
    onto=c("BP", "MF", "CC"), mode=c("GO.db","biomart"), cache=TRUE)
{
    onto <- match.arg(onto)
    gsc.name <- paste(org, "go", tolower(onto), "gs", sep=".") 
    # should a cached version be used?
    if(cache)
    {
        gs <- .getResourceFromCache(gsc.name)
        if(!is.null(gs)) return(gs)
    }

    # from GO.db or from biomart?
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
        useMart <- listDatasets <- useDataset <- getBM <- NULL
        isAvailable("biomaRt", type="software")
        # setting mart
        ensembl <- useMart("ensembl")
        ds <- listDatasets(ensembl)[,"dataset"]
        ds <- grep(paste0("^", org), ds, value=TRUE)
        ensembl <- useDataset(ds, mart=ensembl)

        message("Downloading mapping from BioMart ...")
        message("This may take a few minutes ...")
        
        GO2descr <- getBM(mart=ensembl,
            attributes=c("go_id", "name_1006", "namespace_1003"))
        GO2descr <- GO2descr[GO2descr$go_id != "", ]
        GO2descr <- GO2descr[order(GO2descr[,"go_id"]),]
        ontos <- vapply(GO2descr[,3], 
            function(x)
            {
                spl <- unlist(strsplit(x, "_"))
                onto <- paste(substring(spl,1,1), collapse="")         
                return(toupper(onto))
            }, character(1)) 
        GO2descr <- GO2descr[ontos==onto,1:2]
        
        gene2GO <- getBM(attributes = c("entrezgene", "go_id"), mart=ensembl)
        gene2GO <- gene2GO[apply(gene2GO, 1 , function(r) all(r != "")), ]
        gene2GO <- gene2GO[order(gene2GO[,"go_id"]),]
        gene2GO <- gene2GO[gene2GO[,"go_id"] %in% GO2descr[,"go_id"],]
        gs <- lapply(GO2descr[,"go_id"], 
            function(g) gene2GO[gene2GO[,"go_id"] == g, "entrezgene"])
        gs <- lapply(gs, function(s) as.character(sort(s)))
        names(gs) <- paste(GO2descr[,1], gsub(" ", "_", GO2descr[,2]), sep="_")
    }
    .cacheResource(gs, gsc.name)
    return(gs)
}

#' @export
#' @keywords internal
get.kegg.genesets <- function(pwys, cache=TRUE, gmt.file=NULL)
{    
    .Deprecated("getGenesets")
    .getKEGG(pwys, cache, gmt.file)
}

.getKEGG <- function(pwys, cache=TRUE, gmt.file=NULL)
{
    if(class(pwys) == "character")
    {
        if(length(pwys) == 1 && file.exists(pwys)) pwys <- .extractPwys(pwys)
        else return(.dwnldKeggGS(pwys, cache=cache, gmt.file=gmt.file))
    }
    return(.extractKeggGS(pwys, gmt.file=gmt.file))
}



.dwnldKeggGS <- function(pwys, cache=TRUE, gmt.file=NULL)
{
    # download all gs of organism
    is.org <- length(pwys) == 1 && grepl("^[a-z]{3}$", pwys)
    if(is.org)
    {   
        org <- pwys
        gsc.name <- paste(org, "kegg", "gs", sep=".") 
        # should a cached version be used?
        if(cache)
        {
            gs <- .getResourceFromCache(gsc.name)
            if(!is.null(gs)) return(gs)
        }
        pwys <- KEGGREST::keggList("pathway", org)
        pwy2gene <- KEGGREST::keggLink(org, "pathway")

        gs <- lapply(names(pwys), 
            function(pwy)
            { 
                genes <- pwy2gene[names(pwy2gene) == pwy]
                genes <- sub("^[a-z]{3}:", "", genes)
                genes <- sort(genes)
                names(genes) <- NULL
                return(genes)
            })
    }
    # download selected ids
    else
    {
        gs <- lapply(pwys, 
            function(pwy)
            { 
                info <- KEGGREST::keggLink(paste0("path:", pwy))
                genes <- grep(paste0("^", 
                    substring(pwy, 1, 3), ":"), info, value=TRUE)
                genes <- sub("^[a-z]{3}:", "", genes)
                genes <- sort(genes)
                names(genes) <- NULL
                return(genes)
            })
        .makeTitle <- function(pwy)
        {
            ti <- paste0("map", sub("^[a-z]{3}", "", pwy))
            ti <- KEGGREST::keggList(ti)
            return(ti)
        }
        titles <- vapply(pwys, .makeTitle, character(1))
        names(titles) <- pwys
        pwys <- titles
    }
    names(gs) <- .makeGSNames(names(pwys), pwys)
    if(is.org) .cacheResource(gs, gsc.name)
    if(!is.null(gmt.file)) writeGMT(gs, gmt.file=gmt.file)
    return(gs)
}

# only preferred over '.dwnldKeggGS' 
# when pathway kgmls have already been download
.extractKeggGS <- function(pwys, gmt.file=NULL)
{
    # read in & parse pathways
    if(class(pwys) == "character") pwys <- .extractPwys(pwys)
    
    # get pathway annotations
    nn <- vapply(pwys, getName, character(1))
    tt <- vapply(pwys, getTitle, character(1))
    
    # extract genesets
    gs <- lapply(pwys, 
        function(pwy)
        {
            genes <- .getGenesByPwy(pwy)
            genes <- sub("^[a-z]{3}:", "", genes)
            genes <- sort(genes)
            return(genes)
        })

    names(gs) <- .makeGSNames(nn, tt)
    
    if(!is.null(gmt.file)) writeGMT(gs, gmt.file=gmt.file)
    return(gs)
}

## .extractPwys from zip archive and parse KGML files
.extractPwys <- function(pwy.zip)
{
    pwy.dir <- dirname(pwy.zip)
    unzip(pwy.zip, exdir=pwy.dir, junkpaths=TRUE)
    pwy.files <- list.files(pwy.dir, pattern="*.xml", full.names=TRUE)
    pwys <- sapply(pwy.files, parseKGML)
    ## clean up
    for(f in pwy.files) file.remove(f)
    return(pwys)
}

.getGenesByPwy <- function(pwy)
{
    ts <- vapply(nodes(pwy), getType, character(1))
    genes <- unique(unlist(lapply(nodes(pwy)[ts == "gene"], getName)))
    return(genes)
}

#
# (3) UTILS
#

# build first gmt column: the ID (format: <pwy.nr>_<pwy.title>)
.makeGSNames <- function(ids, titles)
{
    ids <- sub("path:", "", ids)
    titles <- vapply(titles, 
        function(title) unlist(strsplit(title, " - "))[1], 
        character(1))
    titles <- sub("^ +", "", titles)
    titles <- sub(" +$", "", titles)
    titles <- gsub(" ", "_", titles)
    ids <- paste(ids, titles, sep="_")
    return(ids)
}


## parse geneset database
#' @export
#' @keywords internal
parse.genesets.from.GMT <- function(gmt.file)
{
    .Deprecated("getGenesets")
    .parseGMT(gmt.file)
}

.parseGMT <- function(gmt.file)
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

# write genesets to file in GMT format
#' @export
#' @rdname getGenesets
writeGMT <- function(gs, gmt.file)
{
    ## collapse geneset members to one tab separated string 
    gs.strings <- vapply(gs, 
        function(x) paste(x, collapse="\t"),
        character(1))
    
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

# prepare gene sets as gene set collection from a variety of input formats
.prepGS <- function(gsets)
{
    if(!is(gsets, "GeneSetCollection"))
	{	
        if(!is.list(gsets))
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
                    gsets <- getGenesets(gsets)
			    # char vector
			    else gsets <- .createGS(gsets)
		    }
        }
        gsets <- .gsList2Collect(gsets)
	}
    return(gsets)
}

# coerces a list of gene sets (char vecs) into a GeneSetCollection
.gsList2Collect <- function(gs.list)
{
	gs.type <- .detectGSType(names(gs.list)[1])
    ctype <- paste0(gs.type, "Collection")
	gs.list.new <- lapply(names(gs.list), 
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
    names(gs.list.new) <- names(gs.list)
    gsc <- GeneSetCollection(gs.list.new) 
    return(gsc)
}

# create gene sets from a list of gene set ids such as
# 'hsa00010' (KEGG ID) or 'GO:0000002' (GO ID)
.createGS <- function(gs.ids)
{
	gs.type <- .detectGSType(gs.ids[1])
	if(gs.type == "GO") gs <- .getGO(gs.ids)
	else if(gs.type == "KEGG") gs <- .getKEGG(gs.ids)
	else stop(paste("Automatic gene set recognition", 
        "is currently only supported for GO and KEGG"))
	return(gs)
}

# currently supported: GO, KEGG, user.def
.detectGSType <- function(gs.id)
{
	if(substring(gs.id, 1, 3) == "GO:") return("GO")
	else if(grepl("^[a-z]{3}[0-9]{5}", gs.id)) return("KEGG")
	else return("Computed") 
}

