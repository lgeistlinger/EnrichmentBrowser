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
# Update, Dec 2019: including getter for MSigDB and Enrichr       
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
#' either 'go' (default), 'kegg', 'msigdb', or 'enrichr'. 
#' @param cache Logical.  Should a locally cached version used if available?
#' Defaults to \code{TRUE}.
#' @param return.type Character. Determines whether gene sets are returned
#' as a simple list of gene sets (each being a character vector of gene IDs), or
#' as an object of class \code{\linkS4class{GeneSetCollection}}.
#' @param ... Additional arguments for individual gene set databases. 
#' For \code{db = "GO"}: \itemize{ \item onto: Character. Specifies one of the
#' three GO ontologies: 'BP'
#' (biological process), 'MF' (molecular function), 'CC' (cellular component).
#' Defaults to 'BP'. \item mode: Character. Determines in which way the gene 
#' sets are retrieved. This can be either 'GO.db' or 'biomart'. 
#' The 'GO.db' mode creates the gene sets based on BioC annotation packages - 
#' which is fast, but represents not
#' necessarily the most up-to-date mapping. In addition, this option is only
#' available for the currently supported model organisms in BioC.  The
#' 'biomart' mode downloads the mapping from BioMart - which can be time
#' consuming, but allows to select from a larger range of organisms and
#' contains the latest mappings.  Defaults to 'GO.db'.}
#' For \code{db = "msigdb":} \itemize{ \item cat: Character. 
#' MSigDB collection category: 'H' (hallmark), 
#' 'C1' (genomic position), 'C2' (curated databases), 'C3' (binding site motifs),
#' 'C4' (computational cancer), 'C5' (Gene Ontology), 'C6' (oncogenic), 
#' 'C7' (immunologic). See references. 
#' \item subcat: Character. MSigDB collection subcategory. Depends on the
#' chosen MSigDB collection category. For example, 'MIR' to obtain microRNA targets
#' from the 'C3' collection. See references.}
#' For \code{db = "enrichr"}: \itemize{ \item lib: Character. Enrichr gene set 
#' library. For example, 'Genes_Associated_with_NIH_Grants' to obtain gene sets 
#' based on associations with NIH grants. See references.
#' \item show.libs: Logical. Show available gene set libraries? Defaults to 
#' \code{FALSE}.}
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
#'
#' \code{msigdbr::msigdbr} for obtaining gene sets from the MSigDB.
#' @references GO: \url{http://geneontology.org/}
#' 
#' KEGG Organism code: \url{http://www.genome.jp/kegg/catalog/org_list.html}
#'
#' MSigDB: \url{http://software.broadinstitute.org/gsea/msigdb/collections.jsp}
#'
#' Enrichr: \url{https://amp.pharm.mssm.edu/Enrichr/#stats}
#'
#' GMT file format:
#' \url{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}
#' @examples
#' 
#'     # (1) Typical usage for gene set enrichment analysis with GO:
#'     # Biological process terms based on BioC annotation (for human)
#'     go.gs <- getGenesets(org="hsa", db="go")
#'     
#'     # eq.:  
#'     # go.gs <- getGenesets(org="hsa", db="go", onto="BP", mode="GO.db")
#'     \donttest{
#'     # Alternatively:
#'     # downloading from BioMart 
#'     # this may take a few minutes ...
#'     go.gs <- getGenesets(org="hsa", db="go", mode="biomart")
#'     }
#'     # (2) Defining gene sets according to KEGG  
#'     kegg.gs <- getGenesets(org="hsa", db="kegg")
#'     \donttest{
#'     # (3) Obtaining *H*allmark gene sets from MSigDB
#'     hall.gs <- getGenesets(org="hsa", db="msigdb", cat="H")
#'     }
#'     \donttest{
#'     # (4) Obtaining gene sets from Enrichr
#'     tfppi.gs <- getGenesets(org="hsa", db="enrichr", lib="Transcription_Factor_PPIs")
#'     
#'     # displaying available Enrichr gene set libraries
#'     getGenesets(org="hsa", db="enrichr", show.libs=TRUE)        
#'     }
#'     # (6) parsing gene sets from GMT
#'     gmt.file <- system.file("extdata/hsa_kegg_gs.gmt", package="EnrichmentBrowser")
#'     gs <- getGenesets(gmt.file)     
#'     
#'     # (7) writing gene sets to file
#'     writeGMT(gs, gmt.file)
#' 
NULL

#' @export
#' @rdname getGenesets
getGenesets <- function(org, 
    db = c("go", "kegg", "msigdb", "enrichr"), 
    cache = TRUE,
    return.type = c("list", "GeneSetCollection"),
    ...)
{
    stopifnot(length(org) == 1 && is.character(org))
    return.type <- match.arg(return.type)
    if(file.exists(org)) gs <- .parseGMT(org, return.type)
    else
    {
        if(!grepl("^[a-z]{3}$", org))
            stop(paste("\'org\' must be an organism in KEGG three letter code",
                        "or a file in GMT format"))

        db <- match.arg(db)
        if(db == "go") gs <- .getGO(org, cache, return.type, ...)
        else if(db == "kegg") gs <- .getKEGG(org, cache, return.type = return.type)
        else if(db == "msigdb") gs <- .getMSigDb(org, cache, return.type, ...)
        else gs <- .getEnrichr(org, cache, return.type, ...)
    }
    return(gs)
}

.getEnrichr <- function(org, cache, return.type, lib, show.libs = FALSE)
{
    # convert organism
    eorg <- .org2enrichr(org)
    
    # only show available gene set libraries?
    libs <- .enrichrLibs(eorg)
    if(show.libs) return(libs)

    gsc.name <- paste(org, "enrichr", lib, "gs", sep=".") 
    # should a cached version be used?
    if(cache)
    {
        gs <- .getResourceFromCache(gsc.name)
        if(!is.null(gs)) return(gs)
    }
    
    # obtain the gene set library
    if(!(lib %in% libs[,"libraryName"])) stop("Not a valid gene set library") 
    gs <- .getEnrichrLib(eorg, lib)    

    # ID mapping
    from <- ifelse(org == "sce", "GENENAME", "SYMBOL")
    suppressMessages(
        gs <- .idMapGS(gs, org, from = from, to = "ENTREZID")
    )
    if(return.type == "GeneSetCollection")
    {
        ct <- ComputedCollection()
        titles <- vapply(names(gs), .extractTitle, character(1))
        gs <- .makeGSC(gs, titles, org, ct)        
    }
    
    .cacheResource(gs, gsc.name)
    return(gs)
}

.getEnrichrLib <- function(eorg, lib)
{
    getURL <- NULL
    isAvailable("RCurl")               
    
    eurl <- paste0("https://maayanlab.cloud/", eorg, 
                    "/geneSetLibrary?mode=text&libraryName=", lib)
    gs <- getURL(eurl)
    gs <- unlist(strsplit(gs, "\n"))
    gs <- lapply(gs, function(x) unlist(strsplit(x, "\t")))
    names(gs) <- vapply(gs, function(x) x[1], character(1))
    names(gs) <- gsub(" ", "_", names(gs))
    ids <- paste0("ER", seq_along(gs))
    names(gs) <- paste(ids, names(gs), sep = "_")
    gs <- lapply(gs, function(x) x[-c(1,2)])
    .trim <- function(x) vapply(x, 
                                function(y) unlist(strsplit(y, ","))[1], 
                                character(1), USE.NAMES=FALSE)
    lapply(gs, .trim)
}

.org2enrichr <- function(org)
{
    sorgs <- c("cel", "dme", "dre", "hsa", "sce")
    if(!(org %in% sorgs)) stop("Organism not supported") 
    eorgs <- c("Worm", "Fly", "Fish", "", "Yeast")
    names(eorgs) <- sorgs
    eorg <- eorgs[org]
    paste0(eorg, "Enrichr")
}

.enrichrLibs <- function(eorg)
{
    eurl <- paste0("https://maayanlab.cloud/", eorg, "/datasetStatistics")
    GET <- fromJSON <- NULL
    isAvailable("httr", type = "software")
    dbs <- GET(eurl)
    dbs <- dbs$content
    dbs <- intToUtf8(dbs)
    isAvailable("rjson", type = "software")
    dbs <- fromJSON(dbs)
    dfSAF <- getOption("stringsAsFactors")
    options(stringsAsFactors = FALSE)
    dbs <- lapply(dbs$statistics, function(x) do.call(cbind.data.frame, x))
    dbs <- do.call(rbind.data.frame, dbs)
    options(stringsAsFactors = dfSAF)
    dbs
}

.getMSigDb <- function(org, cache, return.type,
                        cat = c("H", paste0("C", 1:7)), 
                        subcat = NA)
{
    cat <- match.arg(cat)
    
    gsc.name <- paste(org, "msigdb", cat, subcat, "gs", sep=".") 
    # should a cached version be used?
    if(cache)
    {
        gs <- .getResourceFromCache(gsc.name)
        if(!is.null(gs)) return(gs)
    }

    msigdbr_show_species <- msigdbr <- NULL
    isAvailable("msigdbr", type = "software")

    ind <- match(org, SPECIES[,"kegg"])
    org <- SPECIES[ind, "tax"]

    if(!(org %in% msigdbr_show_species())) stop("Organism not supported")

    df <- msigdbr(org, cat, subcat)
    gs <- split(as.character(df$entrez_gene), df$gs_id)
    gs.names <- unique(df$gs_name)
    gs.ids <- unique(df$gs_id)
    gs.names <- paste(gs.ids, gs.names, sep = "_")
    names(gs.names) <- gs.ids
    names(gs) <- unname(gs.names[names(gs)])
    
    if(return.type == "GeneSetCollection")
    {
        ct <- BroadCollection(category = tolower(cat), 
                                subCategory = tolower(subcat))
        titles <- vapply(names(gs), .extractTitle, character(1))
        gs <- .makeGSC(gs, titles, org, ct)        
    }
    
    .cacheResource(gs, gsc.name)
    return(gs)
}

.getGO <- function(org, cache, return.type, 
                    onto=c("BP", "MF", "CC"), 
                    mode=c("GO.db","biomart"))
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
    if(mode=="GO.db") gs <- .getGOFromGODB(org, onto, return.type)
    else gs <- .getGOFromBiomart(org, onto, return.type)
    .cacheResource(gs, gsc.name)
    return(gs)
}

.getGOFromGODB <- function(org, onto, return.type)
{
    gs <- topGO::annFUN.org(whichOnto=onto, mapping=.org2pkg(org))
    GO2descr <- AnnotationDbi::as.list(GO.db::GOTERM)
    gs <- gs[intersect(names(gs), names(GO2descr))]
    GO2descr <- GO2descr[names(gs)]
    GO2descr <- vapply(GO2descr, AnnotationDbi::Term, character(1))

    if(return.type == "GeneSetCollection")
    {
        ct <- GOCollection(evidenceCode=NA_character_, ontology=onto)
        gs <- .makeGSC(gs, GO2descr, org, ct)        
    }
    else names(gs) <- paste(names(gs), gsub(" ", "_", GO2descr), sep="_")
    return(gs)
}

.getGOFromBiomart <- function(org, onto, return.type)
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
    .getOnto <- function(x)
    {
        spl <- unlist(strsplit(x, "_"))
        onto <- paste(substring(spl,1,1), collapse="")         
        return(toupper(onto))
    }
    ontos <- vapply(GO2descr[,3], .getOnto, character(1)) 
    GO2descr <- GO2descr[ontos==onto,1:2]
    
    gene2GO <- getBM(attributes = c("entrezgene_id", "go_id"), mart=ensembl)
    gene2GO <- gene2GO[apply(gene2GO, 1 , function(r) all(r != "")), ]
    gene2GO <- gene2GO[order(gene2GO[,"go_id"]),]
    gene2GO <- gene2GO[gene2GO[,"go_id"] %in% GO2descr[,"go_id"],]
    gs <- lapply(GO2descr[,"go_id"], 
        function(g) gene2GO[gene2GO[,"go_id"] == g, "entrezgene_id"])
    gs <- lapply(gs, function(s) sort(as.character(s)))

    ids <- GO2descr[,1]
    GO2descr <- GO2descr[,2] 
    names(gs) <- names(GO2descr) <- ids
    if(return.type == "GeneSetCollection")
    {
        ct <- GOCollection(evidenceCode=NA_character_, ontology=onto)
        gs <- .makeGSC(gs, GO2descr, org, ct)        
    }
    else names(gs) <- paste(ids, gsub(" ", "_", unname(GO2descr)), sep="_")
    return(gs)
}

.makeGSC <- function(gs, titles, org, ctype)
{
    it <- EntrezIdentifier()
    .makeGeneSet <- function(s)
    {
        sname <- gsub("[<>:\\?\\|\"\\*\\/]", "", s)
        is.idcoll <- !is(ctype, "ComputedCollection") && !is(ctype, "BroadCollection")
        if(is.idcoll) ctype@ids <- s 
        gset <- GeneSet(setName=sname,
                    geneIds=gs[[s]],
                    type=it,
                    collectionType=ctype,
                    shortDescription=unname(titles[s]),
                    organism=org)
        return(gset)
    }
    gs <- lapply(names(gs), .makeGeneSet)
    gs <- GeneSetCollection(gs)
}

.getKEGG <- function(pwys, cache, gmt.file=NULL, return.type)
{
    is.org <- length(pwys) == 1 && grepl("^[a-z]{3}$", pwys)
    # download all gs of organism
    if(is.org) gs <- .dwnldAllKeggGS(pwys, cache, return.type)
    # download selected ids
    else gs <- .dwnldSelectedKeggGS(pwys)
    if(!is.null(gmt.file)) writeGMT(gs, gmt.file=gmt.file)
    return(gs)
}


.dwnldAllKeggGS <- function(org, cache, return.type)
{
    gsc.name <- paste(org, "kegg", "gs", sep=".") 
    # should a cached version be used?
    if(cache)
    {
        gs <- .getResourceFromCache(gsc.name)
        if(!is.null(gs)) return(gs)
    }
    pwys <- KEGGREST::keggList("pathway", org)
    pwy2gene <- KEGGREST::keggLink(org, "pathway")

    org.start <- paste0("^", org, ":")
    .getGenes <- function(pwy)
    { 
        genes <- pwy2gene[names(pwy2gene) == pwy]
        genes <- sub(org.start, "", genes)
        genes <- unname(sort(genes))
        return(genes)
    }
    gs <- lapply(names(pwys), .getGenes)
        
    # entrez id?
    id.type <- .detectGeneIdType(gs[[1]][1])
    if(is.na(id.type) || id.type != "entrez") 
    {
        map.k2e <- try(KEGGREST::keggConv("ncbi-geneid", org), silent = TRUE)
        if(is(map.k2e, "try-error"))
        {
            message("NCBI Gene ID mapping is not available for: ", org)
            message("Returning KEGG native gene IDs ...")
        }
        else
        {
            names(map.k2e) <- sub(org.start, "", names(map.k2e))
            map.k2e <- sub("^ncbi-geneid:", "", map.k2e)
            gs <- lapply(gs, function(s) unname(map.k2e[s]))
        }
    }
    if(return.type == "GeneSetCollection") gs <- .makeKeggGSC(gs, pwys, org)
    else names(gs) <- .makeGSNames(names(pwys), pwys)
    .cacheResource(gs, gsc.name)
    return(gs)
}

.makeKeggGSC <- function(gs, pwys, org)
{
    ids <- sub("^path:", "", names(pwys))
    titles <- vapply(unname(pwys), 
        function(title) unlist(strsplit(title, " - "))[1], 
        character(1), USE.NAMES=FALSE)
    names(gs) <- names(titles) <- ids 
    ct <- KEGGCollection()
    gs <- .makeGSC(gs, titles, org, ct)
    return(gs)
}

.dwnldSelectedKeggGS <- function(pwys)
{
    .dwnld <- function(pwy)
    { 
        info <- KEGGREST::keggLink(paste0("path:", pwy))
        genes <- grep(paste0("^", 
            substring(pwy, 1, 3), ":"), info, value=TRUE)
        genes <- sub("^[a-z]{3}:", "", genes)
        genes <- unname(sort(genes))
        return(genes)
    }
    gs <- lapply(pwys, .dwnld) 
    .makeTitle <- function(pwy)
    {
        ti <- paste0("map", sub("^[a-z]{3}", "", pwy))
        ti <- KEGGREST::keggList(ti)
        return(ti)
    }
    titles <- vapply(pwys, .makeTitle, character(1))
    names(titles) <- pwys
    pwys <- titles
    names(gs) <- .makeGSNames(names(pwys), pwys)
    return(gs)
}

# only preferred over '.dwnldKeggGS' 
# when pathway kgmls have already been download
.extractKeggGS <- function(pwys, gmt.file=NULL)
{
    # read in & parse pathways
    if(is.character(pwys)) pwys <- .extractPwys(pwys)
    
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
    ids <- sub("^path:", "", ids)
    titles <- vapply(titles, 
        function(title) unlist(strsplit(title, " - "))[1], 
        character(1))
    titles <- sub("^ +", "", titles)
    titles <- sub(" +$", "", titles)
    titles <- gsub(" ", "_", titles)
    ids <- paste(ids, titles, sep="_")
    return(ids)
}


.extractTitle <- function(n)
{
    spl <- unlist(strsplit(n, "_"))
    paste(spl[-1], collapse = " ")
}

.parseGMT <- function(gmt.file, return.type)
{
    content <- readLines(gmt.file, warn=FALSE)
    le <- length(content)
    gs <- vector("list", length=le)
    gs.names <- vector("character", length=le)
    for(i in seq_len(le))
    {
        line <- content[i]
        spl <- unlist(strsplit(line, "\t"))
        gs.names[i] <- spl[1]
        gs[[i]] <- spl[-c(1,2)]
    }
    if(return.type == "GeneSetCollection")
    {
        spl <- lapply(gs.names, function(s) unlist(strsplit(s, "_"))) 
        ids <- vapply(spl, function(s) s[1], character(1))
        if(length(unique(ids)) < length(ids)) ids <- paste0(ids, seq_along(ids))
        .getTitle <- function(s)
        {
            descr <- "" 
            if(length(s) > 1) descr <- paste(s[2:length(s)], collapse=" ")
            return(descr)
        }
        titles <- vapply(spl, .getTitle, character(1))
        names(gs) <- names(titles) <- ids

        ct <- ComputedCollection()
        gs <- .makeGSC(gs, titles, NA_character_, ct)        
    }
    names(gs) <- gs.names
    return(gs)
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
    # gs type
    first <- names(gs.list)[1]
	gs.type <- .detectGSType(first)
    ctype <- paste0(gs.type, "Collection")
    args <- list()
    if(gs.type == "GO") args <- list(evidenceCode=NA_character_)
    ctype <- do.call(ctype, args)
    
    # org
    spl <- unlist(strsplit(first, "_")) 
    org <- ifelse(gs.type == "KEGG", substring(spl[1], 1, 3), NA_character_)
       
    # ids & titles
    spl <- lapply(names(gs.list), function(s) unlist(strsplit(s, "_"))) 
    ids <- vapply(spl, function(s) s[1], character(1))
    if(length(unique(ids)) < length(ids)) ids <- paste0(ids, seq_along(ids))
    .getTitle <- function(s)
    {
        descr <- "" 
        if(length(s) > 1) descr <- paste(s[2:length(s)], collapse=" ")
        return(descr)
    }
    titles <- vapply(spl, .getTitle, character(1))
    names(gs.list) <- names(titles) <- ids

    gsc <- .makeGSC(gs.list, titles, org, ctype)        
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

