#' @name getMicrobeSets
#'
#' @title Definition of microbe sets according to different sources
#' 
#' @description Functionality for retrieving microbe sets from databases such as
#' BugSigDB and the MicrobeDirectory. 
#'
#' @param db Database from which microbe sets should be retrieved. Currently, 
#' either 'bugsigdb' (default), 'manalyst', 'mpattern', or 'mdirectory'. 
#' @param tax.id.type Character. Taxonomic ID type of the returned microbe sets.
#' Currently either 'ncbi' (default) or 'metaphlan'.
#' @param cache Logical.  Should a locally cached version used if available?
#' Defaults to \code{TRUE}.
#' @param ... Additional arguments for individual microbe set databases.
#' For \code{db = "manalyst"}: \itemize{ \item lib: Character. MicrobiomeAnalyst
#' taxon set library. Options include taxon sets associated with human genetic
#' variations ('gene'), host-intrinsic ('host_int'), host-extrinsic ('host_ext'),
#' environmental ('env'), and microbiome-intrinsic ('mic_int') factors.
#' See references.}
#' @references BugSigDB: \url{https://bugsigdb.org}
#' 
#' MicrobiomeAnalyst: \url{https://www.microbiomeanalyst.ca}
#'
#' @export
getMicrobeSets <- function(db = c("bugsigdb",
                                  "manalyst",
                                  "mpattern",
                                  "mdirectory"),
                           tax.id.type = c("ncbi", "metaphlan", "taxname"),
                           tax.level = "mixed", 
                           cache = TRUE, 
                           ...)
{
    db <- match.arg(db)
    tax.id.type <- match.arg(tax.id.type)
    
    stopifnot(is.character(tax.level))
    if("mixed" %in% tax.level) tax.level <- "mixed" 
    else if(!all(tax.level %in% TAX.LEVELS))
            stop("tax.level must be a subset of { ",
                 paste(TAX.LEVELS, collapse = ", "), " }")

    if(db == "bugsigdb") .getBugSigDB(tax.id.type, tax.level, cache)
    else if(db == "manalyst") .getMAnalyst(tax.id.type, tax.level, cache, ...)
    #else if(db == "mpattern") .getMPattern(tax.id.type)
    #else if(db == "mdirectory") .getMDir(tax.id.type)
}

TAX.LEVELS <- c("kingdom", "phylum", "class", "order",
                "family", "genus", "species", "strain")
MPA.TAX.LEVELS <- c(substring(TAX.LEVELS[1:7], 1, 1), "t")
names(MPA.TAX.LEVELS) <- TAX.LEVELS

.getBugSigDB <- function(id.type, tax.level, cache)
{
    # should a cached version be used?
    msc.name <- paste("bugsigdb", id.type, sep = ".")
    if(cache)
    {
        sigs <- .getResourceFromCache(msc.name)
        if(!is.null(sigs)) return(sigs)
    }
    
    # obtain study, experiment, and signature tables 
    einfo <- .getExpInfo()   
    sigs <- read.csv("https://tinyurl.com/yakgsowm")

    # extract signatures
    snames <- .makeSigNames(sigs, einfo) 
    sigs <- .extractSigs(sigs, id.type, tax.level)
    names(sigs) <- paste(snames$id, snames$titles, sep = "_")
    .cacheResource(sigs, msc.name)
    sigs
}

.extractSigs <- function(sigdf, id.type, tax.level)
{
    id.col <- ifelse(id.type == "ncbi",
                     "NCBI.Taxonomy.IDs",
                     "MetaPhlAn.taxon.names")
    sigs <- sigdf[[id.col]]
    sigs <- strsplit(sigs, ",")
    
    if(tax.level[1] != "mixed")
    {
        if(id.type == "ncbi")
        {
            msigs <- sigdf[["MetaPhlAn.taxon.names"]]
            msigs <- strsplit(msigs, ",")
        }
        else msigs <- sigs
        bugs <- unique(unlist(msigs))

        ind <- lapply(msigs, function(s) match(s, bugs))
        istl <- .isTaxLevel(bugs, tax.level, "metaphlan")
        
        subind <- istl[unlist(ind)]
        subind <- relist(subind, ind)
        for(i in seq_along(sigs)) sigs[[i]] <- sigs[[i]][subind[[i]]]
    }

    if(id.type != "metaphlan")
    {
        sigs <- lapply(sigs, .getTip)
        if(id.type == "taxname")
            sigs <- lapply(sigs, function(s) sub("^[kpcofgst]__", "", s))
    }
    sigs
}

.getTip <- function(n)
{ 
    spl <- strsplit(n, "\\|")
    vapply(spl, function(s) s[length(s)], character(1)) 
}

.isTaxLevel <- function(s, tax.level, id.type)
{
    if(tax.level[1] == "mixed") return(rep(TRUE, length(s)))
    if(id.type == "metaphlan")
    {
        tip <- .getTip(s)
        tip <- substring(tip, 1, 1)
        mtl <- MPA.TAX.LEVELS[tax.level]
        istl <- tip %in% mtl        
    }
    else
    {
        isAvailable("taxize")
        sink(tempfile())
        suppressMessages(
            ranks <- taxize::tax_rank(s, db = "ncbi")
        )
        sink()        
        ranks <- vapply(ranks, function(x) x[1], character(1))
        istl <- ranks %in% tax.level
    }
    return(istl)
}


.getTaxLevel <- function(s, tax.level, id.type)
{
    if(tax.level[1] == "mixed") return(s)
    if(id.type == "metaphlan")
    {
        tip <- .getTip(s)
        tip <- substring(tip, 1, 1)
        mtl <- MPA.TAX.LEVELS[tax.level]
        s <- s[tip %in% mtl]        
    }
    else
    {
        isAvailable("taxize")
        sink(tempfile())
        suppressMessages(
            ranks <- taxize::tax_rank(s, db = "ncbi")
        )
        sink()        
        ranks <- vapply(ranks, function(x) x[1], character(1))
        s <- s[ranks %in% tax.level]
    }
    return(s)
}

.makeSigNames <- function(sigs, einfo)
{
    eid <- sub("^Experiment ", "", sigs[["Experiment"]])
    sid <- sub("^Study ", "", sigs[["Study"]])
    id <- paste(sid, eid, sep = "/")
    sgid <- sub("^Signature ", "", sigs[["Signature.page.name"]])
    up.down <- ifelse(sigs[["Abundance.in.Group.1"]] == "increased", "UP", "DOWN")
    titles <- paste(unname(einfo[id]), up.down, sep = "_")
    id <- paste(id, sgid, sep = "/")
    id <- paste("bsdb", id, sep = ":")
    list(id = id, titles = titles)
}

.getExpInfo <- function()
{
    exps <- read.csv("https://tinyurl.com/yb2fmpa3")
    
    is.study <- grepl("^Study [0-9]+$", exps[["Study"]])
    is.exp <- grepl("^Experiment [0-9]+$", exps[["Experiment.page.name"]])
    exps <- exps[is.study & is.exp,]    

    eid <- sub("^Experiment ", "", exps[["Experiment.page.name"]])
    sid <- sub("^Study ", "", exps[["Study"]])
    id <- paste(sid, eid, sep = "/")
    
    rel.cols <- c("Condition",
                  "Group.1.name",
                  "Group.0.name")        
    exps <- exps[,rel.cols]
    .conc <- function(x) paste(x[1], paste(x[2:3], collapse = "_vs_"), sep = ":")
    exp.str <- apply(exps, 1, .conc)
    exp.str <- gsub(" ", "-", exp.str)
    names(exp.str) <- id
    exp.str
}

.study2pmid <- function()
{
    studs <- read.csv("https://tinyurl.com/ycg8fs9x")

    s2pmid <- studs[["PMID"]]
    names(s2pmid) <- studs[["Study.page.name"]]
    s2pmid
}
 
# TODO: ID mapping
.getMAnalyst <- function(id.type,
                         tax.level,
                         cache, 
                         lib = c("host_int", "host_ext", "env", "mic_int", "gene"))
{
    lib <- match.arg(lib)    

    # cache ?
    msc.name <- paste("mana", lib, id.type, sep = ".")
 
    # should a cached version be used?
    if(cache)
    {
        sigs <- .getResourceFromCache(msc.name)
        if(!is.null(sigs)) return(sigs)
    }
    
    ma.url <- paste0("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/",
                           "resources/lib/tsea/tsea_")

    if(!(lib %in% c("host_int", "host_ext")))
        lib <- switch(lib, 
                        gene = "host_snps_new",
                        mic_int = "microbiome_int",
                        env = "environment")
    ma.url <- paste0(ma.url, lib, ".csv")
    cont <- read.csv(ma.url)
    rel.cols <- c("name", "member", "abund_change")
    cont <- cont[,rel.cols]
    
    sigs <- strsplit(cont[["member"]], "; +")
    up.down <- ifelse(cont[["abund_change"]] == "Increase", "UP", "DOWN")
    titles <- sub(" \\(.+\\)$", "", cont[["name"]]) 
    titles <- gsub(" ", "_", titles)   

    id <- seq_along(sigs)
    id <- paste0("MA", id)

    names(sigs) <- paste(id, titles, up.down, sep = "_")
    .cacheResource(sigs, msc.name)
    sigs
}


