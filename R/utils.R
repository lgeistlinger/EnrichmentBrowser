############################################################
# 
# author: Ludwig Geistlinger
# date: 2015-03-10 13:32:37
# 
# descr: misc utils
# 
############################################################

#' Is a required package available?
#' 
#' Convenience function for checking and installing required packages.
#' 
#' Checks whether a required package is available in the library.  If yes, the
#' package is loaded via \code{\link{require}}.  If not, the package is
#' optionally installed via \code{\link{install}} and,
#' subsequently, loaded via \code{\link{require}}.
#' 
#' @param pkg Character vector of length 1.  A valid name of an existing R
#' package.
#' @param type Character vector of length 1.  What type of package is this?
#' Choose one out of 'annotation', 'software', or 'data' package.
#' @return None. See details.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@sph.cuny.edu>
#' @seealso \code{require}, \code{install}
#' @examples
#' 
#'     \donttest{
#'         isAvailable("EnrichmentBrowser", type="software")
#'     }     
#' 
#' @export isAvailable
isAvailable <- function(pkg, type=c("annotation", "software", "data"))
{
    type <- match.arg(type)
    if(!(pkg %in% .packages(all.available=TRUE)))
    {   
        message(paste0("Corresponding ", type,  " package not found: ", 
            pkg, "\nMake sure that you have it installed."))

        if(interactive()) choice <- readline("Install it now? (y/n): ")
        else choice <- "y"
        
        if(choice == "y") BiocManager::install(pkg, update=FALSE)
        else stop(paste("Package", pkg, "is not available"))
    }   
    require(pkg, character.only=TRUE, quietly=TRUE)
}

SPECIES <- rbind(
        c("anopheles", "Anopheles gambiae", "Ag", "aga", "anoGam", "7165"),
        c("arabidopsis", "Arabidopsis thaliana", "At", "ath", NA, "3702"),
        c("bovine", "Bos taurus", "Bt", "bta", "bosTau", "9913"),
        c("canine", "Canis familiaris", "Cf", "cfa", "canFam", "9615"),
        c("chicken", "Gallus gallus", "Gg", "gga", "galGal", "9031"), 
        c("chimp", "Pan troglodytes", "Pt", "ptr", "PanTro", "9598"),
        c("ecoliK12", "Escherichia coli K12", "EcK12", "eco", NA, "562,83333,511145"), 
        c("ecoliSakai", "Escherichia coli Sakai", "EcSakai", "ecs", NA, "83334"),
        c("fly", "Drosophila melanogaster", "Dm", "dme", "dm", "7227"),
        c("human", "Homo sapiens", "Hs", "hsa", "hg", "9606"),
        c("malaria", "Plasmodium falciparum", "Pf", "pfa", NA, "5833"),
        c("mouse", "Mus musculus", "Mm", "mmu", "mm", "10090"),
        c("pig", "Sus scrofa", "Ss", "ssc", "susScr", "9823"),
        c("rat", "Rattus norvegicus", "Rn", "rno", "rn", "10116"), 
        c("rhesus", "Macaca mulatta", "Mmu", "mcc", "rheMac", "9544"),  
        c("worm", "Caenorhabditis elegans", "Ce", "cel", "ce", "6239"),
        c("xenopus", "Xenopus laevis", "Xl", "xla", "NA", "8355"),
        c("yeast", "Saccharomyces cerevisiae", "Sc", "sce", "sacCer", "4932,559292"),
        c("zebrafish", "Danio rerio", "Dr", "dre", "danRer", "7955")
    )
colnames(SPECIES) <- c("common", "tax", "bioc", "kegg", "ucsc", "ncbi")

.cacheResource <- function(res, rname, ucdir="EnrichmentBrowser")
{
    # are we running interactively?
    cache.dir <- ifelse(interactive(), 
                    tools::R_user_dir(ucdir, which = "cache"),
                    tempdir())

    bfc <- BiocFileCache::BiocFileCache(cache.dir)
    
    # replace existing version if necessary 
    qgsc <-  BiocFileCache::bfcquery(bfc, rname)
    if(BiocFileCache::bfccount(qgsc)) BiocFileCache::bfcremove(bfc, qgsc$rid) 
    
    cache.file <- BiocFileCache::bfcnew(bfc, rname)
    saveRDS(res, file=cache.file)
}

.getResourceFromCache <- function(rname, 
    update.value=7, update.unit="days", ucdir="EnrichmentBrowser")
{
    # are we running interactively?
    cache.dir <- ifelse(interactive(), 
                    tools::R_user_dir(ucdir, which = "cache"),
                    tempdir())

    bfc <- BiocFileCache::BiocFileCache(cache.dir)
    qgsc <-  BiocFileCache::bfcquery(bfc, rname)

    # is there already a cached version?
    res <- NULL
    if(BiocFileCache::bfccount(qgsc))
    {
        # is the cached version outdated?
        dt <- difftime(Sys.time(), qgsc$create_time, units=update.unit)   
        if(is.na(update.value) || dt < update.value)
        {
            if(interactive())
                message(paste("Using cached version from", qgsc$create_time))
            res <- readRDS(BiocFileCache::bfcrpath(bfc, rname))
        }
    }
    return(res)   
}

.detectGeneIdType <- function(id)
{
    type <- NA
    if(grepl("^[Ee][Nn][Ss][A-Za-z]{0,3}[Gg][0-9]+", id)) type <- "ensembl"
    else if(grepl("^[0-9]+$", id)) type <- "entrez"
    else if(grepl("^[Yy][A-Za-z]{2}[0-9]{3}[A-Za-z]", id)) type <- "sgd"
    else if(grepl("^[Aa][Tt][0-9][A-Za-z][0-9]{5}", id)) type <- "tair"
    return(type)
}

.getOrgIdType <- function(org)
{
    it <- "eg"
    if(org == "At") it <- "tair"
    else if(org == "Pf") it <- "plasmo"
    else if(org == "Sc") it <- "sgd"
    return(it)
}   

#.supportedOrganisms <- function() sub(".db0$", "", available.db0pkgs())

.availableOrgPkgs <- function(type=c("OrgDb", "TxDb", "BSgenome"), local=TRUE)
{
    if(local) pkgs <- .packages(all.available=TRUE)
    else
    {
        pver <- substring(packageVersion("EnrichmentBrowser"), 3, 3)
        is.even <- as.integer(pver) %% 2 == 0
        repo.type <- ifelse(is.even, "release", "devel")
        pkgs <- available.packages(paste0("http://bioconductor.org/",
            "packages/", repo.type, "/data/annotation/src/contrib"))[, "Package"]
    }
    type <- match.arg(type)
    org.string <- "^org.[A-Za-z0-9]+.[a-z]+.db$"
    if(type == "TxDb") 
        org.string <- "^TxDb.[A-Z][a-z]+.UCSC.[a-z]{2}[A-Za-z]*[0-9]{1,3}.[a-z]{3,5}Gene$"
    else if(type == "BSgenome") 
        org.string <- "^BSgenome.[A-Z][a-z]+.UCSC.[a-z]{2}[A-Za-z]*[0-9]{1,3}$"
    org.pkgs <- grep(org.string, pkgs, value=TRUE)
    names(org.pkgs) <- NULL 
    return(org.pkgs)
}

.org2pkg <- function(org, type=c("OrgDb", "TxDb", "BSgenome"))
{
    type <- match.arg(type)

        

    # org specification via 
    # (a) 3-letter code, e.g. 'hsa' 
    # (b) genome assembly, e.g. 'hg38'
    is.genome <- sub("[0-9]+$", "", org) %in% SPECIES[,"ucsc"]
    if(is.genome)
    {
        ucsc.id <- org
        i <- grep(sub("[0-9]+$", "", org), SPECIES[,"ucsc"]) 
        bioc.id <- SPECIES[i, "bioc"]
    }
    else
    {
        ind <- apply(SPECIES, 1, function(r) org %in% r)
        if(any(ind)) i <- which(ind)[1]
        else stop(paste0("unrecognized organism ID \'", org, "\'"))
        bioc.id <- SPECIES[i, "bioc"]
        ucsc.id <- SPECIES[i, "ucsc"]
    }

    # TxDB, BSgenome, or OrgDB package?
    if(type %in% c("TxDb", "BSgenome"))
    {
        pkg.string <- paste0("^", type, ".", bioc.id, "[a-z]+.UCSC.", ucsc.id)
        pkg <- grep(pkg.string, .availableOrgPkgs(type), value=TRUE)
        if(length(pkg) == 0)
            pkg <- grep(pkg.string, .availableOrgPkgs(type, local=FALSE), value=TRUE)
        if(length(pkg) == 0)
            stop(paste("No corresponding", type, "package for", org))
        else if(length(pkg) > 1)
        {
            message("Found several genome assemblies")
            for(p in pkg) message(paste(which(pkg==p), p, sep=": "))
            n <- readline(paste0("Choose assembly (1-", length(pkg),") : "))
            pkg <- pkg[as.integer(n)]

            #message("Found several genome assemblies")
            #message(paste("Using latest:", pkg))
            #ver <- sapply(pkg, 
            #    function(p)
            #    {
            #        spl <- unlist(strsplit(p, "\\."))
            #        ind <- length(spl)
            #        if(type == "TxDb") ind <- ind - 1
            #        ass <- spl[ind]
            #        ver <- sub("^[a-zA-Z]+", "", ass)
            #        return(as.integer(ver))
            #    })
            #pkg <- pkg[which.max(ver)]
        }
    }
    else
    {
        id.type <- .getOrgIdType(bioc.id)
        pkg <- paste("org", bioc.id, id.type, "db", sep=".")
    }
    return(pkg)
}

.stdArgs <- function(call, formals) {
    callargs <- as.list(call)[-1]
    toadd <- setdiff(names(formals), names(callargs))
    call[toadd] <- formals[toadd]
    call
}

.matchArgs <- function(fun, call, ...) {
    funfor <- formals(fun)
    exargs <- intersect(names(funfor), names(call))
    c(as.list(call)[-1][exargs], ...)
}

.execArgs <- function(method, mcall, forms, add.args)
{
    call <- .stdArgs(mcall, forms)
    exargs <- .matchArgs(method, call, add.args)
    do.call(method, exargs)
}
