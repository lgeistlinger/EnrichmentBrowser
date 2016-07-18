############################################################
# 
# author: Ludwig Geistlinger
# date: 2015-03-10 13:32:37
# 
# descr: misc utils
# 
############################################################

.isAvailable <- function(pkg, type="annotation")
{
    if(!(pkg %in% .packages(all.available=TRUE)))
    {   
        message(paste0("Corresponding ", type,  " package not found: ", 
            pkg, "\nMake sure that you have it installed."))
        choice <- readline("Install it now? (y/n): ")
        if(choice == "y")
        {   
            biocLite <- NULL
            source("http://bioconductor.org/biocLite.R")
            biocLite(pkg, suppressUpdates=TRUE, suppressAutoUpdate=TRUE)
        }   
        else stop(paste("Package", pkg, "is not available"))
    }   
    require(pkg, character.only = TRUE)
}

.autoDetectGeneIdType <- function(id)
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
    else pkgs <- available.packages(paste0("http://bioconductor.org/",
        "packages/release/data/annotation/src/contrib"))[, "Package"]
    
    type <- match.arg(type)
    org.string <- "^org.[A-z][a-z]+.[a-z]+.db$"
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
            sapply(pkg, function(p) 
                message(paste(which(pkg==p), p, sep=": ")))
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


