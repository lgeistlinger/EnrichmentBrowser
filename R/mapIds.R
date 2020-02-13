############################################################
# 
# author: Ludwig Geistlinger
# date: 2015-12-10 13:32:37
# 
# descr: id mapping
#
# UPDATE 2018-08-01: adding solutions for 1:many and
#   many:1 mappings
# 
############################################################

#' Mapping between gene ID types 
#' 
#' Functionality to map between common gene ID types such as ENSEMBL and ENTREZ
#' for gene expression datasets, gene sets, and gene regulatory networks.
#' 
#' The function 'idTypes' lists the valid values which the arguments 'from'
#' and 'to' can take. This corresponds to the names of the available gene ID
#' types for the mapping.
#' 
#' @aliases map.ids
#' @param obj The object for which gene IDs should be mapped. Supported options
#'  include \itemize{
#' \item Gene expression dataset. An object of class 
#' \code{\linkS4class{SummarizedExperiment}}.
#' Expects the names to be of gene ID type given in argument \code{from}.
#' \item Gene sets. Either a list of gene sets (character vectors of gene
#' IDs) or a \code{\linkS4class{GeneSetCollection}} storing all gene sets.
#' \item Gene regulatory network. A 3-column character matrix;
#' 1st col = IDs of regulating genes; 2nd col = IDs of regulated genes; 
#' 3rd col = regulation effect; Use '+' and '-' for activation / inhibition.}
#' @param org Character. Organism in KEGG three letter code, e.g. \sQuote{hsa}
#' for \sQuote{Homo sapiens}.  See references.
#' @param from Character. Gene ID type from which should be mapped.  Corresponds
#' to the gene ID type of argument \code{obj}. Defaults to \code{ENSEMBL}.
#' @param to Character. Gene ID type to which should be mapped. Corresponds to 
#' the gene ID type the argument \code{obj} should be updated with.
#' If \code{obj} is an expression dataset of class 
#' \code{\linkS4class{SummarizedExperiment}}, \code{to} can also be the name of 
#' a column in the \code{\link{rowData}} 
#' slot to specify user-defined mappings in which conflicts have been 
#' manually resolved. Defaults to \code{ENTREZID}.
#' @param multi.to How to resolve 1:many mappings, i.e. multiple to.IDs for a 
#' single from.ID? This is passed on to the \code{multiVals} argument of 
#' \code{\link{mapIds}} and can thus take several pre-defined values, but also
#' the form of a user-defined function. However, note that this requires that a 
#' single to.ID is returned for each from.ID. Default is \code{"first"},
#' which accordingly returns the first to.ID mapped onto the respective from.ID.
#' @param multi.from How to resolve many:1 mappings, i.e. multiple from.IDs 
#' mapping to the same to.ID? Only applicable if \code{obj} is an expression 
#' dataset of class \code{\linkS4class{SummarizedExperiment}}. 
#' Pre-defined options include: 
#' \itemize{ \item 'first' (Default): returns the first from.ID for each to.ID
#' with multiple from.IDs, 
#' \item 'minp': selects the from.ID with minimum p-value (according to the 
#' \code{\link{rowData}} column \code{PVAL} of \code{obj}),
#' \item 'maxfc': selects the from.ID with maximum absolute log2 fold change 
#' (according to the \code{\link{rowData}} column \code{FC} of \code{obj}).}
#' Note that a user-defined function can also be supplied for custom behaviors.
#' This will be applied for each case where there are multiple from.IDs for a 
#' single to.ID, and accordingly takes the arguments \code{ids} and \code{obj}. 
#' The argument \code{ids} corresponds to the multiple from.IDs from which a 
#' single ID should be chosen, e.g. via information available in argument 
#' \code{obj}. See examples for a case where ids are selected based on a 
#' user-defined \code{\link{rowData}} column.
#' @return idTypes: character vector listing the available gene ID types for
#' the mapping;
#' 
#' idMap: An object of the same class as the input argument \code{obj}, i.e.  
#' a \code{\linkS4class{SummarizedExperiment}} if provided an expression dataset,
#' a list of character vectors or a \code{\linkS4class{GeneSetCollection}} if 
#' provided gene sets, and a character matrix if provided a gene regulatory network.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{\linkS4class{SummarizedExperiment}}, \code{\link{mapIds}},
#' \code{\link{keytypes}}
#' @references KEGG Organism code
#' \url{http://www.genome.jp/kegg/catalog/org_list.html}
#' @examples
#' 
#'     # (1) ID mapping for gene expression datasets 
#'     # create an expression dataset with 3 genes and 3 samples
#'     se <- makeExampleData("SE", nfeat = 3, nsmpl = 3)
#'     names(se) <- paste0("ENSG00000000", c("003", "005", "419"))
#'     idMap(se, org = "hsa")
#'
#'     # user-defined mapping
#'     rowData(se)$MYID <- c("g1", "g1", "g2")
#'     idMap(se, to = "MYID")    
#'
#'     # data-driven resolving of many:1 mappings
#'     
#'     ## e.g. select from.ID with lowest p-value
#'     pcol <- configEBrowser("PVAL.COL")
#'     rowData(se)[[pcol]] <- c(0.001, 0.32, 0.15)
#'     idMap(se, to = "MYID", multi.from = "minp") 
#'    
#'     ## ... or using a customized function
#'     maxScore <- function(ids, se)
#'     {
#'          scores <- rowData(se)[ids, "SCORE"]
#'          ind <- which.max(scores)
#'          return(ids[ind])
#'     }
#'     rowData(se)$SCORE <- c(125.7, 33.4, 58.6)
#'     idMap(se, to = "MYID", multi.from = maxScore) 
#'            
#'     # (2) ID mapping for gene sets 
#'     # create two gene sets containing 3 genes each 
#'     s2 <- paste0("ENSG00000", c("012048", "139618", "141510"))
#'     gs <- list(s1 = names(se), s2 = s2)
#'     idMap(gs, org = "hsa", from = "ENSEMBL", to = "SYMBOL")    
#'
#'     # (3) ID mapping for gene regulatory networks
#'     grn <- cbind(FROM = gs$s1, TO = gs$s2, TYPE = rep("+", 3))
#'     idMap(grn, org = "hsa", from = "ENSEMBL", to = "ENTREZID")  
#' 
#' @export idMap
idMap <- function(obj, org=NA, 
    from="ENSEMBL", to="ENTREZID", 
    multi.to="first", multi.from="first")
{
    isSE <- is(obj, "SummarizedExperiment") || is(obj, "ExpressionSet")
    if(isSE) res <- .idMapSE(obj, org, from, to, multi.to, multi.from) 
    else if(is(obj, "GeneSetCollection")) res <- .idMapGSC(obj, org, from, to)
    else if(is.list(obj)) res <- .idMapGS(obj, org, from, to, multi.to)
    else if(is.matrix(obj)) res <- .idMapGRN(obj, org, from, to, multi.to)
    else stop("Invalid input type for argument \'obj\'") 
    return(res)
}

.idMapSE <- function(se, org=NA, 
    from="ENSEMBL", to="ENTREZID", 
    multi.to="first", multi.from="first")
{
    if(is(se, "ExpressionSet")) se <- as(se, "SummarizedExperiment")
  
    # user-defined mapping in rowData? 
    is.col <- to %in% colnames(rowData(se)) 
    if(!is.col)
    {
        # get annotation
        anno <- org
        if(is.na(anno)) anno <- metadata(se)$annotation
        if(!length(anno)) stop("Organism under investigation not annotated")
        se <- .annotateSE(se, anno)
    }
    else message(paste("Using rowData column", to))

    ids <- names(se)
    # simple ID-based resolving: needs anno, from & to as char
    is.first <- is.character(multi.from) && multi.from == "first" 
    if(is.first && !is.col) x <- .idmap(ids, anno, from, to, multi.to=multi.to) 
    # data-driven resolving 
    else
    { 
	    # 'to' can also be a rowData col  
	    if(!is.col) # needs org, from, to as char
        {
            x <- .idmap(ids, anno, from, to, excl.na=FALSE, 
                            multi.to=multi.to, resolve.multiFrom=FALSE)
            rowData(se)[[to]] <- unname(x)
        }
        x <- .idmapRD(se, to, multi.from) # only needs to as a col
        
        # remove to col
        ind <- colnames(rowData(se)) != to 
        rowData(se) <- rowData(se)[, ind, drop=FALSE]
    }

    se <- se[names(x), ]
    
    # for trace back
    rowData(se)[[from]] <- names(x)
    message(paste("Mapped from.IDs have been added to the rowData column", from))
    
    names(se) <- unname(x)
    return(se)
}

.idMapGRN <- function(grn, org, from, to, multi.to)
{
    sgenes <- unique(as.vector(grn[,1:2]))
    sgenes <- .mapStats(sgenes, org, from, to, multi.to)
    for(i in 1:2)
    {
        grn[,i] <- unname(sgenes[grn[,i]])
        grn <- grn[!is.na(grn[,i]),] 
    }
    return(grn)
}

.idMapGS <- function(gs, org, from, to, multi.to = "first")
{
    sgenes <- unique(unlist(gs))
    sgenes <- .mapStats(sgenes, org, from, to, multi.to)
    gs <- lapply(gs, function(s) unname(sgenes[s]))
    lapply(gs, function(s) s[!is.na(s)])
}

.mapStats <- function(sgenes, org, from, to, multi.to)
{
    orgpkg <- .getAnnoPkg(org)
    sgenes <- sgenes[!is.na(sgenes)]
    sgenes <- sgenes[sgenes != ""]
    suppressMessages(
        sgenes <- AnnotationDbi::mapIds(orgpkg, keys = sgenes, 
                                                column = to, 
                                                keytype = from, 
                                                multiVals = "list")
    )
    sgenes <- .resolveMultiTo(sgenes, orgpkg, from, to, multi.to)
    nr.na <- sum(is.na(sgenes))
	if(nr.na) message(paste("Excluded", nr.na, 
                            "from.IDs without a corresponding to.ID"))
    nr.mf <- sum(table(sgenes) > 1)
    if(nr.mf) message(paste("Encountered", nr.mf, "to.IDs with >1 from.ID"))
    return(sgenes)
}

.idMapGSC <- function(gsc, org, from, to)
{
    if(is.na(org)) org <- GSEABase::organism(gsc[[1]])
    if(!length(org)) stop("Organism under investigation not annotated")

    from <- .createIdentifier(from)
    to <- .createIdentifier(to, org)

    .map <- function(s) GSEABase::mapIdentifiers(s, from=from, to=to)
    mgsc <- lapply(gsc, .map)
    GSEABase::GeneSetCollection(mgsc)
}

# if 
#map <- getAnnMap("ENTREZID", "org.Sc.sgd")
#gs1m <- mapIdentifiers(gs1, from=map, to=EntrezIdentifier())

.createIdentifier <- function(idtype, org=NA)
{
    idtypes <- c("Entrez", "Enzyme", "ENSEMBL", 
                    "Genename", "Refseq", "Symbol", "Unigene", "Uniprot")
    idtypes.uc <- toupper(idtypes)
    if(idtype == "ENTREZID") idtype <- "ENTREZ"

    ind <- match(idtype, idtypes.uc)
    if(is.na(ind))
    { 
        istr <- paste(idtypes.uc, collapse=", ")
        stop(paste("ID type must be one of", istr))
    }
    idobj <- idtypes[ind]
    idobj <- paste0(idobj, "Identifier")

    args <- list()   
    if(!is.na(org))
    { 
        org.pkg <- .org2pkg(org)
        args$annotation <- org.pkg
    }
    idobj <- do.call(idobj, args)
    return(idobj)
}

####

#' @export
#' @keywords internal
map.ids <- function(se, org=NA, from="ENSEMBL", to="ENTREZID")
{
    .Deprecated("idMap")
    idMap(se, org, from, to)
}

#' @rdname idMap
#' @export
idTypes <- function(org)
{
    anno.pkg <- .getAnnoPkg(org)
    return(AnnotationDbi::keytypes(anno.pkg))
}

.idmap <- function(ids, anno, from, to, 
    excl.na=TRUE, multi.to="first", resolve.multiFrom=TRUE)
{
    anno.pkg <- .getAnnoPkg(anno) 
    suppressMessages(
        x <- AnnotationDbi::mapIds(anno.pkg, 
                keys=ids, keytype=from, column=to, multiVals="list")
    )
	# case 1: multiple to.IDs (1:n) -> select one
    x <- .resolveMultiTo(x, anno.pkg, from, to, multi.to)
 
    # case 2: no to.ID -> exclude
    if(excl.na) x <- .exclNaIds(x)
        
	# case 3: multiple from.IDs (n:1) -> select one
    if(resolve.multiFrom) x <- .getFirstToId(x)	    
    
    return(x)
}

.getAnnoPkg <- function(anno)
{
    # org or chip?
    is.org <- grepl("^[a-z]{3}$", anno)
    if(is.org) anno.pkg <- .org2pkg(anno)
    else anno.pkg <- paste0(anno, ".db")
    
    isAvailable(anno.pkg)
    anno.pkg <- get(anno.pkg)
    return(anno.pkg)
}

.annotateSE <- function(se, anno)
{
    is.org <- grepl("^[a-z]{3}$", anno)
    if(!is.org)
    { 
        anno.pkg <- .getAnnoPkg(anno)
        anno <- .annoPkg2Org(anno.pkg)
    }
    metadata(se)$annotation <- anno
    return(se)
}

.resolveMultiTo <- function(x, anno.pkg, from, to, multi.to)
{
    nr.multi <- sum(lengths(x) > 1)
    if(nr.multi)
    { 
        message(paste("Encountered", nr.multi, 
                "from.IDs with >1 corresponding to.ID")) 
        ids <- names(x)
        suppressMessages(
            x <- AnnotationDbi::mapIds(anno.pkg, keys=ids, 
                        keytype=from, column=to, multiVals=multi.to)
        )
        # valid mapping?
        is.valid <- is.vector(x) && is.character(x) && length(x) == length(ids)
        if(!is.valid)
        { 
            e <- paste("\'multi.to\' is required to return",
                        "a single to.ID for each from.ID")
            stop(e)
        }

        # verbose
        prefix <- "(the"
        if(!is.character(multi.to))
        {
            multi.to <- "user-defined" 
            prefix <- "(a"
        }
        m <- paste(prefix, multi.to, "to.ID was chosen for each of them)")   
        message(m)
    }
    else x <- unlist(x)
    return(x)
}

.exclNaIds <- function(x)
{
    isna <- is.na(x)
    nr.na <- sum(isna)
	if(nr.na)
    { 
        message(paste("Excluded", nr.na, "from.IDs without a corresponding to.ID"))
        x <- x[!isna]
    }
    return(x)
}

.getFirstToId <- function(x)
{
    nr.dupl <- sum(table(x) > 1)
	if(nr.dupl)
    { 
        message(paste("Encountered", nr.dupl, 
                "to.IDs with >1 from.ID", 
                "(the first from.ID was chosen for each of them)"))
        x <- x[!duplicated(x)]
    }
    return(x)
}

.idmapRD <- function(se, to.col="ENTREZID", multiFUN="first")
{
	stopifnot(is(se, "SummarizedExperiment"))
	stopifnot(to.col %in% colnames(rowData(se)))
	stopifnot(is.character(multiFUN) || is.function(multiFUN))
    
    fun.opts <- c("first", "minp", "maxfc")
    if(is.character(multiFUN)) stopifnot(multiFUN %in% fun.opts)

	x <- rowData(se)[[to.col]]
    names(x) <- names(se)

	# no to.ID -> exclude
    x <- .exclNaIds(x)

    # to.IDs with > 1 from.ID: select / summarize according to multiFUN
    x <- .resolveDuplicateIDs(x, multiFUN, se)   	
    return(x)
}
    
.resolveDuplicateIDs <- function(x, multiFUN, se)
{
    dupl.tab <- table(x) > 1
	nr.dupl <- sum(dupl.tab)
	if(nr.dupl)
    { 
        message(paste("Encountered", nr.dupl, "to.IDs with >1 from.ID")) 

        # dupls: duplicated entrez 
        dupls <- names(dupl.tab)[dupl.tab]
        
        # for dupls: list entrez -> multiple ensembl        
        dl <- lapply(dupls, function(d) names(x)[x==d])
        names(dl) <- dupls
       
        # unique entrez 
        uto <- unique(x)
        from.grid <- seq_along(x)
        names(from.grid) <- unname(x)
        
        .resolveID <- function(id)
        {
            if(id %in% names(dl)) id <- .chooseID(dl[[id]], se, multiFUN)
            else 
            {
                ind <- from.grid[id]
                id <- names(x)[ind]
                # faster than: id <- names(x)[match(id, x)]
            }
            return(id)
        }

        # resolve dupls
        ufrom <- vapply(uto, .resolveID, character(1))
        x <- names(ufrom)
        names(x) <- unname(ufrom)
    
        prefix <- "(the"
        if(is.function(multiFUN))
        { 
            multiFUN <- "user-defined"
            prefix <- "(a"
        } 
        m <- paste(prefix, multiFUN, "from.ID was chosen for each of them)")   
        message(m)
	}
    return(x)
}


# choose a from.ID data-driven
.chooseID <- function(ids, se, multiFUN)
{
    if(is.function(multiFUN)) return(multiFUN(ids, se)) 
    else if(multiFUN == "first") ind <- 1
    else if(multiFUN == "minp")
	{
        pcol <- configEBrowser("PVAL.COL")
        stopifnot(pcol %in% colnames(rowData(se)))
        ps <- rowData(se)[ids, pcol]
        ind <- which.min(ps)
	}
	else if(multiFUN == "maxfc")
	{
	    fccol <- configEBrowser("FC.COL")
        stopifnot(fccol %in% colnames(rowData(se)))
        fcs <- rowData(se)[ids, fccol]
        ind <- which.max(abs(fcs))
    }
    else stop("Invalid multiFUN")
    return(ids[ind])
}

