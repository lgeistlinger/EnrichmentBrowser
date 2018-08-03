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

#' Mapping between gene ID types for the rownames of a SummarizedExperiment
#' 
#' Functionality to map the rownames of a SummarizedExperiment between common
#' gene ID types such as ENSEMBL and ENTREZ.
#' 
#' The function 'idTypes' lists the valid values which the arguments 'from'
#' and 'to' can take. This corresponds to the names of the available gene ID
#' types for the mapping.
#' 
#' @aliases map.ids
#' @param se An object of class \code{\linkS4class{SummarizedExperiment}}.
#' Expects the names to be of gene ID type given in argument 'from'.
#' @param org Organism in KEGG three letter code, e.g. \sQuote{hsa} for
#' \sQuote{Homo sapiens}.  See references.
#' @param from Gene ID type from which should be mapped.  Corresponds to the
#' gene ID type of the names of argument 'se'. Note that 'from' is ignored if 
#' 'to' is a \code{\link{rowData}} column of 'se'. Defaults to 'ENSEMBL'.
#' @param to Gene ID type to which should be mapped. Corresponds to the gene
#' ID type the rownames of argument 'se' should be updated with.
#' Note that this can also be the name of a column in the \code{\link{rowData}} 
#' slot to specify user-defined mappings in which conflicts have been manually
#' resolved. Defaults to 'ENTREZID'.
#' @param multi.to How to resolve 1:many mappings, i.e. multiple to.IDs for a 
#' single from.ID? This is passed on to the \code{multiVals} argument of 
#' \code{\link{mapIds}} and can thus take several pre-defined values, but also
#' the form of a user-defined function. However, note that this requires that a 
#' single to.ID is returned for each from.ID. Default is 'first', which accordingly
#' returns the first to.ID the respective from.ID has been mapped to. 
#' @param multi.from How to resolve many:1 mappings, i.e. multiple from.IDs map
#' to the same to.ID? Pre-defined options include: 
#' \itemize{ \item 'first' (Default): returns the first from.ID for each to.ID
#' with multiple from.IDs 
#' \item 'minp' selects the from.ID with minimum p-value (according to the 
#' \code{\link{rowData}} column 'PVAL' of 'se')
#' \item 'maxfc' selects the from.ID with maximum absolute log2 fold change 
#' (according to the \code{\link{rowData}} column 'FC' of 'se').}
#' Note that a user-defined function can also be supplied for custom behaviors.
#' This will be applied for each case where there are multiple from.IDs for a 
#' single to.ID, and accordingly takes the arguments 'ids' and 'se'. 
#' The argument 'ids' are the multiple from.IDs from which a single one should 
#' be chosen e.g. data-driven via information available in argument 'se'. 
#' See Examples for a case where ids are selected based on a user-defined  
#' \code{\link{rowData}} column.
#' @return idTypes: character vector listing the available gene ID types for
#' the mapping;
#' 
#' idMap: An object of \code{\linkS4class{SummarizedExperiment}}.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{\linkS4class{SummarizedExperiment}}, \code{\link{mapIds}},
#' \code{\link{keytypes}}
#' @references KEGG Organism code
#' \url{http://www.genome.jp/kegg/catalog/org_list.html}
#' @examples
#' 
#'     # create an expression dataset with 3 genes and 3 samples
#'     se <- makeExampleData("SE", nfeat=3, nsmpl=3)
#'     names(se) <- paste0("ENSG00000000", c("003","005", "419"))
#'     mse <- idMap(se, anno="hsa")
#'
#'     # user-defined mapping
#'     rowData(se)$MYID <- c("g1", "g1", "g2")
#'     mse <- idMap(se, to="MYID")    
#'
#'     # data-driven resolving of many:1 mappings
#'     
#'     ## e.g. select from.ID with lowest p-value
#'     pcol <- configEBrowser("PVAL.COL")
#'     rowData(se)[[pcol]] <- c("0.001", "0.32", "0.15")
#'     mse <- idMap(se, to="MYID", multi.from="minp") 
#'    
#'     ## ... or using a customized function
#'     maxScore <- function(ids, se)
#'     {
#'          scores <- rowData(se, use.names=TRUE)[ids, "SCORE"]
#'          ind <- which.max(scores)
#'          return(ids[ind])
#'     }
#'     rowData(se)$SCORE <- c("125.7", "33.4", "58.6")
#'     mse <- idMap(se, to="MYID", multi.from=maxScore) 
#'            
#'
#' 
#' @export idMap
idMap <- function(se, org=NA, 
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
    }

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
            se <- .annotateSE(se, anno)
        }
        x <- .idmapSE(se, to, multi.from) # only needs to as a col
        
        # remove to col
        ind <- match(to, colnames(rowData(se)))
        rowData(se) <- rowData(se)[,-ind]
    }

    se <- se[names(x), ]
    names(x) <- NULL
    names(se) <- x
    
    
    return(se)
}

#' @export
#' @keywords internal
map.ids <- function(se, org=NA, from="ENSEMBL", to="ENTREZID")
{
    .Deprecated("idMap")
    idMap(se, org, from, to)
}

#' @rdname idMap
#' @export
idTypes <- function(anno)
{
    anno.pkg <- .getAnnoPkg(anno)
    return(keytypes(anno.pkg))
}

.idmap <- function(ids, anno, from, to, 
    excl.na=TRUE, multi.to="first", resolve.multiFrom=TRUE)
{
    anno.pkg <- .getAnnoPkg(anno) 
    suppressMessages(
        x <- mapIds(anno.pkg, keys=ids, keytype=from, column=to, multiVals="list")
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
            x <- mapIds(anno.pkg, keys=ids, 
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

.idmapSE <- function(se, to.col="ENTREZID", multiFUN="first")
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
        ps <- rowData(se, use.names=TRUE)[ids, pcol]
        ind <- which.min(ps)
	}
	else if(multiFUN == "maxfc")
	{
	    fccol <- configEBrowser("FC.COL")
        stopifnot(fccol %in% colnames(rowData(se)))
        fcs <- rowData(se, use.names=TRUE)[ids, fccol]
        ind <- which.max(abs(fcs))
    }
    else stop("Invalid multiFUN")
    return(ids[ind])
}

