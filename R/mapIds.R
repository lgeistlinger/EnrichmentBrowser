############################################################
# 
# author: Ludwig Geistlinger
# date: 2015-12-10 13:32:37
# 
# descr: id mapping
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
#' gene ID type of the names of argument 'se'.  Defaults to 'ENSEMBL'.
#' @param to Gene ID type to which should be mapped.  Corresponds to the gene
#' ID type the featuresNames of argument 'se' should be updated with.  Defaults
#' to 'ENTREZID'.
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
#'     se <- idMap(se, org="hsa")
#' 
#' @export idMap
idMap <- function(se, org=NA, from="ENSEMBL", to="ENTREZID")
{
	# 'to' should also accept a rowData col  

    ### TEMPORARY: will be replaced by as(eSet,SummarizedExperiment)
    if(is(se, "ExpressionSet")) se <- as(se, "RangedSummarizedExperiment")
    ###

    if(is.na(org)) org <- metadata(se)$annotation
    if(!length(org)) stop("Organism under investigation not annotated")
    ids <- rownames(se)
    x <- .idmap(ids, org, from, to) 
    se <- se[names(x), ]
    names(x) <- NULL
    rownames(se) <- x
    if(!length(metadata(se)$annotation)) metadata(se)$annotation <- org
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
idTypes <- function(org)
{
    org.pkg <- .org2pkg(org)
    isAvailable(org.pkg)
    org.pkg <- get(org.pkg) 
    return(keytypes(org.pkg))
}

.idmap <- function(ids, org, from, to)
{
    org.pkg <- .org2pkg(org)
    isAvailable(org.pkg)
    org.pkg <- get(org.pkg) 
    x <- mapIds(org.pkg, keys=ids, keytype=from, column=to, multiVals="list")

	# case 1: multiple to.IDs (1:n) -> take one / first
    nr.multi <- sum(lengths(x) > 1)
    if(nr.multi)
    { 
        message(paste("Encountered", nr.multi, 
            "from.IDs with >1 corresponding to.ID", 
            "(a single to.ID was chosen for each of them)"))
    }
    x <- vapply(x, function(i) i[1], character(1)) 

	# case 2: no to.ID -> exclude
    nr.na <- sum(is.na(x))
    if(nr.na)
    { 
        message(paste("Excluded", nr.na, "genes without a corresponding to.ID"))
        x <- x[!is.na(x)]
    }

	# case 3: multiple from.IDs (n:1) -> select / summarize
	nr.dupl <- sum(table(x) > 1)
    if(nr.dupl)
    { 
        message(paste("Encountered", nr.dupl, 
            "from.IDs that map to the same to.ID", 
            "(a single from.ID was chosen for each of them)"))
        x <- x[!duplicated(x)]
    }
    return(x)
}    
