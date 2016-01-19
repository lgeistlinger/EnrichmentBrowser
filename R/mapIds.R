############################################################
# 
# author: Ludwig Geistlinger
# date: 2015-12-10 13:32:37
# 
# descr: id mapping
# 
############################################################

id.types <- function(org)
{
    org.pkg <- .org2pkg(org)
    .isAvailable(org.pkg)
    org.pkg <- get(org.pkg) 
    return(keytypes(org.pkg))
}

map.ids <- function(eset, org=NA, from="ENSEMBL", to="ENTREZID")
{
    if(is.na(org)) org <- annotation(eset)
    if(!length(org)) stop("Organism under investigation not annotated")
    org.pkg <- .org2pkg(org)
    .isAvailable(org.pkg)
    org.pkg <- get(org.pkg) 
    x <- mapIds(org.pkg, keys=featureNames(eset), keytype=from, column=to)
    nr.na <- sum(is.na(x))
    if(nr.na)
    { 
        message(paste("Excluded", nr.na, "genes without a corresponding to.ID"))
        x <- x[!is.na(x)]
    }
    nr.dupl <- sum(table(x) > 1)
    if(nr.dupl)
    { 
        message(paste("Encountered", nr.dupl, 
            "from.IDs with >1 corresponding to.ID", 
            "(a single to.ID was chosen for each of them)"))
        x <- x[!duplicated(x)]
    }
    eset <- eset[featureNames(eset) %in% names(x), ]
    names(x) <- NULL
    featureNames(eset) <- x
    if(!length(annotation(eset))) annotation(eset) <- org
    return(eset)
}
