############################################################
#
# author: Ludwig Geistlinger
# date: 3 Feb 2011
#
# GGEA - Gene Graph Enrichment Analysis
#
# Update, 09 May 2014: Extension to network-based enrichment
#           analysis 
#
############################################################

nbea.methods <- function() c("ggea", "nea",  "spia")

nbea <- function(
    method=nbea.methods(), 
    eset, 
    gs, 
    grn,
    alpha=0.05, 
    perm=1000, 
    padj.method="BH",
    out.file=NULL,
    browse=FALSE, ...)
{
    GS.MIN.SIZE <- config.ebrowser("GS.MIN.SIZE")
    GS.MAX.SIZE <- config.ebrowser("GS.MAX.SIZE")
    GSP.COL <- config.ebrowser("GSP.COL")

    # restrict eset and gs to intersecting genes
    igenes <- intersect(featureNames(eset), unique(unlist(gs)))
    eset <- eset[igenes,]
    gs <- sapply(gs, function(s) s[s%in% igenes]) 
    lens <- sapply(gs, length)
    gs <- gs[lens >= GS.MIN.SIZE & lens <= GS.MAX.SIZE]

    if(class(method) == "character")
    {
        method <- match.arg(method)
        #if(length(find(method))) res.tbl <- do.call(method, 
        #        list(eset=eset, gs=gs, grn=grn, alpha=alpha, perm=perm, ...))
        #else 
        #if(!(method %in% nbea.methods())) 
        #    stop(paste("\'method\' must be one out of {", 
        #        paste(nbea.methods(), collapse=", "), "}"))
        #else 
        if(method == "nea") 
            res.tbl <- nea.wrapper(eset=eset, 
                gs=gs, grn=grn, alpha=alpha, perm=perm, ...)
        else if(method == "spia") 
            res.tbl <- spia.wrapper(eset=eset, 
                gs=gs, alpha=alpha, perm=perm, ...) 
        else res.tbl <- ggea(eset=eset, gs=gs, 
                grn=grn, alpha=alpha, perm=perm, ...)      
    }
    else if(class(method) == "function") 
        res.tbl <- method(eset=eset, gs=gs, grn=grn, alpha=alpha, perm=perm, ...)
    else stop(paste(method, "is not a valid method for nbea"))

    res.tbl <- data.frame(signif(res.tbl, digits=3))
    sorting.df <- cbind(res.tbl[,ncol(res.tbl)], 
        -res.tbl[,rev(seq_len(ncol(res.tbl)-1))])
    res.tbl <- res.tbl[do.call(order, as.data.frame(sorting.df)),]
    
    res.tbl[,GSP.COL] <- p.adjust(res.tbl[,GSP.COL], method=padj.method)

    res.tbl <- DataFrame(rownames(res.tbl), res.tbl)
    colnames(res.tbl)[1] <- config.ebrowser("GS.COL") 
    rownames(res.tbl) <- NULL
    
    if(!is.null(out.file))
    {
        write.table(res.tbl, 
            file=out.file, quote=FALSE, row.names=FALSE, sep="\t")
        message(paste("Gene set ranking written to", out.file)) 
    }
    else
    {
        res <- list(
            res.tbl=res.tbl, method=method,
            nr.sigs=sum(res.tbl[,GSP.COL] < alpha),
            eset=eset, gs=gs, alpha=alpha)

        if(browse) ea.browse(res)
        else return(res)
    }
}

nea.wrapper <- function(eset, gs, grn, alpha=0.05, perm=100)
{
    ADJP.COL <- config.ebrowser("ADJP.COL")
    GSP.COL <- config.ebrowser("GSP.COL")

    if(perm > 100) perm <- 100
    ags <- featureNames(eset)[fData(eset)[,ADJP.COL] < alpha]
    grn <- unique(grn[,1:2])
    gs.genes <- unique(unlist(gs))
    grn <- grn[(grn[,1] %in% gs.genes) & (grn[,2] %in% gs.genes),]
    network <- apply(grn, 1, function(x) paste(x, collapse=" "))
    message("Computing NEA permutations, this may take a few minutes ...")
    res <- neaGUI::nea(ags=ags, fgs=gs, network=network, nperm=perm)
    res <- res$MainResult
    res <- res[, c("Number_of_Genes", 
        "Number_of_AGS_genes", "Number_links", "Z_score", "P_value")]
    res <- res[order(res[,"Z_score"], decreasing=TRUE), ]
    colnames(res) <- sub("AGS", "de", colnames(res))
    colnames(res) <- sub("Number", "Nr", colnames(res))
    colnames(res) <- sub("_of", "", colnames(res))
    colnames(res) <- gsub("_", ".", colnames(res))
    colnames(res) <- toupper(colnames(res))
    res <- as.matrix(res)
    res <- res[order(res[,GSP.COL]),]
    return(res) 
}

spia.wrapper <- function(eset, gs, alpha=0.05, perm=1000)
{
    FC.COL <- config.ebrowser("FC.COL")
    ADJP.COL <- config.ebrowser("ADJP.COL")
    GSP.COL <- config.ebrowser("GSP.COL")

    de.genes <- fData(eset)[,ADJP.COL] < alpha
    de <- fData(eset)[de.genes, FC.COL]
    names(de) <- featureNames(eset)[de.genes]
    all <- featureNames(eset)
    organism <- substring(names(gs)[1],1,3) 
    res <- SPIA::spia(de=de, all=all, organism=organism, nB=perm)
    res[,"Name"] <- gsub(" ", "_", res[,"Name"])
    rownames(res) <- paste(paste0(organism, res[,"ID"]), res[,"Name"], sep="_")
    res <- res[, c("pSize", "NDE", "tA", "Status", "pG")]
    colnames(res) <- c("SIZE", "NDE", "T.ACT", "STATUS", GSP.COL)
    res[,"STATUS"] <- ifelse(res[,"STATUS"] == "Activated", 1, -1)
    res <- as.matrix(res)
    return(res)
}























