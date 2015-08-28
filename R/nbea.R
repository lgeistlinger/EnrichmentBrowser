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

nbea.methods <- function() c("ggea", "nea",  "spia", "pathnet")

nbea <- function(
    method=nbea.methods(), 
    eset, 
    gs, 
    grn,
    alpha=0.05, 
    perm=1000, 
    padj.method="none",
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
        if(method == "nea") res.tbl <- nea.wrapper(eset, gs, grn, alpha, perm)
        else if(method == "spia") res.tbl <- spia.wrapper(eset, gs, alpha, perm)
        else if(method == "pathnet") 
            res.tbl <- pathnet.wrapper(eset, gs, grn, alpha, perm)
        else res.tbl <- ggea(eset, gs, grn, alpha, perm=perm, ...)      
    }
    else if(class(method) == "function") 
        res.tbl <- method(eset=eset, gs=gs, grn=grn, alpha=alpha, perm=perm, ...)
    else stop(paste(method, "is not a valid method for nbea"))

    res.tbl <- data.frame(signif(res.tbl, digits=3))
    sorting.df <- res.tbl[,ncol(res.tbl)]
    if(ncol(res.tbl) > 1)
        sorting.df <- cbind(sorting.df, -res.tbl[,rev(seq_len(ncol(res.tbl)-1))])
    else colnames(res.tbl)[1] <- GSP.COL
    res.tbl <- res.tbl[do.call(order, as.data.frame(sorting.df)), , drop=FALSE]
    
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
            method=method, res.tbl=res.tbl, 
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

    #if(perm > 100) perm <- 100
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

spia.wrapper <- function(eset, gs, alpha=0.05, perm=1000, beta=1)
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

pathnet.wrapper <- function(eset, gs, grn, alpha=0.05, perm=1000)
{
    ADJP.COL <- config.ebrowser("ADJP.COL")
    GSP.COL <- config.ebrowser("GSP.COL")

    dir.evid <- -log(fData(eset)[,ADJP.COL], base=10)
    dir.evid <- cbind(as.integer(featureNames(eset)), dir.evid)
    colnames(dir.evid) <- c("Gene.ID", "Obs")
    adjm <- grn2adjm(grn)
    pwy.dat <- extr.pwy.dat(gs, grn)
    
    res <- PathNet::PathNet(
            #Enrichment_Analysis = TRUE, 
            #Contextual_Analysis = FALSE, 
            DirectEvidence_info = dir.evid, 
            Column_DirectEvidence = 2,
            Adjacency = adjm, 
            pathway = pwy.dat, 
            n_perm = perm, 
            threshold = alpha)#,
            #use_sig_pathways  = FALSE)

    res <- res$enrichment_results[, 
        c("Name", "No_of_Genes", "Sig_Direct", "Sig_Combi", "p_PathNet")]
    res[,1] <- sapply(res[,1], 
        function(s) grep(unlist(strsplit(s,"_"))[1], names(gs), value=TRUE))
    rownames(res) <- res[,1]
    res <- res[-1]    
    colnames(res) <- c("NR.GENES", "NR.SIG.GENES", "NR.SIG.COMB.GENES", GSP.COL)
    res <- as.matrix(res)
    return(res)
}

# pathnet helper: extract pathway data from gs and grn
extr.pwy.dat <- function(gs, grn)
{
    pwy.dat <- sapply(names(gs), 
        function(n)
        {
            genes <- gs[[n]] 
            sgrn <- query.grn(gs=genes, grn=grn, index=FALSE)
            if(nrow(sgrn))
                dat <- cbind(sgrn[,1:2, drop=FALSE], rep(n, nrow(sgrn)))
            else dat <- NULL
        }
    )
    pwy.dat <- pwy.dat[!sapply(pwy.dat, is.null)]
    pwy.datm <- matrix("", nrow=sum(sapply(pwy.dat, nrow)), ncol=3)
    colnames(pwy.datm) <- c("id1", "id2", "title")
    start <- 1
    for(i in seq_len(length(pwy.dat)))
    {
        end <- start + nrow(pwy.dat[[i]]) - 1
        pwy.datm[start:end,] <- pwy.dat[[i]]
        start <- end + 1
    }
    pwy.dat <- data.frame(id1=as.integer(pwy.datm[,1]), 
        id2=as.integer(pwy.datm[,2]), title=pwy.datm[,3])
    return(pwy.dat)
}

# pathnet helper: converts 3-col grn to adjacency matrix
grn2adjm <- function(grn)
{
    nodes <- sort(unique(as.vector(grn[,1:2])))
    adjm <- sapply(nodes, 
        function(n)
        {
            tgs <- grep(n, grn[,1])
            if(length(tgs))
            {
                tgs <- grn[tgs,2]
                adjv <- as.integer(nodes %in% tgs)
            }
            else adjv <- rep(0, length(nodes))
            return(adjv) 
        }) 
    rownames(adjm) <- nodes
    return(adjm)
}



















