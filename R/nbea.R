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

nbea.methods <- function() 
    c("ggea", "spia", "pathnet", "degraph", 
		"topologygsa", "ganpa", "cepa", "netgsa", "nea")

nbea <- function(
    method=EnrichmentBrowser::nbea.methods(), 
    eset, 
    gs, 
    grn,
    prune.grn=TRUE,
    alpha=0.05, 
    perm=1000, 
    padj.method="none",
    out.file=NULL,
    browse=FALSE, ...)
{
    # get configuration
    GS.MIN.SIZE <- config.ebrowser("GS.MIN.SIZE")
    GS.MAX.SIZE <- config.ebrowser("GS.MAX.SIZE")
    GSP.COL <- config.ebrowser("GSP.COL")
    FC.COL <-  config.ebrowser("FC.COL")
    ADJP.COL <-  config.ebrowser("ADJP.COL")
    
    # dealing with NA's
    eset <- eset[!is.na(fData(eset)[,FC.COL]), ]
    eset <- eset[!is.na(fData(eset)[,ADJP.COL]), ]

    # getting gene sets & grn
    if(!is.list(gs)) gs <- parse.genesets.from.GMT(gs)
    if(!is.matrix(grn)) grn <- read.grn(grn)

    # prune grn
    if(prune.grn) grn <- .pruneGRN(grn)

    # restrict to relevant genes 
    # in the intersection of eset, gs, and grn
    gs.genes <- unique(unlist(gs))
    grn.genes <- unique(c(grn[,1], grn[,2]))
    eset.genes <- featureNames(eset)
    rel.genes <- intersect(intersect(gs.genes, grn.genes), eset.genes)
    eset <- eset[rel.genes,]
    gs <- sapply(gs, function(s) s[s%in% rel.genes])
    lens <- sapply(gs, length)
    gs <- gs[lens >= GS.MIN.SIZE & lens <= GS.MAX.SIZE]
    grn <- grn[grn[,1] %in% rel.genes & grn[,2] %in% rel.genes,] 
    
    # execute ea
    if(class(method) == "character")
    {
        method <- match.arg(method)
        if(method == "spia") res.tbl <- .spia(eset, gs, grn, alpha, perm, ...)
        else if(method == "nea") res.tbl <- .nea(eset, gs, grn, alpha, perm, ...)
        else if(method == "pathnet") res.tbl <- .pathnet(eset, gs, grn, alpha, perm)
        else if(method == "netgsa") res.tbl <- .netgsa(eset, gs, grn)
        else if(method == "ganpa") res.tbl <- .ganpa(eset, gs, grn, perm)
        else if(method == "cepa") res.tbl <- .cepa(eset, gs, grn)
        else if(method == "degraph") res.tbl <- .degraph(eset, gs, grn)
        else if(method == "topologygsa") res.tbl <- .topogsa(eset, gs, grn, alpha, perm, ...)
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

#
# general helpers
#
.pruneGRN <- function(grn)
{
    # remove duplicates
    grn <- unique(grn)

    # rm self edges
    grn <- grn[grn[,1] != grn[,2],]
   
    # rm rev edges 
    genes <- unique(as.vector(grn[,1:2]))
    ggrid <- seq_along(genes)
    names(ggrid) <- genes
    igrn <- cbind(ggrid[grn[,1]], ggrid[grn[,2]])

    n <- nrow(grn)
    grid <- seq_len(n-1)
    ind <- sapply(grid,
        function(i)
        {
            x <- igrn[i,2:1]
            j <- i + 1
            cigrn <- igrn[j:n,,drop=FALSE]
            cigrn <- cigrn[cigrn[,1] == x[1], , drop=FALSE]
            is.rev <- any( cigrn[,2] == x[2] )
            return(is.rev)
        })
    ind <- c(ind, FALSE)
    grn <- grn[!ind,]
}
    
#
# 1 SPIA
#
.spia <- function(eset, gs, grn, 
    alpha=0.05, perm=1000, beta=1, sig.stat=c("xxP", "xxFC", "|", "&")) 
{
    FC.COL <- config.ebrowser("FC.COL")
    ADJP.COL <- config.ebrowser("ADJP.COL")
    GSP.COL <- config.ebrowser("GSP.COL")

    de.genes <- is.sig(fData(eset), alpha, beta, sig.stat)
    de <- fData(eset)[de.genes, FC.COL]
    names(de) <- featureNames(eset)[de.genes]
    all <- featureNames(eset)

    is.kegg <- auto.detect.gs.type(names(gs)[1]) == "KEGG"
    organism <- ""
    data.dir <- NULL
    if(is.kegg) organism <- substring(names(gs)[1],1,3)
    else
    {     
        message("making SPIA data ...")
        path.info <- .make.spia.data(gs, grn)
        data.dir <- system.file("extdata/", package="SPIA")
        save(path.info, file=file.path(data.dir, "SPIA.RData"))
    }
    res <- SPIA::spia(de=de, all=all, organism=organism, data.dir=data.dir, nB=perm)
    res[,"Name"] <- gsub(" ", "_", res[,"Name"])
    rownames(res) <- paste(paste0(organism, res[,"ID"]), res[,"Name"], sep="_")
    res <- res[, c("pSize", "NDE", "tA", "Status", "pG")]
    colnames(res) <- c("SIZE", "NDE", "T.ACT", "STATUS", GSP.COL)
    res[,"STATUS"] <- ifelse(res[,"STATUS"] == "Activated", 1, -1)
    res <- as.matrix(res)
    message("Finished SPIA analysis")
    return(res)
}

# spia helper: create pathway data from gs and grn
.make.spia.data <- function(gs, grn)
{
    rel <- c("activation", "compound", "binding/association", 
            "expression", "inhibition", "activation_phosphorylation", 
            "phosphorylation", "inhibition_phosphorylation", 
            "inhibition_dephosphorylation", "dissociation", "dephosphorylation", 
            "activation_dephosphorylation", "state change", "activation_indirect effect", 
            "inhibition_ubiquination", "ubiquination", "expression_indirect effect", 
            "inhibition_indirect effect", "repression", "dissociation_phosphorylation", 
            "indirect effect_phosphorylation", "activation_binding/association", 
            "indirect effect", "activation_compound", "activation_ubiquination")

    spia.data <- sapply(names(gs), 
        function(s)
        {   
            x <- gs[[s]]
            len <- length(x) 
            nam <- list(x, x)
            m <- matrix(0, nrow=len, ncol=len, dimnames=nam)
            sgrn <- query.grn(gs=x, grn=grn, index=FALSE)
            if(nrow(sgrn) < config.ebrowser("GS.MIN.SIZE")) return(NULL)
            act.grn <- sgrn[sgrn[,3] == "+",,drop=FALSE]
            actm2 <- m
            if(nrow(act.grn))
            {
                if(nrow(act.grn) > 1) actm <- .grn2adjm(act.grn)
                else actm <- matrix(1, nrow=1, ncol=1, dimnames=list(act.grn[1,1], act.grn[1,2]))
                actm2[rownames(actm), colnames(actm)] <- actm
            }
            inh.grn <- sgrn[sgrn[,3] == "-",,drop=FALSE]
            inhm2 <- m
            if(nrow(inh.grn))
            {
                if(nrow(inh.grn) > 1) inhm <- .grn2adjm(inh.grn)
                else inhm <-  matrix(1, nrow=1, ncol=1, dimnames=list(inh.grn[1,1], inh.grn[1,2]))
                inhm2[rownames(inhm), colnames(inhm)] <- inhm
            }
            l <- lapply(rel, 
                function(r)
                {
                    if(r == "activation") return(actm2)
                    else if(r == "inhibition") return(inhm2)
                    else return(m)                    
                } 
            )
            names(l) <- rel    
            l$nodes <- x
            l$title <- s
            l$NumberOfReactions <- 0
            return(l)
        }
    )
    spia.data <- spia.data[!sapply(spia.data, is.null)]
    return(spia.data)
}

#
# 2 NEA
#
.nea <- function(eset, gs, grn, 
    alpha=0.05, perm=100, beta=1, sig.stat=c("xxP", "xxFC", "|", "&"))
{
    nea <- NULL
    .isAvailable("neaGUI", type="software")

    #if(perm > 100) perm <- 100
    isig <- is.sig(fData(eset), alpha, beta, sig.stat)
    ags <- featureNames(eset)[isig]
    grn <- unique(grn[,1:2])
    gs.genes <- unique(unlist(gs))
    grn <- grn[(grn[,1] %in% gs.genes) & (grn[,2] %in% gs.genes),]
    network <- apply(grn, 1, function(x) paste(x, collapse=" "))
    message("Computing NEA permutations, this may take a few minutes ...")
    res <- nea(ags=ags, fgs=gs, network=network, nperm=perm)
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
    GSP.COL <- config.ebrowser("GSP.COL")
    res <- res[order(res[,GSP.COL]),]
    return(res) 
}

#
# 3 Pathnet
#
.pathnet <- function(eset, gs, grn, alpha=0.05, perm=1000)
{
    PathNet <- NULL
    .isAvailable("PathNet", type="software")
    
    ADJP.COL <- config.ebrowser("ADJP.COL")
    GSP.COL <- config.ebrowser("GSP.COL")

    dir.evid <- -log(fData(eset)[,ADJP.COL], base=10)
    dir.evid <- cbind(as.integer(featureNames(eset)), dir.evid)
    colnames(dir.evid) <- c("Gene.ID", "Obs")
    adjm <- .grn2adjm(grn, directed=FALSE)
    pwy.dat <- .extr.pwy.dat(gs, grn)
    
    res <- PathNet(
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
    rownames(res) <- sapply(as.vector(res[,1]), 
        function(s) grep(unlist(strsplit(s,"_"))[1], names(gs), value=TRUE))
    res <- res[-1]    
    colnames(res) <- c("NR.GENES", "NR.SIG.GENES", "NR.SIG.COMB.GENES", GSP.COL)
    res <- as.matrix(res)
    return(res)
}

# pathnet helper: extract pathway data from gs and grn
.extr.pwy.dat <- function(gs, grn)
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
.grn2adjm <- function(grn, directed=TRUE)
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
    adjm <- t(adjm)
  
    if(!directed)
        for(i in seq_along(nodes)) 
            for(j in seq_along(nodes)) 
                if(adjm[i,j]) adjm[j,i] <- 1

    return(adjm)
}

#
# 4 NetGSA
#
.netgsa <- function(eset, gs, grn)
{
     NetGSA <- covsel <- edgelist2adj <- NULL
    .isAvailable("netgsa", type="software")

    x <- exprs(eset)
    y <- pData(eset)[,config.ebrowser("GRP.COL")] + 1

    # prepare gene sets
    #cmat <- .gmt2cmat(gs, rownames(x)) 
    #if(nrow(cmat) < nrow(x)) x <- x[rownames(cmat),]
    f <- file()
    sink(file=f)
    cmat <- safe::getCmatrix(gs, as.matrix=TRUE)
    sink()
    close(f)
    x <- x[rownames(cmat),]

    # prepare network
    out.dir <- config.ebrowser("OUTDIR.DEFAULT")
    if(!file.exists(out.dir)) dir.create(out.dir)
    out.file <- file.path(out.dir, "grn.txt")
    write.table(grn[,1:2], file=out.file, row.names=FALSE)
    adjm <- edgelist2adj(out.file, vertex.names=unique(as.vector(grn[,1:2])))
    file.remove(out.file)
    ind <- intersect(rownames(x), rownames(adjm))
    adjm <- adjm[ind, ind]
    x <- x[rownames(adjm),]
    cmat <- cmat[,rownames(adjm)]
    cmat <- cmat[rowSums(cmat) > 2,]
    
    message("Estimating weighted adjacency matrix for GRN (group 0)")
    A1 <- covsel(t(x[,y==1]), one=adjm, lambda=0.2)
    message("Estimating weighted adjacency matrix for GRN (group 1)")
    A2 <- covsel(t(x[,y==2]), one=adjm, lambda=0.2)

    # execute
    message("Executing NetGSA ...")
    message("This may take a while ...")
    res <- NetGSA(A1$wAdj, A2$wAdj, x, y, B=cmat, directed=TRUE)

    res <- cbind(res$teststat, res$p.value)
    colnames(res) <- c("STAT", config.ebrowser("GSP.COL"))
    rownames(res) <- rownames(cmat)
    return(res)
}

#
# 5 GANPA
#
.ganpa <- function(eset, gs, grn, perm=1000)
{
    GSE.Test.Main <- NULL
    .isAvailable("GANPA", type="software")

    # configure
    GRP.COL <- config.ebrowser("GRP.COL")
    SMPL.COL <- config.ebrowser("SMPL.COL")
    OUT.DIR <- config.ebrowser("OUTDIR.DEFAULT")
    GS.MIN.SIZE <- config.ebrowser("GS.MIN.SIZE")
    GS.MAX.SIZE <- config.ebrowser("GS.MAX.SIZE")
    GSP.COL <- config.ebrowser("GSP.COL")
    
    if(!file.exists(OUT.DIR)) dir.create(OUT.DIR)
    out.prefix <- file.path(OUT.DIR, "ganpa")

    # expression data
    has.scol <- SMPL.COL %in% colnames(pData(eset))
    if(!has.scol) pData(eset)[,SMPL.COL] <- sampleNames(eset)
    sinfo <- pData(eset)[,c(SMPL.COL, GRP.COL)]
    colnames(sinfo) <- c("sampleid", "status")
    expr.obj <- list(gExprs=exprs(eset), sampleinfo=sinfo)
    
    # gene regulatory network
    gnet <- .grn2gnet(grn)    

    # execute
    GSE.Test.Main(gExprs.obj=expr.obj, gsets=gs, gNET=gnet, 
        permN=perm, size.min=GS.MIN.SIZE, size.max=GS.MAX.SIZE,
        msp.correction=FALSE, output.label=out.prefix, permFDR.cutoff=1)

    # read results from output csv
    res <- read.csv(paste(out.prefix, "MeanAbs.OrigW.csv", sep=".")) 
    n <- res[,1]
    res <- as.matrix(res[,c("Size", "S", "NS", "permP")])
    colnames(res)[c(1,4)] <- c("SIZE", GSP.COL)
    rownames(res) <- n
    return(res)
}

.grn2gnet <- function(grn)
{
    ureg <- unique(grn[,1])
    gnet <- sapply(ureg, function(r) grn[grn[,1] == r,2])
    return(gnet)
}

#
# 6 CePa
#
.cepa <- function(eset, gs, grn, perm=1000)
{
    cepa.all <- set.pathway.catalogue <- sampleLabel <- NULL
    .isAvailable("CePa", type="software")

    # define sample groups
    GRP.COL <- config.ebrowser("GRP.COL")
    sl <- sampleLabel(pData(eset)[, GRP.COL], treatment=1, control=0)

    # create pathway catalogue from gs and grn
    # (1) pathway list
    pl <- sapply(gs, function(s) as.character(query.grn(s, grn)))
    pl <- pl[sapply(pl, length) >= config.ebrowser("GS.MIN.SIZE")]

    # (2) interaction list
    il <- data.frame(as.character(seq_len(nrow(grn))), grn[,1:2], stringsAsFactors=FALSE)
    colnames(il) <- c("interaction.id", "input", "output")

    # (3) mapping
    m <- data.frame(node.id=featureNames(eset), symbol=featureNames(eset), stringsAsFactors=FALSE)
    pwy.cat <- set.pathway.catalogue(pathList=pl, interactionList=il, mapping=m)
    
    # executing
    res <- cepa.all(mat=exprs(eset), label=sl, pc=pwy.cat, iter=perm)
    res.mat <- matrix(0.0, nrow=length(res), ncol=7)
    for(i in seq_along(res))
    {
        res.mat[i,1:6] <- sapply(res[[i]], function(x) x$p.value)
        res.mat[i,7] <- min(6 * min(res.mat[i,1:6]), 1) 

    }
    rownames(res.mat) <- names(res)
    n <- paste(toupper(names(res[[1]])), "PVAL" , sep=".")
    colnames(res.mat) <- c(n, config.ebrowser("GSP.COL"))
    return(res.mat)
}

#
# 7 DEGraph
#
.degraph <- function(eset, gs, grn)
{    
    testOneGraph <- NULL
    .isAvailable("DEGraph", type="software")

    grp <- pData(eset)[,config.ebrowser("GRP.COL")]
            
    options(show.error.messages=FALSE) 
    res <- sapply(names(gs),
        function(s)
        {
            gs.grn <- query.grn(gs[[s]], grn, index=FALSE)
            if(nrow(gs.grn) < config.ebrowser("GS.MIN.SIZE")) return(NA)
            am <- .grn2adjm(gs.grn)
            gr <- graphAM(am, "directed")           
            gr <- as(gr, "graphNEL")                 
            r <- try( testOneGraph(graph=gr, data=exprs(eset), 
                       classes=grp, useInteractionSigns=FALSE), silent=TRUE )
            if(is(r, "try-error")) return(NA) else return(r[[1]]$p.value[1])
        })
    options(show.error.messages=TRUE)
    names(res) <- names(gs)
    res <- res[!is.na(res)]
    return(res)
}

#
# 8 topologyGSA 
#
.topogsa <- function(eset, gs, grn, alpha=0.05, perm=1000, test.mean=TRUE)
{    
    # call topologyGSA via clipper's pathQ function
    return(.clipper(eset, gs, grn, alpha, perm))

    # original topologyGSA: deprecated
    # does not terminate on particular gs, eg. hsa04060_Cytokine-cytokine_receptor_interaction
    pathway.mean.test <- pathway.var.test <- NULL
    .isAvailable("topologyGSA", type="software")
  
    is.DAG <- NULL
    .isAvailable("gRbase", type="software")
    
    graph_from_graphnel <- mst <- NULL
     .isAvailable("igraph", type="software")
 
    grp <- pData(eset)[,config.ebrowser("GRP.COL")]
    y1 <- t(exprs(eset)[, grp == 0])
    y2 <- t(exprs(eset)[, grp == 1])

    res <- sapply(names(gs),
        function(s)
        {
            message(s)
            gs.grn <- query.grn(gs[[s]], grn, index=FALSE)
            if(nrow(gs.grn) < config.ebrowser("GS.MIN.SIZE")) return(NA)
            am <- .grn2adjm(gs.grn)
            gr <- graphAM(am, "directed")           
            gr <- as(gr, "graphNEL")                 
            
            if(!is.DAG(gr))
            {
                gr2 <- graph_from_graphnel(gr)
                gr2 <- mst(gr2)
                gr <- as(gr2, "graphNEL")
            }
            if(test.mean) r <- pathway.mean.test(y1, y2, gr, alpha, perm)
            else r <- pathway.var.test(y1, y2, gr, alpha)
            return(r$p.value)
        }
    )
    names(res) <- names(gs)
    res <- res[!is.na(res)]
    return(res)
}

#
# 9 clipper
#
.clipper <- function(eset, gs, grn, alpha=0.05, perm=1000)
{    
    pathQ <- NULL
    .isAvailable("clipper", type="software")
  
    is.DAG <- NULL
    .isAvailable("gRbase", type="software")
    
    graph_from_graphnel <- mst <- NULL
     .isAvailable("igraph", type="software")
 
    grp <- pData(eset)[,config.ebrowser("GRP.COL")] + 1

    a <- Sys.time()
    res <- sapply(names(gs),
        function(s)
        {
            message(s)
            gs.grn <- query.grn(gs[[s]], grn, index=FALSE)
            if(nrow(gs.grn) < config.ebrowser("GS.MIN.SIZE")) return(NA)
            am <- .grn2adjm(gs.grn)
            gr <- graphAM(am, "directed")           
            gr <- as(gr, "graphNEL")                 
            
            if(!is.DAG(gr))
            {
                gr2 <- graph_from_graphnel(gr)
                gr2 <- mst(gr2)
                gr <- as(gr2, "graphNEL")
            }
            r <- try(pathQ(exprs(eset), grp, gr, perm, alpha), silent=TRUE)
            if(is(r, "try-error")) return(NA) else return(r$alphaMean)
        }
    )
    names(res) <- names(gs)
    res <- res[!is.na(res)]
    b <- Sys.time()
    return(res)
}

