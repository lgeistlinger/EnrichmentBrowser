############################################################
#
# author: Ludwig Geistlinger
# date: 06 Dec 2010
#
# Set-based Enrichment Analysis (SBEA)
#
############################################################

# A INPUT FASSADE - wrapping & delegation
sbea.methods <- function() c("ora", "safe", "gsea", "samgs")

sbea <- function(   
    method=sbea.methods(), 
    eset, 
    gs, 
    alpha=0.05, 
    perm=1000, 
    padj.method="none",
    out.file=NULL,
    browse=FALSE)
{   
    GS.MIN.SIZE <- config.ebrowser("GS.MIN.SIZE")
    GS.MAX.SIZE <- config.ebrowser("GS.MAX.SIZE")
    GSP.COL <- config.ebrowser("GSP.COL")

    # restrict eset and gs to intersecting genes
    igenes <- intersect(featureNames(eset), unique(unlist(gs)))
    eset <- eset[igenes,]
    gs <- sapply(gs, function(s) s[s %in% igenes]) 
    lens <- sapply(gs, length)
    gs <- gs[lens >= GS.MIN.SIZE & lens <= GS.MAX.SIZE]

    if(is.character(method))
    { 
        method <- match.arg(method)
        data.type <- experimentData(eset)@other$dataType
        if(is.null(data.type)) data.type <- auto.detect.data.type(exprs(eset))
 
        #if(length(find(method))) gs.ps <- do.call(method, 
        #    list(eset=eset, gs=gs, alpha=alpha, perm=perm))
        #else 
        #if(!(method %in% sbea.methods())) 
        #    stop(paste("\'method\' must be one out of {", 
        #        paste(sbea.methods(), collapse=", "), "}"))
        # else 
        # gsea
        if(data.type != "rseq" & method == "gsea") gs.ps <- gsea(eset, gs, perm)
        else
        {
            cmat <- gmt.2.cmat(gs, featureNames(eset), GS.MIN.SIZE, GS.MAX.SIZE)
            if(nrow(cmat) < nrow(eset)) eset <- eset[rownames(cmat),] 
           
            # hypergeom. ora
            if(method == "ora" & perm==0) gs.ps <- ora(1, eset, cmat, perm, alpha)
            # rseq 
            else if(data.type == "rseq") gs.ps <- rseq.sbea(method, eset, cmat, perm, alpha)
            # resampl ora
            else if(method == "ora") gs.ps <- ora(1, eset, cmat, perm, alpha)
            #safe
            else if(method == "safe") gs.ps <- ora(2, eset, cmat, perm, alpha)
            # samgs
            else if(method == "samgs")
            {
                if(is.null(out.file)) out.dir <- config.ebrowser("OUTDIR.DEFAULT")
                else out.dir <- sub("\\.[a-z]+$","_files", out.file)
                if(!file.exists(out.dir)) dir.create(out.dir)
                samt.file <- file.path(out.dir, "samt.RData")
                GRP.COL <- config.ebrowser("GRP.COL")
                gs.ps <- SAMGS(GS=as.data.frame(cmat), DATA=exprs(eset), 
                        cl=as.factor(as.integer(eset[[GRP.COL]]) + 1), 
                        nbPermutations=perm, 
                        tstat.file=samt.file)
            }
        }
    }
    else if(is.function(method)) 
        gs.ps <- method(eset=eset, gs=gs, alpha=alpha, perm=perm)
    else stop(paste(method, "is not a valid method for sbea"))

    res.tbl <- data.frame(signif(gs.ps, digits=3))
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

gs.ranking <- function(res, signif.only=TRUE)
{
    if(signif.only)
    {
        nr.sigs <- res$nr.sigs
        if(nr.sigs) ranking <- res$res.tbl[seq_len(nr.sigs),]
        else return(NULL)
    }
    else ranking <- res$res.tbl
    return(ranking)
}

# 0 HELPER
gmt.2.cmat <- function(gs, features, min.size=0, max.size=Inf)
{
    if(is.character(gs)) gs <- parse.genesets.from.GMT(gs)
    # transform gs gmt to cmat
    cmat <- sapply(gs, function(x) features %in% x)
    rownames(cmat) <- features

    # restrict to gene sets with valid size
    gs.sizes <- colSums(cmat)
    valid.size <- which((gs.sizes >= min.size) & (gs.sizes <= max.size))
    if(length(valid.size) == 0) stop("No gene set with valid size!")
    cmat <- cmat[, valid.size]
    
    # restrict to genes which are in sets with valid size
    has.set <- which(rowSums(cmat) > 0)
    cmat <- cmat[has.set,]

    return(cmat)
}

###
#
# ENRICHMENT METHODS
#
###
# de.ana as local.stat for safe
local.de.ana <- function (X.mat, y.vec, args.local)
{
    return(function(data, ...) 
    {
        stat <- de.ana(expr=data, grp=y.vec,
            blk=args.local$blk,
            de.method=args.local$de.method, 
            stat.only=TRUE)
        return(stat)
    })
}


rseq.sbea <- function(method, eset, cmat, perm, alpha)
{
    assign("local.de.ana", local.de.ana, envir=.GlobalEnv)
    de.method <- grep(".STAT$", colnames(fData(eset)), value=TRUE)
    de.method <- sub(".STAT$",  "", de.method)
    
    blk <- NULL
    blk.col <- config.ebrowser("BLK.COL") 
    if(blk.col %in% colnames(pData(eset))) blk <- pData(eset)[,blk.col]

    args.local <- list(de.method=de.method, blk=blk)

    args.global <- list(one.sided=FALSE)
    if(method == "ora")
    {
        global <- "Fisher"
        nr.sigs <- sum(fData(eset)[, config.ebrowser("ADJP.COL")] < alpha)
        args.global <- list(one.sided=FALSE, genelist.length=nr.sigs)
    }
    else if(method == "safe") global <- "Wilcoxon" 
    else if(method == "gsea") global <- "Kolmogorov"
    else if(method == "samgs")
    {
        global <- "SAMGS"
        assign("global.SAMGS", global.SAMGS, envir=.GlobalEnv)
    }

    y <- pData(eset)[,config.ebrowser("GRP.COL")]
    gs.ps <- safe::safe(X.mat=exprs(eset), y.vec=y, C.mat=cmat,         
        local="de.ana", args.local=args.local,
        global=global, args.global=args.global, 
        Pi.mat=perm, alpha=alpha, error="none")
 
    res.tbl <- cbind(
            gs.ps@global.stat, 
            gs.ps@global.stat / colSums(cmat), 
            gs.ps@global.pval)
    
    colnames(res.tbl) <- c("GLOB.STAT", "NGLOB.STAT", config.ebrowser("GSP.COL"))
    
    return(res.tbl)
}

# 1 HYPERGEOM ORA
ora.hyperg <- function(ps, cmat, alpha=0.05)
{
    # determine sig. diff. exp. genes of eset, 
    # corresponds to sample size from urn
    is.sig <- ps < alpha
    nr.sigs <- sum(is.sig)
    
    # determine overlap of sig and set genes for each set
    sig.cmat <- cmat & is.sig

    # white balls observed when drawing nr.sigs balls from the urn
    ovlp.sizes <- colSums(sig.cmat)

    # white balls in the urn  (genes in gene set)
    gs.sizes <- colSums(cmat) 
    # black balls in the urn (genes not in gene set)
    uni.sizes <- length(ps) - gs.sizes 

    # determine significance of overlap 
    # based on hypergeom. distribution
    gs.ps <- 1 - phyper(ovlp.sizes, gs.sizes, uni.sizes, nr.sigs)
    
    res.tbl <- cbind(gs.sizes, ovlp.sizes, gs.ps)
    colnames(res.tbl) <- c("NR.GENES", "NR.SIG.GENES", config.ebrowser("GSP.COL"))
    rownames(res.tbl) <- colnames(cmat)

    return(res.tbl)
}

# 2 RESAMPL ORA
# 3 SAFE
#
# wrapper to call safe functionality approriately for
# overrepresentation analysis (ORA)
#
# perm=0 will execute traditional hypergeom. ORA
#
# for perm > 0 use
#   mode=1 ... resampl ORA (fisher)
#   mode=2 ... safe default (wilcoxon)
#
#
ora <- function(mode=2, eset, cmat, perm=1000, alpha=0.05)
{
    GRP.COL <- config.ebrowser("GRP.COL")
    ADJP.COL <- config.ebrowser("ADJP.COL")

    x <- exprs(eset)
    y <- pData(eset)[, GRP.COL]

    # execute hypergeom ORA if no permutations
    ps <- fData(eset)[, ADJP.COL]
    if(perm == 0) res.tbl <- ora.hyperg(ps=ps, cmat=cmat, alpha=alpha)
    # else do resampling using functionality of SAFE
    else{
        # resampl ORA
        if(mode == 1){
            nr.sigs <- sum(fData(eset)[ , ADJP.COL] < alpha)
            args <- list(one.sided=FALSE, genelist.length=nr.sigs)

            gs.ps <- safe::safe(X.mat=x, y.vec=y, global="Fisher", C.mat=cmat, 
                 Pi.mat=perm, alpha=alpha, error="none", args.global=args)
        } 
        # SAFE default                  
        else gs.ps <- safe::safe(X.mat=x, y.vec=y, 
            C.mat=cmat, Pi.mat=perm, alpha=alpha, error="none")
        res.tbl <- cbind(
            gs.ps@global.stat, 
            gs.ps@global.stat / colSums(cmat), 
            gs.ps@global.pval)
        colnames(res.tbl) <- c("GLOB.STAT", "NGLOB.STAT", config.ebrowser("GSP.COL"))
    }
    return(res.tbl)
}

# 4 GSEA
gsea <- function(
    eset, 
    gs.gmt, 
    perm=1000, 
    out.file=NULL)
{        
    GRP.COL <- config.ebrowser("GRP.COL")
    
    # npGSEA
    if(perm==0)
    {
        gsc <- gs.list.2.gs.coll(gs.gmt)
        res <- npGSEA::npGSEA(x=exprs(eset), y=eset[[GRP.COL]], set=gsc)
        ps <- sapply(res, npGSEA::pTwoSided)
        names(ps) <- names(gs.gmt)
        return(ps)
    }

    # build class list
    cls <- list()
    cls$phen <- levels(as.factor(eset[[GRP.COL]]))
    cls$class.v <- ifelse(eset[[GRP.COL]] == cls$phen[1], 0, 1)

    if(is.null(out.file)) 
        out.dir <- config.ebrowser("OUTDIR.DEFAULT") 
    else out.dir <- sub("\\.[a-z]+$", "_files", out.file)
    if(!file.exists(out.dir)) dir.create(out.dir)
        
    # run GSEA
    res <- GSEA(input.ds=as.data.frame(exprs(eset)), 
        input.cls=cls, gs.db=gs.gmt, output.directory=out.dir, nperm=perm)
      
    gs.ps <- S4Vectors::as.matrix(res[,3:5])
    rownames(gs.ps) <- res[,1]

    return(gs.ps)
}



