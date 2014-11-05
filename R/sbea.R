############################################################
#
# author: Ludwig Geistlinger
# date: 06 Dec 2010
#
# Set-based Enrichment Analysis (SEA)
#
# Update, 06 May 2014: extension of file-based SEA 
#           to invokement within R
#
############################################################

# A INPUT FASSADE - wrapping & delegation
sbea.methods <- function() SBEA.METHODS

sbea <- function(   
    method=c("ora", "safe", "gsea", "samgs"), 
    eset, 
    gs, 
    alpha=0.05, 
    perm=1000, 
    out.file=NULL,
    browse=FALSE)
{   
    if(class(eset) == "character") eset <- get(load(eset))
    if(class(method) == "character")
    { 
        method <- method[1]
        if(length(find(method))) gs.ps <- do.call(method, 
            list(eset=eset, gs=gs, alpha=alpha, perm=perm))
        else if(!(method %in% sbea.methods())) 
            stop(paste("\'method\' must be one out of {", 
                paste(sbea.methods(), collapse=", "), "}"))
        # gsea
        else if(method == "gsea") 
            gs.ps <- gsea(eset, gs, perm=perm, alpha=alpha) 
        else
        {
            cmat <- gmt.2.cmat(gs, featureNames(eset), 
                        min.size=GS.MIN.SIZE, max.size=GS.MAX.SIZE)
            if(nrow(cmat) < nrow(eset)) eset <- eset[rownames(cmat),] 
            
            # ora
            if(method == "ora") 
                gs.ps <- ora(mode=1, eset=eset, 
                    cmat=cmat, perm=perm, alpha=alpha)
            # safe
            else if(method == "safe") 
                gs.ps <- ora(eset=eset, cmat=cmat, perm=perm, alpha=alpha) 
            # samgs
            else if(method == "samgs")
            {
                if(is.null(out.file)) 
                    out.dir <- file.path(system.file(
                        package="EnrichmentBrowser"), "results")
                else out.dir <- sub("\\.[a-z]+$","_files", out.file)
                if(!file.exists(out.dir)) dir.create(out.dir)
                samt.file <- file.path(out.dir, "samt.RData")
                gs.ps <- SAMGS(GS=as.data.frame(cmat), DATA=exprs(eset), 
                        cl=as.factor(as.integer(eset[[GRP.COL]]) + 1), 
                        nbPermutations=perm, 
                        tstat.file=samt.file)
            }
        }
    }
    else if(class(method) == "function") 
        gs.ps <- method(eset=eset, gs=gs, alpha=alpha, perm=perm)
    else stop(paste(method, "is not a valid method for sbea"))

    gs.ps <- sort(gs.ps)
    gs.ps <- signif(gs.ps, digits=3)
    res.tbl <- cbind(names(gs.ps), gs.ps)
    colnames(res.tbl) <- c("GENE.SET", "P.VALUE")
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
            nr.sigs=sum(gs.ps < alpha),
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
        if(nr.sigs) ranking <- res$res.tbl[seq_len(nr.sigs),,drop=FALSE]
        else return(NULL)
    }
    else ranking <- res$res.tbl
    return(as.data.frame(ranking))
}

# 0 HELPER
gmt.2.cmat <- function(gs, features, min.size=0, max.size=Inf)
{
    if(class(gs) == "character") gs <- parse.genesets.from.GMT(gs)
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
    names(gs.ps) <- colnames(cmat)

    return(gs.ps)
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
    x <- exprs(eset)
    y <- as.integer(pData(eset)[, GRP.COL]) 

    # execute hypergeom ORA if no permutations
    ps <- fData(eset)[, ADJP.COL]
    if(perm == 0) gs.ps <- ora.hyperg(ps=ps, cmat=cmat, alpha=alpha)
    # else do resampling using functionality of SAFE
    else{
        # resampl ORA
        if(mode == 1){
            nr.sigs <- sum(fData(eset)[ , ADJP.COL] < alpha)
            args <- list(one.sided=FALSE, genelist.length=nr.sigs)

            gs.ps <- safe(  X.mat=x, y.vec=y, 
                    C.mat=cmat, Pi.mat=perm, alpha=alpha, 
                    global="Fisher", args.global=args)
        } 
        # SAFE default                  
        else gs.ps <- safe(X.mat=x, y.vec=y, 
            C.mat=cmat, Pi.mat=perm, alpha=alpha)
        gs.ps <- gs.ps@global.pval
    }
    return(gs.ps)
}

# 4 GSEA
gsea <- function(
    eset, 
    gs.gmt, 
    doc.string="GSEA.Analysis", 
    alpha=0.05, 
    perm=1000, 
    out.file=NULL)
{        
    # build class list
    cls <- list()
    cls$phen <- levels(as.factor(eset[[GRP.COL]]))
    cls$class.v <- ifelse(eset[[GRP.COL]] == cls$phen[1], 0, 1)

    if(is.null(out.file)) 
        out.dir <- file.path(system.file(package="EnrichmentBrowser"), "results")
    else out.dir <- sub("\\.[a-z]+$", "_files", out.file)
    if(!file.exists(out.dir)) dir.create(out.dir)

    if(class(gs.gmt) != "character") 
    {
        gmt.out <- file.path(out.dir, paste(doc.string, "gs.gmt", sep="_"))
        write.gmt(gs.gmt, gmt.file=gmt.out)
        gs.gmt <- gmt.out
    }
    
    # run GSEA
    GSEA(input.ds = as.data.frame(exprs(eset)), 
        input.cls = cls,
        gs.db = gs.gmt, 
        output.directory = file.path(out.dir, ""),
        doc.string            = doc.string,
        nperm                 = perm,
        fdr.q.val.threshold   = -1,
        topgs                 = NROW.TOP.TABLE,
        gs.size.threshold.min = GS.MIN.SIZE,
        gs.size.threshold.max = GS.MAX.SIZE)
      
    res.file <- file.path(out.dir, 
        paste(doc.string, ".SUMMARY.RESULTS.REPORT.1.txt", sep=""))
    res <- as.matrix(read.delim(res.file))
    gs.ps <- as.numeric(res[,6])
    names(gs.ps) <- res[,1]

    tmp.files <- list.files(out.dir, 
        pattern=paste("^", doc.string, sep=""), full.names=TRUE)
    invisible(file.remove(tmp.files))
    
    return(gs.ps)
}

