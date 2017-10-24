############################################################
#
# author: Ludwig Geistlinger
# date: 06 Dec 2010
#
# Set-based Enrichment Analysis (SBEA)
#
############################################################

sbea.methods <- function() 
    c("ora", "safe", "gsea", "samgs", "ebm", "mgsa", 
        "gsa", "padog", "globaltest", "roast", "camera", "gsva")

# INPUT FASSADE - wrapping & delegation
sbea <- function(   
    method=EnrichmentBrowser::sbea.methods(), 
    eset, 
    gs, 
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

    ### TEMPORARY: will be replaced by as(eSet,SummarizedExperiment)
    if(is(eset, "ExpressionSet")) eset <- as(eset, "RangedSummarizedExperiment")
    ###    

    # dealing with NA's
    nr.na <- sum(is.na(rowData(eset)[,FC.COL]))
    if(nr.na) eset <- eset[!is.na(rowData(eset)[,FC.COL]),]
    nr.na <- sum(is.na(rowData(eset)[,ADJP.COL]))
    if(nr.na) eset <- eset[!is.na(rowData(eset)[,ADJP.COL]),]    

    # getting gene sets
    if(!is.list(gs)) gs <- parse.genesets.from.GMT(gs)

    # restrict eset and gs to intersecting genes
    igenes <- intersect(rownames(eset), unique(unlist(gs)))
    eset <- eset[igenes,]
    gs <- lapply(gs, function(s) s[s %in% igenes]) 
    lens <- lengths(gs)
    gs <- gs[lens >= GS.MIN.SIZE & lens <= GS.MAX.SIZE]

    if(is.character(method))
    { 
        method <- match.arg(method)
        data.type <- metadata(eset)$dataType
        if(is.null(data.type)) data.type <- .detectDataType(assay(eset))

        # rseq? 
        if(data.type == "rseq")
        {
            # mgsa
            if(method == "mgsa") gs.ps <- .mgsa(eset, gs, alpha, ...)
            # globaltest
            else if(method == "globaltest") gs.ps <- .globaltest(eset, gs, perm)
            # roast & camera
            else if(method %in% c("roast", "camera"))
                gs.ps <- .roast.camera(method, eset, gs, perm, rseq=TRUE)
		    # gsva
		    else if(method == "gsva") gs.ps <- .gsva(eset, gs, rseq=TRUE)
            else
            {
                # gs2cmat
                #cmat <- .gmt2cmat(gs, rownames(eset), GS.MIN.SIZE, GS.MAX.SIZE)
                #if(nrow(cmat) < nrow(eset)) eset <- eset[rownames(cmat),] 
                f <- file()
                sink(file=f)
                cmat <- safe::getCmatrix(gs, as.matrix=TRUE)
                sink()
                close(f)
                eset <- eset[rownames(cmat),]

                # ora
                if(method == "ora" & perm==0) 
                    gs.ps <- .ora(1, eset, cmat, perm, alpha, ...)
                # ebm
                else if(method == "ebm") gs.ps <- .ebm(eset, cmat)
		        # all others
                else gs.ps <- .rseqSBEA(method, eset, cmat, perm, alpha)
            }
        }
        # microarray
        else
        { 
            # gsea
            if(method == "gsea") gs.ps <- .gsea(eset, gs, perm)
            # gsa
            else if(method == "gsa") gs.ps <- .gsa(eset, gs, perm)
            # padog
	        else if(method == "padog") gs.ps <- .padog(eset, gs, perm)
	        # mgsa
	        else if(method == "mgsa") gs.ps <- .mgsa(eset, gs, alpha, ...)
            # globaltest
            else if(method == "globaltest") gs.ps <- .globaltest(eset, gs, perm)
            # roast & camera
            else if(method %in% c("roast", "camera"))
                gs.ps <- .roast.camera(method, eset, gs, perm)
	        # gsva
	        else if(method == "gsva") gs.ps <- .gsva(eset, gs)
            else
            {
                # gs2cmat
                #cmat <- .gmt2cmat(gs, rownames(eset), GS.MIN.SIZE, GS.MAX.SIZE)
                #if(nrow(cmat) < nrow(eset)) eset <- eset[rownames(cmat),] 
                f <- file()
                sink(file=f)
                cmat <- safe::getCmatrix(gs, as.matrix=TRUE)
                sink()
                close(f)
                eset <- eset[rownames(cmat),]

                # ora
                if(method == "ora") gs.ps <- .ora(1, eset, cmat, perm, alpha, ...)
                #safe
                else if(method == "safe") gs.ps <- .ora(2, eset, cmat, perm, alpha)
                # ebm
                else if(method == "ebm") gs.ps <- .ebm(eset, cmat)
                # samgs
                else if(method == "samgs")
                {
                    if(is.null(out.file)) out.dir <- config.ebrowser("OUTDIR.DEFAULT")
                    else out.dir <- sub("\\.[a-z]+$","_files", out.file)
                    if(!file.exists(out.dir)) dir.create(out.dir)
                    samt.file <- file.path(out.dir, "samt.RData")
                    GRP.COL <- config.ebrowser("GRP.COL")
                    gs.ps <- SAMGS(GS=as.data.frame(cmat), DATA=assay(eset), 
                        cl=as.factor(as.integer(eset[[GRP.COL]]) + 1), 
                        nbPermutations=perm, 
                        tstat.file=samt.file)
                }
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

###
#
# HELPER
# 
###
.gmt2cmat <- function(gs, features, min.size=0, max.size=Inf)
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


###
#
# ENRICHMENT METHODS
#
###

.rseqSBEA <- function(method, eset, cmat, perm, alpha)
{
	assign("eset", eset, envir=.GlobalEnv)
    assign("local.de.ana", local.de.ana, envir=.GlobalEnv)
    de.method <- grep(".STAT$", colnames(rowData(eset)), value=TRUE)
    de.method <- sub(".STAT$",  "", de.method)
    
    blk <- NULL
    blk.col <- config.ebrowser("BLK.COL") 
    if(blk.col %in% colnames(colData(eset))) blk <- colData(eset)[,blk.col]

    args.local <- list(de.method=de.method, blk=blk)

    args.global <- list(one.sided=FALSE)
    if(method == "ora")
    {
        global <- "Fisher"
        nr.sigs <- sum(rowData(eset)[, config.ebrowser("ADJP.COL")] < alpha)
        args.global$genelist.length <- nr.sigs
    }
    else if(method == "safe") global <- "Wilcoxon" 
    else if(method == "gsea") global <- "Kolmogorov"
    else if(method %in% c("samgs", "gsa", "padog"))
    {
        global <- toupper(method)
        global.func <- paste("global", global, sep=".")
        assign(global.func, get(global.func), envir=.GlobalEnv)
        
        if(method == "padog") args.global$gf <- .getGeneFreqWeights(cmat)
    }

	x <- assay(eset)
    y <- colData(eset)[,config.ebrowser("GRP.COL")]
    gs.ps <- safe::safe(X.mat=x, y.vec=y, C.mat=cmat,         
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

.isSig <- function(fdat, alpha=0.05, beta=1, sig.stat=c("xxP", "xxFC", "|", "&"))
{
    FC.COL <- config.ebrowser("FC.COL")
    ADJP.COL <- config.ebrowser("ADJP.COL")

    sig.stat <- sig.stat[1]
    if(grepl("P$", sig.stat))
    {
        if(sig.stat == "xxP") sig <- fdat[, ADJP.COL] < alpha
        else
        {
            perc <- as.integer(substring(sig.stat, 1, 2))
            p <- fdat[,ADJP.COL]
            names(p) <- rownames(fdat) 
            ordp <- sort(p)
            nr.sig <- round( length(p) * (perc / 100) )
            sigs <- names(ordp)[seq_len(nr.sig)]
            sig <- rownames(fdat) %in% sigs
        }
    }
    else if(grepl("FC$", sig.stat))
    { 
        if(sig.stat == "xxFC") sig <- abs(fdat[, FC.COL]) > beta
        else
        {
            perc <- as.integer(substring(sig.stat, 1, 2))
            fc <- fdat[,FC.COL]
            names(fc) <- rownames(fdat) 
            ordfc <- fc[order(abs(fc), decreasing=TRUE)]
            nr.sig <- round( length(fc) * (perc / 100) )
            sigs <- names(ordfc)[seq_len(nr.sig)]
            sig <- rownames(fdat) %in% sigs
        }
    }
    else 
    {
        psig <- fdat[, ADJP.COL] < alpha
        fcsig <- abs(fdat[, FC.COL]) > beta
        sig <- do.call(sig.stat, list(psig, fcsig))
    }
    return(sig)
}

# 1 HYPERGEOM ORA
.oraHypergeom <- function(fdat, cmat, 
    alpha=0.05, beta=1, sig.stat=c("xxP", "xxFC", "|", "&"))
{
    # determine sig. diff. exp. genes of eset, 
    # corresponds to sample size from urn
    isig <- .isSig(fdat, alpha, beta, sig.stat)
    nr.sigs <- sum(isig)
    
    # determine overlap of sig and set genes for each set
    sig.cmat <- cmat & isig

    # white balls observed when drawing nr.sigs balls from the urn
    ovlp.sizes <- colSums(sig.cmat)

    # white balls in the urn  (genes in gene set)
    gs.sizes <- colSums(cmat) 
    # black balls in the urn (genes not in gene set)
    uni.sizes <- nrow(fdat) - gs.sizes 

    # determine significance of overlap 
    # based on hypergeom. distribution
    gs.ps <- phyper(ovlp.sizes-1, gs.sizes, uni.sizes, nr.sigs, lower.tail=FALSE) 
    
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
.ora <- function(mode=2, eset, cmat, perm=1000, 
    alpha=0.05, beta=1, sig.stat=c("xxP", "xxFC", "|", "&"))
{
    GRP.COL <- config.ebrowser("GRP.COL")
    ADJP.COL <- config.ebrowser("ADJP.COL")

    x <- assay(eset)
    y <- colData(eset)[, GRP.COL]

    # execute hypergeom ORA if no permutations
    fdat <- rowData(eset, use.names=TRUE)
    if(perm == 0) res.tbl <- .oraHypergeom(fdat, cmat, alpha, beta, sig.stat)
    # else do resampling using functionality of SAFE
    else{
        # resampl ORA
        if(mode == 1){
            nr.sigs <- sum(.isSig(fdat, alpha, beta, sig.stat))
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
.gsea <- function(
    eset, 
    gs.gmt, 
    perm=1000, 
    out.file=NULL)
{        
    GRP.COL <- config.ebrowser("GRP.COL")
    
    # npGSEA
    if(perm==0)
    {
        npGSEA <- pTwoSided <- NULL
        .isAvailable("npGSEA", type="software")
        gsc <- .gsList2Collect(gs.gmt)
        res <- npGSEA(x=assay(eset), y=eset[[GRP.COL]], set=gsc)
        ps <- sapply(res, pTwoSided)
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
    res <- GSEA(input.ds=as.data.frame(assay(eset)), 
        input.cls=cls, gs.db=gs.gmt, output.directory=out.dir, nperm=perm)
      
    gs.ps <- S4Vectors::as.matrix(res[,3:5])
    rownames(gs.ps) <- res[,1]

    return(gs.ps)
}

# 5 EBM (_E_mpirical _B_rowns _M_ethod)
.ebm <- function(eset, cmat)
{
    empiricalBrownsMethod <- NULL
    .isAvailable("EmpiricalBrownsMethod", type="software")
    pcol <-  rowData(eset, use.names=TRUE)[, config.ebrowser("ADJP.COL")]
    e <- assay(eset)
    gs.ps <- apply(cmat, 2, function(s) empiricalBrownsMethod(e[s,], pcol[s]))
    return(gs.ps)
}


# 6 GSA
.gsa <- function(eset, gs, perm=1000)
{  
    GSA <- NULL
    .isAvailable("GSA", type="software")
 
    minsize <- config.ebrowser("GS.MIN.SIZE")
    maxsize <- config.ebrowser("GS.MAX.SIZE")

    blk <- NULL
    BLK.COL <- config.ebrowser("BLK.COL")
    if(BLK.COL %in% colnames(colData(eset))) blk <- colData(eset)[,BLK.COL] 
    paired <- !is.null(blk)
    resp.type <- ifelse(paired, "Two class paired", "Two class unpaired")

    # prepare input
    x <- assay(eset)
    y <- eset[[config.ebrowser("GRP.COL")]] + 1
    genenames <- rownames(eset)

    # run GSA
    res <- GSA(x=x, y=y, nperms=perm, genesets=gs, resp.type=resp.type,
        genenames=genenames, minsize=minsize, maxsize=maxsize)
   
    # format output
    ps <- cbind(res$pvalues.lo, res$pvalues.hi)
    ps <- 2 * apply(ps, 1, min)
    scores <- res$GSA.scores
    res.tbl <- cbind(scores, ps)
    colnames(res.tbl) <- c("SCORE", config.ebrowser("GSP.COL"))
    rownames(res.tbl) <- names(gs)

    return(res.tbl)
}

# rseq: GSA maxmean stat as global.stat for safe
global.GSA <- function(cmat, u, ...)
{
    # SparseM::as.matrix
    .isAvailable("SparseM", type="software")
    #pos <- grep("SparseM", search())
    am <- getMethod("as.matrix", signature="matrix.csr")#, where=pos)
    tcmat <- t(am(cmat))

    return(
        function(u, cmat2=tcmat) 
        {
            ind.pos <- u > 0 
            
            upos <- u[ind.pos]
            lpos <- rowSums(cmat2[,ind.pos])
            vpos <- as.vector(cmat2[,ind.pos] %*% upos) / lpos
            vpos <- sapply(vpos, function(x) ifelse(is.na(x), 0, x))
            
            uneg <- abs(u[!ind.pos])
            lneg <- rowSums(cmat2[,!ind.pos])
            vneg <- as.vector(cmat2[,!ind.pos] %*% uneg) / lneg
            vneg <- sapply(vneg, function(x) ifelse(is.na(x), 0, x))

            mm <- apply(cbind(vpos, vneg), 1, max)
            return(mm)
        }
    )
}

# 7 PADOG
.padog <- function(eset, gs, perm=1000)
{
    padog <- NULL
    .isAvailable("PADOG", type="software")

    grp <- eset[[config.ebrowser("GRP.COL")]]
    grp <- ifelse(grp == 0, "c", "d") 
  
    blk <- NULL
    BLK.COL <- config.ebrowser("BLK.COL")
    if(BLK.COL %in% colnames(colData(eset))) blk <- colData(eset)[,BLK.COL] 
    paired <- !is.null(blk)
  
    nmin <- config.ebrowser("GS.MIN.SIZE")
  
    res <- padog(esetm=assay(eset), group=grp, 
        paired=paired, block=blk, gslist=gs, Nmin=nmin, NI=perm)
  
    res.tbl <- res[, c("meanAbsT0", "padog0", "PmeanAbsT", "Ppadog")]
    colnames(res.tbl) <- c("MEAN.ABS.T0", "PADOG0", "P.MEAN.ABS.T",  config.ebrowser("GSP.COL"))
    rownames(res.tbl) <- as.vector(res[,"ID"]) 
    return(res.tbl)
}


# compute gene frequencies across genesets
.getGeneFreqWeights <- function(cmat)
{
    gf <- rowSums(cmat)
    if (!all(gf == 1)) 
    {
        q99 <- quantile(gf, 0.99)
        m3sd <- mean(gf) + 3 * sd(gf)
        if(q99 > m3sd) gf[gf > q99] <- q99
        gff <- function(x) 1 + ((max(x) - x)/(max(x) - min(x)))^0.5
        gf <- gff(gf)
    } 
    else 
    {
        gf <- rep(1, nrow(cmat))
        names(gf) <- rownames(cmat)
    }
    return(gf)
}

# rseq: PADOG weighted mean as global.stat for safe
global.PADOG <- function(cmat, u, args.global)
{
    # SparseM::as.matrix
    .isAvailable("SparseM", type="software")
    #pos <- grep("SparseM", search())
    am <- getMethod("as.matrix", signature="matrix.csr")#, where=pos)
    cmat <- t(am(cmat))
    gs.size <- rowSums(cmat) 

    return(
        function(u, cmat2=cmat, gf=args.global$gf, gs.size=rowSums(cmat)) 
        {
            wu <- abs(u) * gf
            return(as.vector(cmat2 %*% wu) / gs.size)
        }
    )
}

# 8a MGSA 
.mgsa <- function(eset, gs, alpha=0.05, beta=1, sig.stat=c("xxP", "xxFC", "|", "&"))
{
    mgsa <- setsResults <- NULL
    .isAvailable("mgsa", type="software")
    
    # extract significant (DE) genes
    isig <- .isSig(rowData(eset, use.names=TRUE), alpha, beta, sig.stat)
    obs <- rownames(eset)[isig]
    pop <- rownames(eset)
  
    # run mgsa
    res <- mgsa(o=obs, sets=gs, population=pop)
    res <- setsResults(res)[,1:3]
    res[,3] <- 1 - res[,3]
    colnames(res)[3] <- config.ebrowser("GSP.COL")
    return(res)
}

# 8b GLOBALTEST
.globaltest <- function(eset, gs, perm=1000)
{
    gt <- NULL
    .isAvailable("globaltest", type="software")

    grp <- colData(eset)[, config.ebrowser("GRP.COL")]
    res <- gt(grp, eset, subsets=gs, permutations=perm)
    res <- res@result[,2:1]
    colnames(res) <- c("STAT", config.ebrowser("GSP.COL"))
    return(res)
}


# 9 ROAST
# 10 CAMERA
.roast.camera <- function(method=c("roast", "camera"), eset, gs, perm=1000, rseq=FALSE)
{
    method <- match.arg(method)

    # design matrix
    grp <- colData(eset)[, config.ebrowser("GRP.COL")]
    blk <- NULL
    BLK.COL <- config.ebrowser("BLK.COL")
    if(BLK.COL %in% colnames(colData(eset))) blk <- colData(eset)[,BLK.COL]
   
    group <- factor(grp)
    paired <- !is.null(blk)
    f <- "~" 
    if(paired) 
    {   
        block <- factor(blk)
        f <- paste0(f, "block + ") 
    }   
    f <- formula(paste0(f, "group"))
    design <- model.matrix(f)

    y <- assay(eset)
    # rseq data
    if(rseq)
    {
        y <- edgeR::DGEList(counts=y,group=grp)
        y <- edgeR::calcNormFactors(y)
        y <- edgeR::estimateDisp(y, design)
    }
    
    # set gene sets
    gs.index <- limma::ids2indices(gs, rownames(eset))
    
    # run roast / camera
    if(method == "roast")
    {
        res <- limma::mroast(y, gs.index, design, 
                                nrot=perm, adjust.method="none", sort="none")
        res <- res[,c("NGenes", "Direction", "PValue")]
        colnames(res) <- c("NR.GENES", "DIR", config.ebrowser("GSP.COL"))
        res[,"DIR"] <- ifelse(res[,"DIR"] == "Up", 1, -1)
    }
    else
    { 
        res <- limma::camera(y, gs.index, design, sort=FALSE)
        res <- res[,1:4]
        colnames(res) <- c("NR.GENES", "COR", "DIR", config.ebrowser("GSP.COL"))
        res[,"DIR"] <- ifelse(res[,"DIR"] == "Up", 1, -1)
    }
    return(res)
}


# 11 GSVA
.gsva <- function(eset, gs, rseq=FALSE)
{
     gsva <- NULL
    .isAvailable("GSVA", type="software")
  
    # compute GSVA per sample enrichment scores
    es <- gsva(expr=assay(eset), gset.idx.list=gs, rnaseq=rseq)
    es <- es$es.obs
  
    # set design matrix
    grp <- colData(eset)[, config.ebrowser("GRP.COL")]
    blk <- NULL
    BLK.COL <- config.ebrowser("BLK.COL")
    if(BLK.COL %in% colnames(colData(eset))) blk <- colData(eset)[,BLK.COL]

    group <- factor(grp)
    paired <- !is.null(blk)
    f <- "~"
    if(paired)
    {
        block <- factor(blk)
        f <- paste0(f, "block + ")
    }
    f <- formula(paste0(f, "group"))
    design <- model.matrix(f)  
   
    # fit the linear model to the GSVA enrichment scores
    fit <- limma::lmFit(es, design)
    fit <- limma::eBayes(fit)
    res <- limma::topTable(fit, number=nrow(es), coef="group1", sort.by="none", adjust.method="none")
    
    # process output
    res <- res[,c("t", "P.Value")]
    colnames(res) <- c("t.SCORE", config.ebrowser("GSP.COL"))
    
    return(res)
}
