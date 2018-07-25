############################################################
#
# author: Ludwig Geistlinger
# date: 06 Dec 2010
#
# Set-based Enrichment Analysis (SBEA)
#
############################################################

#' @rdname sbea
#' @export
sbeaMethods <- function() 
    c("ora", "safe", "gsea", "gsa", "padog", "globaltest", 
        "roast", "camera", "gsva", "samgs", "ebm", "mgsa")

# INPUT FASSADE - wrapping & delegation


#' Set-based enrichment analysis (SBEA)
#' 
#' This is the main function for the enrichment analysis of gene sets.  It
#' implements and wraps existing implementations of several frequently used
#' methods and allows a flexible inspection of resulting gene set rankings.
#' 
#' 'ora': overrepresentation analysis, simple and frequently used test based on
#' the hypergeometric distribution (see Goeman and Buhlmann, 2007, for a
#' critical review).
#' 
#' 'safe': significance analysis of function and expression, generalization of
#' ORA, includes other test statistics, e.g. Wilcoxon's rank sum, and allows to
#' estimate the significance of gene sets by sample permutation; implemented in
#' the safe package (Barry et al., 2005).
#' 
#' 'gsea': gene set enrichment analysis, frequently used and widely accepted,
#' uses a Kolmogorov-Smirnov statistic to test whether the ranks of the
#' p-values of genes in a gene set resemble a uniform distribution (Subramanian
#' et al., 2005).
#' 
#' 'padog': pathway analysis with down-weighting of overlapping genes,
#' incorporates gene weights to favor genes appearing in few pathways versus
#' genes that appear in many pathways; implemented in the PADOG package.
#' 
#' 'roast': rotation gene set test, uses rotation instead of permutation for
#' assessment of gene set significance; implemented in the limma and edgeR
#' packages for microarray and RNA-seq data, respectively.
#' 
#' 'camera': correlation adjusted mean rank gene set test, accounts for
#' inter-gene correlations as implemented in the limma and edgeR packages for
#' microarray and RNA-seq data, respectively.
#' 
#' 'gsa': gene set analysis, differs from GSEA by using the maxmean statistic,
#' i.e. the mean of the positive or negative part of gene scores in the gene
#' set; implemented in the GSA package.
#' 
#' 'gsva': gene set variation analysis, transforms the data from a gene by
#' sample matrix to a gene set by sample matrix, thereby allowing the
#' evaluation of gene set enrichment for each sample; implemented in the GSVA
#' package.
#' 
#' 'globaltest': global testing of groups of genes, general test of groups of
#' genes for association with a response variable; implemented in the
#' globaltest package.
#' 
#' 'samgs': significance analysis of microarrays on gene sets, extends the SAM
#' method for single genes to gene set analysis (Dinu et al., 2007).
#' 
#' 'ebm': empirical Brown's method, combines $p$-values of genes in a gene set
#' using Brown's method to combine $p$-values from dependent tests; implemented
#' in the EmpiricalBrownsMethod package.
#' 
#' 'mgsa': model-based gene set analysis, Bayesian modeling approach taking set
#' overlap into account by working on all sets simultaneously, thereby reducing
#' the number of redundant sets; implemented in the mgsa package.
#' 
#' It is also possible to use additional set-based enrichment methods.  This
#' requires to implement a function that takes 'se', 'gs', 'alpha', and 'perm'
#' as arguments and returns a numeric vector 'ps' storing the resulting p-value
#' for each gene set in 'gs'. This vector must be named accordingly (i.e.
#' names(ps) == names(gs)). See examples.
#' 
#' @aliases ora gsea
#' @param method Set-based enrichment analysis method.  Currently, the
#' following set-based enrichment analysis methods are supported: \sQuote{ora},
#' \sQuote{safe}, \sQuote{gsea}, \sQuote{padog}, \sQuote{roast},
#' \sQuote{camera}, \sQuote{gsa}, \sQuote{gsva}, \sQuote{globaltest},
#' \sQuote{samgs}, \sQuote{ebm}, and \sQuote{mgsa}.  For basic ora also set
#' 'perm=0'. Default is \sQuote{ora}.  This can also be the name of a
#' user-defined function implementing set-based enrichment. See Details.
#' @param se Expression dataset.  An object of class
#' \code{\linkS4class{SummarizedExperiment}}.  Mandatory minimal annotations:
#' \itemize{ \item colData column storing binary group assignment (named
#' "GROUP") \item rowData column storing (log2) fold changes of differential
#' expression between sample groups (named "FC") \item rowData column storing
#' adjusted (corrected for multiple testing) p-values of differential
#' expression between sample groups (named "ADJ.PVAL") } Additional optional
#' annotations: \itemize{ \item colData column defining paired samples or
#' sample blocks (named "BLOCK") \item metadata slot named "annotation" giving
#' the organism under investigation in KEGG three letter code (e.g. "hsa" for
#' Homo sapiens) \item metadata slot named "dataType" indicating the expression
#' data type ("ma" for microarray, "rseq" for RNA-seq) }
#' @param gs Gene sets.  Either a list of gene sets (character vectors of gene
#' IDs) or a text file in GMT format storing all gene sets under investigation.
#' @param alpha Statistical significance level. Defaults to 0.05.
#' @param perm Number of permutations of the sample group assignments.
#' Defaults to 1000. For basic ora set 'perm=0'.  Using
#' method="gsea" and 'perm=0' invokes the permutation approximation from the
#' npGSEA package.
#' @param padj.method Method for adjusting nominal gene set p-values to
#' multiple testing.  For available methods see the man page of the stats
#' function \code{\link{p.adjust}}.  Defaults to'none', i.e. leaves the nominal
#' gene set p-values unadjusted.
#' @param out.file Optional output file the gene set ranking will be written
#' to.
#' @param browse Logical. Should results be displayed in the browser for
#' interactive exploration? Defaults to FALSE.
#' @param ...  Additional arguments passed to individual sbeaMethods.  This
#' includes currently for ORA and MGSA: \itemize{ \item beta: Log2 fold change
#' significance level. Defaults to 1 (2-fold).  \item sig.stat: decides which
#' statistic is used for determining significant DE genes.  Options are:
#' \itemize{ \item 'p' (Default): genes with p-value below alpha.  \item 'fc':
#' genes with abs(log2(fold change)) above beta \item '&': p & fc (logical AND)
#' \item '|': p | fc (logical OR) } }
#' @return sbeaMethods: a character vector of currently supported methods;
#' 
#' sbea: if(is.null(out.file)): an enrichment analysis result object that can
#' be detailedly explored by calling \code{\link{eaBrowse}} and from which a
#' flat gene set ranking can be extracted by calling \code{\link{gsRanking}}.
#' If 'out.file' is given, the ranking is written to the specified file.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso Input: \code{\link{readSE}}, \code{\link{probe2gene}}
#' \code{\link{getGenesets}} to retrieve gene sets from databases such as GO 
#' and KEGG.
#' 
#' Output: \code{\link{gsRanking}} to retrieve the ranked list of gene sets.
#' \code{\link{eaBrowse}} for exploration of resulting gene sets.
#' 
#' Other: \code{\link{nbea}} to perform network-based enrichment analysis.
#' \code{\link{combResults}} to combine results from different methods.
#' @references Goeman and Buhlmann (2007) Analyzing gene expression data in
#' terms of gene sets: methodological issues.  Bioinformatics, 23, 980-7.
#' 
#' Barry et al. (2005) Significance Analysis of Function and Expression.
#' Bioinformatics, 21:1943-9.
#' 
#' Subramanian et al. (2005) Gene Set Enrichment Analysis: a knowledge-based
#' approach for interpreting genome-wide expression profiles.  Proc Natl Acad
#' Sci USA, 102:15545-50.
#' 
#' Dinu et al. (2007) Improving gene set analysis of microarray data by SAM-GS.
#' BMC Bioinformatics, 8:242
#' @examples
#' 
#'     # currently supported methods
#'     sbeaMethods()
#' 
#'     # (1) expression data: 
#'     # simulated expression values of 100 genes
#'     # in two sample groups of 6 samples each
#'     se <- makeExampleData(what="SE")
#'     se <- deAna(se)
#' 
#'     # (2) gene sets:
#'     # draw 10 gene sets with 15-25 genes
#'     gs <- makeExampleData(what="gs", gnames=names(se))
#' 
#'     # (3) make 2 artificially enriched sets:
#'     sig.genes <- names(se)[rowData(se)$ADJ.PVAL < 0.1]
#'     gs[[1]] <- sample(sig.genes, length(gs[[1]])) 
#'     gs[[2]] <- sample(sig.genes, length(gs[[2]]))   
#' 
#'     # (4) performing the enrichment analysis
#'     ea.res <- sbea(method="ora", se=se, gs=gs, perm=0)
#' 
#'     # (5) result visualization and exploration
#'     gsRanking(ea.res)
#' 
#'     # using your own tailored function as enrichment method
#'     dummySBEA <- function(se, gs, alpha, perm)
#'     {
#'         sig.ps <- sample(seq(0, 0.05, length=1000), 5)
#'         nsig.ps <- sample(seq(0.1, 1, length=1000), length(gs)-5)
#'         ps <- sample(c(sig.ps, nsig.ps), length(gs))
#'         names(ps) <- names(gs)
#'         return(ps)
#'     }
#' 
#'     ea.res2 <- sbea(method=dummySBEA, se=se, gs=gs)
#'     gsRanking(ea.res2) 
#' 
#' @export sbea
sbea <- function(   
    method=EnrichmentBrowser::sbeaMethods(), 
    se, 
    gs, 
    alpha=0.05, 
    perm=1000, 
    padj.method="none",
    out.file=NULL,
    browse=FALSE, ...)
{   
    # get configuration
    GS.MIN.SIZE <- configEBrowser("GS.MIN.SIZE")
    GS.MAX.SIZE <- configEBrowser("GS.MAX.SIZE")
    PVAL.COL <- configEBrowser("PVAL.COL")
    FC.COL <-  configEBrowser("FC.COL")
    ADJP.COL <-  configEBrowser("ADJP.COL")

    ### TEMPORARY: will be replaced by as(eSet,SummarizedExperiment)
    if(is(se, "ExpressionSet")) se <- as(se, "RangedSummarizedExperiment")
    ###    

    # dealing with NA's
    nr.na <- sum(is.na(rowData(se)[,FC.COL]))
    if(nr.na) se <- se[!is.na(rowData(se)[,FC.COL]),]
    nr.na <- sum(is.na(rowData(se)[,ADJP.COL]))
    if(nr.na) se <- se[!is.na(rowData(se)[,ADJP.COL]),]    

    # getting gene sets
    if(!is.list(gs)) gs <- getGenesets(gs)

    # restrict se and gs to intersecting genes
    igenes <- intersect(rownames(se), unique(unlist(gs)))
    se <- se[igenes,]
    gs <- lapply(gs, function(s) s[s %in% igenes]) 
    lens <- lengths(gs)
    gs <- gs[lens >= GS.MIN.SIZE & lens <= GS.MAX.SIZE]

    if(is.character(method))
    { 
        method <- match.arg(method)
        data.type <- metadata(se)$dataType
        if(is.null(data.type)) data.type <- .detectDataType(assay(se))

        # rseq? 
        if(data.type == "rseq")
        {
            # mgsa
            if(method == "mgsa") gs.ps <- .mgsa(se, gs, alpha, ...)
            # globaltest
            else if(method == "globaltest") gs.ps <- .globaltest(se, gs, perm)
            # roast & camera
            else if(method %in% c("roast", "camera"))
                gs.ps <- .roast.camera(method, se, gs, perm, rseq=TRUE)
		    # gsva
		    else if(method == "gsva") gs.ps <- .gsva(se, gs, rseq=TRUE)
            else
            {
                # gs2cmat
                #cmat <- .gmt2cmat(gs, rownames(se), GS.MIN.SIZE, GS.MAX.SIZE)
                #if(nrow(cmat) < nrow(se)) se <- se[rownames(cmat),] 
                f <- file()
                sink(file=f)
                cmat <- safe::getCmatrix(gs, as.matrix=TRUE)
                sink()
                close(f)
                se <- se[rownames(cmat),]

                # ora
                if(method == "ora" & perm==0) 
                    gs.ps <- .ora(1, se, cmat, perm, alpha, ...)
                # ebm
                else if(method == "ebm") gs.ps <- .ebm(se, cmat)
		        # all others
                else gs.ps <- .rseqSBEA(method, se, cmat, perm, alpha)
            }
        }
        # microarray
        else
        { 
            # gsea
            if(method == "gsea") gs.ps <- .gsea(se, gs, perm, padj.method)
            # gsa
            else if(method == "gsa") gs.ps <- .gsa(se, gs, perm)
            # padog
	        else if(method == "padog") gs.ps <- .padog(se, gs, perm)
	        # mgsa
	        else if(method == "mgsa") gs.ps <- .mgsa(se, gs, alpha, ...)
            # globaltest
            else if(method == "globaltest") gs.ps <- .globaltest(se, gs, perm)
            # roast & camera
            else if(method %in% c("roast", "camera"))
                gs.ps <- .roast.camera(method, se, gs, perm)
	        # gsva
	        else if(method == "gsva") gs.ps <- .gsva(se, gs)
            else
            {
                # gs2cmat
                #cmat <- .gmt2cmat(gs, rownames(se), GS.MIN.SIZE, GS.MAX.SIZE)
                f <- file()
                sink(file=f)
                cmat <- safe::getCmatrix(gs, as.matrix=TRUE)
                sink()
                close(f)
                se <- se[rownames(cmat),]

                # ora
                if(method == "ora") 
                    gs.ps <- .ora(1, se, cmat, perm, alpha, padj.method, ...)
                #safe
                else if(method == "safe") 
                    gs.ps <- .ora(2, se, cmat, perm, alpha, padj.method)
                # ebm
                else if(method == "ebm") gs.ps <- .ebm(se, cmat)
                # samgs
                else if(method == "samgs")
                {
                    if(is.null(out.file)) out.dir <- configEBrowser("OUTDIR.DEFAULT")
                    else out.dir <- sub("\\.[a-z]+$","_files", out.file)
                    if(!file.exists(out.dir)) dir.create(out.dir, recursive=TRUE)
                    samt.file <- file.path(out.dir, "samt.RData")
                    GRP.COL <- configEBrowser("GRP.COL")
                    gs.ps <- SAMGS(GS=as.data.frame(cmat), DATA=assay(se), 
                        cl=as.factor(as.integer(se[[GRP.COL]])), 
                        nbPermutations=perm, 
                        tstat.file=samt.file)
                }
            }
        }
    }
    else if(is.function(method)) 
        gs.ps <- method(se=se, gs=gs, alpha=alpha, perm=perm)
    else stop(paste(method, "is not a valid method for sbea"))

    res.tbl <- data.frame(signif(gs.ps, digits=3))
    sorting.df <- res.tbl[,ncol(res.tbl)]
    if(ncol(res.tbl) > 1) 
        sorting.df <- cbind(sorting.df, -res.tbl[,rev(seq_len(ncol(res.tbl)-1))])
    else colnames(res.tbl)[1] <- PVAL.COL 
    res.tbl <- res.tbl[do.call(order, as.data.frame(sorting.df)), , drop=FALSE]

    # use built-in FDR control for safe and gsea
    #fdr.methods <- c("gsea", "safe") 
    #if(is.function(method) || !(method %in% fdr.methods))
	
	if(padj.method != "none")
    {   
        adjp <- p.adjust(res.tbl[,PVAL.COL], method=padj.method)
        res.tbl <- cbind(res.tbl, adjp)
        colnames(res.tbl)[ncol(res.tbl)] <- configEBrowser("ADJP.COL")
    }     

    res.tbl <- DataFrame(rownames(res.tbl), res.tbl)
    colnames(res.tbl)[1] <- configEBrowser("GS.COL")
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
            nr.sigs=sum(res.tbl[,PVAL.COL] < alpha),
            se=se, gs=gs, alpha=alpha)
        if(browse) eaBrowse(res)
        else return(res)
    }
}

#' @rdname eaBrowse
#' @export
gsRanking <- function(res, signif.only=TRUE)
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

#' @export
#' @keywords internal
gs.ranking <- function(res, signif.only=TRUE)
{
    .Deprecated("gsRanking")
    gsRanking(res, signif.only)
}

###
#
# HELPER
# 
###
.gmt2cmat <- function(gs, features, min.size=0, max.size=Inf)
{
    if(is.character(gs)) gs <- getGenesets(gs)
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

# deAna as local.stat for safe
local.deAna <- function (X.mat, y.vec, args.local)
{
    return(function(data, ...) 
    {
        stat <- deAna(expr=data, grp=y.vec,
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

.rseqSBEA <- function(method, se, cmat, perm, alpha)
{
	assign("se", se, envir=.GlobalEnv)
    assign("local.deAna", local.deAna, envir=.GlobalEnv)
    de.method <- grep(".STAT$", colnames(rowData(se)), value=TRUE)
    de.method <- sub(".STAT$",  "", de.method)
    
    blk <- NULL
    blk.col <- configEBrowser("BLK.COL") 
    if(blk.col %in% colnames(colData(se))) blk <- colData(se)[,blk.col]

    args.local <- list(de.method=de.method, blk=blk)

    args.global <- list(one.sided=FALSE)
    if(method == "ora")
    {
        global <- "Fisher"
        nr.sigs <- sum(rowData(se)[, configEBrowser("ADJP.COL")] < alpha)
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

	x <- assay(se)
    y <- colData(se)[,configEBrowser("GRP.COL")]
    gs.ps <- safe::safe(X.mat=x, y.vec=y, C.mat=cmat,         
        local="deAna", args.local=args.local,
        global=global, args.global=args.global, 
        Pi.mat=perm, alpha=alpha, error="none")
 
    res.tbl <- cbind(
            gs.ps@global.stat, 
            gs.ps@global.stat / colSums(cmat), 
            gs.ps@global.pval)
    
    colnames(res.tbl) <- c("GLOB.STAT", "NGLOB.STAT", configEBrowser("PVAL.COL"))
    
    return(res.tbl)
}

.isSig <- function(rdat, alpha=0.05, beta=1, sig.stat=c("xxP", "xxFC", "|", "&"))
{
    FC.COL <- configEBrowser("FC.COL")
    ADJP.COL <- configEBrowser("ADJP.COL")

    sig.stat <- sig.stat[1]
    if(grepl("P$", sig.stat))
    {
        if(sig.stat == "xxP") sig <- rdat[, ADJP.COL] < alpha
        else
        {
            perc <- as.integer(substring(sig.stat, 1, 2))
            p <- rdat[,ADJP.COL]
            names(p) <- rownames(rdat) 
            ordp <- sort(p)
            nr.sig <- round( length(p) * (perc / 100) )
            sigs <- names(ordp)[seq_len(nr.sig)]
            sig <- rownames(rdat) %in% sigs
        }
    }
    else if(grepl("FC$", sig.stat))
    { 
        if(sig.stat == "xxFC") sig <- abs(rdat[, FC.COL]) > beta
        else
        {
            perc <- as.integer(substring(sig.stat, 1, 2))
            fc <- rdat[,FC.COL]
            names(fc) <- rownames(rdat) 
            ordfc <- fc[order(abs(fc), decreasing=TRUE)]
            nr.sig <- round( length(fc) * (perc / 100) )
            sigs <- names(ordfc)[seq_len(nr.sig)]
            sig <- rownames(rdat) %in% sigs
        }
    }
    else 
    {
        psig <- rdat[, ADJP.COL] < alpha
        fcsig <- abs(rdat[, FC.COL]) > beta
        sig <- do.call(sig.stat, list(psig, fcsig))
    }
    return(sig)
}

# 1 HYPERGEOM ORA
.oraHypergeom <- function(rdat, cmat, 
    alpha=0.05, beta=1, sig.stat=c("xxP", "xxFC", "|", "&"))
{
    # determine sig. diff. exp. genes of se, 
    # corresponds to sample size from urn
    isig <- .isSig(rdat, alpha, beta, sig.stat)
    nr.sigs <- sum(isig)
    
    # determine overlap of sig and set genes for each set
    sig.cmat <- cmat & isig

    # white balls observed when drawing nr.sigs balls from the urn
    ovlp.sizes <- colSums(sig.cmat)

    # white balls in the urn  (genes in gene set)
    gs.sizes <- colSums(cmat) 
    # black balls in the urn (genes not in gene set)
    uni.sizes <- nrow(rdat) - gs.sizes 

    # determine significance of overlap 
    # based on hypergeom. distribution
    gs.ps <- phyper(ovlp.sizes-1, gs.sizes, uni.sizes, nr.sigs, lower.tail=FALSE) 
    
    res.tbl <- cbind(gs.sizes, ovlp.sizes, gs.ps)
    colnames(res.tbl) <- c("NR.GENES", "NR.SIG.GENES", configEBrowser("PVAL.COL"))
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
.ora <- function(mode=2, se, cmat, perm=1000, alpha=0.05, 
    padj.method="none", beta=1, sig.stat=c("xxP", "xxFC", "|", "&"))
{
    GRP.COL <- configEBrowser("GRP.COL")
    ADJP.COL <- configEBrowser("ADJP.COL")

    x <- assay(se)
    y <- colData(se)[, GRP.COL]

    # execute hypergeom ORA if no permutations
    rdat <- rowData(se, use.names=TRUE)
    if(perm == 0) res.tbl <- .oraHypergeom(rdat, cmat, alpha, beta, sig.stat)
    # else do resampling using functionality of SAFE
    else{
        # use built-in p-adjusting?
        padj <- switch(padj.method,
                        BH = "FDR.BH",
                        fdr = "FDR.BH",
                        bonferroni = "FWER.Bonf",
                        holm = "FWER.Holm",
                        BY = "FDR.YB",
                        "none")

        # resampl ORA
        if(mode == 1){
            nr.sigs <- sum(.isSig(rdat, alpha, beta, sig.stat))
            args <- list(one.sided=FALSE, genelist.length=nr.sigs)

            gs.ps <- safe::safe(X.mat=x, y.vec=y, global="Fisher", C.mat=cmat, 
                 Pi.mat=perm, alpha=alpha, error=padj, args.global=args)
        } 
        # SAFE default                  
        else gs.ps <- safe::safe(X.mat=x, y.vec=y, 
            C.mat=cmat, Pi.mat=perm, alpha=alpha, error=padj)
        pval <- if(padj == "none") gs.ps@global.pval else gs.ps@global.error
        res.tbl <- cbind(
            gs.ps@global.stat, 
            gs.ps@global.stat / colSums(cmat), 
            pval)
        colnames(res.tbl) <- c("GLOB.STAT", "NGLOB.STAT", configEBrowser("PVAL.COL"))
    }
    return(res.tbl)
}

# 4 GSEA
.gsea <- function(
    se, 
    gs.gmt, 
    perm=1000,
    padj.method="none", 
    out.file=NULL)
{        
    GRP.COL <- configEBrowser("GRP.COL")
    
    # npGSEA
    if(perm==0)
    {
        npGSEA <- pTwoSided <- NULL
        isAvailable("npGSEA", type="software")
        gsc <- .gsList2Collect(gs.gmt)
        res <- npGSEA(x=assay(se), y=se[[GRP.COL]], set=gsc)
        ps <- sapply(res, pTwoSided)
        names(ps) <- names(gs.gmt)
        return(ps)
    }

    # build class list
    cls <- list()
    cls$phen <- levels(as.factor(se[[GRP.COL]]))
    cls$class.v <- ifelse(se[[GRP.COL]] == cls$phen[1], 0, 1)

    if(is.null(out.file)) 
        out.dir <- configEBrowser("OUTDIR.DEFAULT") 
    else out.dir <- sub("\\.[a-z]+$", "_files", out.file)
    if(!file.exists(out.dir)) dir.create(out.dir, recursive=TRUE)
     
    # use built-in p-adjusting?
    padj.method <- switch(padj.method,
                        BH = "fdr",
                        fdr = "fdr",
                        bonferroni = "fwer",
                        holm = "fwer",
                        BY = "fdr",
                        "none")
   
    # run GSEA
    res <- GSEA(input.ds=as.data.frame(assay(se)), 
                input.cls=cls, gs.db=gs.gmt, nperm=perm,
                padj.method=padj.method, output.directory=out.dir)
      
    gs.ps <- S4Vectors::as.matrix(res[,3:5])
    rownames(gs.ps) <- res[,1]

    return(gs.ps)
}

# 5 EBM (_E_mpirical _B_rowns _M_ethod)
.ebm <- function(se, cmat)
{
    empiricalBrownsMethod <- NULL
    isAvailable("EmpiricalBrownsMethod", type="software")
    pcol <-  rowData(se, use.names=TRUE)[, configEBrowser("ADJP.COL")]
    e <- assay(se)
    gs.ps <- apply(cmat, 2, function(s) empiricalBrownsMethod(e[s,], pcol[s]))
    return(gs.ps)
}


# 6 GSA
.gsa <- function(se, gs, perm=1000)
{  
    # setup
    GSA <- NULL
    isAvailable("GSA", type="software")
 
    minsize <- configEBrowser("GS.MIN.SIZE")
    maxsize <- configEBrowser("GS.MAX.SIZE")
    GRP.COL <- configEBrowser("GRP.COL")
    BLK.COL <- configEBrowser("BLK.COL")
    
    # prepare input
    x <- assay(se)
    genenames <- names(se)
    y <- se[[GRP.COL]] + 1
   
    # paired?
    blk <- NULL
    if(BLK.COL %in% colnames(colData(se))) blk <- se[[BLK.COL]] 
    paired <- !is.null(blk)
    resp.type <- ifelse(paired, "Two class paired", "Two class unpaired")

    # response vector y need to be differently coded for 2-class paired   
    if(paired)
    {
        y <- blk
        ublk <- unique(blk)
        for(i in seq_along(ublk)) y[blk==ublk[i]] <- c(i,-i)
        y <- as.integer(y)
    }

    # run GSA
    res <- GSA(x=x, y=y, nperms=perm, genesets=gs, resp.type=resp.type,
        genenames=genenames, minsize=minsize, maxsize=maxsize)
   
    # format output
    ps <- cbind(res$pvalues.lo, res$pvalues.hi)
    ps <- 2 * apply(ps, 1, min)
    scores <- res$GSA.scores
    res.tbl <- cbind(scores, ps)
    colnames(res.tbl) <- c("SCORE", configEBrowser("PVAL.COL"))
    rownames(res.tbl) <- names(gs)

    return(res.tbl)
}

# rseq: GSA maxmean stat as global.stat for safe
global.GSA <- function(cmat, u, ...)
{
    # SparseM::as.matrix
    isAvailable("SparseM", type="software")
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
.padog <- function(se, gs, perm=1000)
{
    padog <- NULL
    isAvailable("PADOG", type="software")

    grp <- se[[configEBrowser("GRP.COL")]]
    grp <- ifelse(grp == 0, "c", "d") 
  
    blk <- NULL
    BLK.COL <- configEBrowser("BLK.COL")
    if(BLK.COL %in% colnames(colData(se))) blk <- colData(se)[,BLK.COL] 
    paired <- !is.null(blk)
  
    nmin <- configEBrowser("GS.MIN.SIZE")
    perm <- as.numeric(perm) 
 
    res <- padog(assay(se), group=grp, 
        paired=paired, block=blk, gslist=gs, Nmin=nmin, NI=perm)
  
    res.tbl <- res[, c("meanAbsT0", "padog0", "PmeanAbsT", "Ppadog")]
    colnames(res.tbl) <- c("MEAN.ABS.T0", 
        "PADOG0", "P.MEAN.ABS.T",  configEBrowser("PVAL.COL"))
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
    isAvailable("SparseM", type="software")
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
.mgsa <- function(se, gs, alpha=0.05, beta=1, sig.stat=c("xxP", "xxFC", "|", "&"))
{
    mgsa <- setsResults <- NULL
    isAvailable("mgsa", type="software")
    
    # extract significant (DE) genes
    isig <- .isSig(rowData(se, use.names=TRUE), alpha, beta, sig.stat)
    obs <- rownames(se)[isig]
    pop <- rownames(se)
  
    # run mgsa
    res <- mgsa(o=obs, sets=gs, population=pop)
    res <- setsResults(res)[,1:3]
    res[,3] <- 1 - res[,3]
    colnames(res)[3] <- configEBrowser("PVAL.COL")
    return(res)
}

# 8b GLOBALTEST
.globaltest <- function(se, gs, perm=1000)
{
    gt <- NULL
    isAvailable("globaltest", type="software")

    grp <- colData(se)[, configEBrowser("GRP.COL")]
    se <- as(se, "ExpressionSet")
    res <- gt(grp, se, subsets=gs, permutations=perm)
    res <- res@result[,2:1]
    colnames(res) <- c("STAT", configEBrowser("PVAL.COL"))
    return(res)
}

# 9 ROAST
# 10 CAMERA
.roast.camera <- function(method=c("roast", "camera"), se, gs, perm=1000, rseq=FALSE)
{
    method <- match.arg(method)

    # design matrix
    grp <- colData(se)[, configEBrowser("GRP.COL")]
    blk <- NULL
    BLK.COL <- configEBrowser("BLK.COL")
    if(BLK.COL %in% colnames(colData(se))) blk <- colData(se)[,BLK.COL]
   
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

    y <- assay(se)
    # rseq data
    if(rseq)
    {
        y <- edgeR::DGEList(counts=y,group=grp)
        y <- edgeR::calcNormFactors(y)
        y <- edgeR::estimateDisp(y, design)
    }
    
    # set gene sets
    gs.index <- limma::ids2indices(gs, rownames(se))
    
    # run roast / camera
    if(method == "roast")
        res <- limma::mroast(y, gs.index, design, 
                                nrot=perm, adjust.method="none", sort="none")
    else res <- limma::camera(y, gs.index, design, sort=FALSE)
    res <- res[,c("NGenes", "Direction", "PValue")]
    colnames(res) <- c("NR.GENES", "DIR", configEBrowser("PVAL.COL"))
    res[,"DIR"] <- ifelse(res[,"DIR"] == "Up", 1, -1)

    return(res)
}


# 11 GSVA
.gsva <- function(se, gs, rseq=FALSE)
{
    gsva <- NULL
    isAvailable("GSVA", type="software")
  
    # compute GSVA per sample enrichment scores
    kcdf <- ifelse(rseq, "Poisson", "Gaussian")
    es <- gsva(expr=assay(se), gset.idx.list=gs, kcdf=kcdf)
  
    # set design matrix
    grp <- colData(se)[, configEBrowser("GRP.COL")]
    blk <- NULL
    BLK.COL <- configEBrowser("BLK.COL")
    if(BLK.COL %in% colnames(colData(se))) blk <- colData(se)[,BLK.COL]

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
    colnames(res) <- c("t.SCORE", configEBrowser("PVAL.COL"))
    
    return(res)
}
