############################################################
# 
# author: Ludwig Geistlinger
# date: 2015-06-15 13:32:53
# 
# descr: example datasets for man pages 
# 
############################################################



#' Example data for the EnrichmentBrowser package
#' 
#' Functionality to construct example data sets for demonstration.  This
#' includes expression data, gene sets, gene regulatory networks, and
#' enrichment analysis results.
#' 
#' 
#' @param what Kind of example data set to be constructed.  This should be one
#' out of: \itemize{ \item SE: SummarizedExperiment \item gs: Gene set list
#' \item grn: Gene regulatory network \item ea.res: Enrichment analysis result
#' object as returned by the functions \code{\link{sbea}} and
#' \code{\link{nbea}} }
#' @param ...  Additional arguments to fine-tune the specific example data
#' sets.
#' 
#' For what='SE': \itemize{ \item type: Expression data type. Should be either
#' 'ma' (Microarray intensity measurements) or 'rseq' (RNA-seq read counts).
#' \item nfeat: Number of features/genes. Defaults to 100.  \item nsmpl: Number
#' of samples. Defaults to 12.  \item blk: Create sample blocks. Defaults to
#' TRUE.  \item norm: Should the expression data be normalized? Defaults to
#' FALSE.  \item deAna: Should an differential expression analysis be carried
#' out automatically? Defaults to FALSE.  }
#' 
#' For what='gs': \itemize{ \item gnames: gene names from which the sets will
#' be sampled. Per default the sets will be drawn from c(g1, ..., g100).  \item
#' n: number of sets. Defaults to 10.  \item min.size: minimal set size.
#' Defaults to 15.  \item max.size: maximal set size. Defaults to 25.  }
#' 
#' For what='grn': \itemize{ \item nodes: gene node names for which edges will
#' be drawn.  Per default node names will be c(g1, ..., g100).  \item
#' edge.node.ratio: ratio number of edges / number of nodes.  Defaults to 3,
#' i.e. creates 3 times more edges than nodes.  }
#' 
#' For what='ea.res': \itemize{ \item SE: SummarizedExperiment. Calls
#' makeExampleData(what="SE") per default.  \item gs: Gene sets. Calls
#' makeExampleData(what="gs") per default.  \item method: Enrichment analysis
#' method. Defaults to 'ora'.  \item alpha: Statistical significance level.
#' Defaults to 0.05.  }
#' @return Depends on the 'what' argument.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @examples
#' 
#'     se <- makeExampleData(what="SE")
#' 
#' @export makeExampleData
makeExampleData <- function(
    what=c("SE", "gs", "grn", "ea.res"), ...)
{
    what <- match.arg(what)
    res <- NULL

    if(what=="SE") res <- .makeExmplSE(...)
    else if(what=="gs") res <- .makeExmplGS(...)
    else if(what=="grn") res <- .makeExmplGRN(...)
    else if(what=="ea.res") res <- .makeExmplRes(...)
    else stop("No example data for", what)

    return(res)
}

.makeExmplSE <- function(type=c("ma", "rseq"), 
    nfeat=100, nsmpl=12, blk=TRUE, norm=FALSE, de.ana=FALSE)
{
    type <- match.arg(type)

    nr.cases <- ceiling(nsmpl/2)
    nr.controls <- nsmpl - nr.cases

    # (a) microarray data: intensity measurements
    if(type=="ma")
    { 
        control.grid <- rnorm(nfeat * nr.controls, mean=6, sd=1)
        case.grid <- rnorm(nfeat * nr.cases, mean=7.5, sd=1)
        exmpl.exprs <- cbind(
            matrix(control.grid, nrow=nfeat, ncol=nr.controls),
            matrix(case.grid, nrow=nfeat, ncol=nr.cases)
        )
        
    }
    # (b) RNA-seq data: read counts
    else
    {
        isAvailable("DESeq2", type="software")
        exds <- DESeq2::makeExampleDESeqDataSet(n=nfeat)
        exmpl.exprs <- assay(exds)
    }
    rownames(exmpl.exprs) <- paste0("g", seq_len(nfeat))
    colnames(exmpl.exprs) <- paste0("s", seq_len(nsmpl))
    se <- SummarizedExperiment(assays=list(exprs=exmpl.exprs))
    
    # sample binary group assignment, e.g. two different treatments
    se$GROUP <- c(rep(0, nr.controls), rep(1, nr.cases))

    # sample batch effects: 3 sample blocks
    if(blk) se$BLOCK <- sample(rep(1:3, ceiling(nsmpl/3)))[seq_len(nsmpl)]
    
    if(norm) se <- normalize(se)
    if(de.ana) se <- deAna(se) 
    
    return(se)
}

.makeExmplGS <- function(gnames=NULL, n=10, min.size=15, max.size=25)
{
    if(is.null(gnames)) gnames <- paste0("g", 1:100)
    gs <- replicate(n, sample(gnames, sample(min.size:max.size,1)))
    names(gs) <- paste0("gs", seq_len(n))
    return(gs)
}

.makeExmplGRN <- function(nodes=NULL, edge.node.ratio=3)
{
    if(is.null(nodes)) nodes <- paste0("g", 1:100)
    nr.edges <- edge.node.ratio * length(nodes)           
    grn <- replicate(nr.edges, sample(nodes, 2))
    grn <- cbind(t(grn), sample(c("+","-"), nr.edges, replace=TRUE))
    colnames(grn) <- c("FROM", "TO", "TYPE")
    return(grn)
}

.makeExmplRes <- function(se=NULL, gs=NULL, method="ora", alpha=0.05)
{
    if(is.null(se)) se <- .makeExmplSE()
    if(is.null(gs)) gs <- .makeExmplGS(gnames=rownames(se))
    # (3) make artificial enrichment analysis results:
    # 2 ea methods with 5 significantly enriched gene sets each
    gs.res <- sample(names(gs))
    nr.sigs <- floor(length(gs) / 2)
    p.res <- sort(c(
        runif(nr.sigs, alpha * 10^(-5), alpha), 
        runif(length(gs) - nr.sigs, alpha+0.1, 1)
        ))
    p.res <- round(p.res, digits=3)
    
    res.tbl <- DataFrame(gs.res, p.res)
    colnames(res.tbl) <- c(configEBrowser("GS.COL"), configEBrowser("PVAL.COL"))

    ea.res <- list(method=method,
        res.tbl=res.tbl, nr.sigs=nr.sigs, se=se, gs=gs, alpha=alpha)

}
