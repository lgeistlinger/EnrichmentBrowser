############################################################
# 
# author: Ludwig Geistlinger
# date: 2015-06-15 13:32:53
# 
# descr: example datasets for man pages 
# 
############################################################

make.example.data <- function(
    what=c("eset", "gs", "grn", "ea.res"), ...)
{
    what <- match.arg(what)
    res <- NULL

    if(what=="eset") res <- make.exmpl.eset(...)
    else if(what=="gs") res <- make.exmpl.gs(...)
    else if(what=="grn") res <- make.exmpl.grn(...)
    else if(what=="ea.res") res <- make.exmpl.ea.res(...)
    else stop("No example data for", what)

    return(res)
}

make.exmpl.eset <- function(type=c("ma", "rseq"), 
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
        exds <- DESeq2::makeExampleDESeqDataSet(n=nfeat)
        exmpl.exprs <- SummarizedExperiment::assays(exds)[[1]]
    }
    rownames(exmpl.exprs) <- paste0("g", seq_len(nfeat))
    colnames(exmpl.exprs) <- paste0("s", seq_len(nsmpl))
    eset <- new("ExpressionSet", exprs=exmpl.exprs)
    
    # sample binary group assignment, e.g. two different treatments
    pData(eset)$GROUP <- c(rep(0, nr.controls), rep(1, nr.cases))

    # sample batch effects: 3 sample blocks
    if(blk) pData(eset)$BLOCK <- sample(rep(c(1,2,3), ceiling(nsmpl/3)))[seq_len(nsmpl)]
    
    if(norm) eset <- normalize(eset)
    if(de.ana) eset <- de.ana(eset) 
    
    return(eset)
}

make.exmpl.gs <- function(gnames=NULL, n=10, min.size=15, max.size=25)
{
    if(is.null(gnames)) gnames <- paste0("g", 1:100)
    gs <- replicate(n, sample(gnames, sample(min.size:max.size,1)))
    names(gs) <- paste0("gs", seq_len(n))
    return(gs)
}

make.exmpl.grn <- function(nodes=NULL, edge.node.ratio=3)
{
    if(is.null(nodes)) nodes <- paste0("g", 1:100)
    nr.edges <- edge.node.ratio * length(nodes)           
    grn <- replicate(nr.edges, sample(nodes, 2))
    grn <- cbind(t(grn), sample(c("+","-"), nr.edges, replace=TRUE))
    colnames(grn) <- c("FROM", "TO", "TYPE")
    return(grn)
}

make.exmpl.ea.res <- function(eset=NULL, gs=NULL, method="ora", alpha=0.05)
{
    if(is.null(eset)) eset <- make.exmpl.eset()
    if(is.null(gs)) gs <- make.exmpl.gs(gnames=featureNames(eset))
    # (3) make artificial enrichment analysis results:
    # 2 ea methods with 5 significantly enriched gene sets each
    gs.res <- sample(names(gs))
    nr.sigs <- floor(length(gs) / 2)
    p.res <- sort(c(
        runif(nr.sigs, alpha * 10^(-5), alpha), 
        runif(length(gs) - nr.sigs, alpha+0.1, 1)
        ))
    p.res <- round(p.res, digits=3)

    res.tbl <- cbind(gs.res, p.res)
    colnames(res.tbl) <- c(config.ebrowser("GS.COL"), config.ebrowser("GSP.COL"))

    ea.res <- list(method=method,
        res.tbl=res.tbl, nr.sigs=nr.sigs, eset=eset, gs=gs, alpha=alpha)

}
