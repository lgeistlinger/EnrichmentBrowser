# TODO: Add comment
# 
# Author: ludwig geistlinger
# Date: June 24th 2010
#
# description:  collection of functions that convert outputs of
#               enrichment methods to html, i.e. visualizes them
#               user-friendly
#
###############################################################################

ea.browse <- function(  
    res, nr.show=-1, set.view=TRUE, graph.view=NULL)
{
    html.out <- to.html(res=res, 
        nr.show=nr.show, set.view=set.view, graph.view=graph.view)
    browseURL(paste0("file://", html.out))
}

to.html <- function(    
    res, 
    nr.show=-1, 
    set.view=TRUE, 
    graph.view=NULL,
    html.out=NULL)
{
    method <- res$method
    if(class(method) != "character") method <- NA
    # create out dir
    if(is.null(html.out)){
        res.dir <- file.path(system.file(package="EnrichmentBrowser"), "results")
        if(!file.exists(res.dir)) dir.create(res.dir)
        html.out <- paste(method, Sys.time())
        html.out <- gsub("[^a-z0-9]", "", html.out)
        html.out <- file.path(res.dir, paste(html.out, "html", sep="."))
    }

    out.dir <- dirname(html.out) 

    if(nr.show < 1) nr.show <- res$nr.sigs
    if(nr.show > nrow(res$res.tbl)) nr.show <- nrow(res$res.tbl)

    if(method=="comb") 
        comb2html(res=res, nr.show=nr.show, 
            graph.view=graph.view, html.out=html.out)

    else
    {
        eset <- res$eset
        gs <- res$gs
        alpha <- res$alpha
        res <- res$res.tbl[seq_len(nr.show),,drop=FALSE]

        # create set view if enabled
        sviews <- NULL
        if(set.view) 
            sviews <- t(sapply(res[,1], view.set, 
                eset=eset, gs=gs, out.dir=out.dir, alpha=alpha))

        # graph view
        gviews <- NULL
        if(!is.null(graph.view))
        { 
            if(class(graph.view) == "character") 
                graph.view <- read.grn(graph.view)
            gviews <- t(sapply(res[,1], view.graph, eset=eset, 
                gs=gs, grn=graph.view, alpha=alpha, out.dir=out.dir))
        }
    
        method2html(    
            method=toupper(method), 
            res.tbl=res, 
            set.view=sviews, 
            graph.view=gviews,
            html.out=html.out)
    }
    return(html.out)

}

view.set <- function(s, eset, gs, out.dir, alpha)
{
    set.name <- ifelse(nchar(s) > 8, substring(s,1,8), s)
    is.kegg <- length(grep("^[a-z]{3}[0-9]{5}", set.name))
    if(!is.kegg) set.name <- gsub("[\\W_]", "", set.name, perl=TRUE)

    # get set genes and their expression
    s.genes <- gs[[s]]
    ind <- featureNames(eset) %in% s.genes
    s.genes <- featureNames(eset)[ind]
    expr <- exprs(eset)[ind,]
    gs.eset <- new("ExpressionSet", exprs=expr)
    fData(gs.eset) <- fData(eset)[ind, c(FC.COL, ADJP.COL)]
    pData(gs.eset) <- pData(eset)

    ##
    # 1 PLOT (heatmap, pval distribution, volcano)
    ##
    plot.file <- file.path(out.dir, paste0(set.name, "_setview.pdf"))
    if(!file.exists(plot.file))
    {
        pdf(file=plot.file)
        exprs.heatmap(gs.eset)
        pdistr(gs.eset)
        volcano(gs.eset)
        dev.off()
    }
    views <- basename(plot.file)

    ##
    # 2 BROWSE
    ##
    fDat <- as.matrix(fData(gs.eset))
    if(is.kegg)
    {
        sig.genes <- s.genes[fDat[,2] < alpha]
        url <- get.html.of.marked.pathway.by.objects(pwy=set.name, oids=sig.genes)
        views <- c(url, views)
    }
    
    ##
    # 3 REPORT
    ##
    report.file <- sub("pdf$", "txt", plot.file)
    if(!file.exists(report.file))
    {
        fDat <- signif(fDat[order(fDat[,2]),], digits=2)
        fDat <- cbind(rownames(fDat), fDat)
        colnames(fDat)[1] <- "GENE"
        write.table(fDat, file=report.file, row.names=FALSE, quote=FALSE, sep="\t")
    }
    views <- c(views, basename(report.file))
    return(rev(views))
}
    
view.graph <- function(s, eset, gs, grn, alpha, out.dir)
{
    set.name <- ifelse(nchar(s) > 8, substring(s,1,8), s)
    is.kegg <- length(grep("^[a-z]{3}[0-9]{5}", set.name))
    if(!is.kegg) set.name <- gsub("[\\W_]", "", set.name, perl=TRUE)
    org <- ifelse(is.kegg, substring(set.name, 1, 3), NA)

    genes <- gs[[s]]
    sgrn <- query.grn(genes, grn, index=FALSE)
    if(nrow(sgrn) > 0) 
        sggea.graph <- construct.ggea.graph(
            grn=sgrn, eset=eset, alpha=alpha, org=org)
    else sggea.graph <- NULL

    ##
    # 1 PLOT
    ##
    # open plot device  
    plot.file <- file.path(out.dir, paste0(set.name, "_graphview.pdf"))
    if(!file.exists(plot.file))
    {
        pdf(file=plot.file)
        if(!is.null(sggea.graph)) 
            plot.ggea.graph(sggea.graph, 
                show.scores=(numEdges(sggea.graph) < NROW.TOP.TABLE))
        else plot(NA, axes=FALSE, xlim=c(0,1), ylim=c(0,1), 
            ylab="", xlab="", main="No edges in network for this set!")
        dev.off()
    }
    views <- basename(plot.file)
        
    ##
    # 2 REPORT
    ##
    report.file <- sub("pdf$", "txt", plot.file)
    if(!file.exists(report.file))
    {
        if(!is.null(sggea.graph))
        {
            consistency <- sggea.graph@renderInfo@edges$label
            cons.tbl <- cbind(names(consistency), consistency)
            cons.tbl <- cons.tbl[order(as.numeric(consistency), decreasing=TRUE),]
            colnames(cons.tbl) <- c("EDGE", "CONSISTENCY")
            write.table(cons.tbl, 
                file=report.file, row.names=FALSE, quote=FALSE, sep="\t")
        }
        else cat("No edges in network for this set!", file=report.file)
    }
    views <- c(basename(report.file), views)
    return(views)
}

get.html.of.marked.pathway.by.objects <- function(pwy, oids)
{
    pwy <- sub("^path:", "", pwy)
    oids <- gsub("[a-z]{3}:", "", oids)
    coids <- paste(oids, collapse="+")
    request <- pwy
    if(nchar(coids)) request <- paste(request, coids, sep="+")
    return(paste0(KEGG.SHOW.URL, request))
}

comb2html <- function(res, nr.show, graph.view, html.out)
{
    eset <- res$eset
    gs <- res$gs
    alpha <- res$alpha
    res <- as.matrix(gs.ranking(res))

    # create combinded tables, sorted by each method's ranking 
    rank.cols <- grep("*.RANK", colnames(res), value=TRUE)
    res.tbls <- sapply(rank.cols, 
        function(rc) res[order(as.integer(res[,rc])),], simplify=FALSE)

    res.tbls$COMB.PVAL <- res
    
    out.dir <- dirname(html.out)
    html.outs <- sapply(sub("\\.","",rank.cols), 
        function(rc) sub(".html$", paste0("_", rc, ".html"), html.out))
    html.outs <- c(html.outs, html.out) 

    # produce html
    for(i in seq_along(res.tbls))
    {
        res <- res.tbls[[i]]    
        if(nrow(res) > nr.show) res <- res[seq_len(nr.show),,drop=FALSE]
        if(i == 1)
        { 
            sviews <- t(sapply(res[,1], view.set, 
                eset=eset,gs=gs, out.dir=out.dir, alpha=alpha))
        
            gviews <- NULL
            if(!is.null(graph.view))
            {
                if(class(graph.view) == "character") 
                    graph.view <- read.grn(graph.view)
                gviews <- t(sapply(res[,1], view.graph, eset=eset, 
                    gs=gs, grn=graph.view, alpha=alpha, out.dir=out.dir))
            }
        }
        else
        {
            sviews <- sviews[res[,1],]
            if(!is.null(graph.view)) gviews <- gviews[res[,1],]
        }
        method2html(
            method="eBrowser",
            res.tbl=res,
            set.view=sviews,
            graph.view=gviews,
            html.out=html.outs[i],
            header.links=basename(html.outs))
    } 
}


