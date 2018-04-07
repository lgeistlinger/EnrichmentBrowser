# 
# Author: ludwig geistlinger
# Date: June 24th 2010
#
# description:  collection of functions that convert outputs of
#               enrichment methods to html, i.e. visualizes them
#               user-friendly
#
# Update:
#   27 May 2015: 
#       - major reworking of HTML output based on ReportingTools & hwriter 
#       - interactive graphics based on imageMap from geneplotter & biocGraph
#
###############################################################################

.createIndex <- function(meth, comb)
{
    out.dir <- config.ebrowser("OUTDIR.DEFAULT")
    indexPage <- ReportingTools::HTMLReport(shortName = 'index',
                        title = 'EnrichmentBrowser: Index of Result Files',
                        basePath=out.dir, reportDirectory="reports")
    
    res.files <- list.files(out.dir, pattern="txt|png|html$")
    file.rename(from=file.path(out.dir, res.files), 
                    to=file.path(out.dir, "reports", res.files))

    vcol <- "global_sview.html"
    de.plot <- hwriter::hwriteImage(sub("sview.html", "volc.png", vcol),
        link=vcol, table = FALSE, height=200, width=200, target="_blank")

    de.tag <- hwriter::hwrite(de.plot, border=0, br=TRUE)
    ReportingTools::publish(de.tag, indexPage)
    delink <- ReportingTools::Link("DE Measures for each Gene", "de.txt")
    ReportingTools::publish(delink, indexPage)
    
    if(comb)
    {
        res.tag <- hwriter::hwrite("Combined Results", heading=4)
        ReportingTools::publish(res.tag, indexPage)
        ttlink <- ReportingTools::Link("Top Table", "comb.html")
        ReportingTools::publish(ttlink, indexPage)
        frlink <- ReportingTools::Link("Full Ranking", "comb.txt")
        ReportingTools::publish(frlink, indexPage) 
    }

    for(m in meth)
    {
        res.tag <- hwriter::hwrite(paste(toupper(m), "Results"), heading=4)
        ReportingTools::publish(res.tag, indexPage)
        ttlink <- ReportingTools::Link("Top Table", paste0(m, ".html"))
        ReportingTools::publish(ttlink, indexPage)
        frlink <- ReportingTools::Link("Full Ranking", paste0(m, ".txt"))
        ReportingTools::publish(frlink, indexPage) 
    }
    index <- ReportingTools::finish(indexPage)
    if(interactive())
    {
        message("HTML report: index.html")
        # this can be removed upon release of RStudio v1.2 (Summer 2018)
        if(Sys.getenv("RSTUDIO") == "1") index <- URLencode(index)
        browseURL(index)
    }
}




#' Exploration of enrichment analysis results
#' 
#' Functions to extract a flat gene set ranking from an enrichment analysis
#' result object and to detailedly explore it.
#' 
#' 
#' @aliases ea.browse gs.ranking
#' @param res Enrichment analysis result list (as returned by the functions
#' \code{\link{sbea}} and \code{\link{nbea}}).
#' @param nr.show Number of gene sets to show.  As default all statistically
#' significant gene sets are displayed.
#' @param graph.view Optional.  Should a graph-based summary (reports and
#' visualizes consistency of regulations) be created for the result?  If
#' specified, it needs to be a gene regulatory network, i.e. either an absolute
#' file path to a tabular file or a character matrix with exactly *THREE* cols;
#' 1st col = IDs of regulating genes; 2nd col = corresponding regulated genes;
#' 3rd col = regulation effect; Use '+' and '-' for activation/inhibition.
#' @param html.only Logical.  Should the html file only be written (without
#' opening the browser to view the result page)? Defaults to FALSE.
#' @param signif.only Logical.  Display only those gene sets in the ranking,
#' which satisfy the significance level? Defaults to TRUE.
#' @return gs.ranking: \code{\linkS4class{DataFrame}} with gene sets ranked by
#' the corresponding p-value;
#' 
#' ea.browse: none, opens the browser to explore results.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{\link{sbea}}, \code{\link{nbea}},
#' \code{\link{comb.ea.results}}
#' @examples
#' 
#'     
#'     # real data
#'     # (1) reading the expression data from file
#'     exprs.file <- system.file("extdata/exprs.tab", package="EnrichmentBrowser")
#'     cdat.file <- system.file("extdata/colData.tab", package="EnrichmentBrowser")
#'     rdat.file <- system.file("extdata/rowData.tab", package="EnrichmentBrowser")
#'     probeSE <- readSE(exprs.file, cdat.file, rdat.file)
#'     geneSE <- probe2gene(probeSE) 
#'     geneSE <- deAna(geneSE)
#'     metadata(geneSE)$annotation <- "hsa"
#' 
#'     # artificial enrichment analysis results
#'     gs <- makeExampleData(what="gs", gnames=names(geneSE))
#'     ea.res <- makeExampleData(what="ea.res", method="ora", se=geneSE, gs=gs)
#' 
#'     # (5) result visualization and exploration
#'     gs.ranking(ea.res)
#'     ea.browse(ea.res)
#' 
#' @export ea.browse
ea.browse <- function(res, nr.show=-1, graph.view=NULL, html.only=FALSE)
{
    method <- ifelse( is(res$method, "character"), res$method, NA)
    se <- res$se
    alpha <- res$alpha
    gs <- res$gs
    
    # create out dir
    out.dir <- config.ebrowser("OUTDIR.DEFAULT")
    if(!file.exists(out.dir)) dir.create(out.dir, recursive=TRUE)
    rep.dir <- file.path(out.dir, "reports")
    
    # how many gene sets to show in the output?
    if(nr.show < 1) nr.show <- res$nr.sigs
    if(nr.show > nrow(res$res.tbl)) nr.show <- nrow(res$res.tbl)

    # add description & nr.genes per gene set
    res <- res$res.tbl[seq_len(nr.show),]

    # expecting transition from gene set lists to collections in the near future
    # gsc <- res$gsc
    gsc <- .gsList2Collect(gs[res[,1]])
    res[,1] <- vapply(res[,1], 
        function(s) unlist(strsplit(s, "_"))[1],
        character(1))
    
    is.kegg <- is(collectionType(gsc[[1]]), "KEGGCollection")
    is.go <- is(collectionType(gsc[[1]]), "GOCollection")

    gs.title <- vapply(gsc, description, character(1))
    nr.genes <- vapply(gsc, function(g) length(geneIds(g)), integer(1))

    cnames <- c(colnames(res)[1], "TITLE")
    resn <- DataFrame(res[,1], gs.title)
    if(!("NR.GENES" %in% colnames(res)))
    {
        cnames <- c(cnames, "NR.GENES")
        resn <- DataFrame(resn, nr.genes)
    }
    cnames <- c(cnames, colnames(res)[2:ncol(res)])
    resn <- DataFrame(resn, res[,2:ncol(res)]) 
    colnames(resn) <- cnames
    res <- resn

    # make gene pages
    im <- incidence(gsc)
    org <- organism(gsc[[1]])
    if(org == "") org <- metadata(se)$annotation
    if(!length(org)) stop("Organism annotation not found!\n", 
        "Organism under study must be annotated via metadata(se)$annotation")

    message("Creating gene report ...")
    se <- se[colnames(im),]
    rcols <- sapply(c("FC.COL", "ADJP.COL"), config.ebrowser)
    fDat <- rowData(se, use.names=TRUE)[,rcols]
    fDat <- as.data.frame(fDat)
    gt <- suppressMessages(.geneTable(im, org, fcs=fDat))
    gn.cols <- sapply(c("SYM.COL", "GN.COL"), config.ebrowser)
    rowData(se)[,gn.cols] <- DataFrame(gt[,gn.cols]) 
    gt.reps <- sapply(gsc, function(s) .geneReport(s, gt, out.dir)) 
    
    # make gene set page
    # (1) link .geneReport
    link <- paste0(names(gsc), ".html")
    res[,"NR.GENES"] <- hwriter::hwrite(res[,"NR.GENES"], link=link, table=FALSE)
   
    # set view: volcano, heatmap 
    message("Creating set view ...")
    out.prefix <- file.path(rep.dir, names(gsc))
    names(out.prefix) <- names(gsc)
    vcol <- sapply(gsc, function(s) 
        .viewSet(se[geneIds(s),], out.prefix[setName(s)]))
    vcol <- hwriter::hwriteImage(sub("sview.html", "volc.png", vcol),
        link=vcol, table = FALSE, height=50, width=50, target="_blank")
    res <- DataFrame(res, vcol)
    colnames(res)[ncol(res)] <- "SET.VIEW" 

    # path view: kegg maps
    if(is.kegg)
    { 
        message("Creating kegg view ...")
        isAvailable("pathview", type="software")
        vcol <- sapply(gsc, function(s) 
            .viewPath(setName(s), se[geneIds(s),], out.prefix[setName(s)]))
        vcol <- hwriter::hwriteImage(sub("kview.html", "kpath.png", vcol),
            link=vcol, table = FALSE, height=50, width=50, target="_blank")
        res <- DataFrame(res, vcol)
        colnames(res)[ncol(res)] <- "PATH.VIEW"
    }

    # graph view: ggeaGraphs
    if(!is.null(graph.view)) 
    {
        message("Creating graph view ...")
        vcol <- sapply(gsc, function(s) .viewGraph(se[geneIds(s),], .queryGRN(
            geneIds(s), graph.view, index=FALSE), alpha, out.prefix[setName(s)]))
        vcol <- hwriter::hwriteImage(sub("html$", "png", vcol),
            link=vcol, table = FALSE, height=50, width=50, target="_blank")
        res <- DataFrame(res, vcol)
        colnames(res)[ncol(res)] <- "GRAPH.VIEW"
    }
 
    # (2) link KEGG / GO 
    link <- NULL
    GS.COL <- config.ebrowser("GS.COL")
    if(is.kegg) link <- sapply(gsc, function(s) 
        .getHTMLOfMarkedPathway(setName(s), geneIds(s)[fDat[geneIds(s),2] < alpha]))
    else if(is.go) link <- paste0(config.ebrowser("GO.SHOW.URL"), res[,GS.COL])
    if(!is.null(link)) res[,GS.COL] <- 
        hwriter::hwrite(res[,GS.COL], link=link, table=FALSE)

    htmlRep <- ReportingTools::HTMLReport(shortName=method,
        title=paste(toupper(method), config.ebrowser("RESULT.TITLE"), sep=" - "),
        basePath=out.dir, reportDirectory="reports")
    res <- as.data.frame(res)
    ReportingTools::publish(res, htmlRep) 
        #colClasses = c(rep("sort-string-robust", 3),
        #    rep("sort-num-robust", ncol(gt)-3 )))
    rep <- ReportingTools::finish(htmlRep)
    if(!html.only && interactive())
    { 
        message(paste("Your output files are in", rep.dir))
        message(paste0("HTML report: ", method, ".html"))
        # this can be removed upon release of RStudio v1.2 (Summer 2018)
        if(Sys.getenv("RSTUDIO") == "1") rep <- URLencode(rep)
        browseURL(rep)
    }
}

.makeView <- function(html1, html2, gene.html.pos=c("bottom", "topright"))
{
    gene.html.pos <- match.arg(gene.html.pos)
    s <- unlist(strsplit(html1, "_"))[1]
    
    head <- hwriter::hmakeTag("head", hwriter::hmakeTag("title", s))
    html1.tag <- hwriter::hmakeTag("frame", name="volc", src=html1)
    html2.tag <- hwriter::hmakeTag("frame", name="hmap", src=html2)
    html3.tag <- hwriter::hmakeTag("frame", name="gene", scrolling="auto")
    
    if(gene.html.pos == "topright")
    {
        bkp <- html3.tag
        html3.tag <- html2.tag
        html2.tag <- bkp
    }

    f1.tag <- hwriter::hmakeTag("framse", 
        paste0(sub("</frame>", "", c(html1.tag, html2.tag)), collapse=""), 
        cols=paste0(config.ebrowser("PLOT.WIDTH") + 30, ",*"), border=0)
    html3.tag <- sub("</frame>", "", html3.tag)
    f2.tag <- hwriter::hmakeTag("framse", paste(f1.tag, html3.tag),
        rows=paste0(config.ebrowser("PLOT.HEIGHT") + 30, ",*"), border=0)
    cont <- hwriter::hmakeTag("html", paste0(head,f2.tag))
    return(cont)
}

.viewSet <- function(se, out.prefix)
{
    out.files <- paste(out.prefix, 
        c("sview.html", "volc.png", "hmap.png"), sep="_" )
    
    if(!file.exists(out.files[1]))
    {
        ##
        # 1 PLOT: heatmap & volcano
        ##
        volc.html <- .makeVolcHTML(se, out.files[2])
        hmap.html <- .makeHmapHTML(se, out.files[3])
        cont <- .makeView(volc.html, hmap.html) 
        cat(cont, file=out.files[1])
    }
    views <- basename(out.files[1])
    return(views)

    # flat file set report
    # Do we need this anymore?
#    report.file <- sub("pdf$", "txt", plot.file)
#    if(!file.exists(report.file))
#    {
#        fDat <- signif(fDat[order(fDat[,2]),], digits=2)
#        fDat <- cbind(rownames(fDat), fDat)
#        colnames(fDat)[1] <- "GENE"
#        write.table(fDat, file=report.file, row.names=FALSE, quote=FALSE, sep="\t")
#    }
#    views <- c(views, basename(report.file))
#    return(rev(views))
}

.viewPath <- function(s, se, out.prefix)
{
    org <- substring(s, 1, 3)
    pwy.id <- sub("^[a-z]{3}", "", s)
    fc <- rowData(se)[,config.ebrowser("FC.COL")]
    gn.cols <- sapply(c("SYM.COL", "GN.COL"), config.ebrowser)
    gnam <- apply(rowData(se)[,gn.cols], 1, paste, collapse=": ")
    names(fc) <- names(gnam) <- names(se)

    out.files <- paste(out.prefix, 
        c("kview.html", "kpath.png", "kgraph.png"), sep="_")
    
    if(!file.exists(out.files[1]) && pwy.id != "01100")
    {
        ##
        # 1 PLOT: kegg.native & kegg.graph
        ##
        kpath.html <- .makeKPathHTML(fc, pwy.id, org, out.files[2])
        kgraph.html <- .makeKGraphHTML(fc, gnam, pwy.id, org, out.files[3])
        cont <- .makeView(kgraph.html, kpath.html, gene.html.pos="topright") 
        cat(cont, file=out.files[1])
    }
    
    views <- basename(out.files[1])
    return(views)
} 

.viewGraph <- function(se, sgrn, alpha, out.prefix)
{
    out.files <- paste0(out.prefix, "_gview.", c("html", "png", "txt"))
    
    if(!file.exists(out.files[1]))
    {
        ##
        # 1 PLOT: ggeaGraph
        ##
        ggraph.html <- .makeGGraphHTML(se, sgrn, alpha, out.files[2])
        void.html <- "void.html"
        void.file <- file.path(dirname(out.prefix), void.html)
        cat(hwriter::hmakeTag('html'), file=void.file)
        cont <- .makeView(ggraph.html, void.html)   
        cat(cont, file=out.files[1])
    }
    
    views <- basename(out.files[1])
    return(views)
}

.makeHmapHTML <- function(se, img.file)
{
    width <- config.ebrowser("PLOT.WIDTH") 
    height <- config.ebrowser("PLOT.HEIGHT")

    # 1: make the plot
    # (a) heatmap 1: all features vs all samples
    expr <- assay(se)
    rownames(expr) <- rowData(se)[,config.ebrowser("SYM.COL")]
    png(img.file, width=width, height=height)
    exprs.heatmap(expr=expr, grp=se[[config.ebrowser("GRP.COL")]])
    dev.off()
    img.tag <- hwriter::hwriteImage(basename(img.file))

    # (b) heatmap 2: most signif features
    max.row <- 40
    fc <- abs(rowData(se)[,config.ebrowser("FC.COL")])
    p <- -log(rowData(se)[,config.ebrowser("ADJP.COL")], base=10)
    ind <- (fc >= 1) | (p >= 1)
    se <- se[ind,]
    if(nrow(se) > 1)
    {
        if(nrow(se) > max.row)
        {
            fc <- fc[ind]
            p <- p[ind]
            score <- sqrt(fc^2 + p^2)
            se <- se[order(score, decreasing=TRUE)[seq_len(max.row)],]
            #    # select most variable samples
            #if(ncol(se) > max.col)
            #{
            #    svar <- apply(assay(se), 2, var)
            #    se <- se[,order(svar, decreasing=TRUE)[seq_len(max.col)]]
            #}
        }
        expr <- assay(se)
        rownames(expr) <- rowData(se)[,config.ebrowser("SYM.COL")]
        img2 <- sub(".png$", "2.png", img.file)
        png(img2, width=width, height=height)
        exprs.heatmap(expr=expr, grp=se[[config.ebrowser("GRP.COL")]])
        dev.off()
        img.tag <- paste0(img.tag, hwriter::hwriteImage(basename(img2)))
    }
    
    # 2: make the html
    hmap.html <- sub("png$", "html", img.file)
    cont <- hwriter::hmakeTag('html', hwriter::hmakeTag('body', img.tag))
    cat(cont, file=hmap.html)
    return(basename(hmap.html))
}

.makeVolcHTML <- function(se, img.file)
{    
    width <- config.ebrowser("PLOT.WIDTH") 
    height <- config.ebrowser("PLOT.HEIGHT")

    fc <- rowData(se)[,config.ebrowser("FC.COL")]
    p <- rowData(se)[,config.ebrowser("ADJP.COL")]
    # 1: make the plot
    png(img.file, width=width, height=height)
        volcano(fc, p)
        mai.px <- par('mai') * 75
        usr.grd <- par('usr')
    dev.off()
    p <- -log(p, base=10)
    # 2: make the html

    # INFER LINK COORDS
	#
	# mai: a numerical vector of the form ‘c(bottom, left, top, right)
	#       which gives the margin size specified in inches
	#
	# default plot: (75px == 1 inch)
	#                       btm     left    top     right
	#       par('mai') :    1.02    0.82    0.82    0.42
	#           px:         76.5    61.5    61.5    31.5
	
	# actual plot region:
	#   - origin(x0,y0) on the top left is at (61.5, 61.5)
	#   - width: 500 - 61.5 - 31.5 = 407
	#   - height: 500 - 76.5 - 61.5 = 362
	
	plot.width <- width - sum(mai.px[c(2,4)])
	plot.height <- height - sum(mai.px[c(1,3)])
	
	# re-scale user grid
	min.x <- min(0, usr.grd[1])
	fc <- fc + abs(min.x)
	max.x <- usr.grd[2] + abs(min.x)
	min.y <- min(0, usr.grd[3])
	p <- p + abs(min.y)
	max.y <- usr.grd[4] + abs(min.y)
	
	y.scalef <- plot.height / max.y
	x.scalef <- plot.width / max.x
	
	cx <- fc * x.scalef + mai.px[2]
	cy <- p * y.scalef + mai.px[1]
	
	# turn y-coords around, (0,0) is topleft and not btmleft
	cy <- height - cy
	coord <- cbind(cx, cy, cx+10, cy-10)
        
    # volcano html
    volc.html <- sub("png$", "html", img.file) 
    con <- file(volc.html, open="w")
    
    gn.cols <- sapply(c("SYM.COL", "GN.COL"), config.ebrowser)
    titles <- apply(rowData(se)[,gn.cols], 1, paste, collapse=": ") 
    refs <- paste0(config.ebrowser("GENE.URL"), names(se))
    geneplotter::imageMap(coord, con, 
        list(HREF=refs, TITLE=titles, TARGET=rep("gene", nrow(coord))), 
        basename(img.file))    
    close(con)
    return(basename(volc.html))
}

.makeKPathHTML <- function(fc, pwy.id, org, img.file)
{
    # 1: make the plot
    # run pathview for getting native overplotted image
    out.dir <- dirname(img.file)
    #suppressWarnings(suppressMessages(
    pathview::pathview(gene.data=fc, 
        pathway.id=pwy.id, species=org, kegg.dir=out.dir, out.suffix="kpath")
    # ))
    pv.out <- file.path(getwd(), paste0(org, pwy.id, ".kpath.png"))
    file.rename(from=pv.out, to=img.file)
    
    # 2: make the html
    kpath.html <- sub("png$", "html", img.file)
    img.tag <- hwriter::hwriteImage(basename(img.file))
    cont <- hwriter::hmakeTag('html', hwriter::hmakeTag('body', img.tag))
    cat(cont, file=kpath.html)
    return(basename(kpath.html))
}

.makeKGraphHTML <- function(fc, gname, pwy.id, org, img.file)
{
    width <- config.ebrowser("PLOT.WIDTH") 
    height <- config.ebrowser("PLOT.HEIGHT")

    # 1: make the plot
    # run pathview2 for getting the kegg graph
    out.dir <- dirname(img.file)
    png(img.file, width=width, height=height)
    par(mai=rep(0,4))
    suppressWarnings(suppressMessages(
        gr <- pathview2(gene.data=fc, 
            pathway.id=pwy.id, species=org, kegg.dir=out.dir)
    ))
    dev.off()   
 
    # 2: make the html
    kgraph.html <- sub("png$", "html", img.file)
    if(is(gr, "graph"))
    {
        nd <- nodeRenderInfo(gr)$kegg.ids
        nam <- sapply(names(nd), function(n) 
            ifelse(nd[[n]][1] %in% names(gname), 
                    gname[nd[[n]][1]], 
                    nodeRenderInfo(gr)$label[[n]]))
        names(nam) <- names(nd)
        kstr <- sapply(nd, function(n) 
            paste(paste(org, n, sep=":"), collapse="+"), USE.NAMES=FALSE)
        con <- file(kgraph.html, open="w")
        refs <- paste0(config.ebrowser("KEGG.GENE.URL"), kstr)
        biocGraph::imageMap(gr, con=con,
            tags=list(HREF=refs, TITLE = nam, TARGET = rep("gene", length(nd))),
            imgname=basename(img.file), width=width, height=height)    
        close(con)
    }
    else cat(hwriter::hmakeTag('html'), file=kgraph.html)
    return(basename(kgraph.html))
}

.makeGGraphHTML <- function(se, sgrn, alpha, img.file)
{
    width <- config.ebrowser("PLOT.WIDTH") 
    height <- config.ebrowser("PLOT.HEIGHT")

    sggeaGraph <- NULL
    if(nrow(sgrn) > 0)
        sggeaGraph <- .constructGGEAGraph(grn=sgrn, se=se, alpha=alpha)
    
    # txt report
    report.file <- sub("png$", "txt", img.file)
    if(!is.null(sggeaGraph))
    {
            consistency <- sggeaGraph@renderInfo@edges$label
            cons.tbl <- cbind(names(consistency), consistency)
            cons.tbl <- cons.tbl[order(as.numeric(consistency), decreasing=TRUE),]
            colnames(cons.tbl) <- c("EDGE", "CONSISTENCY")
            write.table(cons.tbl,
                file=report.file, row.names=FALSE, quote=FALSE, sep="\t")
    }
    else cat("No edges in network for this set!", file=report.file)

    # ggeaGraph png
    png(img.file, width=width, height=height)
    par(mai=rep(0,4))
    if(!is.null(sggeaGraph))
        sggeaGraph <- .plotGGEAGraph(sggeaGraph,
            show.scores=(numEdges(sggeaGraph) < config.ebrowser("NR.SHOW")))
    else plot(NA, axes=FALSE, xlim=c(0,1), ylim=c(0,1),
        ylab="", xlab="", main="No edges in network for this set!")
    dev.off()
    
    ggraph.html <- sub("view.png$", "graph.html", img.file)
    if(!is.null(sggeaGraph))
    {
        gn.cols <- sapply(c("SYM.COL", "GN.COL"), config.ebrowser)
        gnam <- apply(rowData(se)[,gn.cols], 1, paste, collapse=": ")
        names(gnam) <- names(se)

        # image map
        nd <- nodes(sggeaGraph)
        nam <- gnam[nd]
        con <- file(ggraph.html, open="w")
        refs <- paste0(config.ebrowser("GENE.URL"),  nd)
        biocGraph::imageMap(sggeaGraph, con=con,
            tags=list(HREF=refs, TITLE = nam, TARGET = rep("gene", length(nd))),
            imgname=basename(img.file), width=width, height=height)    
        close(con)
    }
    else cat(hwriter::hmakeTag('html'), file=ggraph.html)
    return(basename(ggraph.html))
}

.getHTMLOfMarkedPathway <- function(pwy, oids)
{
    pwy <- sub("^path:", "", pwy)
    oids <- gsub("[a-z]{3}:", "", oids)
    coids <- paste(oids, collapse="+")
    request <- pwy
    if(nchar(coids)) request <- paste(request, coids, sep="+")
    return(paste0(config.ebrowser("KEGG.SHOW.URL"), request))
}

.sortGeneTable <- function(gt)
{
    gt <- as.data.frame(gt)
    num.cols <- sapply(gt, is.numeric)
    gt[,num.cols] <- signif(gt[,num.cols], digits=3)
    scols <- sapply(c("ADJP.COL", "FC.COL"), config.ebrowser)
    sorting.df <- gt[, scols]
    sorting.df[,2] <- -abs(sorting.df[,2])
    gt <- gt[do.call(order, sorting.df), , drop=FALSE]
    return(gt)    
}

.geneReport <- function(s, gt, out.dir)
{
    # (0) extract gene information
    gt <- gt[geneIds(s),]
    gt <- .sortGeneTable(gt)
        
    # (1) html table
    htmlRep <- ReportingTools::HTMLReport(basePath=out.dir, reportDirectory="reports",
        shortName=setName(s), title=paste(setName(s), "Gene Report", sep=": "))
    ReportingTools::publish(gt, htmlRep, .modifyDF=list(.ncbiGeneLink))#, .pubmedLink),
        #colClasses = c(rep("sort-string-robust", 3), rep("sort-num-robust", ncol(gt)-3 )))

    # (2) flat file
    fname <- paste0(setName(s), "_genes.txt")
    rep.dir <- file.path(out.dir, "reports")
    ofname <- file.path(rep.dir, fname)
    if(!file.exists(ofname))
        write.table(gt, file=ofname, sep="\t", quote=FALSE, row.names=FALSE)
    dlink <- ReportingTools::Link("Download .txt", fname)
    ReportingTools::publish(dlink, htmlRep)

    rep <- ReportingTools::finish(htmlRep)
    return(rep)
}

.getGeneAnno <- function(ids, org, biotype=FALSE)
{
    # load org pkg
    org.pkg <- .org2pkg(org)
    isAvailable(org.pkg)
    org.pkg <- get(org.pkg)

    # (1) gene identification 
    EZ.COL <- config.ebrowser("EZ.COL")
    gn.cols <- sapply(c("SYM.COL", "GN.COL"), config.ebrowser)
    gt <- suppressMessages(select(org.pkg, keys=ids,
            columns=gn.cols, keytype=EZ.COL))
	
	if(biotype)
	{
		useMart <- listDatasets <- useDataset <- getBM <- NULL
		isAvailable("biomaRt", type="software")
		
		id.type <- "entrezgene"
        message("Connecting to BioMart ...")
        ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")
        ds <- listDatasets(ensembl)[, "dataset"]
        ds <- grep(paste0("^", org), ds, value = TRUE)
        if (length(ds) == 0) 
            stop(paste("Mart not found for:", org))
        else if (length(ds) > 1) {
            message("Found several marts")
            sapply(ds, function(d) message(paste(which(ds == 
                d), d, sep = ": ")))
            n <- readline(paste0("Choose mart (1-", length(ds), 
                ") : "))
            ds <- ds[as.integer(n)]
        }
        ensembl <- useDataset(ds, mart = ensembl)
        attrs <- c(id.type, "gene_biotype")
        biot <- getBM(filters = id.type, attributes = attrs, 
            values = ids, mart = ensembl)
		ind <- match(ids, biot[,1])
		biot <- biot[ind,2]
		gt <- cbind(gt, biot)
		colnames(gt)[ncol(gt)] <- "BIOTYPE"
	}
    return(gt)
}

.geneTable <- function(im, org, fcs=NULL, grn=NULL)#, context="")
{
    # load org pkg
    org.pkg <- .org2pkg(org)
    isAvailable(org.pkg)
    org.pkg <- get(org.pkg)

    # (1) gene identification 
    EZ.COL <- config.ebrowser("EZ.COL")
    gn.cols <- sapply(c("SYM.COL", "GN.COL"), config.ebrowser)
    gt <- select(org.pkg, keys=colnames(im), columns=gn.cols, keytype=EZ.COL)

    # (2) fold change
    if(!is.null(fcs))
    {
        fcs[,1] <- round(fcs[,1], digits=2)
        fcs[,2] <- signif(fcs[,2], digits=2)
        gt <- cbind(gt, fcs) 
    }

    # (3) interactions
    if(!is.null(grn))
    {
        ias.per.gene <- sapply(colnames(im), 
            function(gene) grn[grn[,1] == gene | grn[,2] == gene,,drop=FALSE])
        nr.ias.per.gene <- sapply(ias.per.gene, nrow)
        gt <- cbind(gt, nr.ias.per.gene)
        colnames(gt)[ncol(gt)] <- "INTACTS"
    }

    # (4) nr.sets
#    gene.occ.freq <- colSums(im)
#    gt <- cbind(gt, gene.occ.freq)
#    colnames(gt)[ncol(gt)] <- "SETS"
#
#    # (5) pubmed
#    pmids <- mapIds(org.pkg, keys=colnames(im), 
#        column="PMID", keytype=EZ.COL, multiVals="list")
#
#    # context search ?
#    # if(context != "") pmids <- search.abstracts(pmids, context=context)
#  
#    nr.articles <- sapply(pmids, length)
#    pubmed.urls <- sapply(pmids, function(p) paste0(
#        config.ebrowser("PUBMED.URL"), paste(p, collapse=",")), USE.NAMES=FALSE)
#    articles <- paste(nr.articles, pubmed.urls)
#
#    gt <- cbind(gt, articles)
#    colnames(gt)[ncol(gt)] <- config.ebrowser("PMID.COL")
    return(gt)
}


.ncbiGeneLink <- function(object, ...)
{
    EZ.COL <- config.ebrowser("EZ.COL")
    col <- as.character(object[,EZ.COL])
    link <- paste0(config.ebrowser("GENE.URL"), col)
    object[,EZ.COL] <- hwriter::hwrite(col, link=link, table=FALSE)
    return(object)
}

.pubmedLink <- function(object, ...)
{
    PMID.COL <- config.ebrowser("PMID.COL")
    spl <- sapply(as.character(object[,PMID.COL]), 
        function(s) unlist(strsplit(s, " ")))
    spl <- t(spl)
    spl1 <- as.integer(spl[,1])
    object[,PMID.COL] <- hwriter::hwrite(spl1, link=spl[,2], table=FALSE)
    return(object)    
}
