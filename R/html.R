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

create.index <- function(meth, comb)
{
    out.dir <- config.ebrowser("OUTDIR.DEFAULT")
    indexPage <- HTMLReport(shortName = 'index',
                        title = 'EnrichmentBrowser: Index of Result Files',
                        basePath=out.dir, reportDirectory="reports")
    
    res.files <- list.files(out.dir, pattern="txt|png|html$")
    file.rename(from=res.files, to=file.path("reports", res.files))

    vcol <- "global_sview.html"
    de.plot <- hwriteImage(sub("sview.html", "volc.png", vcol),
        link=vcol, table = FALSE, height=200, width=200, target="_blank")

    publish(hwrite(de.plot, border=0, br=TRUE), indexPage)
    publish(Link("DE Measures for each Gene", "de.txt"), indexPage)
    
    if(comb)
    {
        publish(hwrite("Combined Results", heading=4), indexPage)
        publish(Link("Top Table", "comb.html"), indexPage)
        publish(Link("Full Ranking", "comb.txt"), indexPage) 
    }

    sapply(meth,
        function(m)
        {
            publish(hwrite(paste(toupper(m), "Results"), heading=4), indexPage)
            publish(Link("Top Table", paste0(m, ".html")), indexPage)
            publish(Link("Full Ranking", paste0(m, ".txt")), indexPage) 
        })
    index <- finish(indexPage)
    if (interactive()) browseURL(index)
}


ea.browse <- function(res, nr.show=-1, graph.view=NULL, html.only=FALSE)
{
    method <- ifelse( is(res$method, "character"), res$method, NA)
    eset <- res$eset
    alpha <- res$alpha
    gs <- res$gs

    # create out dir
    out.dir <- config.ebrowser("OUTDIR.DEFAULT")
    if(!file.exists(out.dir)) dir.create(out.dir)
    rep.dir <- file.path(out.dir, "reports")
    
    # how many gene sets to show in the output?
    if(nr.show < 1) nr.show <- res$nr.sigs
    if(nr.show > nrow(res$res.tbl)) nr.show <- nrow(res$res.tbl)

    # add description & nr.genes per gene set
    res <- res$res.tbl[seq_len(nr.show),]

    # expecting transition from gene set lists to collections in the near future
    # gsc <- res$gsc
    gsc <- gs.list.2.gs.coll(gs[res[,1]])
    res[,1] <- sapply(res[,1], function(s) unlist(strsplit(s, "_"))[1])
    
    is.kegg <- is(collectionType(gsc[[1]]), "KEGGCollection")
    is.go <- is(collectionType(gsc[[1]]), "GOCollection")

    gs.title <- sapply(gsc, description)
    nr.genes <- sapply(gsc, function(g) length(geneIds(g)))

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
    # TODO: ensure in sbea and nbea that we are running only
    # on intersecting genes between gs and eset
    im <- incidence(gsc)
    org <- organism(gsc[[1]])
    if(org == "") org <- annotation(eset)
    if(!length(org)) stop("Organism annotation not found!\n", 
        "Organism under study must be annotated via annotation(eset)")

    message("Creating gene report ...")
    eset <- eset[colnames(im),]
    fDat <- fData(eset)[,sapply(c("FC.COL", "ADJP.COL"), config.ebrowser)]
    gt <- suppressMessages(gene.table(im, org, fcs=fDat))
    gn.cols <- sapply(c("SYM.COL", "GN.COL"), config.ebrowser)
    fData(eset)[,gn.cols] <- gt[,gn.cols] 
    gt.reps <- sapply(gsc, function(s) gene.report(s, gt, out.dir)) 
    
    # make gene set page
    # (1) link gene report
    link <- paste0(names(gsc), ".html")
    res[,"NR.GENES"] <- hwrite(res[,"NR.GENES"], link=link, table = FALSE)
   
    # set view: volcano, heatmap 
    message("Creating set view ...")
    out.prefix <- file.path(rep.dir, names(gsc))
    names(out.prefix) <- names(gsc)
    vcol <- sapply(gsc, function(s) 
        view.set(eset[geneIds(s),], out.prefix[setName(s)]))
    vcol <- hwriteImage(sub("sview.html", "volc.png", vcol),
        link=vcol, table = FALSE, height=50, width=50, target="_blank")
    res <- DataFrame(res, vcol)
    colnames(res)[ncol(res)] <- "SET.VIEW" 

    # path view: kegg maps
    if(is.kegg)
    { 
        message("Creating kegg view ...")
        vcol <- sapply(gsc, function(s) 
            view.path(setName(s), eset[geneIds(s),], out.prefix[setName(s)]))
        vcol <- hwriteImage(sub("kview.html", "kpath.png", vcol),
            link=vcol, table = FALSE, height=50, width=50, target="_blank")
        res <- DataFrame(res, vcol)
        colnames(res)[ncol(res)] <- "PATH.VIEW"
    }

    # graph view: ggea graphs
    if(!is.null(graph.view)) 
    {
        message("Creating graph view ...")
        vcol <- sapply(gsc, function(s) view.graph(eset[geneIds(s),], query.grn(
            geneIds(s), graph.view, index=FALSE), alpha, out.prefix[setName(s)]))
        vcol <- hwriteImage(sub("html$", "png", vcol),
            link=vcol, table = FALSE, height=50, width=50, target="_blank")
        res <- DataFrame(res, vcol)
        colnames(res)[ncol(res)] <- "GRAPH.VIEW"
    }
 
    # (2) link KEGG / GO 
    link <- NULL
    GS.COL <- config.ebrowser("GS.COL")
    if(is.kegg) link <- sapply(gsc, function(s) 
        get.html.of.marked.pathway(setName(s), geneIds(s)[fDat[geneIds(s),2] < alpha]))
    else if(is.go) link <- paste0(config.ebrowser("GO.SHOW.URL"), res[,GS.COL])
    if(!is.null(link)) res[,GS.COL] <- 
        hwrite(res[,GS.COL], link=link, table=FALSE)

    htmlRep <- HTMLReport(shortName=method,
        title=paste(toupper(method), config.ebrowser("RESULT.TITLE"), sep=" - "),
        basePath=out.dir, reportDirectory="reports")
    res <- as.data.frame(res)
    publish(res, htmlRep) 
        #colClasses = c(rep("sort-string-robust", 3),
        #    rep("sort-num-robust", ncol(gt)-3 )))
    rep <- finish(htmlRep)
    if(!html.only) if (interactive()) browseURL(rep)
}

make.view <- function(html1, html2, gene.html.pos=c("bottom", "topright"))
{
    gene.html.pos <- match.arg(gene.html.pos)
    s <- unlist(strsplit(html1, "_"))[1]
    
    head <- hmakeTag("head", hmakeTag("title", s))
    html1.tag <- hmakeTag("frame", name="volc", src=html1)
    html2.tag <- hmakeTag("frame", name="hmap", src=html2)
    html3.tag <- hmakeTag("frame", name="gene", scrolling="auto")
    
    if(gene.html.pos == "topright")
    {
        bkp <- html3.tag
        html3.tag <- html2.tag
        html2.tag <- bkp
    }

    f1.tag <- hmakeTag("frameset", 
        paste0(sub("</frame>", "", c(html1.tag, html2.tag)), collapse=""), 
        cols=paste0(config.ebrowser("PLOT.WIDTH") + 30, ",*"), border=0)
    html3.tag <- sub("</frame>", "", html3.tag)
    f2.tag <- hmakeTag("frameset", paste(f1.tag, html3.tag),
        rows=paste0(config.ebrowser("PLOT.HEIGHT") + 30, ",*"), border=0)
    cont <- hmakeTag("html", paste0(head,f2.tag))
    return(cont)
}

view.set <- function(eset, out.prefix)
{
    out.files <- paste(out.prefix, 
        c("sview.html", "volc.png", "hmap.png"), sep="_" )
    
    if(!file.exists(out.files[1]))
    {
        ##
        # 1 PLOT: heatmap & volcano
        ##
        volc.html <- make.volc.html(eset, out.files[2])
        hmap.html <- make.hmap.html(eset, out.files[3])
        cont <- make.view(volc.html, hmap.html) 
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

view.path <- function(s, eset, out.prefix)
{
    org <- substring(s, 1, 3)
    pwy.id <- sub("^[a-z]{3}", "", s)
    fc <- fData(eset)[,config.ebrowser("FC.COL")]
    gn.cols <- sapply(c("SYM.COL", "GN.COL"), config.ebrowser)
    gnam <- apply(fData(eset)[,gn.cols], 1, paste, collapse=": ")
    names(fc) <- names(gnam) <- featureNames(eset)

    out.files <- paste(out.prefix, 
        c("kview.html", "kpath.png", "kgraph.png"), sep="_")
    
    if(!file.exists(out.files[1]) && pwy.id != "01100")
    {
        ##
        # 1 PLOT: kegg.native & kegg.graph
        ##
        kpath.html <- make.kpath.html(fc, pwy.id, org, out.files[2])
        kgraph.html <- make.kgraph.html(fc, gnam, pwy.id, org, out.files[3])
        cont <- make.view(kgraph.html, kpath.html, gene.html.pos="topright") 
        cat(cont, file=out.files[1])
    }
    
    views <- basename(out.files[1])
    return(views)
} 

view.graph <- function(eset, sgrn, alpha, out.prefix)
{
    out.files <- paste0(out.prefix, "_gview.", c("html", "png", "txt"))
    
    if(!file.exists(out.files[1]))
    {
        ##
        # 1 PLOT: ggea.graph
        ##
        ggraph.html <- make.ggraph.html(eset, sgrn, alpha, out.files[2])
        void.html <- "void.html"
        cat(hmakeTag('html'), file=file.path(dirname(out.prefix),void.html))
        cont <- make.view(ggraph.html, void.html)   
        cat(cont, file=out.files[1])
    }
    
    views <- basename(out.files[1])
    return(views)
}

make.hmap.html <- function(eset, img.file)
{
    width <- config.ebrowser("PLOT.WIDTH") 
    height <- config.ebrowser("PLOT.HEIGHT")

    # 1: make the plot
    # (a) heatmap 1: all features vs all samples
    expr <- exprs(eset)
    rownames(expr) <- fData(eset)[,config.ebrowser("SYM.COL")]
    png(img.file, width=width, height=height)
    exprs.heatmap(expr=expr, grp=pData(eset)[,config.ebrowser("GRP.COL")])
    dev.off()
    img.tag <- hwriteImage(basename(img.file))

    # (b) heatmap 2: most signif features
    max.row <- 40
    fc <- abs(fData(eset)[,config.ebrowser("FC.COL")])
    p <- -log(fData(eset)[,config.ebrowser("ADJP.COL")], base=10)
    ind <- (fc >= 1) | (p >= 1)
    eset <- eset[ind,]
    if(nrow(eset) > 1)
    {
        if(nrow(eset) > max.row)
        {
            fc <- fc[ind]
            p <- p[ind]
            score <- sqrt(fc^2 + p^2)
            eset <- eset[order(score, decreasing=TRUE)[seq_len(max.row)],]
            #    # select most variable samples
            #if(ncol(eset) > max.col)
            #{
            #    svar <- apply(exprs(eset), 2, var)
            #    eset <- eset[,order(svar, decreasing=TRUE)[seq_len(max.col)]]
            #}
        }
        expr <- exprs(eset)
        rownames(expr) <- fData(eset)[,config.ebrowser("SYM.COL")]
        img2 <- sub(".png$", "2.png", img.file)
        png(img2, width=width, height=height)
        exprs.heatmap(expr=expr, grp=pData(eset)[,config.ebrowser("GRP.COL")])
        dev.off()
        img.tag <- paste0(img.tag, hwriteImage(basename(img2)))
    }
    
    # 2: make the html
    hmap.html <- sub("png$", "html", img.file)
    cont <- hmakeTag('html', hmakeTag('body', img.tag))
    cat(cont, file=hmap.html)
    return(basename(hmap.html))
}

make.volc.html <- function(eset, img.file)
{    
    width <- config.ebrowser("PLOT.WIDTH") 
    height <- config.ebrowser("PLOT.HEIGHT")

    fc <- fData(eset)[,config.ebrowser("FC.COL")]
    p <- fData(eset)[,config.ebrowser("ADJP.COL")]
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
    titles <- apply(fData(eset)[,gn.cols], 1, paste, collapse=": ") 
    refs <- paste0(config.ebrowser("GENE.URL"), featureNames(eset))
    geneplotter::imageMap(coord, con, 
        list(HREF=refs, TITLE=titles, TARGET=rep("gene", nrow(coord))), 
        basename(img.file))    
    close(con)
    return(basename(volc.html))
}

make.kpath.html <- function(fc, pwy.id, org, img.file)
{
    # 1: make the plot
    # run pathview for getting native overplotted image
    out.dir <- dirname(img.file)
    #suppressWarnings(suppressMessages(
    pathview(gene.data=fc, 
        pathway.id=pwy.id, species=org, kegg.dir=out.dir, out.suffix="kpath")
    # ))
    pv.out <- file.path(getwd(), paste0(org, pwy.id, ".kpath.png"))
    file.rename(from=pv.out, to=img.file)
    
    # 2: make the html
    kpath.html <- sub("png$", "html", img.file)
    cont <- hmakeTag('html', hmakeTag('body', hwriteImage(basename(img.file))))
    cat(cont, file=kpath.html)
    return(basename(kpath.html))
}

make.kgraph.html <- function(fc, gname, pwy.id, org, img.file)
{
    width <- config.ebrowser("PLOT.WIDTH") 
    height <- config.ebrowser("PLOT.HEIGHT")

    # 1: make the plot
    # run pathview2 for getting the kegg graph
    out.dir <- dirname(img.file)
    png(img.file, width=width, height=height)
    par(mai=rep(0,4))
    gr <- 
        suppressWarnings(suppressMessages(
        pathview2(gene.data=fc, 
        pathway.id=pwy.id, species=org, kegg.dir=out.dir)
         ))
    dev.off()   
 
    # 2: make the html
    kgraph.html <- sub("png$", "html", img.file)
    if(is(gr, "graph"))
    {
        nd <- nodeRenderInfo(gr)$kegg.ids
        nam <- sapply(names(nd), function(n) 
            ifelse(nd[[n]][1] %in% names(gname), gname[nd[[n]][1]], nodeRenderInfo(gr)$label[[n]]))
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
    else cat(hmakeTag('html'), file=kgraph.html)
    return(basename(kgraph.html))
}

make.ggraph.html <- function(eset, sgrn, alpha, img.file)
{
    width <- config.ebrowser("PLOT.WIDTH") 
    height <- config.ebrowser("PLOT.HEIGHT")

    sggea.graph <- NULL
    if(nrow(sgrn) > 0)
        sggea.graph <- construct.ggea.graph(grn=sgrn, eset=eset, alpha=alpha)
    
    # txt report
    report.file <- sub("png$", "txt", img.file)
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

    # ggea graph png
    png(img.file, width=width, height=height)
    par(mai=rep(0,4))
    if(!is.null(sggea.graph))
        sggea.graph <- plot.ggea.graph(sggea.graph,
            show.scores=(numEdges(sggea.graph) < config.ebrowser("NR.SHOW")))
    else plot(NA, axes=FALSE, xlim=c(0,1), ylim=c(0,1),
        ylab="", xlab="", main="No edges in network for this set!")
    dev.off()
    
    ggraph.html <- sub("view.png$", "graph.html", img.file)
    if(!is.null(sggea.graph))
    {
        gn.cols <- sapply(c("SYM.COL", "GN.COL"), config.ebrowser)
        gnam <- apply(fData(eset)[,gn.cols], 1, paste, collapse=": ")
        names(gnam) <- featureNames(eset)

        # image map
        nd <- nodes(sggea.graph)
        nam <- gnam[nd]
        con <- file(ggraph.html, open="w")
        refs <- paste0(config.ebrowser("GENE.URL"),  nd)
        biocGraph::imageMap(sggea.graph, con=con,
            tags=list(HREF=refs, TITLE = nam, TARGET = rep("gene", length(nd))),
            imgname=basename(img.file), width=width, height=height)    
        close(con)
    }
    else cat(hmakeTag('html'), file=ggraph.html)
    return(basename(ggraph.html))
}

get.html.of.marked.pathway <- function(pwy, oids)
{
    pwy <- sub("^path:", "", pwy)
    oids <- gsub("[a-z]{3}:", "", oids)
    coids <- paste(oids, collapse="+")
    request <- pwy
    if(nchar(coids)) request <- paste(request, coids, sep="+")
    return(paste0(config.ebrowser("KEGG.SHOW.URL"), request))
}


gene.report <- function(s, gt, out.dir)
{
    htmlRep <- HTMLReport(basePath=out.dir, reportDirectory="reports",
        shortName=setName(s), title=paste(setName(s), "Gene Report", sep=": "))
    publish(gt[geneIds(s),], htmlRep, .modifyDF=list(ncbi.gene.link))#, pubmed.link),
        #colClasses = c(rep("sort-string-robust", 3), rep("sort-num-robust", ncol(gt)-3 )))
    rep <- finish(htmlRep)
    return(rep)
}

get.gene.annotation <- function(ids, org, biotype=TRUE)
{
    # load org pkg
    org.pkg <- .org2pkg(org)
    .isAvailable(org.pkg)
    org.pkg <- get(org.pkg)

    # (1) gene identification 
    EZ.COL <- config.ebrowser("EZ.COL")
    gn.cols <- sapply(c("SYM.COL", "GN.COL"), config.ebrowser)
    gt <- suppressMessages(select(org.pkg, keys=ids,
            columns=gn.cols, keytype=EZ.COL))
	
	if(biotype)
	{
		useMart <- listDatasets <- useDataset <- getBM <- NULL
		.isAvailable("biomaRt", type="software")
		
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

gene.table <- function(im, org, fcs=NULL, grn=NULL)#, context="")
{
    # load org pkg
    org.pkg <- .org2pkg(org)
    .isAvailable(org.pkg)
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


ncbi.gene.link <- function(object, ...)
{
    EZ.COL <- config.ebrowser("EZ.COL")
    col <- as.character(object[,EZ.COL])
    link <- paste0(config.ebrowser("GENE.URL"), col)
    object[,EZ.COL] <- hwrite(col, link=link, table=FALSE)
    return(object)
}

pubmed.link <- function(object, ...)
{
    PMID.COL <- config.ebrowser("PMID.COL")
    spl <- sapply(as.character(object[,PMID.COL]), 
        function(s) unlist(strsplit(s, " ")))
    spl <- t(spl)
    object[,PMID.COL] <- hwrite(as.integer(spl[,1]), link=spl[,2], table=FALSE)
    return(object)    
}




