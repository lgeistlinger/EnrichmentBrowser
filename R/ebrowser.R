######################################################################
#
# author: Ludwig Geistlinger
# date: 25 Feb 2011
#
#
# THE ENRICHMENT BROWSER
#
# update: 05 June 2014, eBrowser2.0 all-in-one wrapper function
#
######################################################################

# CONSTANTS
ADJ.METH <- "BH"
ADJP.COL <- "ADJ.PVAL"
FC.COL <- "FC"
GN.COL <- "GENE"
GRP.COL <- "GROUP"
GRPS <- c("1", "0")
PRB.COL <- "PROBE"
RAWP.COL <- "RAW.PVAL"
SMPL.COL <- "SAMPLE"
BLK.COL <- "BLOCK"
TABLE.OF.RESULTS <- "Table of Results"

CSS <- "eBro.css"
OUTDIR.DEFAULT <- "./"
KEGG.SHOW.URL <- "http://www.genome.jp/dbget-bin/show_pathway?"
VIEW.TAGS <- paste0("<img src=\"file://IMG.DIR/",
            c("report_icon.jpg", "plot_icon.jpg", "browse_icon.gif"),
            "\" width=\"30\" height=\"30\" border=\"0\" alt=\"report\">")
NODE.LWD <- 3
NROW.TOP.TABLE <- 20
NR.SHOW.DEFAULT <- 10

NBEA.METHODS <- c("ggea", "nea",  "spia")
SBEA.METHODS <- c("ora", "safe", "gsea", "samgs")
METHODS <- c(SBEA.METHODS, NBEA.METHODS)
GS.MIN.SIZE <- 5
GS.MAX.SIZE <- 500
ALPHA.DEFAULT <- 0.05
PERM.DEFAULT <- 1000
PERM.BLOCK.LENGTH <- 100

##
# set the output directory
##
set.out.dir <- function(out.dir)
{
    if(!file.exists(dirname(out.dir)))
        stop(paste0("Not a valid output directory path \'",out.dir,"\'"))
    if(!file.exists(out.dir)){
        message(paste("Output directory", out.dir, 
            "does not exist.\nThus, the directory is going to be created."))
        dir.create(out.dir)
    }
}


##
##
# eBrowser functionality
##
ebrowser <- function(meth, exprs, pdat, fdat, gs, grn=NULL, perm=1000, 
    alpha=0.05, beta=1, comb=FALSE, browse=TRUE, nr.show=-1, out.dir=NULL)
{
    sapply(meth, function(m)
        if(!(m %in% METHODS))
            stop(paste("No such method:\'",m,
                    "\'!\nChoose method out of {", 
                    paste(METHODS, collapse=", "), "}.")))
    
    if(any(NBEA.METHODS %in% meth)) 
        if(is.null(grn))
            stop("\'grn\' must be not null")

    if(is.null(out.dir)) 
        out.dir <- 
            file.path(system.file(package="EnrichmentBrowser"), "results")
    set.out.dir(out.dir)
    
    # execution
    message("Read expression data ...")
    eset <- read.eset( exprs.file=exprs, pdat.file=pdat, fdat.file=fdat )
    
    message("Transform probe expression to gene expression ...")    
    gene.eset <- probe.2.gene.eset( probe.eset=eset, 
                    heatm.file=file.path(out.dir,"heatmap.png"),
                    distr.file=file.path(out.dir,"pdistr.png"),
                    volc.file=file.path(out.dir, "volcano.png"))

    nr.meth <- length(meth)
    if(comb) res.list <- vector("list", length=nr.meth)
    for(i in seq_len(nr.meth))
    {
        m <- meth[i]
        message(paste("Execute", toupper(m), "..."))
        out.file <- file.path(out.dir, paste0("eBrowser_", m, "_RESULT.txt"))

        if(m %in% NBEA.METHODS) 
            res <- nbea( method=m, eset=gene.eset, gs=gs, 
                    grn=grn, alpha=alpha, beta=beta, perm=perm )

        else res <- sbea( method=m, eset=gene.eset, 
                    gs=gs, alpha=alpha, perm=perm )

        write.table(res$res.tbl, file=out.file, 
            quote=FALSE, row.names=FALSE, sep="\t")

        # produce html reports, if desired
        if(browse)
        {
            to.html( res=res, nr.show=nr.show,
                graph.view=grn, html.out=sub("txt$", "html", out.file) )
        }
        
        # link gene statistics
        if(m == "samgs")
            fData(gene.eset)$SAM.T <- 
                round(get(load(file.path(out.dir,"samt.RData"))), digits=2)

        if(m == "gsea")
            fData(gene.eset)$GSEA.S2N <- 
                round(get(load(file.path(out.dir,"gsea_s2n.RData"))), digits=2)

        if(comb) res.list[[i]] <- res
    }

    # write genewise differential expression
    gene.diffexp.file <- file.path(out.dir, "gene_diffexp.txt")
    gene.diffexp <- apply(as.matrix(fData(gene.eset)), 2, as.numeric) 
    ord <- order(gene.diffexp[,ADJP.COL])
    gene.diffexp <- signif(gene.diffexp[ord,], digits=2)
    gene.diffexp <- cbind(featureNames(gene.eset)[ord], gene.diffexp)
    colnames(gene.diffexp)[1] <- GN.COL

    write.table(gene.diffexp, 
        file=gene.diffexp.file, row.names=FALSE, quote=FALSE, sep="\t")

    message(paste("Genewise differential expression written to", 
        gene.diffexp.file))    

    if(comb)
    {
        #combine results (average ranks, ztrans p-vals)
        message("Combine results ...")
        if(nr.meth == 1) 
            message(paste("Only one method given, \'comb\' ignored."))
        else
        {
            # combine results in a specified out file
            out.file <- file.path(out.dir, "eBrowser_comb_RESULT.txt")
            res <- comb.ea.results(res.list=res.list)
            write.table(res$res.tbl, 
                file=out.file, quote=FALSE, row.names=FALSE, sep="\t")

            # produce html report, if desired
            if(browse)
                to.html(res=res, 
                    nr.show=nr.show,
                    graph.view=grn,
                    html.out=sub("txt$", "html", out.file))
        }
    }
    
    message(paste("Your output files are in",out.dir,"!"))
    
    # create INDEX.html
    if(browse)
    { 
        message("Produce html report ...")
        create.index(out.dir=out.dir, meth=meth, comb=comb)
        browseURL(file.path(out.dir, "index.html"))
    }
}

