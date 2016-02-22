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


# CONSTANTS: 
# TODO resolve:
# options: 
# 1) par-like style --> current solution
#
# 2) methods/functions for getting/setting, e.g. group(obj), adjp(obj), ... 
#
# 3) convert them all to function parameters 

.ebrowser_config_cache <- new.env(parent=emptyenv())

config.ebrowser <- function(key, value=NULL) 
{
    .key_readonly <- c(
        "PRB.COL", "EZ.COL", "GN.COL", "SYM.COL", "PMID.COL", 
        "NCBI.URL", "PUBMED.URL", "GENE.URL", "KEGG.URL", "KEGG.GENE.URL",
        "KEGG.SHOW.URL", "GO.SHOW.URL")
 
    if(is.null(value)) .ebrowser_config_cache[[key]]
    else if(!(key %in% .key_readonly)) .ebrowser_config_cache[[key]] <- value
}
    
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
    setwd(out.dir)
}


##
##
# eBrowser functionality
##
ebrowser <- function(
    meth, exprs, pdat, fdat, org, data.type=c(NA, "ma", "rseq"),
    norm.method="quantile", de.method="limma",
    gs, grn=NULL, perm=1000, alpha=0.05, beta=1, 
    comb=FALSE, browse=TRUE, nr.show=-1)
{
    GRP.COL <- config.ebrowser("GRP.COL")
    FC.COL <- config.ebrowser("FC.COL")
    ADJP.COL <- config.ebrowser("ADJP.COL")
    EZ.COL <- config.ebrowser("EZ.COL")    

    METHODS <- c(sbea.methods(), nbea.methods())
    is.method <- (meth %in% METHODS) | is.function(meth)
    if(!all(is.method)) stop("No such method: ", meth[!is.method])
    
    if(any(nbea.methods() %in% meth)) 
        if(is.null(grn))
            stop("\'grn\' must be not null")
  
    out.dir <- config.ebrowser("OUTDIR.DEFAULT")
    set.out.dir(out.dir)
    
    # execution
    # read expression data
    data.type <- match.arg(data.type)
    if(!is(exprs, "ExpressionSet"))
    {
        message("Read expression data ...")
        eset <- read.eset( exprs.file=exprs, 
            pdat.file=pdat, fdat.file=fdat, data.type=data.type)
    }
    else
    { 
        eset <- exprs
        if(is.na(data.type))
            data.type <- auto.detect.data.type(exprs(eset))
        experimentData(eset)@other$dataType <- data.type
    }
    
    # normalize?
    if(norm.method != 'none')
    {
        message("Normalize ...")
        eset <- normalize(eset, norm.method=norm.method)
    }
    
    # probe 2 gene if data.type = ma
    # ... and it's not already a gene level eset
    if(experimentData(eset)@other$dataType == "ma")
    {
        has.pcol <- config.ebrowser("PRB.COL") %in% colnames(fData(eset))
        has.anno <- ifelse(length(annotation(eset)), nchar(annotation(eset)) > 3, FALSE)
        if(has.pcol || has.anno)
        {
            message("Transform probe expression to gene expression ...")    
            gene.eset <- probe.2.gene.eset(eset)
        }
        else gene.eset <- eset
    }
    else gene.eset <- eset

    message("DE analysis ...")    
    gene.eset <- de.ana(gene.eset, de.method=de.method)
    if(missing(org)) org <- annotation(gene.eset) 
    else annotation(gene.eset) <- org
        
    nr.meth <- length(meth)
    if(comb) res.list <- vector("list", length=nr.meth)
    for(i in seq_len(nr.meth))
    {
        m <- meth[i]
        message(paste("Execute", toupper(m), "..."))
        out.file <- paste0(m, ".txt")

        if(m %in% nbea.methods()) 
            res <- nbea( method=m, eset=gene.eset, gs=gs, 
                    grn=grn, alpha=alpha, beta=beta, perm=perm )

        else res <- sbea( method=m, eset=gene.eset, 
                    gs=gs, alpha=alpha, perm=perm )

        write.table(res$res.tbl, file=out.file, 
            quote=FALSE, row.names=FALSE, sep="\t")

        # produce html reports, if desired
        if(browse) ea.browse( res, nr.show, graph.view=grn, html.only=TRUE )
        
        # link gene statistics
        if(m == "samgs" && file.exists("samt.RData"))
            fData(gene.eset)$SAM.T <- 
                round(get(load("samt.RData")), digits=2)
       

        if(m == "gsea" && file.exists("gsea_s2n.RData"))
            fData(gene.eset)$GSEA.S2N <- 
                round(get(load("gsea_s2n.RData")), digits=2)
        
        if(comb) res.list[[i]] <- res
    }

    # write genewise differential expression
    gene.diffexp.file <- "de.txt"
    ind <- which(colnames(fData(gene.eset)) == FC.COL)
    ind <- ind:ncol(fData(gene.eset))
    gene.diffexp <- apply(as.matrix(fData(gene.eset)[,ind]), 2, as.numeric) 
    ord <- order(gene.diffexp[,ADJP.COL])
    gene.diffexp <- signif(gene.diffexp[ord,], digits=2)
    gene.diffexp <- cbind(featureNames(gene.eset)[ord], gene.diffexp)
    colnames(gene.diffexp)[1] <- EZ.COL

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
            out.file <- "comb.txt"
            res <- comb.ea.results(res.list=res.list)
            write.table(res$res.tbl, 
                file=out.file, quote=FALSE, row.names=FALSE, sep="\t")

            # produce html report, if desired
            if(browse)
                ea.browse(res, nr.show, graph.view=grn, html.only=TRUE)
        }
    }
    
    message(paste("Your output files are in", out.dir, "!"))
    
    # create INDEX.html
    if(browse)
    { 
        # plot DE
        if(length(featureNames(gene.eset)) > 1000)
        {
            message("Restricting global view to the 1000 most significant genes ...")
            gene.eset <- gene.eset[order(fData(gene.eset)[,ADJP.COL])[1:1000],]
        }
        gt <- get.gene.symbol.and.name(featureNames(gene.eset), org)
        fData(gene.eset)[,colnames(gt)[2:3]] <- gt[,2:3]
        vs <- view.set(gene.eset, out.prefix=file.path(out.dir,"global"))
        
        message("Produce html report ...")
        create.index(meth, comb)
    }
}

