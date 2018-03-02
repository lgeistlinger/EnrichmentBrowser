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
        "KEGG.SHOW.URL", "GO.SHOW.URL", "SBEA.PKGS", "NBEA.PKGS")
 
    if(is.null(value)) .ebrowser_config_cache[[key]]
    else if(!(key %in% .key_readonly)) .ebrowser_config_cache[[key]] <- value
}
    
##
# set the output directory
##
.checkOutDir <- function(out.dir)
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
    if(!file.exists(out.dir)) dir.create(out.dir, recursive=TRUE)    

    # execution
    # read expression data
    data.type <- match.arg(data.type)
    if(is.character(exprs))
    {
        message("Read expression data ...")
        eset <- readSE( assay.file=exprs, 
            cdat.file=pdat, rdat.file=fdat, data.type=data.type )
    }
    else
    { 
        ### TEMPORARY: will be replaced by as(eSet,SummarizedExperiment)
        if(is(eset, "ExpressionSet")) 
            eset <- as(eset, "RangedSummarizedExperiment")
        ###  
        eset <- exprs
        if(is.na(data.type))
            data.type <- .detectDataType(assay(eset))
        metadata(eset)$dataType <- data.type
    }
    
    # normalize?
    if(norm.method != 'none')
    {
        message("Normalize ...")
        eset <- normalize(eset, norm.method=norm.method)
    }
    
    # probe 2 gene if data.type = ma
    # ... and it's not already a gene level eset
    if(metadata(eset)$dataType == "ma")
    {
        has.pcol <- config.ebrowser("PRB.COL") %in% colnames(rowData(eset))
        anno <- metadata(eset)$annotation
        has.anno <- ifelse(length(anno), nchar(anno) > 3, FALSE)
        if(has.pcol || has.anno)
        {
            message("Transform probe expression to gene expression ...")    
            gene.eset <- probe2gene(eset)
        }
        else gene.eset <- eset
    }
    else gene.eset <- eset

    message("DE analysis ...")    
    gene.eset <- de.ana(gene.eset, de.method=de.method)
    if(missing(org)) org <- metadata(gene.eset)$annotation 
    else metadata(gene.eset)$annotation <- org
        
    nr.meth <- length(meth)
	if(length(perm) != nr.meth) perm <- rep(perm[1], nr.meth)
    if(comb) res.list <- vector("list", length=nr.meth)
    for(i in seq_len(nr.meth))
    {
        m <- meth[i]
        message(paste("Execute", toupper(m), "..."))
        out.file <- file.path(out.dir, paste0(m, ".txt"))

        if(m %in% nbea.methods()) 
            res <- nbea( method=m, eset=gene.eset, gs=gs, 
                    grn=grn, alpha=alpha, beta=beta, perm=perm[i] )

        else res <- sbea( method=m, eset=gene.eset, 
                    gs=gs, alpha=alpha, perm=perm[i] )

        write.table(res$res.tbl, file=out.file, 
            quote=FALSE, row.names=FALSE, sep="\t")

        # produce html reports, if desired
        if(browse) ea.browse( res, nr.show, graph.view=grn, html.only=TRUE )
        
        # link gene statistics
        sam.file <- file.path(out.dir, "samt.RData")    
        if(m == "samgs" && file.exists(sam.file))
            rowData(gene.eset)$SAM.T <- 
                round(get(load(sam.file)), digits=2)
       
        s2n.file <- file.path(out.dir, "gsea_s2n.RData")
        if(m == "gsea" && file.exists(s2n.file))
            rowData(gene.eset)$GSEA.S2N <- 
                round(get(load(s2n.file)), digits=2)
        
        if(comb) res.list[[i]] <- res
    }

    # write genewise differential expression
    gt.file <- file.path(out.dir, "de.txt")
	message("Annotating genes ...")
    gt <- .getGeneAnno(names(gene.eset), org)
    gt <- cbind(gt, rowData(gene.eset, use.names=TRUE))
    gt <- .sortGeneTable(gt)
    ind <- gt[,config.ebrowser("EZ.COL")]
    gene.eset <- gene.eset[ind, ]
    rowData(gene.eset) <- gt

    message(paste("Genewise differential expression written to", gt.file))    
    write.table(gt, file=gt.file, row.names=FALSE, quote=FALSE, sep="\t")

    if(comb)
    {
        #combine results (average ranks, ztrans p-vals)
        message("Combine results ...")
        if(nr.meth == 1) 
            message(paste("Only one method given, \'comb\' ignored."))
        else
        {
            # combine results in a specified out file
            out.file <- file.path(out.dir, "comb.txt")
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
        if(length(gene.eset) > 1000)
        {
            message("Restricting global view to the 1000 most significant genes")
            gene.eset <- gene.eset[1:1000,]
        }
        vs <- .viewSet(gene.eset, out.prefix=file.path(out.dir, "global"))
        
        message("Produce html report ...")
        .createIndex(meth, comb)
    }
}

