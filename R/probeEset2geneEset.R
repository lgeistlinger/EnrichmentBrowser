# 
# Author: ludwig geistlinger
# Date: May 27th, 2010
#
# transforms a probe eset to a gene eset
###############################################################################

### FUNCTIONS

# fast conversion of probe 2 gene expression
probe.2.gene.exprs <- function(eset)
{
    if(!(GN.COL %in% colnames(fData(eset)))) eset <- anno.p2g(eset)
    
    # remove probes without gene annotation
    not.na <- !is.na(fData(eset)[, GN.COL])
    if(sum(not.na) < nrow(eset)) eset <- eset[not.na,]

    probe.exprs <- as.matrix(exprs(eset))
    probe2gene <- as.vector(fData(eset)[, GN.COL])
    
    # determine unique genes
    genes <- unique(probe2gene)
    
    # compute gene expression
    gene.grid <- seq_len(length(genes))
    names(gene.grid) <- genes
    gene.int.map <- gene.grid[probe2gene]
    gene.exprs <- t(sapply(gene.grid, 
        function(g)
        {
            curr.probes <- which(gene.int.map == g)
            curr.exprs <- probe.exprs[curr.probes,]
            if(is.matrix(curr.exprs)) 
                curr.exprs <- colMeans(curr.exprs, na.rm=TRUE)
            return(curr.exprs)
        })) 
    colnames(gene.exprs) <- sampleNames(eset)
    rownames(gene.exprs) <- genes

    # create new eset
    gene.eset <- new("ExpressionSet", exprs=gene.exprs)
    pData(gene.eset) <- pData(eset)
    return(gene.eset)
}

anno.p2g <- function(eset) 
{
    anno.pkg <- paste0(annotation(eset), ".db")
    if(!(anno.pkg %in% .packages(all.available=TRUE))) 
        stop(paste("Gene annotation not found!", 
            "Make sure that \'annotation(eset)\' is correct and that you",
            "have you have the corresponding annotation package (.db)",
            "installed.\n",
            "Alternatively, you can add a \'GENE\' column to the fData slot."))
    require(anno.pkg, character.only = TRUE)
    p2g.map <- get(paste0(annotation(eset),"ENTREZID"))
    p2g.map <- get("as.list", envir=as.environment("package:AnnotationDbi"))(p2g.map)
    p2g.map <- p2g.map[featureNames(eset)]
    fData(eset)[,GN.COL] <- unlist(p2g.map)
    return(eset)
}

probe.2.gene.eset <- function(
    probe.eset, 
    value.type=c("log2count", "log2ratio"),
    gene.eset.file=NULL,
    heatm.file=NULL,
    distr.file=NULL,
    volc.file=NULL)
{
    value.type <- value.type[1]

    # load the probe eset ...
    eset <- probe.eset
    if(class(eset) == "character") eset <- get(load(probe.eset)) 

    # ... transform it to a gene eset ...
    gene.eset <- probe.2.gene.exprs(eset)

    # ... perform de analysis ...
    gene.eset <- de.ana(gene.eset, value.type=value.type)
    
    # ... and put out the gene eset
    if(!is.null(gene.eset.file)) save(gene.eset, file=gene.eset.file)

    # plot de
    # (a) pval distribution
    if(!is.null(distr.file)) pdistr(gene.eset, distr.file)
    
    # (b) volcano plot
    if(!is.null(volc.file)) volcano(gene.eset, volc.file)
    
    # (c) heatmap
    if(!is.null(heatm.file)) exprs.heatmap(gene.eset, heatm.file)

    return(gene.eset)
}
