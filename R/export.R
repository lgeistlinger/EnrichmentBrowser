#' Export results from set- and network-based enrichment analysis 
#' 
#' This function exports results of differential expression (DE) analysis such
#' as enrichplot and GOPlot. 
#' @param res Enrichment analysis results as returned by \code{\link{sbea}}
#' or \code{\link{nbea}}.
#' @param to Character. Downstream package to which export enrichment analysis
#' results to. Defaults to \code{"enrichplot"}, currently the only supported
#' export option.
export <- function(res,
                   to = c("enrichplot", "GOPlot"))
{
    to <- match.arg(to)
    if(to == "enrichplot") .exportToEnrichplot(res)
    else .exportToGOPlot(res)
}

.exportToGOPlot <- function(res) NULL

.exportToEnrichplot <- function(res)
{
    isAvailable("DOSE", type = "software")
    if(res$method == "ora")
    { 
        new("enrichResult",
            result = .res2enrichResult(res),
            pvalueCutoff = res$alpha,
            pAdjustMethod = "BH", 
            qvalueCutoff = 0.2,
            organism = "UNKNOWN", # metadata(res$se)$annotation,
            ontology = "UNKOWN", # .detectGSType(names(res$gs)[1]),
            gene = rownames(res$se)[.isSig(rowData(res$se))],
            keytype = "UNKNOWN", #.detectGeneIdType(names(res$se)[1]),
            universe = rownames(res$se),
            geneSets = res$gs[res$res.tbl$GENE.SET],
            readable = TRUE)
    }
    else if(res$method == "gsea") new("gseaResult")
}

.res2enrichResult <- function(res)
{
    PVAL.COL <- configEBrowser("PVAL.COL")    
    ADJP.COL <- configEBrowser("ADJP.COL")    

    res.tbl <- res$res.tbl 
    spl <- strsplit(res.tbl$GENE.SET, "_") 
    ids <- vapply(spl, `[`, character(1), x = 1)
    .makeDesc <- function(x) paste(x[2:max(2, length(x))], collapse = " ")
    desc <- vapply(spl, .makeDesc, character(1))

    is.sig <- .isSig(rowData(res$se))
    nr.sigs <- sum(is.sig)
    msigs <- rownames(res$se)[is.sig]
    gs <- res$gs[res.tbl$GENE.SET]
    gs <- lapply(gs, function(x) intersect(msigs, x))
    .collapseGenes <- function(s) paste(s, collapse = "/")
    gids <- vapply(gs, .collapseGenes, character(1), USE.NAMES = FALSE)

    if(ADJP.COL %in% colnames(res.tbl)) adjp <- res.tbl[[ADJP.COL]]
    else adjp <- p.adjust(res.tbl[[PVAL.COL]], method = "BH") 
    
    data.frame(ID = ids,
              Description = desc, 
              GeneRatio = paste(res.tbl$NR.SIG.GENES, nr.sigs, sep = "/"),
              BgRatio = paste(res.tbl$NR.GENES, nrow(res$se), sep = "/"), 
              pvalue = res.tbl[[PVAL.COL]],
              p.adjust = adjp,
              qvalue = adjp,
              geneID = gids,
              Count = res.tbl$NR.SIG.GENES) 
}
