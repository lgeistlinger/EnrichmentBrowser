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
            result = as.data.frame(res$res.tbl),
            pvalueCutoff = res$alpha,
            pAdjustMethod = "BH", 
            qvalueCutoff = 0.2,
            organism = metadata(se)$annotation,
            ontology = .detectGSType(names(res$gs)[1]),
            gene = NULL, # the DE genes,
            keytype = .detectGeneIdType(names(se)[1]),
            universe = names(se),
            geneSets = res$gs,
            readable = TRUE)
    }
    else if(res$method == "gsea") new("gseaResult")
}


.res2result <- function(res)
{
    res.tbl <- res$res.tbl 
    spl <- lapply(res.tbl$GENE.SET, function(x) unlist(strsplit(x, "_"))) 
    ids <- vapply(spl, function(x) x[1], character(1))
    .makeDesc <- function(x) paste(x[2:length(x)], collapse = " ")
    desc <- vapply(spl, .makeDesc, character(1))

    .collapseGenes <- function(s) paste(s, collapse = "/")
    gids <- vapply(res$gs[res.tbl$GENE.SET], .collapseGenes,
                   character(1), USE.NAMES = FALSE)
    
    data.frame(ID = ids,
              Description = desc, 
              GeneRatio = paste(res.tbl$NR.SIG.GENES,
                                sum(rowData(se)$ADJ.PVAL < res$alpha), 
                                sep = "/"),
              BgRatio = paste(res.tbl$NR.GENES, nrow(res$se), sep = "/"), 
              pvalue = res.tbl$PVAL,
              p.adjust = NA, #res.tbl$ADJ.PVAL,
              qvalue = NA,
              geneID = gids,
              Count = res.tbl$NR.SIG.GENES) 
}


