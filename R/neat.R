# NEAT, 30 July 2020
#

###
#
# NEAT - main function
#
# @param:   
#   se        ... SummarizedExperiment R object
#   gs        ... Gene sets
#   grn       ... Gene regulatory network
#             (2 or 3 cols. Col 1: Regulator, col 2: Target. If present, col3 (effect) is ignored)
#   alpha     ... Significance level. Defaults to 0.05.
#   beta      ... Significant log2 fold change level. Defaults to 1 (two-fold).
#   sig.stat  ... criterion to determine which genes in se are DE
#   mtc       ... Multiple testing correction for the NEAT p-values. Default is 'fdr',
#             which corresponds to the Benjamini-Hockberg method. To know the shortcuts 
#             for other multiple testing correction methods, see ?p.adjust
#
# @returns: the NEAT enrichment table, reformatted following the EnrichmentBrowser format
#
###   

.neat = function(se, gs, grn, alpha, beta=1, sig.stat=c("p", "fc", "|", "&")) {
  isAvailable('neat', type='software')
  # derive genes for alist and blist
  isig <- .isSig(rowData(se), alpha, beta, sig.stat)
  de.list <- list('DE gene set' = rownames(se)[isig])
  fgs.list <- gs
  # get network inputs
  network <- grn[ , 1:2]
  all.nodes <- unique(as.vector(network))
  # double-check that all sets in fgs.list have at least one gene in network
  fgs.list <- .pathcheck(pathway.list = fgs.list, nodes = all.nodes)
  # compute neat
  res <- neat(alist = de.list, blist = fgs.list, network = network, 
             nettype = 'directed', nodes = all.nodes, mtc.type = 'fdr')
  res <- as.data.frame(res)
  # restructure output of neat
  rownames(res) <- res$B
  res = res[ , -c(1:2)]
  names(res) = c('n_AB', 'E(N_AB|H_0)', 'PVAL.COL', 'P.ADJUSTED')
  res <- as.matrix(res)
  return(res) 
}

.pathcheck <- function(pathway.list, nodes) {
  get.n.genes = function(path) length(which(path %in% nodes))
  n.nodes = sapply(pathway.list, get.n.genes, simplify = T)
  path.ids = which(n.nodes == 0)
  if (length(path.ids) > 0) {
    warn1 = paste('The following pathways were removed because they are not present in the network:',
                  c(names(pathway.list)[path.ids]))
    warning(warn1)
    pathway.list = pathway.list[-path.ids]
  }
  return(pathway.list)
}
