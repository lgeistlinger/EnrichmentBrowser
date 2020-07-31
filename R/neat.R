# NEAT
# v1: July 30, 2020
# v2: July 31, 2020

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
#   sig.stat  ... criterion to determine which genes in se are DE. Default is 'p'
#   directed  ... Logical. It specifies if grn corresponds to a directed (T) or 
#             an undirected (F) network. Default is T
# @returns: the NEAT enrichment table, reformatted following the EnrichmentBrowser format
#
###   

.neat = function(se, gs, grn, alpha, beta=1, sig.stat=c("p", "fc", "|", "&"), directed = T) {
  isAvailable('neat', type='software')
  # derive genes for alist and blist
  isig <- .isSig(rowData(se), alpha, beta, sig.stat)
  de.list <- list('DE gene set' = rownames(se)[isig])
  fgs.list <- gs
  # get network inputs
  network <- grn[ , 1:2]
  all.nodes <- unique(as.vector(network))
  # compute neat
  nettype <- ifelse(directed, 'directed', 'undirected')
  res <- neat(alist = de.list, blist = fgs.list, network = network, 
              nettype = nettype, nodes = all.nodes, mtc.type = 'fdr')
  res <- as.data.frame(res)
  # restructure output of neat
  rownames(res) <- res$B
  res <- res[ , c('nab', 'expected_nab', 'pvalue')]
  res <- as.matrix(res)
  colnames(res) <- c('n_AB', 'E(N_AB|H_0)', configEBrowser("PVAL.COL"))
  return(res) 
}
