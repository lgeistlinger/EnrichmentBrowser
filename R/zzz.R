############################################################
# 
# author: Ludwig Geistlinger
# date: 2015-06-23 10:25:54
# 
# descr: 
# 
############################################################

.onLoad <- function(libname, pkgname) 
{
    # (1) important pData, fData, and result column names
    .ebrowser_config_cache[["ADJP.COL"]] <- "ADJ.PVAL"
    .ebrowser_config_cache[["FC.COL"]] <- "FC"
    .ebrowser_config_cache[["PRB.COL"]] <- "PROBEID"
    .ebrowser_config_cache[["EZ.COL"]] <- "ENTREZID"
    .ebrowser_config_cache[["GN.COL"]] <- "GENENAME"
    .ebrowser_config_cache[["SYM.COL"]] <- "SYMBOL"
    .ebrowser_config_cache[["GRP.COL"]] <- "GROUP"
    .ebrowser_config_cache[["SMPL.COL"]] <- "SAMPLE"
    .ebrowser_config_cache[["BLK.COL"]] <- "BLOCK"
    .ebrowser_config_cache[["GS.COL"]] <- "GENE.SET"
    .ebrowser_config_cache[["GSP.COL"]] <- "P.VALUE"
    .ebrowser_config_cache[["PMID.COL"]] <- "PUBMED"

    # (2) URLs
    .ebrowser_config_cache[["NCBI.URL"]] <- "http://www.ncbi.nlm.nih.gov/"
    .ebrowser_config_cache[["PUBMED.URL"]] <-
        paste0(.ebrowser_config_cache[["NCBI.URL"]], "pubmed/")
    .ebrowser_config_cache[["GENE.URL"]] <-
        paste0(.ebrowser_config_cache[["NCBI.URL"]], "gene/")
    .ebrowser_config_cache[["KEGG.URL"]] <- "http://www.genome.jp/dbget-bin/"
    .ebrowser_config_cache[["KEGG.GENE.URL"]] <-
        paste0(.ebrowser_config_cache[["KEGG.URL"]], "www_bget?")
    .ebrowser_config_cache[["KEGG.SHOW.URL"]] <-
        paste0(.ebrowser_config_cache[["KEGG.URL"]], "show_pathway?")
    .ebrowser_config_cache[["GO.SHOW.URL"]] <- 
        "http://amigo.geneontology.org/amigo/term/"
    
    # (3) file paths
    .ebrowser_config_cache[["EBROwSER.HOME"]] <- system.file(package="EnrichmentBrowser")
    .ebrowser_config_cache[["OUTDIR.DEFAULT"]] <- 
        file.path(.ebrowser_config_cache[["EBROwSER.HOME"]], "results")
    
    # (4) methodological defaults 
    .ebrowser_config_cache[["GS.MIN.SIZE"]] <- 5
    .ebrowser_config_cache[["GS.MAX.SIZE"]] <- 500

    # (5) output appearance
    .ebrowser_config_cache[["RESULT.TITLE"]] <- "Table of Results" 
    .ebrowser_config_cache[["NR.SHOW"]] <- 20
    .ebrowser_config_cache[["PLOT.WIDTH"]] <- 500
    .ebrowser_config_cache[["PLOT.HEIGHT"]] <- 500
        
    # (6) misc
}



