############################################################
# 
# author: Ludwig Geistlinger
# last update: 2018-04-16 
# 
# descr: differential expression analysis of RNA-seq data
# 
# call: Rscript de_rseq.R   <exprs.file> 
#                           <cdat.file> 
#                           <rdat.file> 
#                           <de.method> 
#                           <out.file>
#
# choose <de.method> out of {'limma', 'edgeR', 'DESeq'}
############################################################

if(length(commandArgs()) != 10) 
{
    message(paste("usage: Rscript de_rseq.R <exprs.file> <cdat.file>", 
                                    "<rdat.file> <de.method> <out.file>"))
    quit(save="no")
}

message("Loading EnrichmentBrowser")
suppressPackageStartupMessages(library(EnrichmentBrowser))

exprs.file <- commandArgs()[6]
cdat.file <- commandArgs()[7]
rdat.file <- commandArgs()[8]
de.method <- commandArgs()[9]
out.file <- commandArgs()[10]

message("Reading data ...")
se <- readSE(exprs.file, cdat.file, rdat.file)

message("DE analysis ...")
se <- deAna(se, de.method=de.method, padj.method="none")

de.tbl <- rowData(se)[,sapply(c("FC.COL","ADJP.COL"), configEBrowser)]
de.tbl <- cbind(de.tbl, p.adjust(de.tbl[,2], method="BH"))
de.tbl <- cbind(names(se), de.tbl)
colnames(de.tbl) <- c("GENE.ID", "log2FC", "RAW.PVAL", "ADJ.PVAL")

write.table(de.tbl, file=out.file, row.names=FALSE, quote=FALSE, sep="\t")
message("DE table written to ", out.file)
