++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
THE ENRICHMENT BROWSER
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Author: Ludwig Geistlinger
Last update: 25 Apr 2018

Contents:	
    1. What's the EnrichmentBrowser?
    2. Usage
    3. Examples
    4. Troubleshooting

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
1. What's the EnrichmentBrowser?
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The EnrichmentBrowser is an R/Bioconductor package for the enrichment analyis of
gene expression data.
See http://bioconductor.org/packages/EnrichmentBrowser for more information.

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
2. Usage 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

CALL: Rscript eBrowserCMD.R

ARGS:   --meth=<meth>   Enrichment method(s)
        --exprs=<exprs> Expression matrix; tab-separated text file
        --cdat=<cdat>   Column data (sample annotation); tab-separated text file
        --rdat=<rdat>   Row data (gene annotation); tab-separated text file
        --org=<org>     Organism under investigation; 
                        3-letter code (e.g. "hsa" for "Homo sapiens")
        --gs=<gs>       Gene sets; tab-separated text file (GMT format)
                        Alternatively, use "GO" or "KEGG" to retrieve 
                        corresponding predefined gene sets

OPTS:   --grn=<grn>     Gene regulatory network; tab separated text file; 
                        alternatively, use "KEGG" to compile a GRN from KEGG
        --dtype=<dtype> Expression data type (Default: NA); 
                        Use "rseq" for RNA-seq read count data
                        or "ma" for microarray intensity measurements 
        --nmeth=<nmeth> Normalization method (Default: none)
        --dmeth=<dmeth> Differential expression method (Default: limma)
        --alpha=<alpha> Statistical significance level (Default: 0.05)
        --beta=<alpha>  Log2 fold change significance level (Default: 1)
        --perm=<perm>   Nr of permutations (Default: 1000)
        --outd=<outd>   Output directory (Default: ./)
        --show=<show>   Nr of gene sets to include in HTML output (Default: 10)
            
FLAGS:  --help, -h      Calls this usage message
        --html          Produce html report in addition to plain text output
        --comb          Combination of results of multiple enrichment methods

See http://www.bioconductor.org/packages/release/bioc/manuals/EnrichmentBrowser/man/EnrichmentBrowser.pdf#page.10
	
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
3. Examples
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Example files are in the 'extdata' folder of the installation directory.

> system.file("extdata", package="EnrichmentBrowser")

# Basic ORA
Rscript eBrowserCMD.R --meth=ora --perm=0 --org=hsa \
    --exprs=extdata/exprs.tab \
    --cdat=extdata/colData.tab \
    --rdat=extdata/rowData.tab \
    --gs=extdata/hsa_kegg_gs.gmt

# Basic ORA with HTML report
Rscript eBrowserCMD.R --meth=ora --perm=0 --org=hsa \
    --exprs=extdata/exprs.tab \
    --cdat=extdata/colData.tab \
    --rdat=extdata/rowData.tab \
    --gs=extdata/hsa_kegg_gs.gmt --html

# Using KEGG genesets
Rscript eBrowserCMD.R --meth=ora --perm=0 --org=hsa \
    --exprs=extdata/exprs.tab \
    --cdat=extdata/colData.tab \
    --rdat=extdata/rowData.tab \
    --gs=KEGG --html

# GGEA
Rscript eBrowserCMD.R --meth=ggea --org=hsa \
    --exprs=extdata/exprs.tab \
    --cdat=extdata/colData.tab \
    --rdat=extdata/rowData.tab \
    --gs=extdata/hsa_kegg_gs.gmt \
    --grn=KEGG

# Using GO gene sets and KEGG GRN
Rscript eBrowserCMD.R --meth=ggea --org=hsa \ 
    --exprs=extdata/exprs.tab \
    --cdat=extdata/colData.tab \ 
    --rdat=extdata/rowData.tab \
    --gs=GO --grn=KEGG --html
 
# Multiple Methods 
Rscript eBrowserCMD.R --meth=ora,gsea --org=hsa \
    --exprs=extdata/exprs.tab \
    --cdat=extdata/colData.tab \
    --rdat=extdata/rowData.tab \
    --gs=extdata/hsa_kegg_gs.gmt

# Multiple Methods - Combined Output
Rscript eBrowserCMD.R --meth=ora,gsea --comb --org=hsa \
    --exprs=extdata/exprs.tab \
    --cdat=extdata/colData.tab \
    --rdat=extdata/rowData.tab \
    --gs=KEGG --html

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
4. Troubleshooting
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Report issues on https://github.com/lgeistlinger/EnrichmentBrowser/issues
Contact: Ludwig.Geistlinger@sph.cuny.edu

