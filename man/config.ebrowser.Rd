\name{config.ebrowser}
\alias{config.ebrowser}
\title{
    Configuring the EnrichmentBrowser    
}
\description{
    Function to get and set configuration parameters determining
    the default behavior of the EnrichmentBrowser
}
\usage{
    config.ebrowser( key, value = NULL )
}
\arguments{
    \item{key}{
        Configuration parameter.
    }
    \item{value}{
        Value to overwrite the current value of key.
    }
}   
\details{
    Important pData, fData, and result column names:
    \itemize{
        \item SMPL.COL: pData column storing the sample IDs (default: "SAMPLE")
        \item GRP.COL: pData column storing binary group assignment (default: "GROUP")
        \item BLK.COL: pData column defining paired samples or sample blocks (default: "BLOCK")

        \item PRB.COL: fData column storing probe/feature IDs ("PROBEID", read-only)
        \item EZ.COL: fData column storing gene ENTREZ IDs ("ENTREZID", read-only)
        \item SYM.COL: fData column storing gene symbols ("SYMBOL", read-only)
        \item GN.COL: fData column storing gene names ("GENENAME", read-only)

        \item FC.COL: fData column storing (log2) fold changes of differential 
            expression between sample groups (default: "FC")
        \item ADJP.COL: fData column storing adjusted (corrected for multiple testing) 
            p-values of differential expression between sample groups (default: "ADJ.PVAL")

        \item GS.COL: result table column storing gene set IDs (default: "GENE.SET")
        \item GSP.COL: result table column storing gene set significance (default: "P.VALUE")
        \item PMID.COL: gene table column storing PUBMED IDs ("PUBMED", read-only)
    }
    
    Important URLs (all read-only):
    \itemize{
        \item NCBI.URL: http://www.ncbi.nlm.nih.gov/
        \item PUBMED.URL: http://www.ncbi.nlm.nih.gov/pubmed/
        \item GENE.URL: http://www.ncbi.nlm.nih.gov/gene/
        \item KEGG.URL: http://www.genome.jp/dbget-bin/
        \item KEGG.GENE.URL: http://www.genome.jp/dbget-bin/www_bget?
        \item KEGG.SHOW.URL: http://www.genome.jp/dbget-bin/show_pathway?
        \item GO.SHOW.URL: http://amigo.geneontology.org/amigo/term/
    }

    Default output directory:
    \itemize{
        \item EBROWSER.HOME: system.file(package="EnrichmentBrowser")
        \item OUTDIR.DEFAULT: file.path(EBROWSER.HOME, "results")    
    }

    Gene set size:
    \itemize{
        \item GS.MIN.SIZE: minimum number of genes per gene set (default: 5)
        \item GS.MAX.SIZE: maximum number of genes per gene set (default: 500)
    }

    Result appearance:
    \itemize{
        \item RESULT.TITLE: (default: "Table of Results")
        \item NR.SHOW: maximum number of entries to show (default: 20)
    } 
}      
\value{
    If is.null(value) this returns the value of the selected configuration parameter.
    Otherwise, it updates the selected parameter with the given value.
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@bio.ifi.lmu.de>
}
\examples{
    # getting config information
    config.ebrowser("GS.MIN.SIZE") 

    # setting config information
    # WARNING: this is for advanced users only!
    # inappropriate settings will impair EnrichmentBrowser's functionality
    config.ebrowser(key="GS.MIN.SIZE", value=3)  
}

