\name{make.example.data}
\alias{make.example.data}
\title{
    Example data for the EnrichmentBrowser package
}
\description{
    Functionality to construct example data sets for demonstration.
    This includes expression data, gene sets, gene regulatory networks, 
    and enrichment analysis results.
}
\usage{
    make.example.data( what = c("eset", "gs", "grn", "ea.res"), ... )
}
\arguments{
    \item{what}{
        Kind of example data set to be constructed.
        This should be one out of: 
        \itemize{
            \item eset: Expression set
            \item gs: Gene set list
            \item grn: Gene regulatory network
            \item ea.res: Enrichment analysis result object as returned 
                by the functions \code{\link{sbea}} and \code{\link{nbea}}
        }
    } 
    \item{...}{
        Additional arguments to fine-tune the specific example data sets.
        
        For what='eset':
        \itemize{
        \item type: Expression data type. Should be either 
            'ma' (Microarray intensity measurements) or 
            'rseq' (RNA-seq read counts). 
        \item nfeat: Number of features/genes. Defaults to 100.
        \item nsmpl: Number of samples. Defaults to 12.
        \item blk: Create sample blocks. Defaults to TRUE.
        \item norm: Should the expression data be normalized? Defaults to FALSE.
        \item de.ana: Should an differential expression analysis be carried out 
            automatically? Defaults to FALSE.
        }
        
        For what='gs':
        \itemize{
        \item gnames: gene names from which the sets will be sampled. Per default
            the sets will be drawn from c(g1, ..., g100).
        \item n: number of sets. Defaults to 10.
        \item min.size: minimal set size. Defaults to 15.
        \item max.size: maximal set size. Defaults to 25.
        }

        For what='grn':
        \itemize{
        \item nodes: gene node names for which edges will be drawn.
            Per default node names will be c(g1, ..., g100).
        \item edge.node.ratio: ratio number of edges / number of nodes. 
            Defaults to 3, i.e. creates 3 times more edges than nodes.
        }

        For what='ea.res':
        \itemize{
        \item eset: Expression set. Calls make.example.data(what="eset") per default.
        \item gs: Gene sets. Calls make.example.data(what="gs") per default.
        \item method: Enrichment analysis method. Defaults to 'ora'.
        \item alpha: Statistical significance level. Defaults to 0.05. 
        }
    }
}
\value{
   Depends on the 'what' argument. 
} 
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@bio.ifi.lmu.de>
}
\examples{
    eset <- make.example.data(what="eset")
}

