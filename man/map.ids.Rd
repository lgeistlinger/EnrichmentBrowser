\name{map.ids}
\alias{map.ids}
\alias{id.types}
\title{ 
    Mapping between gene ID types for feature names of an expression set 
}
\description{
    Functionality to map feature names of an expression set between common gene ID
    types such as ENSEMBL and ENTREZ.
}
\usage{
    id.types( org )

    map.ids( eset, org=NA, from="ENSEMBL", to="ENTREZID" )
}
\arguments{
    \item{eset}{
        An object of class \code{\linkS4class{SummarizedExperiment}}.
        Expects the names to be of gene ID type given in argument 'from'.
    }
    \item{org}{
        Organism in KEGG three letter code, 
        e.g. \sQuote{hsa} for \sQuote{Homo sapiens}.
        See references.
    }
    \item{from}{
        Gene ID type from which should be mapped.
        Corresponds to the gene ID type of the featureNames of argument 'eset'.
        Defaults to 'ENSEMBL'.
    }
    \item{to}{
        Gene ID type to which should be mapped.
        Corresponds to the gene ID type the featuresNames of argument 'eset' should
        be updated with.
        Defaults to 'ENTREZID'. 
    }
}
\value{
    id.types: character vector listing the available gene ID types for the mapping;

    map.ids: An object of \code{\linkS4class{SummarizedExperiment}}.
}
\details{
    The function 'id.types' lists the valid values which the arguments 'from' 
    and 'to' can take. This corresponds to the names of the available gene ID 
    types for the mapping.
}
\references{
    KEGG Organism code
    \url{http://www.genome.jp/kegg/catalog/org_list.html}
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@sph.cuny.edu>
}
\seealso{
    \code{\linkS4class{SummarizedExperiment}},
    \code{\link{mapIds}}, \code{\link{keytypes}}
}
\examples{
    # create an expression dataset with 3 genes and 3 samples
    eset <- make.example.data("eset", nfeat=3, nsmpl=3)
    names(eset) <- paste0("ENSG00000000", c("003","005", "419"))
    eset <- map.ids(eset, org="hsa")
}

