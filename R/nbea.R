############################################################
#
# author: Ludwig Geistlinger
# date: 3 Feb 2011
#
# GGEA - Gene Graph Enrichment Analysis
#
# Update, 09 May 2014: Extension to network-based enrichment
#           analysis 
#
############################################################

#' @rdname nbea
#' @export
nbeaMethods <- function() 
    c("ggea", "spia", "pathnet", "degraph", 
		"ganpa", "cepa", "topologygsa", "netgsa")


#' Network-based enrichment analysis (NBEA)
#' 
#' This is the main function for network-based enrichment analysis.  It
#' implements and wraps existing implementations of several frequently used
#' methods and allows a flexible inspection of resulting gene set rankings.
#' 
#' 'ggea': gene graph enrichment analysis, scores gene sets according to
#' consistency within the given gene regulatory network, i.e. checks activating
#' regulations for positive correlation and repressing regulations for negative
#' correlation of regulator and target gene expression (Geistlinger et al.,
#' 2011). When using 'ggea' it is possible to estimate the statistical
#' significance of the consistency score of each gene set in two different
#' ways: (1) based on sample permutation as described in the original
#' publication (Geistlinger et al., 2011) or (2) using an approximation in the
#' spirit of Bioconductor's npGSEA package that is much faster.
#' 
#' 'spia': signaling pathway impact analysis, combines ORA with the probability
#' that expression changes are propagated across the pathway topology;
#' implemented in Bioconductor's SPIA package (Tarca et al., 2009).
#' 
#' 'pathnet': pathway analysis using network information, applies ORA on
#' combined evidence for the observed signal for gene nodes and the signal
#' implied by connected neighbors in the network; implemented in Bioconductor's
#' PathNet package.
#' 
#' 'degraph': differential expression testing for gene graphs, multivariate
#' testing of differences in mean incorporating underlying graph structure;
#' implemented in Bioconductor's DEGraph package.
#' 
#' 'topologygsa': topology-based gene set analysis, uses Gaussian graphical
#' models to incorporate the dependence structure among genes as implied by
#' pathway topology; implemented in CRAN's topologyGSA package.
#' 
#' 'ganpa': gene association network-based pathway analysis, incorporates
#' network-derived gene weights in the enrichment analysis; implemented in
#' CRAN's GANPA package.
#' 
#' 'cepa': centrality-based pathway enrichment, incorporates network
#' centralities as node weights mapped from differentially expressed genes in
#' pathways; implemented in CRAN's CePa package.
#' 
#' 'netgsa': network-based gene set analysis, incorporates external information
#' about interactions among genes as well as novel interactions learned from
#' data; implemented in CRAN's NetGSA package.
#' 
#' It is also possible to use additional network-based enrichment methods.
#' This requires to implement a function that takes 'se', 'gs', 'grn', 'alpha',
#' and 'perm' as arguments and returns a numeric matrix 'res.tbl' with a
#' mandatory column named 'PVAL' storing the resulting p-value for each gene
#' set in 'gs'. The rows of this matrix must be named accordingly (i.e.
#' rownames(res.tbl) == names(gs)). See examples.
#' 
#' @aliases ggea spia
#' @param method Network-based enrichment analysis method.  Currently, the
#' following network-based enrichment analysis methods are supported:
#' \sQuote{ggea}, \sQuote{spia}, \sQuote{pathnet}, \sQuote{degraph},
#' \sQuote{topologygsa}, \sQuote{ganpa}, \sQuote{cepa}, and \sQuote{netgsa}.
#' Default is 'ggea'.  This can also be the name of a
#' user-defined function implementing network-based enrichment. See Details.
#' @param se Expression dataset.  An object of class
#' \code{\linkS4class{SummarizedExperiment}}.  Mandatory minimal annotations:
#' \itemize{ \item colData column storing binary group assignment (named
#' "GROUP") \item rowData column storing (log2) fold changes of differential
#' expression between sample groups (named "FC") \item rowData column storing
#' adjusted (corrected for multiple testing) p-values of differential
#' expression between sample groups (named "ADJ.PVAL") } Additional optional
#' annotations: \itemize{ \item colData column defining paired samples or
#' sample blocks (named "BLOCK") \item metadata slot named "annotation" giving
#' the organism under investigation in KEGG three letter code (e.g. "hsa" for
#' Homo sapiens) \item metadata slot named "dataType" indicating the expression
#' data type ("ma" for microarray, "rseq" for RNA-seq) }
#' @param gs Gene sets.  Either a list of gene sets (character vectors of gene
#' IDs) or a text file in GMT format storing all gene sets under investigation.
#' @param grn Gene regulatory network.  Either an absolute file path to a
#' tabular file or a character matrix with exactly *THREE* cols; 1st col = IDs
#' of regulating genes; 2nd col = corresponding regulated genes; 3rd col =
#' regulation effect; Use '+' and '-' for activation/inhibition.
#' @param prune.grn Logical.  Should the GRN be pruned? This removes
#' duplicated, self, and reversed edges. Defaults to TRUE.
#' @param alpha Statistical significance level. Defaults to 0.05.
#' @param perm Number of permutations of the expression matrix to estimate the
#' null distribution. Defaults to 1000. If using method=\sQuote{ggea}, it is
#' possible to set 'perm=0' to use a fast approximation of gene set
#' significance to avoid permutation testing. See Details.
#' @param padj.method Method for adjusting nominal gene set p-values to
#' multiple testing.  For available methods see the man page of the stats
#' function \code{\link{p.adjust}}.  Defaults to 'none', i.e. leaves the
#' nominal gene set p-values unadjusted.
#' @param out.file Optional output file the gene set ranking will be written
#' to.
#' @param browse Logical. Should results be displayed in the browser for
#' interactive exploration? Defaults to FALSE.
#' @param ...  Additional arguments passed to individual nbeaMethods.  This
#' includes currently: \itemize{ \item beta: Log2 fold change significance
#' level. Defaults to 1 (2-fold).  } For SPIA and NEA: \itemize{ \item
#' sig.stat: decides which statistic is used for determining significant DE
#' genes.  Options are: \itemize{ \item 'p' (Default): genes with p-value below
#' alpha.  \item 'fc': genes with abs(log2(fold change)) above beta \item '&':
#' p & fc (logical AND) \item '|': p | fc (logical OR) } } For GGEA: \itemize{
#' \item cons.thresh: edge consistency threshold between -1 and 1.  Defaults to
#' 0.2, i.e. only edges of the GRN with consistency >= 0.2 are included in the
#' analysis. Evaluation on real datasets has shown that this works best to
#' distinguish relevant gene sets.  Use consistency of -1 to include all edges.
#' \item gs.edges: decides which edges of the grn are considered for a gene set
#' under investigation. Should be one out of c('&', '|'), denoting logical AND
#' and OR. respectively. Accordingly, this either includes edges for which
#' regulator AND / OR target gene are members of the investigated gene set.  }
#' @return nbeaMethods: a character vector of currently supported methods;
#' 
#' nbea: if(is.null(out.file)): an enrichment analysis result object that can
#' be detailedly explored by calling \code{\link{eaBrowse}} and from which a
#' flat gene set ranking can be extracted by calling \code{\link{gsRanking}}.
#' If 'out.file' is given, the ranking is written to the specified file.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso Input: \code{\link{readSE}}, \code{\link{probe2gene}},
#' \code{\link{getGenesets}} to retrieve gene set definitions from databases 
#' such as GO and KEGG.
#' \code{\link{compileGRN}} to construct a GRN from pathway databases.
#' 
#' Output: \code{\link{gsRanking}} to rank the list of gene sets.
#' \code{\link{eaBrowse}} for exploration of resulting gene sets.
#' 
#' Other: \code{\link{sbea}} to perform set-based enrichment analysis.
#' \code{\link{combResults}} to combine results from different methods.
#' @references Geistlinger et al. (2011) From sets to graphs: towards a
#' realistic enrichment analysis of transcriptomic systems.  Bioinformatics,
#' 27(13), i366-73.
#'
#' Tarca et al. (2009) A novel signaling pathway impact analysis. 
#' Bioinformatics, 25(1):75-82.
#'
#' @examples
#' 
#'     # currently supported methods
#'     nbeaMethods()
#' 
#'     # (1) expression data: 
#'     # simulated expression values of 100 genes
#'     # in two sample groups of 6 samples each
#'     se <- makeExampleData(what="SE")
#'     se <- deAna(se)
#' 
#'     # (2) gene sets:
#'     # draw 10 gene sets with 15-25 genes
#'     gs <- makeExampleData(what="gs", gnames=names(se))
#' 
#'     # (3) make 2 artificially enriched sets:
#'     sig.genes <- names(se)[rowData(se)$ADJ.PVAL < 0.1]
#'     gs[[1]] <- sample(sig.genes, length(gs[[1]])) 
#'     gs[[2]] <- sample(sig.genes, length(gs[[2]]))   
#'    
#'     # (4) gene regulatory network 
#'     grn <- makeExampleData(what="grn", nodes=names(se))
#'     
#'     # (5) performing the enrichment analysis
#'     ea.res <- nbea(method="ggea", se=se, gs=gs, grn=grn)
#' 
#'     # (6) result visualization and exploration
#'     gsRanking(ea.res, signif.only=FALSE)
#' 
#'     # using your own tailored function as enrichment method
#'     dummyNBEA <- function(se, gs, grn, alpha, perm)
#'     {
#'         sig.ps <- sample(seq(0,0.05, length=1000),5)
#'         insig.ps <- sample(seq(0.1,1, length=1000), length(gs)-5)
#'         ps <- sample(c(sig.ps, insig.ps), length(gs))
#'         score <- sample(1:100, length(gs), replace=TRUE)
#'         res.tbl <- cbind(score, ps)
#'         colnames(res.tbl) <- c("SCORE", "PVAL")
#'         rownames(res.tbl) <- names(gs)
#'         return(res.tbl[order(ps),])
#'     }
#' 
#'     ea.res2 <- nbea(method=dummyNBEA, se=se, gs=gs, grn=grn)
#'     gsRanking(ea.res2) 
#' 
#' @export nbea
nbea <- function(
    method=EnrichmentBrowser::nbeaMethods(), 
    se, 
    gs, 
    grn,
    prune.grn=TRUE,
    alpha=0.05, 
    perm=1000, 
    padj.method="none",
    out.file=NULL,
    browse=FALSE, ...)
{
    # get configuration
    GS.MIN.SIZE <- configEBrowser("GS.MIN.SIZE")
    GS.MAX.SIZE <- configEBrowser("GS.MAX.SIZE")
    PVAL.COL <- configEBrowser("PVAL.COL")
    FC.COL <-  configEBrowser("FC.COL")
    ADJP.COL <-  configEBrowser("ADJP.COL")

    if(is(se, "ExpressionSet")) se <- as(se, "SummarizedExperiment")
    
    # TODO: disentangle DE and EA analysis
    if(!(FC.COL %in% colnames(rowData(se))))
        stop(paste("Required rowData column", FC.COL, "not found"))   
    if(!(ADJP.COL %in% colnames(rowData(se))))
        stop(paste("Required rowData column", ADJP.COL, "not found"))   

    # dealing with NA's
    se <- se[!is.na(rowData(se)[,FC.COL]), ]
    se <- se[!is.na(rowData(se)[,ADJP.COL]), ]

    # getting gene sets & grn
    if(!is.list(gs)) gs <- getGenesets(gs)
    if(!is.matrix(grn)) grn <- .readGRN(grn)

    # prune grn
    if(prune.grn) grn <- .pruneGRN(grn)

    # restrict to relevant genes 
    # in the intersection of se, gs, and grn
    gs.genes <- unique(unlist(gs))
    grn.genes <- unique(c(grn[,1], grn[,2]))
    se.genes <- rownames(se)
    rel.genes <- intersect(intersect(gs.genes, grn.genes), se.genes)
    se <- se[rel.genes,]
    gs <- lapply(gs, function(s) s[s%in% rel.genes])
    lens <- lengths(gs)
    gs <- gs[lens >= GS.MIN.SIZE & lens <= GS.MAX.SIZE]
    grn <- grn[grn[,1] %in% rel.genes & grn[,2] %in% rel.genes,] 
    
    # execute ea
    if(class(method) == "character")
    {
        method <- match.arg(method)
        if(method == "spia") res.tbl <- .spia(se, gs, grn, alpha, perm, ...)
        else if(method == "pathnet") res.tbl <- .pathnet(se, gs, grn, alpha, perm)
        else if(method == "netgsa") res.tbl <- .netgsa(se, gs, grn)
        else if(method == "ganpa") res.tbl <- .ganpa(se, gs, grn, perm)
        else if(method == "cepa") res.tbl <- .cepa(se, gs, grn)
        else if(method == "degraph") res.tbl <- .degraph(se, gs, grn)
        else if(method == "topologygsa") res.tbl <- .topogsa(se, gs, grn, alpha, perm, ...)
        else res.tbl <- .ggea(se, gs, grn, alpha, perm=perm, ...)      
    }
    else if(class(method) == "function") 
        res.tbl <- method(se=se, gs=gs, grn=grn, alpha=alpha, perm=perm, ...)
    else stop(paste(method, "is not a valid method for nbea"))

    res.tbl <- data.frame(signif(res.tbl, digits=3))
    sorting.df <- res.tbl[,ncol(res.tbl)]
    if(ncol(res.tbl) > 1)
        sorting.df <- cbind(sorting.df, -res.tbl[,rev(seq_len(ncol(res.tbl)-1))])
    else colnames(res.tbl)[1] <- PVAL.COL
    res.tbl <- res.tbl[do.call(order, as.data.frame(sorting.df)), , drop=FALSE]
    
    res.tbl[,PVAL.COL] <- p.adjust(res.tbl[,PVAL.COL], method=padj.method)

    res.tbl <- DataFrame(rownames(res.tbl), res.tbl)
    colnames(res.tbl)[1] <- configEBrowser("GS.COL") 
    rownames(res.tbl) <- NULL
    
    if(!is.null(out.file))
    {
        write.table(res.tbl, 
            file=out.file, quote=FALSE, row.names=FALSE, sep="\t")
        message(paste("Gene set ranking written to", out.file)) 
    }
    
    res <- list(
            method=method, res.tbl=res.tbl, 
            nr.sigs=sum(res.tbl[,PVAL.COL] < alpha),
            se=se, gs=gs, alpha=alpha)

    if(browse) eaBrowse(res)
    return(res)
}

#
# general helpers
#
.pruneGRN <- function(grn)
{
    # remove duplicates
    grn <- grn[!duplicated(grn[,1:2]),]

    # rm self edges
    grn <- grn[grn[,1] != grn[,2],]
   
    # rm rev edges 
    genes <- unique(as.vector(grn[,1:2]))
    ggrid <- seq_along(genes)
    names(ggrid) <- genes
    igrn <- cbind(ggrid[grn[,1]], ggrid[grn[,2]])

    n <- nrow(grn)
    grid <- seq_len(n-1)
    ind <- vapply(grid,
        function(i)
        {
            x <- igrn[i,2:1]
            j <- i + 1
            cigrn <- igrn[j:n,,drop=FALSE]
            cigrn <- cigrn[cigrn[,1] == x[1], , drop=FALSE]
            is.rev <- any( cigrn[,2] == x[2] )
            return(is.rev)
        }, logical(1))
    ind <- c(ind, FALSE)
    grn <- grn[!ind,]
}
    
#
# 1 SPIA
#
.spia <- function(se, gs, grn, 
    alpha=0.05, perm=1000, beta=1, sig.stat=c("xxP", "xxFC", "|", "&")) 
{
    FC.COL <- configEBrowser("FC.COL")
    ADJP.COL <- configEBrowser("ADJP.COL")
    PVAL.COL <- configEBrowser("PVAL.COL")

    de.genes <- .isSig(rowData(se, use.names=TRUE), alpha, beta, sig.stat)
    de <- rowData(se, use.names=TRUE)[de.genes, FC.COL]
    names(de) <- rownames(se)[de.genes]
    all <- rownames(se)

    is.kegg <- .detectGSType(names(gs)[1]) == "KEGG"
    organism <- ""
    data.dir <- NULL
    if(is.kegg) organism <- substring(names(gs)[1],1,3)
    else
    {     
        message("making SPIA data ...")
        path.info <- .makeSPIAData(gs, grn)
        data.dir <- system.file("extdata", package="SPIA")
        save(path.info, file=file.path(data.dir, "SPIA.RData"))
        data.dir <- paste0(data.dir, "/")
    }
    res <- SPIA::spia(de=de, all=all, organism=organism, data.dir=data.dir, nB=perm)
    res[,"Name"] <- gsub(" ", "_", res[,"Name"])
    rownames(res) <- paste(paste0(organism, res[,"ID"]), res[,"Name"], sep="_")
    res <- res[, c("pSize", "NDE", "tA", "Status", "pG")]
    colnames(res) <- c("SIZE", "NDE", "T.ACT", "STATUS", PVAL.COL)
    res[,"STATUS"] <- ifelse(res[,"STATUS"] == "Activated", 1, -1)
    res <- as.matrix(res)
    message("Finished SPIA analysis")
    return(res)
}

# spia helper: create pathway data from gs and grn
.makeSPIAData <- function(gs, grn)
{
    rel <- c("activation", "compound", "binding/association", 
            "expression", "inhibition", "activation_phosphorylation", 
            "phosphorylation", "inhibition_phosphorylation", 
            "inhibition_dephosphorylation", "dissociation", "dephosphorylation", 
            "activation_dephosphorylation", "state change", "activation_indirect effect", 
            "inhibition_ubiquination", "ubiquination", "expression_indirect effect", 
            "inhibition_indirect effect", "repression", "dissociation_phosphorylation", 
            "indirect effect_phosphorylation", "activation_binding/association", 
            "indirect effect", "activation_compound", "activation_ubiquination")

    spia.data <- sapply(names(gs), 
        function(s)
        {   
            x <- gs[[s]]
            len <- length(x) 
            nam <- list(x, x)
            m <- matrix(0, nrow=len, ncol=len, dimnames=nam)
            sgrn <- .queryGRN(gs=x, grn=grn, index=FALSE)
            if(nrow(sgrn) < configEBrowser("GS.MIN.SIZE")) return(NULL)
            act.grn <- sgrn[sgrn[,3] == "+",,drop=FALSE]
            actm2 <- m
            if(nrow(act.grn))
            {
                if(nrow(act.grn) > 1) actm <- .grn2adjm(act.grn)
                else actm <- matrix(1, nrow=1, ncol=1, dimnames=list(act.grn[1,1], act.grn[1,2]))
                actm2[rownames(actm), colnames(actm)] <- actm
            }
            inh.grn <- sgrn[sgrn[,3] == "-",,drop=FALSE]
            inhm2 <- m
            if(nrow(inh.grn))
            {
                if(nrow(inh.grn) > 1) inhm <- .grn2adjm(inh.grn)
                else inhm <-  matrix(1, nrow=1, ncol=1, dimnames=list(inh.grn[1,1], inh.grn[1,2]))
                inhm2[rownames(inhm), colnames(inhm)] <- inhm
            }
            l <- lapply(rel, 
                function(r)
                {
                    if(r == "activation") return(actm2)
                    else if(r == "inhibition") return(inhm2)
                    else return(m)                    
                } 
            )
            names(l) <- rel    
            l$nodes <- x
            l$title <- s
            l$NumberOfReactions <- 0
            return(l)
        }
    )
    spia.data <- spia.data[!sapply(spia.data, is.null)]
    return(spia.data)
}

#
# 2 NEA
#
.nea <- function(se, gs, grn, 
    alpha=0.05, perm=100, beta=1, sig.stat=c("xxP", "xxFC", "|", "&"))
{
    nea <- NULL
   isAvailable("neaGUI", type="software")

    #if(perm > 100) perm <- 100
    isig <- .isSig(rowData(se, use.names=TRUE), alpha, beta, sig.stat)
    ags <- rownames(se)[isig]
    grn <- unique(grn[,1:2])
    gs.genes <- unique(unlist(gs))
    grn <- grn[(grn[,1] %in% gs.genes) & (grn[,2] %in% gs.genes),]
    network <- apply(grn, 1, function(x) paste(x, collapse=" "))
    message("Computing NEA permutations, this may take a few minutes ...")
    res <- nea(ags=ags, fgs=gs, network=network, nperm=perm)
    res <- res$MainResult
    res <- res[, c("Number_of_Genes", 
        "Number_of_AGS_genes", "Number_links", "Z_score", "P_value")]
    res <- res[order(res[,"Z_score"], decreasing=TRUE), ]
    colnames(res) <- sub("AGS", "de", colnames(res))
    colnames(res) <- sub("Number", "Nr", colnames(res))
    colnames(res) <- sub("_of", "", colnames(res))
    colnames(res) <- gsub("_", ".", colnames(res))
    colnames(res) <- toupper(colnames(res))
    res <- as.matrix(res)
    PVAL.COL <- configEBrowser("PVAL.COL")
    res <- res[order(res[,PVAL.COL]),]
    return(res) 
}

#
# 3 Pathnet
#
.pathnet <- function(se, gs, grn, alpha=0.05, perm=1000)
{
    PathNet <- NULL
    isAvailable("PathNet", type="software")
    
    ADJP.COL <- configEBrowser("ADJP.COL")
    PVAL.COL <- configEBrowser("PVAL.COL")

    dir.evid <- -log(rowData(se, use.names=TRUE)[,ADJP.COL], base=10)
    dir.evid <- cbind(as.integer(rownames(se)), dir.evid)
    colnames(dir.evid) <- c("Gene.ID", "Obs")
    adjm <- .grn2adjm(grn, directed=FALSE)
    pwy.dat <- .extrPwyDat(gs, grn)
    
    res <- PathNet(
            #Enrichment_Analysis = TRUE, 
            #Contextual_Analysis = FALSE, 
            DirectEvidence_info = dir.evid, 
            Column_DirectEvidence = 2,
            Adjacency = adjm, 
            pathway = pwy.dat, 
            n_perm = perm, 
            threshold = alpha)#,
            #use_sig_pathways  = FALSE)

    res <- res$enrichment_results[, 
        c("Name", "No_of_Genes", "Sig_Direct", "Sig_Combi", "p_PathNet")]
    rownames(res) <- vapply(as.vector(res[,1]), 
        function(s) grep(unlist(strsplit(s,"_"))[1], names(gs), value=TRUE),
        character(1))
    res <- res[-1]    
    colnames(res) <- c("NR.GENES", "NR.SIG.GENES", "NR.SIG.COMB.GENES", PVAL.COL)
    res <- as.matrix(res)
    return(res)
}

# pathnet helper: extract pathway data from gs and grn
.extrPwyDat <- function(gs, grn)
{
    pwy.dat <- lapply(names(gs), 
        function(n)
        {
            genes <- gs[[n]] 
            sgrn <- .queryGRN(gs=genes, grn=grn, index=FALSE)
            if(nrow(sgrn))
                dat <- cbind(sgrn[,1:2, drop=FALSE], rep(n, nrow(sgrn)))
            else dat <- NULL
        }
    )
    ind <- vapply(pwy.dat, is.null, logical(1))
    pwy.dat <- pwy.dat[!ind]
    nr <- sum( vapply(pwy.dat, nrow, integer(1)) )
    pwy.datm <- matrix("", nrow=nr, ncol=3)
    colnames(pwy.datm) <- c("id1", "id2", "title")
    start <- 1
    for(i in seq_len(length(pwy.dat)))
    {
        end <- start + nrow(pwy.dat[[i]]) - 1
        pwy.datm[start:end,] <- pwy.dat[[i]]
        start <- end + 1
    }
    pwy.dat <- data.frame(id1=as.integer(pwy.datm[,1]), 
        id2=as.integer(pwy.datm[,2]), title=pwy.datm[,3])
    return(pwy.dat)
}

# pathnet helper: converts 3-col grn to adjacency matrix
.grn2adjm <- function(grn, directed=TRUE)
{
    nodes <- sort(unique(as.vector(grn[,1:2])))
    adjm <- sapply(nodes, 
        function(n)
        {
            tgs <- grep(n, grn[,1])
            if(length(tgs))
            {
                tgs <- grn[tgs,2]
                adjv <- as.integer(nodes %in% tgs)
            }
            else adjv <- rep(0, length(nodes))
            return(adjv) 
        }) 
    rownames(adjm) <- nodes
    adjm <- t(adjm)
  
    if(!directed)
        for(i in seq_along(nodes)) 
            for(j in seq_along(nodes)) 
                if(adjm[i,j]) adjm[j,i] <- 1

    return(adjm)
}

#
# 4 NetGSA
#
.netgsa <- function(se, gs, grn)
{
     NetGSA <- covsel <- edgelist2adj <- NULL
    isAvailable("netgsa", type="software")

    x <- assay(se)
    y <- colData(se)[,configEBrowser("GRP.COL")] + 1

    # prepare gene sets
    #cmat <- .gmt2cmat(gs, rownames(x)) 
    #if(nrow(cmat) < nrow(x)) x <- x[rownames(cmat),]
    f <- file()
    sink(file=f)
    cmat <- safe::getCmatrix(gs, as.matrix=TRUE)
    sink()
    close(f)
    x <- x[rownames(cmat),]

    # prepare network
    out.dir <- configEBrowser("OUTDIR.DEFAULT")
    if(!file.exists(out.dir)) dir.create(out.dir, recursive=TRUE)
    out.file <- file.path(out.dir, "grn.txt")
    write.table(grn[,1:2], file=out.file, row.names=FALSE)
    vnames <- sort(unique(as.vector(grn[,1:2])))
    adjm <- edgelist2adj(out.file, vertex.names=vnames)
    file.remove(out.file)
    ind <- intersect(rownames(x), rownames(adjm))
    adjm <- adjm[ind, ind]
    x <- x[rownames(adjm),]
    cmat <- cmat[rownames(adjm),]
    cmat <- cmat[rowSums(cmat) > 2,]
    
    message("Estimating weighted adjacency matrix for GRN (group 0)")
    A1 <- covsel(t(x[,y==1]), one=adjm, lambda=0.2)
    message("Estimating weighted adjacency matrix for GRN (group 1)")
    A2 <- covsel(t(x[,y==2]), one=adjm, lambda=0.2)

    # execute
    message("Executing NetGSA ...")
    message("This may take a while ...")
    res <- NetGSA(A1$wAdj, A2$wAdj, x, y, B=cmat, directed=TRUE)

    res <- cbind(res$teststat, res$p.value)
    colnames(res) <- c("STAT", configEBrowser("PVAL.COL"))
    rownames(res) <- rownames(cmat)
    return(res)
}

#
# 5 GANPA
#
.ganpa <- function(se, gs, grn, perm=1000)
{
    GSE.Test.Main <- NULL
    isAvailable("GANPA", type="software")

    # configure
    GRP.COL <- configEBrowser("GRP.COL")
    SMPL.COL <- configEBrowser("SMPL.COL")
    OUT.DIR <- configEBrowser("OUTDIR.DEFAULT")
    GS.MIN.SIZE <- configEBrowser("GS.MIN.SIZE")
    GS.MAX.SIZE <- configEBrowser("GS.MAX.SIZE")
    PVAL.COL <- configEBrowser("PVAL.COL")
    
    if(!file.exists(OUT.DIR)) dir.create(OUT.DIR, recursive=TRUE)
    out.prefix <- file.path(OUT.DIR, "ganpa")

    # expression data
    has.scol <- SMPL.COL %in% colnames(colData(se))
    if(!has.scol) colData(se)[,SMPL.COL] <- colnames(se)
    sinfo <- colData(se)[,c(SMPL.COL, GRP.COL)]
    colnames(sinfo) <- c("sampleid", "status")
    expr.obj <- list(gExprs=assay(se), sampleinfo=sinfo)
    
    # gene regulatory network
    gnet <- .grn2gnet(grn)    

    # execute
    GSE.Test.Main(gExprs.obj=expr.obj, gsets=gs, gNET=gnet, 
        permN=perm, size.min=GS.MIN.SIZE, size.max=GS.MAX.SIZE,
        msp.correction=FALSE, output.label=out.prefix, permFDR.cutoff=1)

    # read results from output csv
    res <- read.csv(paste(out.prefix, "MeanAbs.OrigW.csv", sep=".")) 
    n <- res[,1]
    res <- as.matrix(res[,c("Size", "S", "NS", "permP")])
    colnames(res)[c(1,4)] <- c("SIZE", PVAL.COL)
    rownames(res) <- n
    return(res)
}

.grn2gnet <- function(grn)
{
    ureg <- unique(grn[,1])
    gnet <- sapply(ureg, function(r) grn[grn[,1] == r,2])
    return(gnet)
}

#
# 6 CePa
#
.cepa <- function(se, gs, grn, perm=1000)
{
    cepa.all <- set.pathway.catalogue <- sampleLabel <- NULL
    isAvailable("CePa", type="software")

    # define sample groups
    GRP.COL <- configEBrowser("GRP.COL")
    sl <- sampleLabel(colData(se)[, GRP.COL], treatment=1, control=0)

    # create pathway catalogue from gs and grn
    # (1) pathway list
    pl <- sapply(gs, function(s) as.character(.queryGRN(s, grn)))
    pl <- pl[sapply(pl, length) >= configEBrowser("GS.MIN.SIZE")]

    # (2) interaction list
    il <- data.frame(as.character(seq_len(nrow(grn))), grn[,1:2], stringsAsFactors=FALSE)
    colnames(il) <- c("interaction.id", "input", "output")

    # (3) mapping
    m <- data.frame(node.id=rownames(se), symbol=rownames(se), stringsAsFactors=FALSE)
    pwy.cat <- set.pathway.catalogue(pathList=pl, interactionList=il, mapping=m)
    
    # executing
    res <- cepa.all(mat=assay(se), label=sl, pc=pwy.cat, iter=perm)
    res.mat <- matrix(0.0, nrow=length(res), ncol=7)
    for(i in seq_along(res))
    {
        res.mat[i,1:6] <- sapply(res[[i]], function(x) x$p.value)
        res.mat[i,7] <- min(6 * min(res.mat[i,1:6]), 1) 

    }
    rownames(res.mat) <- names(res)
    n <- paste(toupper(names(res[[1]])), "PVAL" , sep=".")
    colnames(res.mat) <- c(n, configEBrowser("PVAL.COL"))
    return(res.mat)
}

#
# 7 DEGraph
#
.degraph <- function(se, gs, grn)
{    
    testOneGraph <- NULL
    isAvailable("DEGraph", type="software")

    grp <- colData(se)[,configEBrowser("GRP.COL")]
            
    options(show.error.messages=FALSE) 
    res <- sapply(names(gs),
        function(s)
        {
            gs.grn <- .queryGRN(gs[[s]], grn, index=FALSE)
            if(nrow(gs.grn) < configEBrowser("GS.MIN.SIZE")) return(NA)
            am <- .grn2adjm(gs.grn)
            gr <- graphAM(am, "directed")           
            gr <- as(gr, "graphNEL")                 
            r <- try( testOneGraph(graph=gr, data=assay(se), 
                       classes=grp, useInteractionSigns=FALSE), silent=TRUE )
            if(is(r, "try-error")) return(NA) else return(r[[1]]$p.value[1])
        })
    options(show.error.messages=TRUE)
    names(res) <- names(gs)
    res <- res[!is.na(res)]
    return(res)
}

#
# 8 topologyGSA 
#
.topogsa <- function(se, gs, grn, alpha=0.05, perm=1000, test.mean=TRUE)
{    
    # call topologyGSA via clipper's pathQ function
    return(.clipper(se, gs, grn, alpha, perm))

    # original topologyGSA: deprecated
    # does not terminate on particular gs, eg. hsa04060_Cytokine-cytokine_receptor_interaction
    pathway.mean.test <- pathway.var.test <- NULL
    isAvailable("topologyGSA", type="software")
  
    is.DAG <- NULL
    isAvailable("gRbase", type="software")
    
    graph_from_graphnel <- mst <- NULL
     isAvailable("igraph", type="software")
 
    grp <- colData(se)[,configEBrowser("GRP.COL")]
    y1 <- t(assay(se)[, grp == 0])
    y2 <- t(assay(se)[, grp == 1])

    res <- sapply(names(gs),
        function(s)
        {
            message(s)
            gs.grn <- .queryGRN(gs[[s]], grn, index=FALSE)
            if(nrow(gs.grn) < configEBrowser("GS.MIN.SIZE")) return(NA)
            am <- .grn2adjm(gs.grn)
            gr <- graphAM(am, "directed")           
            gr <- as(gr, "graphNEL")                 
            
            if(!is.DAG(gr))
            {
                gr2 <- graph_from_graphnel(gr)
                gr2 <- mst(gr2)
                gr <- as(gr2, "graphNEL")
            }
            if(test.mean) r <- pathway.mean.test(y1, y2, gr, alpha, perm)
            else r <- pathway.var.test(y1, y2, gr, alpha)
            return(r$p.value)
        }
    )
    names(res) <- names(gs)
    res <- res[!is.na(res)]
    return(res)
}

#
# 9 clipper
#
.clipper <- function(se, gs, grn, alpha=0.05, perm=1000)
{    
    pathQ <- NULL
    isAvailable("clipper", type="software")
  
    is.DAG <- NULL
    isAvailable("gRbase", type="software")
    
    graph_from_graphnel <- mst <- NULL
     isAvailable("igraph", type="software")
 
    grp <- colData(se)[,configEBrowser("GRP.COL")] + 1

    res <- sapply(names(gs),
        function(s)
        {
            message(s)
            gs.grn <- .queryGRN(gs[[s]], grn, index=FALSE)
            if(nrow(gs.grn) < configEBrowser("GS.MIN.SIZE")) return(NA)
            am <- .grn2adjm(gs.grn)
            gr <- graphAM(am, "directed")           
            gr <- as(gr, "graphNEL")                 
            
            if(!is.DAG(gr))
            {
                gr2 <- graph_from_graphnel(gr)
                gr2 <- mst(gr2)
                gr <- as(gr2, "graphNEL")
            }
            r <- try(pathQ(assay(se), grp, gr, perm, alpha), silent=TRUE)
            if(is(r, "try-error")) return(NA) else return(r$alphaMean)
        }
    )
    names(res) <- names(gs)
    res <- res[!is.na(res)]
    return(res)
}

