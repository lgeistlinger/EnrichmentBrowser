############################################################
#
# author: Ludwig Geistlinger
# date: 8 Feb 2011
#
# visualization utility for GRNs and expression data
#
############################################################


#' GGEA graphs of consistency between regulation and expression
#' 
#' Gene graph enrichment analysis (GGEA) is a network-based enrichment analysis
#' method implemented in the EnrichmentBrowser package.  The idea of GGEA is to
#' evaluate the consistency of known regulatory interactions with the observed
#' gene expression data.  A GGEA graph for a gene set of interest displays the
#' consistency of each interaction in the network that involves a gene set
#' member.  Nodes (genes) are colored according to expression
#' (up-/down-regulated) and edges (interactions) are colored according to
#' consistency, i.e. how well the interaction type (activation/inhibition) is
#' reflected in the correlation of the expression of both interaction partners.
#' 
#' 
#' @aliases ggea.graph ggea.graph.legend
#' @param gs Gene set under investigation.  This should be a character vector
#' of gene IDs.
#' @param grn Gene regulatory network.  Character matrix with exactly *THREE*
#' cols; 1st col = IDs of regulating genes; 2nd col = corresponding regulated
#' genes; 3rd col = regulation effect; Use '+' and '-' for
#' activation/inhibition.
#' @param se Expression data given as an object of class
#' \code{\linkS4class{SummarizedExperiment}}.
#' @param alpha Statistical significance level. Defaults to 0.05.
#' @param beta Log2 fold change significance level. Defaults to 1 (2-fold).
#' @param max.edges Maximum number of edges that should be displayed.  Defaults
#' to 50.
#' @param cons.thresh Consistency threshold.  Graphical parameter that
#' correspondingly increases line width of edges with a consistency above the
#' chosen threshold (defaults to 0.7).
#' @param show.scores Logical. Should consistency scores of the edges be 
#' displayed? Defaults to FALSE.
#' @return None, plots to a graphics device.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{\link{nbea}} to perform network-based enrichment analysis.
#' \code{\link{eaBrowse}} for exploration of resulting gene sets.
#' @examples
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
#'     # (3) compiling artificial regulatory network 
#'     grn <- makeExampleData(what="grn", nodes=names(se))
#' 
#'     # (4) plot consistency graph
#'     ggeaGraph(gs=gs[[1]], grn=grn, se=se)
#' 
#'     # (5) get legend
#'     ggeaGraphLegend()
#' 
#' @export ggeaGraph
ggeaGraph <- function(gs, grn, se, 
    alpha=0.05, beta=1, max.edges=50, cons.thresh=0.7, show.scores=FALSE)
{
    sgrn <- .queryGRN(gs, grn, index=FALSE)
    g <- .constructGGEAGraph(grn=sgrn, se=se, alpha=alpha, beta=beta, 
                                max.edges=max.edges, cons.thresh=cons.thresh)
    .plotGGEAGraph(g, show.scores=show.scores)
}

#' @export
#' @keywords internal
ggea.graph <- function(gs, grn, se, 
    alpha=0.05, beta=1, max.edges=50, cons.thresh=0.7, show.scores=FALSE)
{
    .Deprecated("ggeaGraph")
    ggeaGraph(gs, grn, se, alpha=alpha, beta=beta, 
        max.edges=max.edges, cons.thresh=cons.thresh, show.scores=show.scores)
}

##
## .plotGGEAGraph
##
## function plots a graph of a gene regulatory network
## that has been previously constructed via .constructGGEAGraph
##
.plotGGEAGraph <- function(graph, show.scores=FALSE, title="GGEA Graph")
{
    graphRenderInfo(graph) <- list(main=title)

    if(!show.scores)
    {
        newLabels <- rep("", numEdges(graph))   
        names(newLabels) <- names(edgeRenderInfo(graph)$label)
        edgeRenderInfo(graph)$label <- newLabels 
    }

    # layout, render & plot graph
    ncolor <- nodeRenderInfo(graph)$col
    ecolor <- edgeRenderInfo(graph)$col
    elwd <- edgeRenderInfo(graph)$lwd
    graph <- Rgraphviz::layoutGraph(graph)

    nodeRenderInfo(graph) <- list(col=ncolor)
    edgeRenderInfo(graph) <- list(col=ecolor, lwd=elwd)
    graph <- Rgraphviz::renderGraph(graph)
    return(invisible(graph))
}

##
## .constructGGEAGraph
##
## function construct and renders a graph of a gene regulatory network
##
.constructGGEAGraph <- function(grn, se, 
    alpha=0.05, beta=1, max.edges=50, cons.thresh=0.7)
{
    if(is(se, "ExpressionSet")) se <- as(se, "SummarizedExperiment")

    # consistency
    nodes <- intersect(names(se), grn[,1:2]) 
    nr.nodes <- length(nodes)
    
    node.grid <- seq_len(nr.nodes)
    names(node.grid) <- nodes

    grn <- .transformGRN(grn, node.grid)
	if(nrow(grn) < 2) return(NULL)

    fDat <- as.matrix(rowData(se, use.names=TRUE)[nodes, 
        sapply(c("FC.COL", "ADJP.COL"), configEBrowser)])
    de <- .compDE(fDat, alpha=alpha, beta=beta)
    grn.de <- cbind(de[grn[,1]], de[grn[,2]])
    if(ncol(grn) > 2) grn.de <- cbind(grn.de, grn[,3])
    edge.cons <- apply(grn.de, MARGIN=1, FUN=.isConsistent)
    
    ord <- order(abs(edge.cons), decreasing=TRUE)
    edge.cons <- edge.cons[ord]
    grn <- grn[ord,]

    # restrict to c most consistent and i inconsistent edges,
    # s.t. i + c = nr.edges	
    if(nrow(grn) > max.edges)
    {
        edge.cons <- head(edge.cons, max.edges)
        grn <- head(grn, max.edges)
    }

    # rse to restricted representation
    ind <- unique(as.vector(grn[,1:2]))
    node.grid <- node.grid[ind]
    nodes <- names(node.grid)
    nr.nodes <- length(nodes)
    fDat <- fDat[nodes,]

    grn[,1] <- vapply(grn[,1], 
        function(x) which(node.grid == x),
        integer(1))
    grn[,2] <- vapply(grn[,2], 
        function(x) which(node.grid == x),
        integer(1)) 
    ord <- do.call(order, as.data.frame(grn))
    grn <- grn[ord,]
    edge.cons <- edge.cons[ord]
    node.grid <- seq_len(nr.nodes)
    names(node.grid) <- nodes
    
    # init graph: nodes and edge list
    edgeL <- lapply(node.grid, 
        function(i) list(edges=as.integer(grn[which(grn[,1]==i), 2])))
    names(edgeL) <- nodes

    gr <- new("graphNEL", nodes=nodes, edgeL=edgeL, edgemode="directed")

    # render graph
    # (a) render nodes
    nColor <- apply(fDat, 1, .determineNodeColor)
    nLwd <- rep(3, nr.nodes)
    names(nLwd) <- nodes
    org <- metadata(se)$annotation
    if(!length(org)) nLabel <- nodes
    else nLabel <- .getKEGGDisplayName(nodes, org=org)
    names(nLabel) <- nodes

    nodeRenderInfo(gr) <- list(label=nLabel, col=nColor, lwd=nLwd)

    # (b) render edges
    nr.edges <- numEdges(gr)
    edges <- sub("\\|", "~", names(edgeData(gr)))
    
    # colors & lwds, etc
    eLty <- rep("solid", nr.edges)
    names(eLty) <- edges
    eCol <- vapply(edge.cons, .determineEdgeColor, character(1))
    names(eCol) <- edges
    eLabel <- round(edge.cons, digits=1)
    names(eLabel) <- edges
    eLwd <- vapply(edge.cons, 
        function(e) .determineEdgeLwd(e, cons.thresh=cons.thresh),
        numeric(1))
    names(eLwd) <- edges

    edgeRenderInfo(gr) <-
        list(lty=eLty, col=eCol, textCol=eCol, label=eLabel, lwd=eLwd)

    if(ncol(grn) > 2) 
    {
        eArrowhead <- ifelse(grn[,3] == 1, "open", "tee")
        names(eArrowhead) <- edges
        edgeRenderInfo(gr)$arrowhead=eArrowhead
    }

    return(gr)
}

.getKEGGDisplayName <- function(gene.id, org)
{
    entry <- KEGGREST::keggList(paste(org, gene.id, sep=":"))
    dnames <- vapply(entry, 
        function(e) unlist(strsplit(e, "[,;] "))[1],
        character(1))
    return(dnames)
}

##
## .determineEdgeLwd
##
## function returns line width depending on edge consistency
##
.determineEdgeLwd <- function(edge.cons, cons.thresh=0.7) 
    ifelse(abs(edge.cons) > cons.thresh, 
        1 + (abs(edge.cons) - cons.thresh) * 10, 1)

##
## .determineEdgeColor
##
## function returns a color depending on edge consistency
##
.determineEdgeColor <- function(edge.cons) 
    ifelse(edge.cons < 0, rgb(0,0,abs(edge.cons)), rgb(abs(edge.cons),0,0))

##
## .determineNodeColor
##
## function returns a color depending on fc sign and significance
##
.determineNodeColor <- function(node.info)
{
    color <- rgb(0,0,0)
    if(any(is.na(node.info))) return(color)

    # intensity of color is 1 - P,
    # i.e. the more significant, the more intense
    intens <- 1 - node.info[2]

    # fold change positive, i.e. upregulated -> green
    if(sign(node.info[1]) == 1) color <- rgb(0, intens, 0)
    # fold change negative, i.e. downregulated -> red
    else if(sign(node.info[1]) == -1) color <- rgb(intens, 0 , 0)

    return(color)
}


##
## ggeaGraphLegend
##
## function plots the legend for a ggeaGraph
#
#' @rdname ggeaGraph
#' @export
#
ggeaGraphLegend <- function() 
{
    opar <- par(mar=c(0,0,3,0), mgp=c(0,0,0))
    on.exit(par(opar))

    ## init plot
    plot(1,1, type="n", xlim=c(0.2,2), ylim=c(0,12),
        axes=FALSE, xlab="", ylab="", main="GGEA graph legend")

    text(0.8, 1, "inhibition", pos=2, cex=1.2)
    segments(1, 1, 2, 1, col="black", lty=1)
    text(1.97, 1, "|", pos=4, col="black")
    text(0.8, 2, "activation", pos=2, cex=1.2)
    segments(1,2, 2, 2, col="black", lty=1)
    text(1.95, 2, ">", pos=4, col="black")

    text(0.8, 3, "EDGE TYPES", pos=2, cex=1.5)  

    ## EDGE COLORS
    text(1.5, 4, "(the clearer the color appears, the more significant it is)")

    ## inconsistent edges
    i = 5
    text(0.8, i, "inconsistent (blue)", pos=2, cex=1.2)
    segments(1,i,1.45,i,col="blue", lty=1)
    text(1.4, i, ">", pos=4, col="blue")
    segments(1.5,i,1.95,i, col="blue", lty=1)
    text(1.915, i, "|", pos=4, col="blue")

    ## consistent edges
    i = 6
    text(0.8, i, "consistent (red)", pos=2, cex=1.2)
    segments(1,i,1.45,i, col="red", lty=1)
    text(1.4, i, ">", pos=4, col="red")
    segments(1.5,i,1.95,i, col="red", lty=1)
    text(1.915, i, "|", pos=4, col="red")

    text(0.8, 7, "EDGE COLORS", pos=2, cex=1.5)

    ## NODES
    text(1.5, 8, "(the clearer the color appears, the more significant it is)")

    i = 9
    text(0.8, i, "down-regulated", pos=2, cex=1.2)
    symbols(x=1.5, y=i, circles=0.025, inches=FALSE, add=TRUE, fg="red", lwd=2)

    i = 10
    text(0.8, i, "up-regulated", pos=2, cex=1.2)
    symbols(x=1.5, y=i, circles=0.025, inches=FALSE, add=TRUE, fg="green", lwd=2)

    text(0.8, 11, "NODE COLORS", pos=2, cex=1.5)
}

#' @export
#' @keywords internal
ggea.graph.legend <- function() 
{
    .Deprecated("ggeaGraphLegend")
    ggeaGraphLegend()
}
