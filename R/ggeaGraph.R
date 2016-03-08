############################################################
#
# author: Ludwig Geistlinger
# date: 8 Feb 2011
#
# visualization utility for GRNs and expression data
#
############################################################

ggea.graph <- function(gs, grn, eset, 
    alpha=0.05, beta=1, max.edges=50, cons.thresh=0.7)
{
    sgrn <- query.grn(gs, grn, index=FALSE)
    g <- construct.ggea.graph(grn=sgrn, eset=eset,
        alpha=alpha, beta=beta, max.edges=max.edges, cons.thresh=cons.thresh)
    plot.ggea.graph(g)
}

##
## plot.ggea.graph
##
## function plots a graph of a gene regulatory network
## that has been previously constructed via construct.ggea.graph
##
plot.ggea.graph <- function(graph, show.scores=FALSE, title="GGEA Graph")
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
## construct.ggea.graph
##
## function construct and renders a graph of a gene regulatory network
##
construct.ggea.graph <- function(grn, eset, 
    alpha=0.05, beta=1, max.edges=50, cons.thresh=0.7)
{
    # consistency
    nodes <- intersect(featureNames(eset), grn[,1:2]) 
    nr.nodes <- length(nodes)
    
    node.grid <- seq_len(nr.nodes)
    names(node.grid) <- nodes

    grn <- transform.grn(grn, node.grid)

    fDat <- as.matrix(fData(eset)[nodes, 
        sapply(c("FC.COL", "ADJP.COL"), config.ebrowser)])
    de <- comp.de(fDat, alpha=alpha, beta=beta)
    grn.de <- cbind(de[grn[,1]], de[grn[,2]])
    if(ncol(grn) > 2) grn.de <- cbind(grn.de, grn[,3])
    edge.cons <- apply(grn.de, MARGIN=1, FUN=is.consistent)
    
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

    # reset to restricted representation
    ind <- unique(as.vector(grn[,1:2]))
    node.grid <- node.grid[ind]
    nodes <- names(node.grid)
    nr.nodes <- length(nodes)
    fDat <- fDat[nodes,]

    grn[,1] <- sapply(grn[,1], function(x) which(node.grid == x))
    grn[,2] <- sapply(grn[,2], function(x) which(node.grid == x)) 
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
    nColor <- apply(fDat, 1, determine.node.color)
    nLwd <- rep(3, nr.nodes)
    names(nLwd) <- nodes
    org <- annotation(eset)
    if(!length(org)) nLabel <- nodes
    else nLabel <- get.kegg.display.name(nodes, org=org)
    names(nLabel) <- nodes

    nodeRenderInfo(gr) <- list(label=nLabel, col=nColor, lwd=nLwd)

    # (b) render edges
    nr.edges <- numEdges(gr)
    edges <- sub("\\|", "~", names(edgeData(gr)))
    
    # colors & lwds, etc
    eLty <- rep("solid", nr.edges)
    names(eLty) <- edges
    eCol <- sapply(edge.cons, determine.edge.color)
    names(eCol) <- edges
    eLabel <- round(edge.cons, digits=1)
    names(eLabel) <- edges
    eLwd <- sapply(edge.cons, function(e) 
        determine.edge.lwd(e, cons.thresh=cons.thresh))
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

get.kegg.display.name <- function(gene.id, org)
{
    entry <- keggList(paste(org, gene.id, sep=":"))
    dnames <- sapply(entry, function(e) unlist(strsplit(e, "[,;] "))[1])
    return(dnames)
}

##
## determine.edge.lwd
##
## function returns line width depending on edge consistency
##
determine.edge.lwd <- function(edge.cons, cons.thresh=0.7) 
    ifelse(abs(edge.cons) > cons.thresh, 
        1 + (abs(edge.cons) - cons.thresh) * 10, 1)

##
## determine.edge.color
##
## function returns a color depending on edge consistency
##
determine.edge.color <- function(edge.cons) 
    ifelse(edge.cons < 0, rgb(0,0,abs(edge.cons)), rgb(abs(edge.cons),0,0))

##
## determine.node.color
##
## function returns a color depending on fc sign and significance
##
determine.node.color <- function(node.info)
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
## ggea.graph.legend
##
## function plots the legend for a ggea.graph
##
ggea.graph.legend <- function() 
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

