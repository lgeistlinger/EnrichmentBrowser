############################################################
# 
# author: Ludwig Geistlinger
# date: 2015-05-21 11:34:11
# 
# descr: adapting pathview to EnrichmentBrowser needs 
# 
############################################################

pathview2 <- function(
    gene.data, pathway.id, species = "hsa", kegg.dir=".",
    gene.idtype="entrez", gene.annotpkg=NULL, min.nnodes=3, 
    map.null=TRUE, map.symbol=TRUE, node.sum="sum", limit=1, bins=10,
    both.dirs=TRUE, low="green", mid="gray", high="red",
    na.col="transparent", afactor=1, text.width = 15
)
{
    gd.names=names(gene.data)
    ng=length(gene.data)
    nsamp.g=1
    gene.idtype=toupper(gene.idtype)

    bods <- get(data("bods", package="pathview"))
    gene.idtype.list <- get(data("gene.idtype.list", package="pathview"))
    species.data=kegg.species.code(species, na.rm=T, code.only=FALSE)
  
    species=species.data["kegg.code"]
    entrez.gnodes=species.data["entrez.gnodes"]==1
    if(is.na(species.data["ncbi.geneid"])){
        if(!is.na(species.data["kegg.geneid"])){
            msg.fmt="Only native KEGG gene ID is supported for this species,\nmake sure it looks like \"%s\"!"
            msg=sprintf(msg.fmt, species.data["kegg.geneid"])
            message("Note: ", msg)
        } else{
        stop("This species is not annotated in KEGG!")
        }
    }
    if(is.null(gene.annotpkg)) gene.annotpkg=bods[match(species, bods[,3]),1]
    if(length(grep("ENTREZ|KEGG", gene.idtype))<1){
        if(is.na(gene.annotpkg)) stop("No proper gene annotation package available!")
        if(!gene.idtype %in% gene.idtype.list) stop("Wrong input gene ID type!")
        gene.idmap=id2eg(gd.names, category=gene.idtype, pkg.name=gene.annotpkg)
        gene.data=mol.sum(gene.data, gene.idmap)
        gene.idtype="ENTREZ"
    }
    if(gene.idtype=="ENTREZ" & !entrez.gnodes){
        message("Info: Getting gene ID data from KEGG...")
        gene.idmap=keggConv("ncbi-geneid", species)
        message("Info: Done with data retrieval!")
        kegg.ids=gsub(paste(species, ":", sep=""), "", names(gene.idmap))
        ncbi.ids=gsub("ncbi-geneid:", "", gene.idmap)
        gene.idmap=cbind(ncbi.ids, kegg.ids)
        gene.data=mol.sum(gene.data, gene.idmap)
        gene.idtype="KEGG"
    }

    #parse  
    warn.fmt="Parsing %s file failed, please check the file!"
    pathway.name = paste(species, pathway.id, sep = "")
    kfiles=list.files(path=kegg.dir, pattern="[.]xml|[.]png")
    tfiles=paste(pathway.name, c("xml","png"), sep=".")
    if(!all(tfiles %in% kfiles)){
        dstatus=download.kegg(pathway.id = pathway.id, species = species, kegg.dir=kegg.dir)
        if(dstatus=="failed") {
        warn.fmt="Failed to download KEGG xml/png files, %s skipped!"
        warn.msg=sprintf(warn.fmt, pathway.name)
        message("Warning: ", warn.msg)
        return(invisible(0))
        }
    }
    
    xml.file <- file.path(kegg.dir, paste(pathway.name, "xml", sep="."))
    gR1=try(pathview:::parseKGML2Graph2(xml.file, genes=F, expand=FALSE, split.group=FALSE), silent=TRUE)
    node.data=try(node.info(gR1), silent=T)
    if(class(node.data)=="try-error"){
      warn.msg=sprintf(warn.fmt, xml.file)
      message("Warning: ", warn.msg)
      return(invisible(0))
    }

    gene.node.type="gene"
    plot.data.gene=node.map(gene.data, node.data, 
        node.types=gene.node.type, node.sum=node.sum, entrez.gnodes=entrez.gnodes)
    kng=plot.data.gene$kegg.names
    kng.char=gsub("[0-9]", "", unlist(kng))
    if(any(kng.char>"")) entrez.gnodes=FALSE
    if(map.symbol & entrez.gnodes) {
              if(is.na(gene.annotpkg)) {
                warn.fmt="No annotation package for the species %s, gene symbols not mapped!"
                warn.msg=sprintf(warn.fmt, species)
                message("Warning: ", warn.msg)
              } else {
                plot.data.gene$labels <- eg2id(as.character(
                    plot.data.gene$kegg.names), category="SYMBOL", pkg.name=gene.annotpkg)[,2]
                mapped.gnodes <- rownames(plot.data.gene)
                node.data$labels[mapped.gnodes] <- plot.data.gene$labels
              }
    }
    cols.ts.gene <- node.color(plot.data.gene, limit, bins, both.dirs=both.dirs,
        discrete=FALSE, low=low, mid=mid, high=high, na.col=na.col)
           
    #group nodes mapping and merge
    grp.idx=node.data$size>1
    if(sum(grp.idx)>0){
        sub2grp=cbind(unlist(node.data$component[grp.idx], use.names=F), rep(names(grp.idx)[grp.idx], node.data$size[grp.idx]))
        du.idx=duplicated(sub2grp[,1])
        if(sum(du.idx)>0){
            du.rn=sub2grp[,1] %in% sub2grp[du.idx,1]
            message("Warning: reconcile groups sharing member nodes!")
            print(sub2grp[du.rn,])
            du.grps=sub2grp[du.idx,]
            rn=which(du.idx)
            for(r in rn){
                comps=node.data$component[[sub2grp[r,2]]]
                comps=comps[comps!=sub2grp[r,1]]
                node.data$component[[sub2grp[r,2]]]=comps
                node.data$size[sub2grp[r,2]]=node.data$size[sub2grp[r,2]]-1
            }
            sub2grp=sub2grp[!du.idx,]
        }
        rownames(sub2grp)=sub2grp[,1]
        for(gn in names(grp.idx)[grp.idx])
            gR1=combineKEGGnodes(node.data$component[[gn]], gR1, gn)
    } else sub2grp=NULL
    nNames=nodes(gR1)
    nSizes=node.data$size[nNames]

    # GENES ONLY!
    nNames <- nNames[node.data$type[nNames] == "gene"]
    gR1 <- subKEGGgraph(nNames, gR1)

    #unconnected nodes processing
    deg=degree(gR1)
    deg=deg$inDegree+deg$outDegree
    if(sum(deg<1)>0){
        gR2=subKEGGgraph(nNames[deg>0], gR1)
        nNames=nNames[deg>0]
        nSizes=nSizes[deg>0]
        if(!is.null(sub2grp))
            sub.idx=sub2grp[,1] %in% nNames |sub2grp[,2] %in% nNames
        else sub.idx=0
    } else {
        gR2=gR1
        if(!is.null(sub2grp)) sub.idx=rep(T, nrow(sub2grp))
        else sub.idx=0
    }

    if(length(nNames)<2){
        msg=sprintf("%s not rendered, 0 or 1 connected nodes!\nTry \"kegg.native=T\" instead!", pathway.name)
        message("Note: ", msg)
        return(list())
    }
 
     
    #give up the KEGG positions, use graphviz layout
    #general attributes
    attrs=list()
    attrs$graph$rankdir="LR"
    attrs$node <- list(fixedsize=FALSE)

    #node attributes
    ntype=node.data$type[nNames]
    cpd.idx=ntype=="compound"
    map.idx=ntype=="map"
    rect.idx=!(cpd.idx|map.idx)
    nAttr=list()
    nAttr$label=rep('', length(nNames))
    shapes=node.data$shape[nNames]
    if(any(cpd.idx)) shapes[cpd.idx]="ellipse"
    if(any(map.idx)) shapes[map.idx]="plaintext"
    nAttr$shape=shapes
    nAttr$height=.75*17/46*nSizes*afactor
    nAttr$width=rep(.75, length(nNames))*afactor
    
    if(any(cpd.idx)){
        nAttr$height[cpd.idx]=nAttr$height[cpd.idx]*1.5
        nAttr$width[cpd.idx]=nAttr$width[cpd.idx]*1.5
    }
    if(any(map.idx)){
        nAttr$height[map.idx]=nAttr$height[map.idx]*1.5
        nAttr$width[map.idx]=nAttr$width[map.idx]*2
    }
    nAttr<- lapply(nAttr, function(x){ names(x) <- nNames; return(x) })

    na.col=rgb(t(col2rgb(na.col)), maxColorValue=255)
    fillcol=rep(na.col, length(nNames))
    names(fillcol)=nNames

    #edge attributes
    subdisplay <- subtypeDisplay(gR2)
    if(length(subdisplay)<1) eAttrs=list() else{
        KEGGEdgeSubtype <- get(data("KEGGEdgeSubtype", package="pathview"))
        na.rn=apply(subdisplay, 2, function(x) sum(is.na(x))==7)
        if(sum(na.rn)>0) subdisplay[,na.rn]=KEGGEdgeSubtype[KEGGEdgeSubtype[,1]=="others",rownames(subdisplay)]
        eLabel <- subdisplay["label", ]
        eCol <- subdisplay["color", ]
        eTextCol <- subdisplay["fontcolor", ]
        eLty <- subdisplay["style", ]
        eArrowhead <- subdisplay["arrowhead", ]
        if (ncol(subdisplay) == 1) {
            tmp <- colnames(subdisplay)[1]
            names(eLabel) <- names(eCol) <- names(eTextCol) <- tmp
            names(eLty) <- names(eArrowhead) <- tmp
        }
        eAttrs <- list(lty = eLty, col = eCol, textCol = eTextCol, 
                 label = eLabel, arrowhead = eArrowhead)
    }
  
    gR2.layout=gR2
    edgeRenderInfo(gR2.layout)=eAttrs
    layoutType= "dot"
    gR2.layout <- Rgraphviz::layoutGraph(gR2.layout, attrs = attrs, nodeAttrs=nAttr, layoutType=layoutType)
    edgeRenderInfo(gR2.layout)=eAttrs
    nri=nodeRenderInfo(gR2.layout)
    loc=list(x=nri$nodeX, y=nri$nodeY)
    if(sum(rect.idx)>0){
        w.unit=min(nri$lWidth[rect.idx])
        h.unit=min(nri$height[rect.idx])
    }
    cni=nSizes>1
    if(sum(cni)>0){
        xloc=rep(loc[[1]][cni], nSizes[cni])
        sn.y=unlist(sapply(nSizes[cni], function(x) seq(-(x-1)/2, (x-1)/2,1)),use.names =F)
        yloc=rep(loc[[2]][cni], nSizes[cni])+h.unit*sn.y
    } else xloc=yloc=NULL
    xloc.nd=c(xloc,loc[[1]][nSizes==1 & rect.idx])
    yloc.nd=c(yloc,loc[[2]][nSizes==1 & rect.idx])
    labs=node.data$labels
    labs[nNames[map.idx]]=sapply(labs[nNames[map.idx]],wordwrap,width=text.width, break.word=F)
    labs[nNames[cpd.idx]]=sapply(labs[nNames[cpd.idx]],wordwrap,width=text.width, break.word=T)

    cols.ts.gene=cbind(cols.ts.gene)
    nc.gene=max(ncol(cols.ts.gene),0)
    pn.suffix=colnames(cols.ts.gene)

    nidx.gene=which(nNames %in% rownames(cols.ts.gene))
    cidx.gene=match(nNames[nidx.gene], rownames(cols.ts.gene))
    sci.gene=match(sub2grp[sub.idx,1], rownames(cols.ts.gene))
    sci.node=match(sub2grp[sub.idx,1], nNames)

    #initialize node colors, independent of user data
    if(sum(rect.idx)>0){
        cn.col=rep(NA, sum(sub.idx))
        cn.col=fillcol[sci.node]
        names(cn.col)=sub2grp[sub.idx,1]
        rect.col=c(cn.col,fillcol[nSizes==1 & rect.idx])
        rect.col[rect.col==na.col]=NA
        rect.col=matrix(rect.col, ncol=1)
    }
    if(sum(cpd.idx)>0){
        ell.col=fillcol[cpd.idx]
        ell.col[ell.col==na.col]=NA
        ell.col=matrix(ell.col, ncol=1)
        w.e=min(nri$lWidth[cpd.idx])
        h.e=min(nri$height[cpd.idx])
        xloc.e=loc[[1]][cpd.idx]
        yloc.e=loc[[2]][cpd.idx]
    }

    fillcol=rep(na.col, length(nNames))
    names(fillcol)=nNames

    if(!is.null(cols.ts.gene) & sum(rect.idx)>0){
        fillcol[nidx.gene]=cols.ts.gene[cidx.gene,]
        cn.col=matrix(NA, nrow=sum(sub.idx), ncol=nc.gene)
        cn.col[]=cols.ts.gene[sci.gene,]
        rownames(cn.col)=sub2grp[sub.idx,1]
        if(nc.gene>1) rect.col=rbind(cn.col,fillcol[nSizes==1 & rect.idx])
        else rect.col=c(cn.col,fillcol[nSizes==1 & rect.idx])
        rect.col[rect.col==na.col]=NA
    }
  
    lab.names <- labs[nNames[!cpd.idx]]
    names(lab.names) <- nodes(gR2.layout)
    nodeRenderInfo(gR2.layout)$label <- lab.names
    nodeRenderInfo(gR2.layout)$fill <- fillcol
    
    # make edge colors
    eri <- edgeRenderInfo(gR2.layout)
    etype <- ifelse(eri$arrowhead == "normal", 
        1, ifelse(eri$arrowhead == "tee", -1, 0))
    ind <- etype %in% c(1,-1)
    nd <- unique(plot.data.gene[,c("labels","mol.data")])
    mol.data <- as.vector(nd[,"mol.data"])
    mol.data <- sapply(mol.data, function(i) ifelse(i > 1, 1, ifelse(i < -1, -1, i)))
    names(mol.data) <- as.vector(nd[,"labels"])
    nmol.data <- sapply(nNames, function(n) mol.data[lab.names[[n]]], USE.NAMES=FALSE)
    names(nmol.data) <- nNames

    from <- nmol.data[eri$enamesFrom[ind]]
    to <- nmol.data[eri$enamesTo[ind]]
    grn <- cbind(from, to, etype[ind])
    edge.cons <- apply(grn, 1, 
        function(x) if(any(is.na(x))) return(0) else return(is.consistent(x)))
    ecol <- determine.edge.color(edge.cons)
    edgeRenderInfo(gR2.layout)$col[ind] <- ecol 

    rg <- Rgraphviz::renderGraph(gR2.layout)

    nodeRenderInfo(rg)$kegg.ids <- 
        node.data$kegg.names[nodes(gR2.layout)]

    return(invisible(rg))
}

#parseKGML2Graph2 <-function (file, ...)
#{
#    pathway <- parseKGML2(file)
#    gR <- KEGGpathway2Graph2(pathway, ...)
#    return(gR)
#}
#
#parseKGML2 <- function (file)
#{
#    doc <- XML::xmlTreeParse(file, getDTD = FALSE)
#    r <- XML::xmlRoot(doc)
#    childnames <- sapply(XML::xmlChildren(r), XML::xmlName)
#    isEntry <- childnames == "entry"
#    isRelation <- childnames == "relation"
#    isReaction <- childnames == "reaction"
#    kegg.pathwayinfo <- parsePathwayInfo(r)
#    kegg.nodes <- sapply(r[isEntry], parseEntry)
#    kegg.edges <- sapply(r[isRelation], parseRelation)
#    kegg.reactions <- sapply(r[isReaction], parseReaction2)
#    names(kegg.nodes) <- sapply(kegg.nodes, getEntryID)
#    pathway <- new("KEGGPathway", pathwayInfo = kegg.pathwayinfo,
#                   nodes = kegg.nodes, edges = kegg.edges, reactions = kegg.reactions)
#    return(pathway)
#}
#
#parseReaction2 <- function (reaction)
#{
#    attrs <- XML::xmlAttrs(reaction)
#    name <- attrs[["name"]]
#    type <- attrs[["type"]]
#    children <- XML::xmlChildren(reaction)
#    childrenNames <- names(children)
#    substrateIndices <- grep("^substrate$", childrenNames)
#    productIndices <- grep("^product$", childrenNames)
#    substrateName <- substrateAltName <- vector("character",
#                                                length(substrateIndices))
#    productName <- productAltName <- vector("character", length(productIndices))
#    for (i in seq(along = substrateIndices)) {
#      ind <- substrateIndices[i]
#      substrate <- children[[ind]]
#      substrateName[i] <- XML::xmlAttrs(substrate)[["id"]]
#      substrateChildren <- XML::xmlChildren(substrate)
#      if (length(substrateChildren) > 0) {
#        substrateAlt <- substrateChildren$alt
#        substrateAltName[i] <- XML::xmlAttrs(substrateAlt)[["name"]]
#      }
#      else {
#        substrateAlt <- as.character(NA)
#        substrateAltName[i] <- as.character(NA)
#      }
#    }
#    for (i in seq(along = productIndices)) {
#      ind <- productIndices[i]
#      product <- children[[ind]]
#      productName[i] <- XML::xmlAttrs(product)[["id"]]
#      productChildren <- XML::xmlChildren(product)
#      if (length(productChildren) > 0) {
#        productAlt <- productChildren$alt
#        productAltName[i] <- XML::xmlAttrs(productAlt)[["name"]]
#      }
#      else {
#        productAlt <- as.character(NA)
#        productAltName[i] <- as.character(NA)
#      }
#    }
#    new("KEGGReaction", name = name, type = type, substrateName = substrateName,
#        substrateAltName = substrateAltName, productName = productName,
#        productAltName = productAltName)
#}
#
#KEGGpathway2Graph2 <- function (pathway, genesOnly = TRUE, 
#    expandGenes = TRUE, split.group=FALSE, check.reaction=TRUE)
#{
#    stopifnot(is(pathway, "KEGGPathway"))
#    if(split.group) pathway <- splitKEGGgroup(pathway)
#    rdata=(pathway@reactions)
#    if (expandGenes){
#      if(check.reaction & length(rdata)>0) message("Note: ", "Gene nodes not expanded when reactions are converted to edges!")
#      else pathway <- expandKEGGPathway(pathway)
#    }
#    knodes <- nodes(pathway)
#    kedges <- edges(pathway)
#    node.entryIDs <- getEntryID(knodes)
#    edge.entryIDs <- getEntryID(kedges)
#    V <- node.entryIDs
#    edL <- vector("list", length = length(V))
#    names(edL) <- V
#    if (is.null(nrow(edge.entryIDs))) {
#      for (i in seq(along = edL)) {
#        edL[[i]] <- list()
#      }
#    }
#    else {
#      for (i in 1:length(V)) {
#        id <- node.entryIDs[i]
#        hasRelation <- id == edge.entryIDs[, "Entry1ID"]
#        if (!any(hasRelation)) {
#          edL[[i]] <- list(edges = NULL)
#        }
#        else {
#          entry2 <- unname(unique(edge.entryIDs[hasRelation,
#                                                "Entry2ID"]))
#          edL[[i]] <- list(edges = entry2)
#        }
#      }
#    }
#    gR <- new("graphNEL", nodes = V, edgeL = edL, edgemode = "directed")
#
#    if(check.reaction & length(rdata)>0){
#      r2e.res=reaction2edge(pathway, gR)
#      gR=r2e.res[[1]]
#      kedges=r2e.res[[2]]
#      knodes=r2e.res[[3]]
#    }
#
#    names(kedges) <- sapply(kedges, function(x) paste(getEntryID(x),
#                                                      collapse = "~"))
#    env.node <- new.env()
#    env.edge <- new.env()
#    assign("nodes", knodes, envir = env.node)
#    assign("edges", kedges, envir = env.edge)
#    nodeDataDefaults(gR, "KEGGNode") <- env.node
#    edgeDataDefaults(gR, "KEGGEdge") <- env.edge
#    if (genesOnly) {
#      gR <- subGraphByNodeType(gR, "gene")
#    }
#    return(gR)
#}
#

