# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

# EDIT, 17 Sep 2014, Ludwig Geistlinger
#   adapting for use in the EnrichmentBrowser package
#
# EDIT, 03 Aug 2015, Ludwig Geistlinger
#   adapting for use in SAFE framework (local and global stat)

# gsea signal2noise ratio as local.stat for safe
local.s2n <- function (X.mat, y.vec, ...)
{

    stopifnot(length(unique(y.vec)) == 2)    
    if (!all(sort(unique(y.vec)) == c(0,1))) {
        warning("y.vec is not (0,1), thus Group 1 == ", y.vec[1])
        y.vec <- (y.vec == y.vec[1]) * 1
    }
    return(function(data, vec = y.vec, ...) 
    {

        A <- data + 0.00000001

        ind1 <- which(vec==0)
        n1 <- length(ind1)    
        M1 <- rowMeans(A[,ind1])
        A2 <- A*A    
        S1 <- rowMeans(A2[,ind1])   
        S1 <- S1 - M1*M1    
        S1 <- sqrt(abs((n1/(n1-1)) * S1))   
            
        ind2 <- which(vec==1)
        n2 <- length(ind2)
        M2 <- rowMeans(A[,ind2])
        S2 <- rowMeans(A2[,ind2])   
        S2 <- S2 - M2*M2    
        S2 <- sqrt(abs((n2/(n2-1)) * S2))   
        
        # small sigma "fix" as used in GeneCluster
        S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
        S2 <- ifelse(S2 == 0, 0.2, S2) 
        S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
        S1 <- ifelse(S1 == 0, 0.2, S1) 
        M1 <- M1 - M2
        S1 <- S1 + S2
        s2n <- M1/S1

        return(s2n)
    })
}

# TODO: gsea KS stat as global.stat for safe
global.GSEA <-
function (C.mat, u, args.global)
{
    m2 <- length(u)
    size2 <- (rep(1, m2) %*% C.mat)[1, ]
    if (!args.global$one.sided) {
        return(function(u, C.mat2 = as.matrix(C.mat), m = m2,
            g.vec = size2) {
            G <- rep(1, m) %*% t(g.vec)
            ranked.Cmatrix <- C.mat2[order(-abs(u)), ] * sqrt((m -
                G)/G) - (1 - C.mat2[order(-abs(u)), ]) * sqrt(G/(m -
                G))
            return(apply(apply(ranked.Cmatrix, 2, cumsum), 2,
                max))
        })
    }
    else {
        return(function(u, C.mat2 = as.matrix(C.mat), m = m2,
            g.vec = size2) {
            G <- rep(1, m) %*% t(g.vec)
            ranked.Cmatrix <- C.mat2[order(-u), ] * sqrt((m -
                G)/G) - (1 - C.mat2[order(-u), ]) * sqrt(G/(m -
                G))
            return(apply(apply(ranked.Cmatrix, 2, cumsum), 2,
                max))
        })
    }
}

# G S E A -- Gene Set Enrichment Analysis
# This is a methodology for the analysis of global molecular profiles called 
# Gene Set Enrichment Analysis (GSEA). It determines states (e.g. phenotypes). 
# GSEA operates on all genes from an experiment, rank ordered by the signal to 
# noise ratio and determines whether members of an a priori defined gene set are
# nonrandomly distributed towards the top or bottom of the list and thus may 
# correspond to an important biological process. To assess significance the 
# program uses an empirical permutation procedure to test deviation from random
# that preserves correlations between genes. 
#
# For details see Subramanian et al 2005
#
# Inputs:
#   input.ds: Input gene expression Affymetrix dataset file in RES or GCT format
#   input.cls:  Input class vector (phenotype) file in CLS format 
#   gs.file: Gene set database in GMT format 
#   output.directory: Directory where to store output and results (default: .) 
#   reshuffling.type: Type of permutation reshuffling: 
#       "sample.labels" or "gene.labels" (default: "sample.labels") 
#   nperm: Number of random permutations (default: 1000) 
#   weighted.score.type: Enrichment correlation-based weighting: 
#       0=no weight (KS), 1=standard weigth, 2 = over-weigth (default: 1) 
#   gs.size.threshold.min: Minimum size (in genes) 
#       for database gene sets to be considered (default: 25) 
#   gs.size.threshold.max: Maximum size (in genes) 
#       for database gene sets to be considered (default: 500) 
#   random.seed: Random number generator seed. (default: 123456) 
#
#   Output:
#    The results of the method are stored in the "output.directory" 
#       specified by the user as part of the input parameters. 
#      The results files are:
#    - Two tab-separated global result text files (one for each phenotype). 
#       These files are labeled according to the doc string prefix and the 
#       phenotype name from the CLS file: 
#       <doc.string>.SUMMARY.RESULTS.REPORT.<phenotype>.txt
#    - One set of global plots. They include 
#       a.- gene list correlation profile, 
#       b.- global observed and null densities, 
#       c.- heat map for the entire sorted dataset, and 
#       d.- p-values vs. NES plot. 
#       These plots are in a single JPEG file named 
#       <doc.string>.global.plots.<phenotype>.jpg. When the program is run 
#       interactively these plots appear on a window in the R GUI.
#    - A variable number of tab-separated gene result text files according to 
#       how many sets pass any of the significance thresholds 
#      ("nom.p.val.threshold," "fwer.p.val.threshold," and "fdr.q.val.threshold")
#       and how many are specified in the "topgs" parameter. These files are 
#       named: <doc.string>.<gene set name>.report.txt. 
#   - A variable number of gene set plots (one for each gene set report file). 
#       These plots include a.- Gene set running enrichment
#      "mountain" plot, b.- gene set null distribution and c.- heat map for 
#       genes in the gene set. These plots are stored in a 
#      single JPEG file named <doc.string>.<gene set name>.jpg.
# The format (columns) for the global result files is as follows.
# GS : Gene set name.
# SIZE : Size of the set in genes.
# SOURCE : Set definition or source.
# ES : Enrichment score.
# NES : Normalized (multiplicative rescaling) normalized enrichment score.
# NOM p-val : Nominal p-value (from the null distribution of the gene set).
# FDR q-val: False discovery rate q-values
# FWER p-val: Family wise error rate p-values.
# Tag %: Percent of gene set before running enrichment peak.
# Gene %: Percent of gene list before running enrichment peak.
# Signal : enrichment signal strength.
# FDR (median): FDR q-values from the median of the null distributions.
# glob.p.val: P-value using a global statistic 
#   (number of sets above the set's NES).
# 
# The rows are sorted by the NES values 
#   (from maximum positive or negative NES to minimum)
# 
# The format (columns) for the gene set result files is as follows.
# 
# #: Gene number in the (sorted) gene set
# GENE : gene name. For example the probe accession number, 
#   gene symbol or the gene identifier gin the dataset.
# SYMBOL : gene symbol from the gene annotation file.
# DESC : gene description (title) from the gene annotation file.
# LIST LOC : location of the gene in the sorted gene list.
# S2N : signal to noise ratio (correlation) of the gene in the gene list.
# RES : value of the running enrichment score at the gene location.
# CORE_ENRICHMENT: is this gene is the "core enrichment" section of the list? 
#   Yes or No variable specifying in the gene location is before (positive ES) 
#   or after (negative ES) the running enrichment peak.
# 
# The rows are sorted by the gene location in the gene list.
# The function call to GSEA returns a  two element list containing the two 
#   global result reports as data frames ($report1, $report2).
# 
# results1: Global output report for first phenotype 
# result2:  Global putput report for second phenotype
# Start of GSEA methodology 

 #######ebrowser usage
#GSEA(input.ds = as.data.frame(exprs(eset)),
#        input.cls = cls,
#        gs.db = gs.gmt,
#        output.directory = out.dir,
#        nperm                 = perm,
#        gs.size.threshold.min = GS.MIN.SIZE,
#        gs.size.threshold.max = GS.MAX.SIZE)
################

# Main GSEA Analysis Function that implements the entire methodology

GSEA <- function(
input.ds, 
input.cls, 
gs.db, 
output.directory = "", 
reshuffling.type = "sample.labels", 
nperm = 1000, 
weighted.score.type = 1, 
#gs.size.threshold.min = 25, 
#gs.size.threshold.max = 500, 
random.seed = 123456) 
{
     
    if (.Platform$OS.type == "windows") 
    {
        memory.limit(6000000000)
        memory.limit()
    }
  
    # Read input data matrix
    set.seed(seed=random.seed, kind = NULL)
    adjust.param <- 0.5
    dataset <- input.ds
    gene.labels <- row.names(dataset)
    sample.names <- names(dataset)
    A <- data.matrix(dataset)
    cols <- ncol(A)
    rows <- nrow(A)
  
    # Read input class vector
    CLS <- input.cls
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
    phen1 <- class.phen[1]
    phen2 <- class.phen[2]
  
    # sort samples according to phenotype
    col.index <- order(class.labels, decreasing=FALSE)
    class.labels <- class.labels[col.index]
    sample.names <- sample.names[col.index]
    A <- A[, col.index]
    colnames(A) <- sample.names
  
    # Read input gene set database
#    temp <- readLines(gs.db, warn=FALSE)
#    max.Ng <- length(temp)
#    temp.size.G <- sapply(temp, 
#        function(t) length(unlist(strsplit(t, "\t"))) - 2)

#   max.size.G <- max(temp.size.G)      
#    gs <- matrix(rep("null", max.Ng*max.size.G), nrow=max.Ng, ncol= max.size.G)
#    temp.names <- temp.desc <- vector(length = max.Ng, mode = "character")
#    gs.count <- 1
#    for (i in seq_len(max.Ng)) 
#    {
#        spl <- unlist(strsplit(temp[[i]], "\t"))
#        gene.set.size <- length(spl) - 2
#        gs.line <- noquote(spl)
#        gene.set.name <- gs.line[1] 
#        gene.set.desc <- gs.line[2] 
#        gene.set.tags <- 
#            sapply(seq_len(gene.set.size), function(j) gs.line[j + 2])
#        existing.set <- is.element(gene.set.tags, gene.labels)
#        set.size <- sum(existing.set)
#        if ((set.size >= gs.size.threshold.min) && 
#              (set.size <= gs.size.threshold.max))
#        {
#            temp.size.G[gs.count] <- set.size
#            gs[gs.count,] <- c(gene.set.tags[existing.set], 
#                rep(NA, max.size.G - temp.size.G[gs.count]))
#            temp.names[gs.count] <- gene.set.name
#            temp.desc[gs.count] <- gene.set.desc
#            gs.count <- gs.count + 1
#        }
#    } 

    Ng <- length(gs.db)
    gs.names <- names(gs.db)
    size.G <- sapply(gs.db, length) 
    gs <- matrix(NA, nrow=Ng, ncol=max(size.G))
    for(i in seq_len(Ng)) gs[i, seq_len(size.G[i])] <- gs.db[[i]]

    N <- nrow(A)
    Ns <- ncol(A)
    #all.gene.descs <- all.gene.symbols <-  gene.labels[i]
    
    Obs.indicator <- Obs.RES <- matrix(nrow= Ng, ncol=N)
    Obs.ES <- Obs.arg.ES <- Obs.ES.norm <- vector(length = Ng, mode = "numeric")
  
    # Compute observed and random permutation gene rankings
    obs.s2n <- vector(length=N, mode="numeric")
    signal.strength <- tag.frac <- gene.frac <- 
        coherence.ratio <- vector(length=Ng, mode="numeric")
    obs.phi.norm <- matrix(nrow = Ng, ncol = nperm)
    correl.matrix <- obs.correl.matrix <- order.matrix <- 
        obs.order.matrix <- matrix(nrow = N, ncol = nperm)
  
    nperm.per.call <- 100
    n.groups <- nperm %/% nperm.per.call
    n.rem <- nperm %% nperm.per.call
    n.perms <- c(rep(nperm.per.call, n.groups), n.rem)
    n.ends <- cumsum(n.perms)
    n.starts <- n.ends - n.perms + 1
    n.tot <- ifelse(n.rem == 0, n.groups, n.groups + 1)
    
    for (nk in seq_len(n.tot)) 
    {
        call.nperm <- n.perms[nk]
        message(paste("Permutations:", n.starts[nk], "--", n.ends[nk]))
        O <- GSEA.GeneRanking(A, class.labels, gene.labels, call.nperm)
        order.matrix[,n.starts[nk]:n.ends[nk]] <- O$order.matrix
        obs.order.matrix[,n.starts[nk]:n.ends[nk]] <- O$obs.order.matrix
        correl.matrix[,n.starts[nk]:n.ends[nk]] <- O$s2n.matrix
        obs.correl.matrix[,n.starts[nk]:n.ends[nk]] <- O$obs.s2n.matrix
        rm(O)
    }

    message("Processing ...")
    # using median to assign enrichment scores
    obs.s2n <- apply(obs.correl.matrix, 1, median)    
    names(obs.s2n) <- gene.labels
    save(obs.s2n, file=file.path(output.directory, "gsea_s2n.RData"))  
          
    obs.index <- order(obs.s2n, decreasing=TRUE)            
    obs.s2n <- obs.s2n[obs.index]       
    #obs.gene.labels <- gene.labels[obs.index]       
    #obs.gene.descs <- all.gene.descs[obs.index]       
    #obs.gene.symbols <- all.gene.symbols[obs.index]       
 
    for (r in seq_len(nperm)) 
    {
        correl.matrix[, r] <- correl.matrix[order.matrix[,r], r]
        obs.correl.matrix[, r] <- obs.correl.matrix[obs.order.matrix[,r], r]
    }

    gene.list2 <- obs.index
    for (i in seq_len(Ng)) 
    {
        gene.set <- gs[i,!is.na(gs[i,])]
        gene.set2 <- match(gene.set, gene.labels)
        GSEA.results <- GSEA.EnrichmentScore(
            gene.list=gene.list2, gene.set=gene.set2, 
            weighted.score.type=weighted.score.type, correl.vector = obs.s2n)
        Obs.ES[i] <- GSEA.results$ES
        Obs.arg.ES[i] <- GSEA.results$arg.ES
        Obs.RES[i,] <- GSEA.results$RES
        Obs.indicator[i,] <- GSEA.results$indicator
        if (Obs.ES[i] >= 0) 
        {  
            # compute signal strength
            tag.frac[i] <- sum(Obs.indicator[i,seq_len(Obs.arg.ES[i])])/size.G[i]
            gene.frac[i] <- Obs.arg.ES[i]/N
        } 
        else 
        {
            tag.frac[i] <- sum(Obs.indicator[i, Obs.arg.ES[i]:N])/size.G[i]
            gene.frac[i] <- (N - Obs.arg.ES[i] + 1)/N
        }
        signal.strength[i] <- tag.frac[i] * 
            (1 - gene.frac[i]) * (N / (N - size.G[i]))
    }

    # Compute enrichment for random permutations 
    phi <- phi.norm <- obs.phi <- matrix(nrow = Ng, ncol = nperm)
    if (reshuffling.type == "sample.labels") 
    { 
        # reshuffling phenotype labels
        for (i in seq_len(Ng)) 
        {
            gene.set <- gs[i,!is.na(gs[i,])]
            gene.set2 <- match(gene.set, gene.labels)
            for (r in seq_len(nperm)) 
            {
                gene.list2 <- order.matrix[,r]
                GSEA.results <- GSEA.EnrichmentScore2(
                    gene.list=gene.list2, gene.set=gene.set2, 
                    weighted.score.type=weighted.score.type, 
                    correl.vector=correl.matrix[, r])   
                phi[i, r] <- GSEA.results$ES
            }
            obs.gene.list2 <- obs.order.matrix[,1]
            GSEA.results <- GSEA.EnrichmentScore2(gene.list=obs.gene.list2, 
            gene.set=gene.set2, weighted.score.type=weighted.score.type, 
            correl.vector=obs.correl.matrix[, nperm])
            obs.phi[i, ] <- GSEA.results$ES
        }
    }
    else if (reshuffling.type == "gene.labels") 
    { 
        # reshuffling gene labels
        for (i in seq_len(Ng)) 
        {
            gene.set <- gs[i,!is.na(gs[i,])]
            gene.set2 <- match(gene.set, gene.labels)
            for (r in seq_len(nperm)) 
            {
                reshuffled.gene.labels <- sample(1:rows)
                GSEA.results <- GSEA.EnrichmentScore2(
                    gene.list=reshuffled.gene.labels, gene.set=gene.set2, 
                    weighted.score.type=weighted.score.type, 
                    correl.vector=obs.s2n)   
                phi[i, r] <- GSEA.results$ES
            }
            obs.gene.list2 <- obs.order.matrix[,1]
            GSEA.results <- GSEA.EnrichmentScore2(gene.list=obs.gene.list2, 
            gene.set=gene.set2, weighted.score.type=weighted.score.type, 
            correl.vector=obs.correl.matrix[, nperm])   
            obs.phi[i, ] <- GSEA.results$ES
        }
    }
    # Compute 3 types of p-values
    # Find nominal p-values       
    p.vals <- matrix(0, nrow = Ng, ncol = 2)
    for (i in seq_len(Ng)) 
    {
        ind <- phi[i,] >= 0
        pos.phi <- phi[i, ind]
        neg.phi <- phi[i, !ind] 
        ES.value <- Obs.ES[i]
        p.vals[i, 1] <- signif(ifelse(ES.value >= 0,
                        sum(pos.phi >= ES.value)/length(pos.phi),
                        sum(neg.phi <= ES.value)/length(neg.phi)), digits=5)
        
        # Rescaling normalization for each gene set null
        pos.m <- mean(pos.phi)
        neg.m <- mean(abs(neg.phi))
        pos.phi <- pos.phi/pos.m
        neg.phi <- neg.phi/neg.m
        for (j in seq_len(nperm))
        { 
            phi.norm[i, j] <- 
                phi[i, j] / ifelse(phi[i, j] >= 0, pos.m, neg.m)
            obs.phi.norm[i, j] <-
               obs.phi[i, j] / ifelse(obs.phi[i, j] >= 0, pos.m, neg.m)
        }
        Obs.ES.norm[i] <- Obs.ES[i] / ifelse(Obs.ES[i] >= 0, pos.m, neg.m)
    }
    
    # Compute FWER p-vals
#    max.ES.vals.p <- NULL
#    max.ES.vals.n <- NULL
#    for (j in seq_len(nperm)) 
#    {
#        ind <- phi.norm[,j] >= 0
#        pos.phi <- phi.norm[ind, j]
#        neg.phi <- phi[!ind, j] 
#
#        if (length(pos.phi)) max.ES.vals.p <- c(max.ES.vals.p, max(pos.phi))
#        if (length(neg.phi)) max.ES.vals.n <- c(max.ES.vals.n, min(neg.phi))
#    }
#
#    for (i in seq_len(Ng)) 
#    {
#        ES.value <- Obs.ES.norm[i]
#        p.vals[i, 2] <- signif(ifelse(ES.value >= 0, 
#                            sum(max.ES.vals.p >= ES.value),
#                            sum(max.ES.vals.n <= ES.value)) / 
#                            length(max.ES.vals.p), digits=5)
#    }
#
#    # Compute FDRs 
#    NES <- phi.norm.mean <- obs.phi.norm.mean <- phi.norm.median <- 
#        obs.phi.norm.median <- phi.norm.mean <- obs.phi.mean <- 
#        FDR.mean <- FDR.median <- phi.norm.median.d <- 
#        obs.phi.norm.median.d <- vector(length=Ng, mode="numeric")
#
#    Obs.ES.index <- order(Obs.ES.norm, decreasing=TRUE)
#    Orig.index <- seq(1, Ng)
#    Orig.index <- Orig.index[Obs.ES.index]
#    Orig.index <- order(Orig.index, decreasing=FALSE)
#    Obs.ES.norm.sorted <- Obs.ES.norm[Obs.ES.index]
#    gs.names.sorted <- gs.names[Obs.ES.index]
#
#    NES <- Obs.ES.norm.sorted
#    for (k in seq_len(Ng)) 
#    {
#        ES.value <- NES[k]
#        count.col <- obs.count.col <- vector(length=nperm, mode="numeric")
#        for (i in seq_len(nperm)) 
#        {
#            phi.vec <- phi.norm[,i]
#            obs.phi.vec <- obs.phi.norm[,i]
#            if (ES.value >= 0) 
#            {
#                count.col.norm <- sum(phi.vec >= 0)
#                obs.count.col.norm <- sum(obs.phi.vec >= 0)
#                count.col[i] <- ifelse(count.col.norm > 0, 
#                    sum(phi.vec >= ES.value)/count.col.norm, 0)
#                obs.count.col[i] <- ifelse(obs.count.col.norm > 0, 
#                    sum(obs.phi.vec >= ES.value)/obs.count.col.norm, 0)
#            } 
#            else 
#            {
#                count.col.norm <- sum(phi.vec < 0)
#                obs.count.col.norm <- sum(obs.phi.vec < 0)
#                count.col[i] <- ifelse(count.col.norm > 0, 
#                    sum(phi.vec <= ES.value)/count.col.norm, 0)
#                obs.count.col[i] <- ifelse(obs.count.col.norm > 0, 
#                    sum(obs.phi.vec <= ES.value)/obs.count.col.norm, 0)
#            }
#        }
#        phi.norm.mean[k] <- mean(count.col)
#        obs.phi.norm.mean[k] <- mean(obs.count.col)
#        phi.norm.median[k] <- median(count.col)
#        obs.phi.norm.median[k] <- median(obs.count.col)
#        FDR.mean[k] <- ifelse(phi.norm.mean[k]/obs.phi.norm.mean[k] < 1, 
#            phi.norm.mean[k]/obs.phi.norm.mean[k], 1)
#        FDR.median[k] <- ifelse(phi.norm.median[k]/obs.phi.norm.median[k] < 1, 
#            phi.norm.median[k]/obs.phi.norm.median[k], 1)
#    }
#
#    # adjust q-values
#    if (adjust.FDR.q.val ==TRUE) 
#    {
#        pos.nes <- sum(NES >= 0)
#        min.FDR.mean <- FDR.mean[pos.nes]
#        min.FDR.median <- FDR.median[pos.nes]
#        for(k in seq(pos.nes - 1, 1, -1)) 
#        {
#            if(FDR.mean[k] < min.FDR.mean) min.FDR.mean <- FDR.mean[k]
#            if(min.FDR.mean < FDR.mean[k]) FDR.mean[k] <- min.FDR.mean
#        }
#        neg.nes <- pos.nes + 1
#        min.FDR.mean <- FDR.mean[neg.nes]
#        min.FDR.median <- FDR.median[neg.nes]
#        for (k in seq(neg.nes + 1, Ng)) 
#        {
#             if(FDR.mean[k] < min.FDR.mean) min.FDR.mean <- FDR.mean[k]
#             if (min.FDR.mean < FDR.mean[k]) FDR.mean[k] <- min.FDR.mean
#        }
#    }   
#
#    obs.phi.norm.mean.sorted <- obs.phi.norm.mean[Orig.index]
#    phi.norm.mean.sorted <- phi.norm.mean[Orig.index]
#    FDR.mean.sorted <- FDR.mean[Orig.index]
#    FDR.median.sorted <- FDR.median[Orig.index]
#    
#    #   Compute global statistic
#    glob.p.vals <- vector(length=Ng, mode="numeric")
#    NULL.pass <- OBS.pass <- vector(length=nperm, mode="numeric")
#    for (k in seq_len(Ng)) 
#    {
#        if (NES[k] >= 0) 
#        {
#            for (i in seq_len(nperm)) 
#            {
#                NULL.pos <- sum(phi.norm[,i] >= 0)            
#                NULL.pass[i] <- ifelse(NULL.pos > 0, 
#                    sum(phi.norm[,i] >= NES[k])/NULL.pos, 0)
#                OBS.pos <- sum(obs.phi.norm[,i] >= 0)
#                OBS.pass[i] <- ifelse(OBS.pos > 0, 
#                    sum(obs.phi.norm[,i] >= NES[k])/OBS.pos, 0)
#            }
#        } 
#        else 
#        {
#            for (i in seq_len(nperm)) 
#            {
#                NULL.neg <- sum(phi.norm[,i] < 0)
#                NULL.pass[i] <- ifelse(NULL.neg > 0, 
#                    sum(phi.norm[,i] <= NES[k])/NULL.neg, 0)
#                OBS.neg <- sum(obs.phi.norm[,i] < 0)
#                OBS.pass[i] <- ifelse(OBS.neg > 0, 
#                    sum(obs.phi.norm[,i] <= NES[k])/OBS.neg, 0)
#            }
#        }
#        glob.p.vals[k] <- sum(NULL.pass >= mean(OBS.pass))/nperm
#    }
#    glob.p.vals.sorted <- glob.p.vals[Orig.index]

    # Produce results report
    Obs.ES <- signif(Obs.ES, digits=5)
    Obs.ES.norm <- signif(Obs.ES.norm, digits=5)
    p.vals <- signif(p.vals, digits=4)
#    signal.strength <- signif(signal.strength, digits=3)
#    tag.frac <- signif(tag.frac, digits=3)
#    gene.frac <- signif(gene.frac, digits=3)
#    FDR.mean.sorted <- signif(FDR.mean.sorted, digits=5)
#    FDR.median.sorted <-  signif(FDR.median.sorted, digits=5)
#    glob.p.vals.sorted <- signif(glob.p.vals.sorted, digits=5)

    report <- DataFrame(gs.names, size.G, Obs.ES, Obs.ES.norm, p.vals[,1])
#            , FDR.mean.sorted, p.vals[,2], tag.frac, 
#         gene.frac, signal.strength, FDR.median.sorted, glob.p.vals.sorted))
    colnames(report) <- c("GS", "SIZE", "ES", "NES", config.ebrowser("GSP.COL"))#, 
    rownames(report) <- NULL
    return(report)
#         "FDR q-val", "FWER p-val", "Tag \\%", "Gene \\%", "Signal", 
#         "FDR (median)", "glob.p.val")
#    report2 <- report[order(Obs.ES.norm, decreasing=TRUE),]
#    report3 <- report[order(Obs.ES.norm, decreasing=FALSE),]
#    phen1.rows <- sum(Obs.ES.norm >= 0)
#    phen2.rows <- length(Obs.ES.norm) - phen1.rows
#    report.phen1 <- report2[seq_len(phen1.rows),]
#    report.phen2 <- report3[seq_len(phen2.rows),]

#    if (output.directory != "")  
#    {
#        if (phen1.rows > 0) 
#        {
#            filename <- paste(output.directory, doc.string, 
#                ".SUMMARY.RESULTS.REPORT.", phen1,".txt", sep="", collapse="")
#            write.table(report.phen1, 
#                file = filename, quote=FALSE, row.names=FALSE, sep = "\t")
#        }
#        if (phen2.rows > 0) {
#            filename <- paste(output.directory, doc.string, 
#                ".SUMMARY.RESULTS.REPORT.", phen2,".txt", sep="", collapse="")
#            write.table(report.phen2, 
#                file = filename, quote=FALSE, row.names=FALSE, sep = "\t")
#        }
#    }
#
#    return(list(report1 = report.phen1, report2 = report.phen2))

}  # end of definition of GSEA.analysis

# Auxiliary functions and definitions 

GSEA.GeneRanking <- function(A, class.labels, gene.labels, nperm) 
{ 

# This function ranks the genes according to the signal to noise ratio for the 
# actual phenotype and also random permutations and bootstrap subsamples of both
# the observed and random phenotypes. It uses matrix operations to implement the
# signal to noise calculation in stages and achieves fast execution speed. It 
# supports two types of permutations: random (unbalanced) and balanced. It also
# supports subsampling and bootstrap by using masking and multiple-count
# variables.  When "fraction" is set to 1 (default) the there is no subsampling
# or boostrapping and the matrix of observed signal to noise ratios will have 
# the same value for all permutations. This is wasteful but allows to support 
# all the multiple options with the same code. Notice that the second matrix for
# the null distribution will still have the values for the random permutations 
# (null distribution). This mode (fraction = 1.0) is the defaults, the 
# recommended one and the one used in the examples.It is also the one that has 
# be tested more thoroughly. The resampling and boostrapping options are 
# intersting to obtain smooth estimates of the observed distribution but its is
# left for the expert user who may want to perform some sanity checks before 
# trusting the code.
#
# Inputs:
#   A: Matrix of gene expression values (rows are genes, columns are samples) 
#   class.labels: Phenotype of class disticntion of interest. 
#       A vector of binary labels having first the 1's and then the 0's 
#   gene.labels: gene labels. Vector of probe ids or 
#       accession numbers for the rows of the expression matrix 
#   nperm: Number of random permutations/bootstraps to perform 
#
# Outputs:
#   s2n.matrix: Matrix with random permuted or bootstraps signal to noise ratios
#       (rows are genes, columns are permutations or bootstrap subsamplings)
#   obs.s2n.matrix: Matrix with observed signal to noise ratios 
#       (rows are genes, columns are boostraps subsamplings. 
#           If fraction is set to 1.0 then all the columns have the same values
#   order.matrix: Matrix with the orderings that will 
#       sort the columns of the obs.s2n.matrix in decreasing s2n order
#   obs.order.matrix: Matrix with the orderings that will 
#       sort the columns of the s2n.matrix in decreasing s2n order
#
    A <- A + 0.00000001

    N <- nrow(A)
    Ns <- ncol(A)

    subset.mask <- reshuffled.class.labels1 <- reshuffled.class.labels2 <- 
        class.labels1 <- class.labels2 <- matrix(0, nrow=Ns, ncol=nperm)
    order.matrix <- obs.order.matrix <- s2n.matrix <- 
        obs.s2n.matrix <- matrix(0, nrow = N, ncol = nperm)
    #obs.gene.labels <- obs.gene.descs <- 
    #    obs.gene.symbols <- vector(length = N, mode="character")
    M1 <- M2 <- S1 <- S2 <- matrix(0, nrow = N, ncol = nperm)
    C <- split(class.labels, class.labels)
    class1.size <- length(C[[1]])
    class2.size <- length(C[[2]])
    class1.index <- seq_len(class1.size)
    class2.index <- (class1.size + 1):(class1.size + class2.size)

    for (r in seq_len(nperm)) 
    {
        class1.subset <- sample(class1.index, size = ceiling(class1.size))
        class2.subset <- sample(class2.index, size = ceiling(class2.size))
        subset.class1 <- as.integer(class1.index %in% class1.subset)
        subset.class2 <- as.integer(class2.index %in% class2.subset)
       
        subset.mask[, r] <- as.numeric(c(subset.class1, subset.class2))
        fraction.class1 <- class1.size/Ns
        fraction.class2 <- class2.size/Ns
        # random (unbalanced) permutation
        full.subset <- c(class1.subset, class2.subset)
        label1.subset <- sample(full.subset, size = Ns * fraction.class1)
        for (i in seq_len(Ns)) 
        {
            m1 <- sum(!is.na(match(label1.subset, i)))
            m2 <- sum(!is.na(match(full.subset, i)))
            reshuffled.class.labels1[i, r] <- m1
            reshuffled.class.labels2[i, r] <- m2 - m1
            if (i <= class1.size) 
            {
                class.labels1[i, r] <- m2
                class.labels2[i, r] <- 0
            } 
            else 
            {
             class.labels1[i, r] <- 0
             class.labels2[i, r] <- m2
            }
        }
    }

    # compute S2N for the random permutation matrix
    P <- reshuffled.class.labels1 * subset.mask
    n1 <- sum(P[,1])         
    M1 <- A %*% P
    M1 <- M1/n1      
    A2 <- A*A        
    S1 <- A2 %*% P   
    S1 <- S1/n1 - M1*M1    
    S1 <- sqrt(abs((n1/(n1-1)) * S1))   
    P <- reshuffled.class.labels2 * subset.mask
    n2 <- sum(P[,1])           
    M2 <- A %*% P           
    M2 <- M2/n2          
    A2 <- A*A           
    S2 <- A2 %*% P      
    S2 <- S2/n2 - M2*M2 
    S2 <- sqrt(abs((n2/(n2-1)) * S2))
    rm(P)
    rm(A2)
   
    # small sigma "fix" as used in GeneCluster
    S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
    S2 <- ifelse(S2 == 0, 0.2, S2)
    S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
    S1 <- ifelse(S1 == 0, 0.2, S1)
    M1 <- M1 - M2
    rm(M2)
    S1 <- S1 + S2
    rm(S2)
    s2n.matrix <- M1/S1
    order.matrix <- apply(s2n.matrix, 2, order, decreasing=TRUE)
    
    # compute S2N for the "observed" permutation matrix
    P <- class.labels1 * subset.mask
    n1 <- sum(P[,1])         
    M1 <- A %*% P
    M1 <- M1/n1      
    A2 <- A*A        
    S1 <- A2 %*% P   
    S1 <- S1/n1 - M1*M1    
    S1 <- sqrt(abs((n1/(n1-1)) * S1))   
    P <- class.labels2 * subset.mask
    n2 <- sum(P[,1])           
    M2 <- A %*% P           
    M2 <- M2/n2          
    A2 <- A*A           
    S2 <- A2 %*% P      
    S2 <- S2/n2 - M2*M2 
    S2 <- sqrt(abs((n2/(n2-1)) * S2))
    rm(P)
    rm(A2)

    # small sigma "fix" as used in GeneCluster
    S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
    S2 <- ifelse(S2 == 0, 0.2, S2)
    S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
    S1 <- ifelse(S1 == 0, 0.2, S1)
    M1 <- M1 - M2
    rm(M2)
    S1 <- S1 + S2
    rm(S2)
    obs.s2n.matrix <- M1/S1
    obs.order.matrix <- apply(obs.s2n.matrix, 2, order, decreasing=TRUE)

    return(list(s2n.matrix = s2n.matrix, 
                obs.s2n.matrix = obs.s2n.matrix, 
                order.matrix = order.matrix,
                obs.order.matrix = obs.order.matrix))
}

GSEA.EnrichmentScore <- function(
    gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) 
{  
#
# Computes the weighted GSEA score of gene.set in gene.list. 
# The weighted score type is the exponent of the correlation 
# weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). 
# When the score type is 1 or 2 it is necessary to input the correlation vector
# with the values in the same order as in the gene list.
#
# Inputs:
#   gene.list: The ordered gene list 
#       (e.g. integers indicating the original position in the input dataset)  
#   gene.set: A gene set (e.g. integers indicating 
#       the location of those genes in the input dataset) 
#   weighted.score.type: Type of score: 
#       weight: 0 (unweighted = Kolmogorov-Smirnov), 
#       1 (weighted), and 2 (over-weighted)  
#   correl.vector: A vector with the coorelations (e.g. signal to noise scores)
#       corresponding to the genes in the gene list 
#
# Outputs:
#   ES: Enrichment score (real number between -1 and +1) 
#   arg.ES: Location in gene.list where the peak 
#       running enrichment occurs (peak of the "mountain") 
#   RES: Numerical vector containing the running 
#       enrichment score for all locations in the gene list 
#   tag.indicator: Binary vector indicating 
#       the location of the gene sets (1's) in the gene list 

    # notice that the sign is 0 (no tag) or 1 (tag) 
    tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))
    no.tag.indicator <- 1 - tag.indicator 
    N <- length(gene.list) 
    Nh <- length(gene.set) 
    Nm <-  N - Nh 
    if (weighted.score.type == 0) correl.vector <- rep(1, N)
    alpha <- weighted.score.type
    correl.vector <- abs(correl.vector**alpha)
    sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
    norm.tag    <- 1.0 / sum.correl.tag
    norm.no.tag <- 1.0 / Nm
    RES <- cumsum(tag.indicator * correl.vector * 
        norm.tag - no.tag.indicator * norm.no.tag)      
    max.ES <- max(RES)
    min.ES <- min(RES)
    if (max.ES > -min.ES) 
    {
        ES <- signif(max.ES, digits = 5)
        arg.ES <- which.max(RES)
    } 
    else 
    {
        ES <- signif(min.ES, digits=5)
        arg.ES <- which.min(RES)
    }
    return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
}


GSEA.EnrichmentScore2 <- function(
    gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) 
{  
#
# Computes the weighted GSEA score of gene.set in gene.list. It is the same 
# calculation as in GSEA.EnrichmentScore but faster (x8) without producing the 
# RES, arg.RES and tag.indicator outputs.
# This call is intended to be used to asses the enrichment of random 
# permutations rather than the observed one.
# The weighted score type is the exponent of the correlation 
# weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). 
# When the score type is 1 or 2 it is necessary to input the correlation vector
# with the values in the same order as in the gene list.
#
# Inputs:
#   gene.list: The ordered gene list 
#       (e.g. integers indicating the original position in the input dataset)  
#   gene.set: A gene set 
#       (e.g. integers indicating the location of those genes in the input dataset) 
#   weighted.score.type: Type of score: 
#       weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)  
#   correl.vector: A vector with the coorelations (e.g. signal to noise scores)
#       corresponding to the genes in the gene list 
#
# Outputs:
#   ES: Enrichment score (real number between -1 and +1) 
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.
    
    N <- length(gene.list) 
    Nh <- length(gene.set) 
    Nm <-  N - Nh 
    
    peak.res.vector <- valley.res.vector <- 
    tag.diff.vector <- vector(length=Nh, mode="numeric")
    
    loc.vector <- vector(length=N, mode="numeric")
    loc.vector[gene.list] <- seq_len(N)
    tag.loc.vector <- loc.vector[gene.set]
    tag.loc.vector <- sort(tag.loc.vector, decreasing =FALSE)
 
    if (weighted.score.type == 0) tag.correl.vector <- rep(1, Nh)
    else if (weighted.score.type == 1) 
    {
        tag.correl.vector <- correl.vector[tag.loc.vector]
        tag.correl.vector <- abs(tag.correl.vector)
    } 
    else if (weighted.score.type == 2) 
    {
        tag.correl.vector <- 
            correl.vector[tag.loc.vector] * correl.vector[tag.loc.vector]
        tag.correl.vector <- abs(tag.correl.vector)
    } 
    else 
    {
        tag.correl.vector <- correl.vector[tag.loc.vector] * weighted.score.type
        tag.correl.vector <- abs(tag.correl.vector)
    }
 
    norm.tag <- 1.0/sum(tag.correl.vector)
    tag.correl.vector <- tag.correl.vector * norm.tag
    norm.no.tag <- 1.0/Nm
    tag.diff.vector[1] <- (tag.loc.vector[1] - 1) 
    tag.diff.vector[2:Nh] <- 
        tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
    tag.diff.vector <- tag.diff.vector * norm.no.tag
    peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
    valley.res.vector <- peak.res.vector - tag.correl.vector
    max.ES <- max(peak.res.vector)
    min.ES <- min(valley.res.vector)
    ES <- signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)
 
    return(list(ES = ES))
}

