# SAMGS : original method from Dinu et al - BMC Bioinformatics - 2007
#
# EDIT 17 Sep 2014, Ludwig Geistlinger: 
#   adapting for use in EnrichmentBrowser package
#
# EDIT 02 Aug 2015, Ludwig Geistlinger
#   adapting for use in SAFE framework (local and global stat)

# SAMGS stat as global.stat for safe
global.SAMGS <-
function(C.mat, u, ...)
{
    return(
        function(u, C.mat2 = C.mat) 
        {
            return(as.numeric(t(SparseM::as.matrix(C.mat2)) %*% u^2))
        }
    )
}

# SAM t-like stat as local.stat for safe
local.t.SAM <- function (X.mat, y.vec, ...)
{
    rMV <- function(data)
    {
        m <- rowMeans(data)
        dif <- data - m 
        ssd <- rowSums(dif^2)
        l <- list("means"=m,
        "sumsquaredif" =ssd,
        "vars" =ssd/(ncol(data)-1),
        "centered rows"=dif)
        return(l)
    }

    nb.Groups <- 100
    mad.Const <- .64

    stopifnot(length(unique(y.vec)) == 2)    
    if (!all(sort(unique(y.vec)) == c(0,1))) {
        warning("y.vec is not (0,1), thus Group 1 == ", y.vec[1])
        y.vec <- (y.vec == y.vec[1]) * 1
    }
    return(function(data, vec = y.vec, ...) 
    {
        C1 <- which(vec==0)
        C2 <- which(vec==1)
        nb.Genes    <- nrow(data)
        nb.Samples  <- ncol(data)
        C1.size     <- length(C1)
        C2.size     <- length(C2)

        stat.C1            <- rMV(data[,C1])
        stat.C2            <- rMV(data[,C2])
        diffmean.C1C2      <- stat.C1$means - stat.C2$means
        pooledSqrtVar.C1C2 <- sqrt( (1/C1.size+1/C2.size) * 
            (stat.C1$sumsquaredif + stat.C2$sumsquaredif) / (nb.Samples-2) )
    
        tmp <- as.data.frame(cbind(pooledSqrtVar.C1C2,diffmean.C1C2))
        tmp <- tmp[order(tmp[,1]),]
        group.Size     <- as.integer(nb.Genes/nb.Groups)
        percentiles    <- seq(0,1,.05)
        nb.Percentiles <- length(percentiles)
        s0.quantiles   <- quantile(pooledSqrtVar.C1C2,percentiles)

        tt <- matrix(NA,nb.Groups,nb.Percentiles)
        coeffvar <- as.data.frame(cbind(s0.quantiles,rep(NA,nb.Percentiles)))
        group.grid <- seq_len(group.Size*nb.Groups)
        denom <- tmp[group.grid,1]
        for(j in seq_len(nb.Percentiles))
        {
            nom <- tmp[group.grid,2] + s0.quantiles[j]
            x <- matrix(denom/nom, group.Size, nb.Groups)
            tt[,j] <- apply(x, 2, mad, constant=mad.Const)
            coeffvar[j,2] <- sd(tt[,j]) / mean(tt[,j])
        }
        s0 <-  min(s0.quantiles[coeffvar[,2] == min(coeffvar[,2])])
        TlikeStat <-  diffmean.C1C2 / (pooledSqrtVar.C1C2+s0)            
        return(TlikeStat)
    })
}

## START: original code
rowMeansVars <- function(DATA,margin=1)
{
 if(margin==2) DATA <- t(DATA)
 m <- rowMeans(DATA)
 dif <- DATA - m
 ssd <- rowSums(dif^2)
 list("means"=m,
      "sumsquaredif" =ssd,
      "vars" =ssd/(ncol(DATA)-1),
      "centered rows"=dif )
}

sam.TlikeStat <- function(DATA,
                         cl=NULL,
                         s0=NULL,
                         s0.param=list(nb.Groups=100,mad.Const=.64) )
{
# DATA : expression data
#     -> dataframe with  rows=genes,
#                        columns=samples,

# cl : factor defining a bipartition of the samples
#      IN THE SAME ORDER AS IN DATA
# NB : use     'C1' + 'C2'     XOR   'cl'  alone
   
   if(!is.null(cl))
   {
        cl <- as.factor(cl)
        C1 <- which(as.numeric(cl)==1)
        C2 <- which(as.numeric(cl)==2)
   }
   if(is.null(C1) | is.null(C2))
       stop("Error -  sam.TlikeStat : classes 1 and 2 are undefined.")

    nb.Genes    <- nrow(DATA)
    nb.Samples  <- ncol(DATA)
    C1.size     <- length(C1)
    C2.size     <- length(C2)

    stat.C1            <- rowMeansVars(DATA[,C1])
    stat.C2            <- rowMeansVars(DATA[,C2])
    diffmean.C1C2      <- stat.C1$means - stat.C2$means
    pooledSqrtVar.C1C2 <- sqrt( (1/C1.size+1/C2.size) * 
       (stat.C1$sumsquaredif + stat.C2$sumsquaredif) / (nb.Samples-2) )

    if(is.null(s0)){
          nb.Groups <- s0.param$nb.Groups
          mad.Const <- s0.param$mad.Const

          tmp <- as.data.frame(cbind(pooledSqrtVar.C1C2,diffmean.C1C2))
          tmp <- tmp[order(tmp[,1]),]

          group.Size     <- as.integer(nb.Genes/nb.Groups)
          percentiles    <- seq(0,1,.05)
          nb.Percentiles <- length(percentiles)
          s0.quantiles   <- quantile(pooledSqrtVar.C1C2,percentiles)

          tt <- matrix(NA,nb.Groups,nb.Percentiles)
          coeffvar <- as.data.frame(cbind(s0.quantiles,rep(NA,nb.Percentiles)))
          for(j in seq_len(nb.Percentiles)){
             x <- matrix( tmp[seq_len(group.Size*nb.Groups),1] / 
                           (tmp[seq_len(group.Size*nb.Groups),2] + 
                               s0.quantiles[j]), group.Size, nb.Groups)
             tt[,j] <- apply(x, 2, mad, constant=mad.Const)
             coeffvar[j,2] <- sd(tt[,j]) / mean(tt[,j])
          }

          s0 <-  min(s0.quantiles[coeffvar[,2] == min(coeffvar[,2])])
    }

    list(s0            = s0,
         diffmean      = diffmean.C1C2,
         pooledSqrtVar = pooledSqrtVar.C1C2,
         TlikeStat     = diffmean.C1C2/(pooledSqrtVar.C1C2+s0))
}


# GS.format.dataframe.to.list was called below in the function SAMGS

GS.format.dataframe.to.list <- function(GS)
{
    if(is.data.frame(GS)){
        genes <- rownames(GS)
        L <- NULL
        for(ags in names(GS)){
           w <- which(GS[,ags]==1)
           if(length(w)>0)  {
              L <- c(L,list(genes[w]))
              names(L)[length(L)] <- ags
           }
        }
        L
    }else{
        GS
    }
}



SAMGS <- function(GS,
                 DATA,
                 cl,
                 nbPermutations=1000,
                 silent=FALSE,
                 tstat.file=NULL)
{

# GS : gene sets
#      -> a dataframe with rows=genes,
#                        columns= gene sets,
#                        GS[i,j]=1 if gene i in gene set j
#                        GS[i,j]=0 otherwise
#         OR
#         a list with each element corresponding 
#           to a gene set = a vector of strings (genes identifiers)
#
# DATA : expression data
#     -> a dataframe with  rows=genes,
#                        columns=samples
#
# cl : a factor defining a bipartition of the 
#           samples IN THE SAME ORDER AS IN DATA
#
    genes <- rownames(DATA)
    GS <-  GS.format.dataframe.to.list(GS)
    GS <-  lapply(GS,function(z) which(genes %in% z))

    nb.Samples  <- ncol(DATA)
    nb.GeneSets <- length(GS)   # nb of gene sets
    GeneSets.sizes <- sapply(GS,length) # size of each gene set
    C1.size <- table(cl)[1]  # nb of samples in class 1

    # finding constant s0 for SAM-like test
    tmp     <- sam.TlikeStat(DATA,cl=cl)
    s0      <- tmp$s0

    # stats obtained on 'true' data
    # SAM T-like statistic for each gene

    samT.ok <- tmp$TlikeStat                        
    if(!is.null(tstat.file)) save(samT.ok, file=tstat.file)    
    
    # SAMGS statitic for each gene set
    sam.sumsquareT.ok <- sapply(GS,function(z) sum(samT.ok[z]^2))  
    # stats obtained on 'permuted' data
    permut.C1 <- matrix(NA,nbPermutations,C1.size)
    sam.sumsquareT.permut <- matrix(NA,nbPermutations,nb.GeneSets)
    cl.permut=cl #initial group values
    
    for(i in seq_len(nbPermutations)) 
    {
        C1.permut <-  permut.C1[i,] <- sample(nb.Samples,C1.size)
        C2.permut <-  seq_len(nb.Samples)[-C1.permut]
        cl.permut[C1.permut] <- 1
        cl.permut[C2.permut] <- 2
        samT.permut <- sam.TlikeStat(DATA,cl=cl.permut,s0=s0)$TlikeStat         
        sam.sumsquareT.permut[i,] <- 
            sapply(GS,function(z) sum(samT.permut[z]^2)) 
        # SAMGS statitic for each gene set  - for current permutation
        if(!silent && (i%%50 == 0)) message(paste(i," permutations done."))
    }

    GeneSets.pval <- apply(
        t(sam.sumsquareT.permut) >= sam.sumsquareT.ok, 1, sum) / nbPermutations
    #GeneSets.qval <- qvalue::qvalue(GeneSets.pval)$qvalues
    #GeneSets.pval <- GeneSets.qval
    norm.stat <- sam.sumsquareT.ok / sapply(GS, length)
    res.tbl <- cbind(sam.sumsquareT.ok, norm.stat, GeneSets.pval)
    colnames(res.tbl) <- c("SUMSQ.STAT", "NSUMSQ.STAT", config.ebrowser("GSP.COL")) 
    rownames(res.tbl) <- names(GS)
    
    return(res.tbl)
}

