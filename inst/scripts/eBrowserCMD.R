usage <- function() {
    use.message <- paste(
        "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++", 
        "\tTHE ENRICHMENT BROWSER\n", 
        "\tCALL:\tRscript eBrowserCMD.R\n", 
        "\tARGS:\t--meth=<meth>\tEnrichment method(s)", 
        "\t\t--exprs=<exprs>\tExpression matrix; tab-separated text file", 
        "\t\t--cdat=<cdat>\tColumn data (sample annotation); tab-separated text file", 
        "\t\t--rdat=<rdat>\tRow data (gene annotation); tab-separated text file", 
        "\t\t--org=<org>\tOrganism under investigation; 3-letter code", 
        "\t\t--gs=<gs>\tGene sets; tab-separated text file (GMT format)\n", 
        "\tOPTS:\t--grn=<grn>\tGene regulatory network; tab-separated text file", 
        "\t\t--dtype=<dtype>\tExpression data type (Default: NA)", 
        "\t\t--nmeth=<nmeth>\tNormalization method (Default: none)", 
        "\t\t--dmeth=<dmeth>\tDifferential expression method (Default: limma)", 
        "\t\t--alpha=<alpha>\tStatistical significance level (Default: 0.05)", 
        "\t\t--beta=<alpha>\tLog2 fold change significance level (Default: 1)", 
        "\t\t--perm=<perm>\tNr of permutations (Default: 1000)", 
        "\t\t--outd=<outd>\tOutput directory (Default: ./)", 
        "\t\t--show=<show>\tNr of gene sets to include in HTML output (Default: 10)", 
        "\t\t\t\t", "\tFLAGS:\t--help, -h\tCalls this usage message", 
        "\t\t--html\t\tProduce html report in addition to plain text output", 
        "\t\t--comb\t\tCombination of results of multiple enrichment methods", 
        "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n", 
        sep = "\n")
    
    message(use.message)
    message(paste0("See http://www.bioconductor.org/packages/release/bioc/manuals",
                    "/EnrichmentBrowser/man/EnrichmentBrowser.pdf#page.10\n"))
    quit(save = "no")
}


## parse command line arguments
parseCmdArgs <- function(cmd.args) {
    spl.args <- sapply(cmd.args, function(arg) unlist(strsplit(arg, "=")))
    ebrowser.args <- spl.args[sapply(spl.args, function(spl) length(spl) == 2)]
    flags <- sapply(ebrowser.args, function(arg) arg[1])
    args <- sapply(ebrowser.args, function(arg) arg[2])
    names(args) <- flags
    return(args)
}

## get a particular argument from the command line parameters
getArg <- function(arg, args) {
    if (!(arg %in% names(args))) 
        stop(paste("The required argument", arg, "is missing!"))
    return(args[[arg]])
}

## clean up output structure
cleanUp <- function(out.dir) {
    dirs <- list.files(out.dir, pattern = "_RESULT_files$", full.names = TRUE)
    
    for (d in dirs) {
        files <- list.files(d, full.names = TRUE)
        for (f in files) file.rename(from = f, to = paste(out.dir, basename(f), sep = ""))
        file.remove(d)
    }
}

autoDetectDataType <- function(expr) 
    ifelse(all(isWholenumber(expr), na.rm = TRUE), "rseq", "ma")

isWholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol


## eBrowser functionality
main <- function(cmd.args) {
    if (length(cmd.args) < 6 || any(c("--help", "-h") %in% cmd.args)) usage()
    
    args <- parseCmdArgs(cmd.args)
    
    # required params
    exprs <- getArg("--exprs", args)
    cdat <- getArg("--cdat", args)
    rdat <- getArg("--rdat", args)
    org <- getArg("--org", args)
    gs <- getArg("--gs", args)
    
    meth <- NULL
    if ("--meth" %in% names(args)) 
        meth <- tolower(unlist(strsplit(getArg("--meth", args), ",")))
    
    message("Loading EnrichmentBrowser")
    suppressPackageStartupMessages(library(EnrichmentBrowser))
    
    # optional params
    grn <- NULL
    if ("--grn" %in% names(args)) {
        grn <- args[["--grn"]]
        if (grn == "KEGG") {
            message("Compiling GRN from KEGG ...")
            grn <- compileGRN(org, "kegg")
        } else if (file.exists(grn)) 
            grn <- gsub(" ", "", as.matrix(read.delim(grn)))
    }
    
    data.type <- NULL
    if ("--dtype" %in% names(args)) 
        data.type <- args[["--dtype"]]
    
    norm.method <- "none"
    if ("--nmeth" %in% names(args)) 
        norm.method <- args[["--nmeth"]]
    
    de.method <- "limma"
    if ("--dmeth" %in% names(args)) 
        de.method <- args[["--dmeth"]]
    
    alpha <- 0.05
    if ("--alpha" %in% names(args)) 
        alpha <- as.numeric(args[["--alpha"]])
    
    beta <- 1
    if ("--beta" %in% names(args)) 
        alpha <- as.numeric(args[["--beta"]])
    
    perm <- 1000
    if ("--perm" %in% names(args)) 
        perm <- as.integer(args[["--perm"]])
    
    out.dir <- file.path(getwd(), "eBrowserOut")
    if (file.exists(out.dir)) 
        system(command = paste("rm -rf", file.path(out.dir, "*")))
    if ("--outd" %in% names(args)) {
        outd <- args[["--outd"]]
        if (dirname(outd) %in% list.files(".")) 
            outd <- file.path(getwd(), outd)
        if (file.exists(outd)) 
            system(command = paste("rm -rf", file.path(outd, "*")))
    }
    
    nr.show <- 10
    if ("--show" %in% names(args)) 
        nr.show <- as.integer(args[["--show"]])
    
    # flags
    html <- "--html" %in% cmd.args
    comb <- "--comb" %in% cmd.args
    
    # execution
    message("Executing EnrichmentBrowser")
    message("Retrieving gene sets")
    if (file.exists(gs)) 
        gs <- getGenesets(gs) else if (gs == "GO") 
        gs <- getGenesets(org, db = "go") else if (gs == "KEGG") 
        gs <- getGenesets(org, db = "kegg")
    configEBrowser("OUTDIR.DEFAULT", out.dir)
    se <- readSE(exprs, cdat, rdat)
    if (substring(names(se)[1], 1, 3) == "ENS") {
        message("Detected gene ID type: ENSEMBL")
        message("Map IDs: ENSEMBL -> ENTREZ")
        se <- idMap(se, org = org, from = "ENSEMBL", to = "ENTREZID")
    }
    
    if (is.null(data.type)) 
        data.type <- autoDetectDataType(assay(se))
    if (is.null(meth)) {
        message("No enrichment method selected")
        meth <- sbeaMethods()
        if (!is.null(grn)) 
            meth <- c(meth, nbeaMethods())
        message(paste("Executing", paste(meth, collapse = ",")))
    }
    
    ebrowser(meth = meth, exprs = se, org = org, data.type = data.type, norm.method = norm.method, 
        de.method = de.method, gs = gs, grn = grn, perm = perm, alpha = alpha, beta = beta, 
        comb = comb, browse = html, nr.show = nr.show)
    
    if ("--outd" %in% names(args)) {
        setwd(dirname(out.dir))
        message(paste("Moving files to", outd))
        system(command = paste("cp -r", file.path(out.dir, "reports", "*"), outd))
        message("Cleaning up ...")
        system(command = paste("rm -rf", out.dir))
    }
    
}

## MAIN
main(commandArgs())
