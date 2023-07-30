#' @name plots
#' 
#' @title Visualization of gene expression
#' 
#' @description Visualization of differential gene expression via heatmap, p-value histogram
#' and volcano plot (fold change vs. p-value).
#' 
#' 
#' @param p Numeric vector of p-values for each gene.
#' @param fc Numeric vector of fold changes (typically on log2 scale).
#' @param expr Expression matrix. Rows correspond to genes, columns to samples.
#' @param grp *BINARY* group assignment for the samples.  Use '0' and '1' for
#' unaffected (controls) and affected (cases) samples, respectively.
#' @param scale.rows Should rows of the expression matrix be scaled for better
#' visibility of expression differences between sample groups? Defaults to
#' TRUE.
#' @return None, plots to a graphics device.
#' @author Ludwig Geistlinger
#' @seealso \code{\link{deAna}} for differential expression analysis,
#' \code{ComplexHeatmap::Heatmap}, and \code{\link{hist}} for generic plotting.
#' @examples
#' 
#'     # (1) simulating expression data: 100 genes, 12 samples
#'     se <- makeExampleData(what="SE") 
#'     
#'     # plot heatmap
#'     exprsHeatmap(expr=assay(se), grp=as.factor(se$GROUP))
#' 
#'     # (2) DE analysis
#'     se <- deAna(se)
#'     pdistr(rowData(se)$ADJ.PVAL)
#'     volcano(fc=rowData(se)$FC, p=rowData(se)$ADJ.PVAL)
#' 
NULL


# P-Value Distribution
#' @rdname plots
#' @export
pdistr <- function(p)
{
    hist(p, breaks=100, prob=TRUE, col="cyan",
            main="P-Value Distribution", xlab="P-Value", ylab="Frequency")
}

# Volcano Plot (fold change vs. p-value)
#' @export
#' @rdname plots
volcano <- function(fc, p)
{
    plot(x=fc, 
        y=-log(p, base=10), col="red", 
            main="Volcano Plot", xlab="log2(foldChange)", ylab="-log10(p)")
}

# Heatmap: based on ComplexHeatmap
#' @export
#' @rdname plots
exprsHeatmap <- function(expr, grp, scale.rows = TRUE)
{
    isAvailable("ComplexHeatmap", type = "software")

    dtype <- .detectDataType(expr)
    if(dtype == "rseq") edgeR::cpm(expr, log = TRUE) 
  
    # scale?
    expr <- t(scale(t(expr)))
    expr[is.nan(expr)] <- 0

	# group colors
	grp <- as.factor(grp)
	coll <- c("#B62A84", "#2AB68C")
	names(coll) <- levels(grp)
	coll <- list(Group=coll)

	# annotation
	df <- data.frame(Group = grp)
    ha <- ComplexHeatmap::HeatmapAnnotation(df = df, col=coll)

	# plot
    print(ComplexHeatmap::Heatmap(expr, name="Expression", top_annotation = ha, 
        show_row_names=(nrow(expr) < 41), show_column_names=(ncol(expr) < 41),
        column_title="Samples", row_title="Features")) 
}

