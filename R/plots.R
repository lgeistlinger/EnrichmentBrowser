# P-Value Distribution
pdistr <- function(p)
{
    MASS::truehist(p, nbins=100, prob=TRUE,
            main="P-Value Distribution", xlab="P-Value", ylab="Frequency")
}

# Volcano Plot (fold change vs. p-value)
volcano <- function(fc, p)
{
    plot(x=fc, 
        y=-log(p, base=10), col="red", 
            main="Volcano Plot", xlab="log2(foldChange)", ylab="-log10(p)")
}

# Heatmap: based on ComplexHeatmap
exprs.heatmap <- function(expr, grp)
{
    #message("Transforming expression data to log2 scale for heatmap visualization")
    if(max(expr, na.rm=TRUE) - min(expr, na.rm=TRUE) > 100) expr <- log(expr + 1, base=2)
    df <- data.frame(Group = as.factor(grp))
    ha <- ComplexHeatmap::HeatmapAnnotation(df = df)
    print(ComplexHeatmap::Heatmap(expr, name="Expression", top_annotation = ha, 
        show_row_names=(nrow(expr) < 41), show_column_names=(ncol(expr) < 41),
        column_title="Samples", row_title="Features")) 
}

