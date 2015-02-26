###
#
# HTML Writers
#
###

## 
## ## ## method2html
## 
## Generic visualization of results in table format

method2html <- function(    method, 
                res.tbl, 
                set.view=NULL, 
                graph.view=NULL, 
                html.out, 
                header.links=NULL)
{
    # SET VIEW
    if(!is.null(set.view))
    {
        res.tbl <- cbind(res.tbl, make.view.col(set.view))
        colnames(res.tbl)[ncol(res.tbl)] <- "SET.VIEW"
    }

    # GRAPH VIEW
    if(!is.null(graph.view))
    {
        res.tbl <- cbind(res.tbl, make.view.col(graph.view))
        colnames(res.tbl)[ncol(res.tbl)] <- "GRAPH.VIEW"
    }

    ## create and write the header to html
    if(file.exists(html.out)) file.remove(html.out)
    html.out <- file(html.out, open='a')
    header <- create.html.header(method=method, 
        title=TABLE.OF.RESULTS, html.out=html.out)

    ## create and write the result table to html
    tbl <- r.table.2.html.table(res.tbl, html.out, header.links=header.links)
    writeLines("</ul></body></html>", html.out)
    close(html.out)
}

create.index <- function(out.dir, meth, comb=FALSE)
{
    index <- file.path(out.dir, "index.html")
    if(file.exists(index)) file.remove (index)
    index <- file(index, open='a')
    header <- create.html.header(method="eBrowser", 
        title="Index of Result Files", html.out=index)
    
    res.files <- list.files(out.dir, full.names=TRUE)
    has.subdirs <- length(grep("RESULT_files$", res.files))

    # report combined results
    if(comb)
    {
        comb.res <- "eBrowser_comb_RESULT."
        elems <- c("<h4 class=\"title30\">Combined Results</h4>",
            "<ul class=\"simple\">",
            paste0("<li><a href=\"", comb.res, 
                "html\">", "Top Table", "</a></li>"),
            paste0("<li><a href=\"", comb.res, 
                "txt\">", "Full Ranking", "</a></li>"),
            "</ul>")
        invisible(sapply(elems, function(e) writeLines(e, index)))    
    }

    # report results of each method
    sapply(meth, 
        function(m)
        {
            m.resdir <- ""
            if(has.subdirs) 
                m.resdir <- paste("eBrowser_", m, "_RESULT_files/", sep="") 
            m.res <- paste(m.resdir, "eBrowser_", m, "_RESULT.", sep="")

            elems <- c(
                paste0("<h4 class=\"title30\">", toupper(m), " Results</h4>"),
                "<ul class=\"simple\">",
                paste0("<li><a href=\"", m.res, 
                    "html\">", "Top Table", "</a></li>"),
                paste0("<li><a href=\"", m.res, 
                    "txt\">", "Full Ranking", "</a></li>"),
                "</ul>")
    
            invisible(sapply(elems, function(e) writeLines(e, index)))    
        })

    # link genewise differential expression
    elems <- c(
        "<h4 class=\"title30\">Genewise Differential Expression</h4>",
        "<ul class=\"simple\">",
        "<li><a href=\"gene_diffexp.txt\">DE Measures for each Gene</a></li>",
        "<li><a href=\"heatmap.png\">Heatmap</a></li>",
        "<li><a href=\"pdistr.png\">P-Value Distribution</a></li>",
        "<li><a href=\"volcano.png\">Volcano Plot</a></li>",
        "</ul>", "</ul>", "</body>", "</html>")
    invisible(sapply(elems, function(e) writeLines(e, index)))
    close(index)
}

## converts an R data table into an HTML table
## returns it as a string
r.table.2.html.table <- function(r.table, html.out, header.links=NULL)
{
    r.table <- as.matrix(r.table)
    writeLines("<table border=\"1\">", html.out)
    ## write the table header
    if(is.null(header.links)) header.parts <- sapply(colnames(r.table), 
                            function(x) paste("<th>",x,"</th>", sep=""))
    else
    {
        hll <- length(header.links)
    
        header.parts <- c(
            # the normal header: GENE.SET
            paste0("<th>", colnames(r.table)[1], "</th>"),
            # the rank headers: ORA.RANK, ..., AV.RANK
            sapply(seq_len(hll-1),
                function(i) paste0("<th><a href=\"", header.links[i], 
                "\">", colnames(r.table)[i+1], "</a></th>")),
            # the pvalue headers:
            # (a) the methods: ORA.PVAL, ...
            sapply(seq_len(hll-2),
                function(i) paste0("<th><a href=\"", header.links[i],
                "\">", colnames(r.table)[i+hll], "</a></th>")),
            # (b) the comb: MERGED.PVAL
            paste0("<th><a href=\"", header.links[hll], 
            "\">", colnames(r.table)[2*hll-1], "</a></th>"),
            # (c) the views 
            sapply(colnames(r.table)[(2*hll):ncol(r.table)],
                function(x) paste0("<th>",x,"</th>")))
    }

    header <- c("<tr>", header.parts, "</tr>")
    invisible(sapply(header, function(e) writeLines(e, html.out)))
    
    ## write the table body
    rows <- apply(r.table, 1, function(r) write.html.row(r, html.out))
    writeLines("</table>", html.out)
}

## create the header of an html output file
create.html.header <- function(method, title, html.out)
{
    images <- list.files(system.file("images/", package="EnrichmentBrowser"))
    method.icon <- grep(method, images, value=TRUE)
    method.icon <- system.file( 
        paste0("images/", method.icon), package="EnrichmentBrowser")

    ## (1) create the html head
    header.strs <-
        c("<?xml version=\"1.0\" encoding=\"UTF-8\"?>",
        "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\"",
        " \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">",
        "<html xmlns=\"http://www.w3.org/1999/xhtml\">", "<head>",
        paste0("<title>", method, " Result</title>"),
        paste0("<link rel=\"stylesheet\" href=\"file://",
            file.path(system.file("images", package="EnrichmentBrowser"), CSS), 
            "\" type=\"text/css\" />"), "</head>",
        "<body bgcolor=\"#ffffff\" link=\"#003399\" vlink=\"#003399\">",
        "<table width=\"700\">", "<tr>", "<td width=\"180\" align=\"center\">",
        paste0("<img src=\"", method.icon, "\" alt=\"", method, 
        "\" height=\"70\" width=\"140\" border=\"0\" />"), "</td>", "<td>",
        paste0("<h2 class=\"title30\">", method, " - ", title, "</h2>"),
        "</td>", "</tr>", "</table>", "<br />")
    
    invisible(sapply(header.strs, function(s) writeLines(s, html.out)))
}

## write a row of a R table in Html format
write.html.row <- function(row, html.out)
{
    row.elems <- c("<tr>", 
        sapply(row, function(x) paste("<td>",x,"</td>", sep="")), "</tr>")
    writeLines(paste(row.elems, collapse=""), html.out)
}

## make a href col 
make.href.col <- function(ref, tag)
    paste0("<a href=\"", ref, "\">", tag, "</a>")    

## make a number of href cols
make.href.cols <- function(refs, tag)
    sapply(refs, make.href.col, tag=tag)

# collapse all hrefs of a view together
make.view.col <- function(view)
{
    VIEW.TAGS <- sub("IMG.DIR", 
        system.file("images", package="EnrichmentBrowser"), VIEW.TAGS)
    view.hrefs <- sapply(seq_len(ncol(view)), 
        function(i) make.href.cols( refs=view[,i], tag=VIEW.TAGS[i]))
    if(is.matrix(view.hrefs)) 
        view.col <- apply(view.hrefs, 1, paste, collapse=" ")
    else view.col <- paste(view.hrefs, collapse=" ")
    return(view.col)
}

## add a
add.href.cols.2.res.table <- function(method, res.tbl, refs, tag)
{
    cols <- make.href.cols(refs, tag=tag[2])
    res.tbl <- cbind(res.tbl, cols)
    colnames(res.tbl)[ncol(res.tbl)] <- paste(method, tag[1], sep = ".")
    return(res.tbl)
}

