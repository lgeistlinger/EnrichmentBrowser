# 
# Author: Ludwig Geistlinger
# Date: 5th August 2010
#
# download all pathways of a particular organism from the KEGG FTP
#
# UPDATE: 16 Jan 2014 - switch to KEGGREST functionality as KEGG FTP 
#           is no longer accessible
#
###############################################################################

#
# USING KEGGREST
#

download.kegg.pathways <- function(org, out.dir=NULL, zip=FALSE)
{
    pwys <- keggList("pathway", org)
    message("download:")
    kgmls <- vapply(names(pwys), 
        function(pwy)
        {
            message(pwy)
            pwy <- sub("^path:", "", pwy)
            kgml <- keggGet(pwy, "kgml")
            return(kgml)
        }, character(1))
    names(kgmls) <- sub("^path:", "", names(kgmls))

    if(is.null(out.dir))
    {
        pwys <- sapply(kgmls, parseKGML)
        return(pwys)
    }

    message("write kgml files ...")
    if(!file.exists(out.dir)) dir.create(out.dir)
    for(n in names(kgmls)) 
        cat(kgmls[n], file=file.path(out.dir, paste(n, ".xml", sep="")))

    if(zip) 
    {
        message("zip ...")
        pattern <- file.path(out.dir, "*.xml")
        zip.file <- file.path(out.dir, paste(org, "zip", sep="."))
        system(paste("zip", zip.file, pattern))
        for(f in list.files(out.dir, pattern=pattern)) file.remove(f)
    }
}

