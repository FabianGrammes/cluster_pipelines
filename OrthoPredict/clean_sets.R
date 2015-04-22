
args <- commandArgs(trailingOnly = TRUE)
print(args)

cmd.in <- args[1]

ls <- list.files( pattern = 'csv', path = cmd.in, full.names = TRUE )

library(biomaRt)

do.mart <- function( file ){
    bc <- read.csv( file )
    # run BioMart
    ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
    attr <- c('affy_hg_u133a', 'ensembl_gene_id', 'hgnc_symbol', 'description')
    bm <- getBM(attributes = attr,  
                filters = 'affy_hg_u133a',
                values = bc$Human.Probesetid,
                mart = ensembl)
    # remove genes with one 2 many hits (for affy probes.)
    tb <- table( bm$affy_hg_u133a )
    bm <- bm[ bm$affy_hg_u133a %in% names(tb[ tb == 1]), ]
    # write
    out.file1 <- gsub('.csv', '.processed.csv', file )
    out.file2 <- gsub('.csv', '.processed.ids', file )
    write.csv( bm, file = out.file1 )
    write.table( bm$ensembl_gene_id, file = out.file2, row.names = FALSE,
                col.names = FALSE, quote = FALSE)
    cat('FINISHED', file, 'Nrow:', nrow(bc), '| processed Nrow:', nrow(bm),'\n' )
}

sapply( ls, do.mart )
