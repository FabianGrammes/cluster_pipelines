
# Collects all gene IDs from all *.ssa.id files in a given path and writes them
# to table

args <- commandArgs(trailingOnly = TRUE)
print(args)

cmd.path <- args[1]
cmd.set <- args[2]
cmd.file <- args[3]

#-------------------------------------------------------------------------------
require( CigSsa.db )
ls <- list.files( pattern = '.ssa.id', path = cmd.path, full.names = TRUE )

read.ids <- function( x, gs ){
    t <- read.table( x, header = F )[,1]
    cn <- unlist( strsplit( x, '\\/' ) )
    cn <- gsub('.ssa.ids', '', cn[ length(cn) ])
    df <- data.frame( set = gs,
                      tree.id = cn,
                      gene.id = gsub('.t[0123456789]+', '', t )
                     )
    return( df )
}

id.list <- lapply( ls , read.ids, gs = cmd.set )
id.df <- do.call( rbind, id.list )

gi <- get.genes( id.df$gene.id, mode = 'full' )

id.df <- merge( id.df, gi, by.x = 'gene.id', by.y = 'gene_id', all.x = TRUE, sort = FALSE)

write.table( id.df, file = cmd.file, row.names = FALSE, sep = '\t', quote = FALSE )
