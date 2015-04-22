
# Rscript ortho_pipeline-ReRun.R ACO2 test



#===============================================================================
#======= FUNCTIONS =============================================================
#===============================================================================

# Function to remove sequnces from .fasta file after they have been identified to not belog in the tree
# ids is a vector with the ids to keep
remove.seqs <- function( ids, fa.in, fa.out ){
    require( ape )
    fa <- read.FASTA( fa.in )
    s.name <- sapply(names(fa), function(x) unlist(strsplit(x, '\\|'))[2], USE.NAMES = F)
    fa <- fa[ (!grepl('CIGSSA', s.name)) | (s.name %in% ids) ]
    write.dna( fa, file = fa.out, format = 'fasta'  )
    t = table(s.name %in% ids)
    cat( 'Seqs remaining: ',t[1], ' removed: ', t[2], '\n' )
}

mafft <- function( fa.in, out, show = FALSE ){
    cmd <- paste('mafft --auto --thread 5 --quiet', fa.in, '>', out, sep = ' ')
    system( cmd )
}

# function to set negative edge.length to 0
zero.el <- function(x){
    x$edge.length[x$edge.length < 0 ] <- 0
    return( x)
}

parse.tips <- function(x, part = 2){
    out <- sapply(x, function(y) unlist(strsplit(y, '\\|'))[part], USE.NAMES = F)
    return( out )
}

# function from phangorn includes a small fix
prop.clades <- function (phy, obj, part = NULL, rooted = FALSE) 
{
    if (is.null(part)) {
        if (length(obj) == 1 && class(obj[[1]]) != "phylo") 
            obj <- unlist(obj, recursive = FALSE)
        part <- prop.part(obj, check.labels = TRUE)
    }
    if (!identical(phy$tip.label, attr(part, "labels"))) {
        i <- match(phy$tip.label, attr(part, "labels"))
        j <- match(seq_len(Ntip(phy)), phy$edge[, 2])
        phy$edge[j, 2] <- i
        phy$tip.label <- attr(part, "labels")
    }
    bp <- prop.part(phy)
    if (!rooted) {
        bp <- postprocess.prop.part(bp)
        part <- postprocess.prop.part(part)
    }
    n <- numeric(phy$Nnode)
    for (i in seq_along(bp)) {
        for (j in seq_along(part)) {
            if (identical(bp[[i]], part[[j]])) {
                n[i] <- attr(part, "number")[j]
                done <- TRUE
                break
            }
        }
    }
    n
}

#===============================================================================
#======= SCRIPT ================================================================
#===============================================================================

# Read command line options

options(echo=TRUE) 
args <- commandArgs(trailingOnly = TRUE)
print(args)

cmd.ids <- args[1]

cmd.out <- paste(args[2], '/', args[1],
                 c( '.spike.fa', '.spike.algn', '.tree1.pdf', '.tree1.tree', '.ssa.ids'),
                 sep ='')

system( paste('mkdir', args[4]) )

#if( grepl('.txt', cmd.ids) ){
#    cmd.ids <- read.table( file = cmd.ids, header = F )
#}

library( biomaRt )
library( reshape )


#---- Part A -------------------------------------------------------------------
# read in the Ssa gene ids to keep...

ids <- read.table(cmd.out[5], header = F )[,1]

# remove the sequences from *.spike.fa
remove.seqs( ids = ids, fa.in = cmd.out[1], fa.out = cmd.out[1])

# ReRun MAFFT 
mafft( cmd.out[1], cmd.out[2] )

#---- Part D -------------------------------------------------------------------

# remake the tree
library(phangorn)
algn <- read.phyDat( cmd.out[2], type = 'AA', format = 'fasta')
dm = dist.ml(algn, model = 'WAG')
tree1 <- nj(dm)
tree1 <- zero.el(tree1)
if( any( grepl('pmarinus', tree1$tip.label) ) ){
    rootseq <-  tree1$tip.label[grepl('pmarinus', tree1$tip.label)]
    if( length(rootseq) > 1)
        tree1 <- root(tree1, rootseq[1] )
    else
        tree1 <- root(tree1, rootseq )
}

m0 = pml(tree1, algn)
m.opt <- optim.pml(m0,  model = "WAG", optNni = TRUE, optInv = TRUE)

bs.log <- capture.output({
    bs = bootstrap.pml(m.opt, bs = 100, optNni = TRUE, multicore = TRUE  )
})

bs.tree <- m.opt$tree
bs.tree <- zero.el(bs.tree)
bs.tree$node.label <- prop.clades(m.opt$tree, bs)
#write.tree(bs.tree, file = cmd.out[4])

col <- ifelse(grepl("CIGSSA", bs.tree$tip.label), 'white', 'lightgrey')

pdf( file = cmd.out[3], width = 6, height = 8)
plot( bs.tree , cex = 0.7, no.margin = TRUE, show.node.label = TRUE,
     tip.color = "white")
tiplabels( bs.tree$tip.label, adj = c(-0.05, 0.5), cex = 0.7,
           bg = col )
dev.off()
