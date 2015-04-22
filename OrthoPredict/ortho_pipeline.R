# USAGE EKSAMPLE

# Rscript ortho_pipeline.R ENSGALG00000008900 ensembl_gene_id ggallus_gene_ensembl test Cores

# The script will produce the following files (* as specified by argument 4 e.g. test/vitello )
# *.fa         = .fasta
# *.blast.xml  = blastp results (will be deleted)
# *.blast.tab  = blastp results
# *.spike.fa   = .fasta file with Ssa sequences idientified by blast
# *.spike.algn = mafft amino acid alignment file
# *.tree1.pdf  = phylogenetic tree figure (rooted, if a Lamprey sequence was found among the orthologous)
# *.tree1.tree = tree file
# *.ssa.ids    = 


# IMPORTANT load modules:
#  - blast+
#  - mafft
#  - python 

# Rscript needs (at least) 5 CPUs

#===============================================================================
#======= FUNCTIONS =============================================================
#===============================================================================

# Input:
# - ids = gene ids, as specified by filter (hgnc_symbols)...
# - filter = biomart accepted
# - dataset = Start biomart dataset 
get.ortho <- function(ids, filter = 'hgnc_symbol',  dataset = 'hsapiens_gene_ensembl' ){
    require( biomaRt )
    require( reshape )
    ensembl = useMart( "ensembl", dataset = dataset )
# it's only possible to retrieve max 6 orthologous, so split into 2
    attr.1 <- c('external_gene_name',
                'pformosa_homolog_ensembl_gene','pformosa_homolog_ensembl_peptide',    # Amazon Molly
                'gmorhua_homolog_ensembl_gene', 'gmorhua_homolog_ensembl_protein' ,    # Cod
                'amexicanus_homolog_ensembl_gene','amexicanus_homolog_ensembl_peptide',# Cave fish
                'trubripes_homolog_ensembl_gene','trubripes_homolog_ensembl_peptide',  # Fugu
                'olatipes_homolog_ensembl_gene','olatipes_homolog_ensembl_peptide',    # Medaka
                'oniloticus_homolog_ensembl_gene','oniloticus_homolog_ensembl_peptide' # Nile tillapia 
                )
    attr.2 <- c('external_gene_name',
                'xmaculatus_homolog_ensembl_gene','xmaculatus_homolog_ensembl_peptide',# Platyfish
                'loculatus_homolog_ensembl_gene','loculatus_homolog_ensembl_peptide',  # Spotted gar
                'gaculeatus_homolog_ensembl_gene','gaculeatus_homolog_ensembl_peptide',# Stickleback
                'tnigroviridis_homolog_ensembl_gene','tnigroviridis_homolog_ensembl_peptide', # Tetraodon
                'drerio_homolog_ensembl_gene','drerio_homolog_ensembl_peptide',         # Zebrafish
                'pmarinus_homolog_ensembl_gene','pmarinus_homolog_ensembl_peptide'      # Lamprey (OUTGROUP)
    )
    # Run BioMArt
    mart.1 <-getBM(attributes = attr.1,
                   filters = filter,
                   mart = ensembl,
                   values = ids )
    mart.2 <-getBM(attributes = attr.2,
                   filters = filter,
                   mart = ensembl,
                   values = ids )
    # shape data
    mm.1 <- melt(mart.1, id.vars = 'external_gene_name')
    mm.2 <- melt(mart.2, id.vars = 'external_gene_name')
    mm <- do.call(rbind, list(mm.1, mm.2))
    # FILTER
    mm <- mm[ ! (mm$value =='' | is.na(mm$value)), ]
    mm <- mm[ ! grepl('_gene', mm$variable), ]
    mm <- mm[ ! duplicated(mm), ]
    colnames(mm)[ grepl( 'external_gene_name', colnames(mm))] <- 'cluster_name'
    mm$species = sapply( as.character(mm$variable), function(x) unlist(strsplit(x, "_"))[1])
    # Finished
    return( mm )
}


# Function to get all sequences takes as input the output of get.ortho
get.fasta <- function( df ){
    # split DF to list of ensembl gene/protein IDs per species
    id.list <- split( df$value, df$species )
    # internal BioMart function to get sequences for Ids
    fetch.seq <- function( spec, ids ){
        ens.species <- paste( spec, 'gene', 'ensembl', sep ='_' )
        filter <- 'ensembl_peptide_id'
        ensembl = useMart( "ensembl", dataset = ens.species )
        # Query to BioMart
        martA <-getBM(attributes = c( 'peptide', 'ensembl_peptide_id' ),
                     filters = filter,
                     mart = ensembl,
                     values = ids )
        martB <-getBM(attributes = c( 'ensembl_peptide_id', 'hgnc_symbol', 'external_gene_name' ),
                     filters = filter,
                     mart = ensembl,
                     values = ids )
        mart <- merge( martA, martB, by = 'ensembl_peptide_id', all = TRUE, sort = FALSE )
        return( mart )
    }
    # Loop over species
    N = length( id.list )
    fasta.ls = vector( mode = 'list', length =N )
    names(fasta.ls) <- names( id.list )[1:2]
    for( i in 1:N ){
        fasta.ls[[i]] <- fetch.seq( names(id.list)[i], id.list[[i]] )
        cat(names(id.list)[i], " OK\n" )
    }
    fasta.ls <- do.call(rbind, fasta.ls)
    # merge with the input data frame
    out <- merge( df, fasta.ls, by.x = 'value', by.y = 'ensembl_peptide_id',
                 all = TRUE, sort = FALSE )
    return( out )
}


write.fasta <- function(x, file.n ){
    sink( file.n )
    for( i in 1:nrow(x)){
        header <- paste(">", x[i,4], "|", x[i,1], "|", x[i,7], "\n", sep ="")
        cat( header )
        cat( paste(x[i,5],"\n",sep="") )
    }
    sink()
}

# The all.transcripts function is thought for results when you blast against all
# possible peptides. 
process.blast <- function( file, all.transcripts = FALSE ){
    require( CigSsa.db )
    blast <- read.table( file = file, header = T, sep ='\t')
    blast$coverage <-  round(blast$align_length/blast$q.length, 4)
    # discard hits with less than 30% coverage
    blast <- blast[ blast$coverage > 0.3, ]
    if( all.transcripts ){
        idf <- function(x) unlist(strsplit(unlist(strsplit(x, '\\|'))[3], ' '))[1]
        blast$s.id2 <- sapply( blast$s.description, idf, USE.NAMES = FALSE)
    }else{
        blast$s.id2 <- gsub('.t[0123456789]+','', blast$s.id)
    }
    blast <- blast[! duplicated(blast[,which(colnames(blast) %in% c('q.id','s.id2'))]), ]
    gn <- get.genes( blast$s.id2 )
    out <-  merge(blast, gn, by.x = 's.id2', by.y = 'gene_id', all.x = T, sort = F)
    out <- out[, ! grepl( 'start|end|description', colnames(out))]
    # write output (Supressed)
    # write.table( out, file = gsub('.tab', '.filt.tab', file), row.names = FALSE, sep = '\t')
    return( unique(blast$s.id) )
}


# Function to add selected sequence information (ids) from .fasta file 2 into .fasta file 1
spike.seqs <- function( ids, fa.out, fa.f1, fa.f2 = '/mnt/users/fabig/Ssa_transcriptome/cigene3.6_chrom.blastDB/CIGSSA_Longest_HomologySupport.pep' ){
    require( ape )
    require( CigSsa.db )
    fa <- read.FASTA( fa.f1)
    db <- read.FASTA( fa.f2)
    db.sel <- db[ names(db) %in% ids ]
    names(db.sel) <- paste( 'ssalar|', names(db.sel), '|',
                           get.genes(gsub('.t[0123456789]+','', names(db.sel)), all=T)[,2],
                           sep ='')
    out <- c(fa, db.sel)
    write.dna( out, file = fa.out, format = 'fasta' )
}

mafft <- function( fa.in, out, show = FALSE, cores = 1 ){
    cmd <- paste('mafft --auto --thread', cores,
                 '--quiet', fa.in, '>', out, sep = ' ')
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

#options(echo=TRUE) 
args <- commandArgs(trailingOnly = TRUE)
print(args)

cmd.ids <- args[1]
cmd.filt <- args[2]
cmd.ds <- args[3]

cmd.out <- paste(args[4], '/', args[1],
                 c('.fa','.blast.xml','.blast.tab', '.spike.fa',
                   '.spike.algn', '.tree1.pdf', '.tree1.tree', '.ssa.ids'),
                 sep ='')

cmd.cores <- args[5]


system( paste('mkdir', args[4]) )

#if( grepl('.txt', cmd.ids) ){
#    cmd.ids <- read.table( file = cmd.ids, header = F )
#}

library( biomaRt )
library( reshape )




#---- Part A -------------------------------------------------------------------
# get all orthologous from fish
ortho <- get.ortho( cmd.ids, filter = cmd.filt,  dataset = cmd.ds)
# fetch protein sequences from ensembl
ortho.fa <- get.fasta( ortho )
# write .fasta
write.fasta( ortho.fa, cmd.out[1] )
 
#---- Part B -------------------------------------------------------------------
# BlastP
blastcmd <- paste('blastp -query', cmd.out[1],
                  '-db /mnt/users/fabig/Ssa_transcriptome/cigene3.6_chrom.blastDB/CIGSSA_Longest_HomologySupport.pep',
                  '-out', cmd.out[2],
                  '-num_threads', cmd.cores,
                  '-max_target_seqs 5 -evalue 1e-30 -outfmt 5',
                  sep = ' ')
cat( blastcmd )
system( blastcmd )

pythoncmd <- paste('python ~/script_python/blast_xml2tab.py --i', cmd.out[2],
                   '--o', cmd.out[3])
system( pythoncmd )
system( paste( 'rm', cmd.out[2] ))
#---- Part C -------------------------------------------------------------------

# spike in ssa sequences
sel <- process.blast( cmd.out[3] )
spike.seqs( ids = sel,
            fa.out = cmd.out[4],
            fa.f1 = cmd.out[1] )

# remove the original .fa file
system( paste( 'rm', cmd.out[1]) )

#---- Part D -------------------------------------------------------------------
# run MAFFT and produce the 1st 'try and error' alignment: Need to inspect and
# remove outliers manually !!!
mafft( cmd.out[4], cmd.out[5], cores = cmd.cores )

# make the tree
library(phangorn)
algn <- read.phyDat( cmd.out[5], type = 'AA', format = 'fasta')
if( length(algn) <= 5 ){
    file.remove( cmd.out[1:5])
    stop( 'Too few Fish Sequences to make a tree')
}else{
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
    mo.log <- capture.output({
        m.opt <- optim.pml(m0,  model = "WAG", optNni = TRUE, optInv = TRUE)
    })

    # NO bootstrapping in the first round...
    #bs.log <- capture.output({
    #    bs = bootstrap.pml(m.opt, bs = 100, optNni = TRUE, multicore = TRUE  )
    #})
    
    bs.tree <- m.opt$tree
    bs.tree <- zero.el(bs.tree)
    #bs.tree$node.label <- prop.clades(m.opt$tree, bs)
    #write.tree(bs.tree, file = cmd.out[7])

    col <- ifelse(grepl("CIGSSA", bs.tree$tip.label), 'white', 'lightgrey')

    pdf( file = cmd.out[6], width = 6, height = 8)
    plot( bs.tree , cex = 0.7, no.margin = TRUE, show.node.label = FALSE,
         tip.color = "white")
    tiplabels( bs.tree$tip.label, adj = c(-0.05, 0.5), cex = 0.7,
              bg = col )
    dev.off()

    # write *.ssa.id
    ids <- parse.tips(tree1$tip.label)
    ids <- ids[ grepl('CIGSSA', ids) ]
    write.table( ids, file = cmd.out[8], row.names = F, col.names = F)
}
