#options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)


infile <- args[1]     # input file (BLASTP hits)
outfile <- args[2]    # name of the output file
dat.file <- args[3]   # name of the "MASTER DB" file
ens.mart <- args[4]   # MART name
ens.arg1 <- as.character(args[5])
ens.arg2 <- as.character(args[6])

library(biomaRt)
library(plyr)
#-------------------------------------------------------------------------------
# read in the BLAST hit file
dat <- read.delim( file = infile, header = F, stringsAsFactors = FALSE)
# BLAST header
header <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart','send','eval', 'bitscore', 'qlen',
            'slen','nident' )
colnames(dat) <- header

#-------------------------------------------------------------------------------
# identify HIT source
ens.inf <- read.table(dat.file, header = TRUE, sep ='\t',stringsAsFactors = FALSE )
# add row to identify swissprot hits
ens.inf <- rbind(ens.inf, c(NA,NA,NA, 'swiss','gi' ))

ident.db <- function(x){
    ri <- which(sapply(ens.inf$prefix, function(y) grepl(paste('^',y,'\\|',sep =''), x),
                       USE.NAMES = FALSE)
                )
    return(ens.inf[ri,'mart'])
}

# setup new column indicating hit source
dat$db <- sapply(dat$sseqid, ident.db, USE.NAMES = FALSE)
# CHECK if all sources could be identified
if( any(dat$db == "character(0)")){
    stop("Some blast hits could not be identifyed\nCheck the input data.table!")
}

#-------------------------------------------------------------------------------
# fun to correct the uniprot protein ids
parse.id <- function(x){
    if(grepl('^gi\\|',x)){
        x <- unlist(strsplit(x, split = '\\|' ))
        if( length(x) == 3 ){
            out <- x[2] 
        }
        else{
            out <- gsub('\\.[0123456789]+', '', x[4])
        }        
        return( out )
    }else{
        x <- unlist(strsplit(x, split = '\\|' ))
        return(x[2])
    }
}

dat$sid <- sapply(dat$sseqid, parse.id, USE.NAMES = FALSE)
dat$coverage <- dat$length / dat$slen

id.list <- split(dat$sid, dat$db)
cat('-----------------------------------------------------------\n')
cat('Sending IDs to BioMart...\n')

mart.fun <- function(id, ds, ens.mart){
    uni.id <- unique(id)
    uni.seq <- seq(0, length(uni.id), by = 10000 )
    uni.seq <- c(uni.seq, length(uni.id))
    mart <- list()
    # Loop in chunks of 10000
    for(i in 1:(length(uni.seq)-1)){
        st = uni.seq[i]+1
        en = uni.seq[i+1]
        vals = uni.id[st:en]
        if( ds == 'swiss'){   
            uniprot = useMart("unimart",dataset="uniprot")
            mart[[i]] <- getBM(attributes=c('accession',
                               'gene_name','protein_name','protein_evidence'),
                          filters = 'accession',
                          mart = uniprot,
                          values = vals)
        }else{    
            # step 1: retrieve ensembl hits
            ensembl <- useMart(ens.mart ,dataset = ds)
            hits <- getBM(attributes=c('ensembl_peptide_id', ens.arg1, ens.arg2,
                              'uniprot_sptrembl'),
                          filters = 'ensembl_peptide_id',
                          mart = ensembl,
                          values = vals)
            # step 2: get evidence from uniprot
            uniprot = useMart("unimart",dataset="uniprot")
            hits2 <- getBM(attributes=c('accession','protein_evidence'),
                           filters = 'accession',
                           mart = uniprot,
                           values = unique(hits$uniprot_sptrembl))
            # step 3: merge information
            hits.f <- merge(hits, hits2, by.x = 'uniprot_sptrembl', by.y = 'accession',
                            all.x = TRUE, sort = FALSE) 
            # order so that in case of duplicated protein hits, the one with better supprort is on top
            hits.f <- hits.f[ order(hits.f$ensembl_peptide_id, hits.f$protein_evidence),]
            # remove duplicates
            hits.f <- hits.f[ ! duplicated( hits.f$ensembl_peptide_id ),]
            # clean the descriptions; remove everything in brakets
            hits.f$description <- gsub(" \\[.*\\]", "", hits.f$description)
            mart[[i]] <- hits.f
        }
        #cut(paste('chunk', i))
    }
    # join list
    mart <- ldply(mart)
    mart$origin <- ds
    cat('Found', (nrow(mart)/uni.seq) * 100, '%\n')
    # arrange an name the columns
    if( ds == 'swiss'){
        colnames(mart) <- c('uniprot_acc', 'gene_symbol', 'gene_name', 'protein_evidence','origin')
        mart$ensembl_peptide_id <- ''
    }
    else{
        colnames(mart) <- c('uniprot_acc', 'ensembl_peptide_id', 'gene_name', 'gene_symbol','protein_evidence','origin')
    }
    mart <- mart[,c('uniprot_acc','ensembl_peptide_id','gene_name','gene_symbol','protein_evidence','origin')]
    return(mart)   
}


mart.list <- list()
for( i in 1:length(id.list)){
    mart.list[[i]] <- mart.fun(id.list[[i]], names(id.list)[i], ens.mart = ens.mart)
}

names(mart.list) <- names(id.list)
mart.df <- ldply(mart.list)

mart.df$sid <- ''
mart.df$sid <- mart.df$ensembl_peptide_id
mart.df$sid[mart.df$ensembl_peptide_id == ''] <- mart.df$uniprot_acc[mart.df$ensembl_peptide_id == '']
# Set missing protein_evidences to:level 4: Predicted
mart.df$protein_evidence[is.na(mart.df$protein_evidence)] <- '4: Predicted'

#-------------------------------------------------------------------------------
dat <- dat[,c('qseqid','sid','pident','nident','coverage', 'eval' )]
dat2 <- merge(dat, mart.df, by = 'sid', all.x = TRUE, sort = F)
dat2 <- dat2[order(dat2$qseqid, dat2$eval),]

dat2.split <- dlply( dat2, .variable = 'qseqid' )
#-------------------------------------------------------------------------------


# fun to filter out the [best+most informative+highest evidence] hit per gene 
filter.fun <- function(x){
    x <- x[! is.na(x$protein_evidence),]
    upe <- levels(as.factor(x$protein_evidence))
    p.rank <- list()
    for(i in seq_along(upe)){
        p.rank[[i]] <- which(grepl( upe[i], x$protein_evidence))
    }
    # identify uniformative names
    rank.l1 <- lapply(p.rank, function(y) y[! y %in% which(
        grepl('Uncharacterized protein|genomic scaffold|like|whole genome shotgun sequence',
              x$protein_name))] )
    rank.l2 <- lapply(p.rank, function(y) y[! y %in% which(
        grepl('Uncharacterized protein|genomic scaffold|whole genome shotgun sequence', x$protein_name))] ) 
    # best informative hit information
    out <- NULL
    co = x$coverage
    i1 = 1
    for( i in rev(rank.l1)){
        if( length(i) == 0){
            next
        }
        else{
            i2 = min(i)
            if( (co[i2]+0.2 >= co[i1]) | (co[i1] >= 0.7 & co[i2] >= 0.7) ){
                out <- cbind(x[i2,],i2)
            }
            else{
                out <- cbind(x[i2,],i2)
            }
        }
    }
    if(is.null(out)){
        for( i in rev(rank.l2)){
            if( length(i) == 0){
                next
            }
            else{
                i2 = min(i)
                if( (co[i2]+0.2 >= co[i1]) | (co[i1] >= 0.7 & co[i2] >= 0.7) ){
                    out <- cbind(x[i2,],i2)
                }
                else{
                    out <- cbind(x[i2,],i2)
                }
            }
        }       
    }
    if(is.null(out)){
         for( i in rev(p.rank)){
            if( length(i) == 0){
                next
            }
            else{
                i2 = min(i)
                if( (co[i2]+0.2 >= co[i1]) | (co[i1] >= 0.7 & co[i2] >= 0.7) ){
                    out <- cbind(x[i2,],i2)
                }
                else{
                    out <- cbind(x[i2,],i2)
                }
            }
        }
     }
     return(out)
}

# Run the function
tt <- ldply( dat2.split, .fun = filter.fun )

#-------------------------------------------------------------------------------
prefix.fun <- function(x){
    if(x['coverage'] >= 0.7 & grepl('1',x['protein_evidence']) & x['pident'] >= 50 ){
       out = 'conserved'
    }else if(x['coverage'] >= 0.7 & grepl('2|3|4',x['protein_evidence'])){
       out = 'conserved hypothetical'
    }else if(x['coverage'] < 0.7 & x['coverage'] >= 0.2  ){
       out = 'hypothetical'
    }else if(x['coverage'] < 0.2 ){
       out = 'hypothetical, weakly similar'
    }else{
       out = 'predicted protein'
    }
    return(out)
}

tt$prefix <- apply( tt, 1, FUN = prefix.fun)
#-------------------------------------------------------------------------------
tt$qseqid <- gsub('\\.t[0123456789]+', '', tt$qseqid )
tt$gene_name <- tolower(tt$gene_name)

tt$gene_name[grepl('genomic scaffold', tt$gene_name) | nchar(tt$gene_name) == 0 ] <- 'uncharacterized protein'
tt$gene_name <- paste(tt$gene_name, tt$prefix, sep = '|')

tt.final <- tt[,c('qseqid','uniprot_acc','ensembl_peptide_id','gene_name','gene_symbol', 'eval','coverage','origin')]
tt.final$origin <- gsub('_eg_gene','', tt.final$origin )

write.table( tt.final, file = outfile, row.names = FALSE, sep = '\t')

cat('-----------------------------------------------------------\n')
l1 <- length(tt.final$qseqid)
l2 <- length(unique(dat$qseqid))
cat(paste('Output: ',round((l1/l2)*100,1),'% : ',l1, ' of ', l2,'\n'))
cat('-----------------------------------------------------------\n')
