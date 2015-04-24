#Rshell.sh Transdecoder_filtering.R -fasta.file -bed.file -blast.file -pfam.file -namebase.transdecoder.files -out.file

# variables

args = commandArgs(TRUE)
print(args)

#help options
if(length(grep("-h", args))>=1|length(args)<6||length(args)>6) {    
  cat("\nHelp for TD-transcript-filtering...\n")
  cat("\n*seqinr and data.table packages must be installed*\n")
  cat("USAGE:  Rsript scriptname <Trinity_fasta.file> <TD_bed.file> <BLASTP_concat.file> <TD_pfam.file> <Output_Dir>\n")
  cat("        gtf.path: <full path to gtf file>\n")
  cat("        workdir: <full path to directory containing TD-files>\n")
  cat("        gene_prefix: <Prefix to name the genes>\n\n") 
  quit(save="no")} 

if(length(args)==6) {
    require(seqinr);   
    require(data.table)
    cat('Arguments:\n')
    for(i in args){
        cat(i, '\n')
    }
}else{
    cat('Error in input arguments...')
    break()
}

#################################################################################
## READ VARIABLES ##
#################################################################################
  
fasta.file  <- args[1]
bed.file  <- args[2]
b.file <- args[3]
p.file <- args[4]
trin.file <- args[5]
out.file <- args[6]


#################################################################################
## FUNCTIONS ##
#################################################################################
cat('\n\nReading functions....\n')

output.table=function(table, table.name){
  write.table(table, quote=FALSE, col.names=TRUE, row.names=FALSE,
              sep='\t', file= table.name )
}

filter.longest=function(data, support=c(TRUE, FALSE), complete=c(TRUE, FALSE), only.longest=c(TRUE, FALSE)){
  dat  <-  data
  if(support==TRUE)      { dat  <- dat[! grepl('NO', dat$evidence), ] } 
  if(complete==TRUE)     { dat  <- dat[dat$TDclass %in% 'complete', ] }   
  if(only.longest==TRUE) { dat  <- dat[!duplicated(dat$LOC), ] } 
  cat(dim(dat))
  dat
}

require(data.table); require(seqinr)

##############################################
###### Starting to compute ###################
##############################################
fasta <- read.fasta( fasta.file, as.string = TRUE )

fasta <- sapply( fasta, function(x) getAnnot(x), USE.NAMES = FALSE)
fasta <- sapply( fasta, function(x) unlist(strsplit(x, ' '))[1])
fasta <- gsub('>', '', fasta)

trans <- data.frame('LOC' = sapply(fasta, function(x)
                        paste(unlist(strsplit(x, '_'))[1:2], collapse = '_') ) ,
                    'TRA' = fasta,
                    stringsAsFactors = FALSE
                    )
cat('FASTA-- genes:', length(unique(trans$LOC)),
    'transcripts:', length(unique(trans$TRA)), '\n')

#========================================================


b  <- read.table(b.file, sep='\t', stringsAsFactors = FALSE)
cat('BLAST:', dim(b), '\n')

p  <- read.table(p.file,sep ='\t', head=F, skip =3, quote = "", 
                 row.names = NULL, stringsAsFactors = FALSE)
cat('PFAM:', dim(p), '\n')

bed  <- read.table(bed.file, sep='\t', header=F,
                   stringsAsFactors=FALSE, fill=TRUE)
bed  <- bed[-grep('^track', bed$V1),]
cat('BED:', dim(bed), '\n')

# getting bed info:
splt.bed  <- strsplit(bed$V4, ';')
bed$TDmodel  <-  gsub('ID=', '', sapply(splt.bed, '[', 1))
bed$TDclass  <-  sapply(strsplit(sapply(splt.bed, '[', 3), 'type:|_len'), '[', 2)
bed$CDS.length  <- abs(bed$V8-bed$V7)
bed$TDframe  <- rep( NA, nrow(bed))
f1  <- seq(1, 1000000, by=3); bed$TDframe[bed$V7 %in% f1]  <- '1'
f2  <- seq(2, 1000000, by=3); bed$TDframe[bed$V7 %in% f2]  <- '2'
f3  <- seq(3, 1000000, by=3); bed$TDframe[bed$V7 %in% f3]  <- '3'

# sort bed and trim collumns
bed.sort  <- bed[order(bed$CDS.length, decreasing=T),c(1, 6, 13, 14, 15)]
colnames(bed.sort)  <- c('TRA', 'TDstrand', 'TDmodel', 'TDclass', 'CDS.length')

# merge gtf and bed
trans.bed  <- merge(x=trans, y=bed.sort, by='TRA', all.x=TRUE)
cat('Merger 1:', dim(trans.bed),'\n')
#############################################
### ADDING evidence COLLUMNS TO TD-MATRIX ###
#############################################

# adding blast results 
cat('Adding blast results...\n')
b  <- b[!duplicated(b$V1), ] # keeping only best hit per query

trans.bed$target  <- b$V2[match(trans.bed$TDmodel, b$V1)]
trans.bed$evalue  <- b$V13[match(trans.bed$TDmodel, b$V1)]
trans.bed$BLAST.hitlength  <- b$V4[match(trans.bed$TDmodel, b$V1)]

# adding pfam match
cat('Adding pfam results...\n')
trans.bed$PFAM  <- rep(NA, nrow(trans.bed))
trans.bed$PFAM[!is.na(match(trans.bed$TDmodel, p$V4))]  <- 'yes' 

## add new evidence collumn
trans.bed$evidence  <- rep(NA, nrow(trans.bed))
trans.bed$evidence[!is.na(trans.bed$target) & !is.na(trans.bed$PFAM)]  <- 'BLAST+PFAM'
trans.bed$evidence[is.na(trans.bed$target)  & !is.na(trans.bed$PFAM)]  <- 'PFAM'
trans.bed$evidence[!is.na(trans.bed$target) & is.na(trans.bed$PFAM)]  <- 'BLAST'
trans.bed$evidence[is.na(trans.bed$target)  & is.na(trans.bed$PFAM)]  <- 'NO BLAST OR PFAM'
trans.bed$evidence[is.na(trans.bed$TDmodel)]  <- 'NO TDmodel'
cat('Merger 2:', dim(trans.bed),'\n')
#output.table(trans.bed.filt, 'All_TD_models.tab')


# keeping only longest ORF model per trans
t <- trans.bed
t$evidenceN <- ifelse( grepl('NO',t$evidence ), 1, 2)
t <- t[!is.na(t$TDmodel), ] # remove seqs that have no TDmodel at all
#t      <- t[t$TDclass=='complete', ]
t   <- t[order(t$LOC, t$evidenceN, t$CDS.length, decreasing=T),]
t   <- t[!duplicated(t$TRA),]
cat( 'Merger 3:',dim(t),'\n' )

final.transcripts <- t[,1:10]

#final.transcripts  <- final.transcripts[, -c(3,4)]

##################################################################
# Her skal vi    fitrere                                         #
##################################################################

longest.support.dat  <- filter.longest(data=final.transcripts, support=TRUE,
                                       complete=FALSE, only.longest=TRUE)
longest.all.dat  <- filter.longest(data=final.transcripts, support=FALSE,
                                   complete=FALSE, only.longest=TRUE)

cat(paste('Locus number all:', nrow(longest.all.dat), '\n'),
    paste('Locus number with homology support:', nrow(longest.support.dat), '\n'))
   

output.table(longest.all.dat, 'Longest_All.annot.tab')
output.table(longest.support.dat, 'Longest_Support.annot.tab')


fasta.writer <- function(infile, outfile, filter){
    cat('Writing:', outfile ,'\n')
    fasta <- read.fasta( infile, as.string=TRUE )
    fasta <- fasta[names(fasta) %in% longest.all.dat[,filter]]
    write.fasta(fasta, names=names(fasta),
                file = outfile)
}

fasta.writer(infile = fasta.file, outfile = paste(out.file, '.fasta', sep=''),
             filter = 'TRA')
fasta.writer(infile = paste(trin.file, '.pep', sep=''),
             outfile = paste(out.file, '.pep', sep=''),
             filter = 'TDmodel')
fasta.writer(infile = paste(trin.file, '.cds', sep=''),
             outfile = paste(out.file, '.cds', sep=''),
             filter = 'TDmodel')







