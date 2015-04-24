#Rshell.sh TD_transcript_filtering.R -gtf.file -work.dir -gene_prefix

# variables

#gtf.file  <- '/mnt/SALMON-SEQDATA/aquagenome/Ssa/Synteny/Harr/CUFF/Harr/merged_asm/merged.gtf'
#out.dir  <- '/mnt/users/srsand/Salmon/Gene_models/harr/'
#gene.prefix  <- 'CIGTt'


args = commandArgs(TRUE)
print(args)

#help options
if(length(grep("-h", args))>=1|length(args)<3||length(args)>3) {    
  cat("\nHelp for TD-transcript-filtering...\n")
  cat("\n*seqinr and data.table packages must be installed*\n")
  cat("USAGE:  Rsript scriptname <Cufflink_gtf.file> <TD_gff.file> <BLASTP_concat.file> <TD_pfam.file> <Output_Dir>\n")
  cat("        gtf.path: <full path to gtf file>\n")
  cat("        workdir: <full path to directory containing TD-files>\n")
  cat("        gene_prefix: <Prefix to name the genes>\n\n") 
  quit(save="no")} 

if(length(args)==3) {
  require(seqinr);   
  require(data.table)
  cat('Arguments:\n')
  cat(args)
  } else { cat('Error in input arguments...'); break()}

###############
## FUNCTIONS ##
###############
cat('\n\nReading functions....\n')

output.table=function(table, table.name){
  write.table(table, quote=F, col.names=T, row.names=F, sep='\t', file=paste(out.dir, paste(gene.prefix, table.name, sep='_'), sep=''))
}

filter.longest=function(data, support=c(TRUE, FALSE), complete=c(TRUE, FALSE), only.longest=c(TRUE, FALSE)){
  dat  <-  data
  if(support==TRUE)      { dat  <- dat[!(dat$evidence=='No'), ] } 
  if(complete==TRUE)     { dat  <- dat[dat$TDclass %in% 'complete', ] }   
  if(only.longest==TRUE) { dat  <- dat[!duplicated(dat$LOC), ] } 
  cat(dim(dat))
  dat
}

# Longest support
# Longest all --> filtered on 

output.fasta=function(filter.table=longest.support.dat, file.name.tag){
  
  fa.pep.all  <- td.pep
  names(fa.pep.all)  <- gsub('cds.', '', names(fa.pep.all))
  fa.pep  <- fa.pep.all[names(fa.pep.all) %in% filter.table$TDmodel]
  if(!length(fa.pep)==nrow(filter.table)) { cat('sequence file error'); break()}
  
  fa.cds.all  <- td.cds
  names(fa.cds.all)  <- gsub('cds.', '', names(fa.cds.all))
  fa.cds  <- fa.cds.all[names(fa.cds.all) %in% filter.table$TDmodel]
  if(!length(fa.cds)==nrow(filter.table)) { cat('sequence file error'); break()}
  
  names(fa.pep)  <-  filter.table$CIGnames[match(names(fa.pep), filter.table$TDmodel)]
  names(fa.cds)  <-  filter.table$CIGnames[match(names(fa.cds), filter.table$TDmodel)]
  write.fasta(fa.pep, names=names(fa.pep), file=paste(out.dir, paste(gene.prefix, paste(file.name.tag, '.pep', sep=''), sep='_'), sep=''))
  write.fasta(fa.cds, names=names(fa.cds), file=paste(out.dir, paste(gene.prefix, paste(file.name.tag, '.cds', sep=''), sep='_'), sep=''))
  
}


require(data.table); require(seqinr)

##############################################
###### Starting to compute ###################
##############################################

cat('\n...Making evidence table...\n')
  
#Run script
gtf.file  <- args[1]
out.dir  <- args[2]
gene.prefix  <- args[3]

cat('making file names....')
b.file    <- paste(out.dir, 'BLASTP_concatenated.tab', sep='')
p.file    <- paste(out.dir,'TD_concatenated_FIXCOLSEP.pfam', sep='')
bed.file  <- paste(out.dir,'TD_concatenated.bed' , sep='')

# read in data
cat('reading in datafiles\n *reading GTF file')
gtf  <- fread(gtf.file)
gtf  <-  as.data.frame(gtf)
split.col9  <-  strsplit(as.character(gtf$V9), '\\;')
gtf$LOC    <- gsub('gene_id |\\"', '', sapply(split.col9 , '[', 1))
gtf$TRA    <- gsub(' transcript_id |\\"', '', sapply(split.col9 , '[', 2))
#head(gtf);dim(gtf)

# reducing gtf file...
trans  <- gtf[!duplicated(gtf$TRA), c(1, 7, 10, 11)]
names(trans)  <- c('scaff', 'ref.strand', 'LOC', 'TRA')
#head(trans); dim(trans)


cat('*reading BLAST-output\n')
b  <- read.table(b.file, sep='\t')

cat('*reading PFAM-output\n')
p  <- read.csv(p.file, sep='\t', head=F)

cat('*reading bed-output\n')
bed  <- read.table(bed.file, sep='\t', header=F, stringsAsFactors=F, fill=T)
bed  <- bed[-grep('^track', bed$V1),]

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
trans.bed  <- merge(x=trans, y=bed.sort, by='TRA', all=FALSE)
trans.bed$CIG.LOC  <- gsub('XLOC', gene.prefix, trans.bed$LOC)
trans.bed$ref.strand  <- as.character(trans.bed$ref.strand) 

# check this:CIGSSA_119114.t1, XLOC_119114
#head(trans.bed)
#trans.bed[trans.bed$LOC %in% 'XLOC_119114' , ]
# dim(trans.bed[trans.bed$LOC %in% 'XLOC_119114' , ])
#gtf[gtf$LOC %in% 'XLOC_119114', ]
# bed.sort[bed.sort$TRA %in% 'TCONS_00287837', ]
#  TCONS  <-  trans.bed$TDmodel[trans.bed$LOC %in% 'XLOC_119114']
#bcmo.pep = read.fasta('/mnt/users/srsand/Salmon/Gene_models/Ssa_megamerge_v3/TD_concatenated.pep', as.string=T)
#getLength(bcmo.pep[paste('cds.', TCONS, sep='')])

# remove strand eror: keeps + +/- AND - + AND .
trans.bed.filt  <- trans.bed[trans.bed$ref.strand=='+' &  trans.bed$TDstrand == '+' | trans.bed$ref.strand=='-' &  trans.bed$TDstrand == '+' | trans.bed$ref.strand=='.', ]

#trans.bed.filt[trans.bed.filt$LOC %in% 'XLOC_119114' , ]
# dim(trans.bed.filt[trans.bed.filt$LOC %in% 'XLOC_119114' , ])

#############################################
### ADDING evidence COLLUMNS TO TD-MATRIX ###
#############################################

# adding blast results 
cat('Adding blast results...')
b  <- b[!duplicated(b$V1), ] # keeping only best hit per query
b$V1  <- gsub('cds.', '', b$V1)
trans.bed.filt$target  <- b$V2[match(trans.bed.filt$TDmodel, b$V1)]
trans.bed.filt$evalue  <- b$V13[match(trans.bed.filt$TDmodel, b$V1)]
trans.bed.filt$BLAST.hitlength  <- b$V4[match(trans.bed.filt$TDmodel, b$V1)]

# adding pfam match
cat('Adding pfam results...')
trans.bed.filt$PFAM  <- rep(NA, nrow(trans.bed.filt))
trans.bed.filt$PFAM[!is.na(match(trans.bed.filt$TDmodel, p$V4))]  <- 'yes' 

## add new evidence collumn
trans.bed.filt$evidence  <- rep(NA, nrow(trans.bed.filt))
trans.bed.filt$evidence[!is.na(trans.bed.filt$target) & !is.na(trans.bed.filt$PFAM)]  <- 'BLAST+PFAM'
trans.bed.filt$evidence[is.na(trans.bed.filt$target)  & !is.na(trans.bed.filt$PFAM)]  <- 'PFAM'
trans.bed.filt$evidence[!is.na(trans.bed.filt$target) & is.na(trans.bed.filt$PFAM)]  <- 'BLAST'
trans.bed.filt$evidence[is.na(trans.bed.filt$target)  & is.na(trans.bed.filt$PFAM)]  <- 'No'

output.table(trans.bed.filt, 'All_TD_models.tab')


# keeping only longest ORF model per trans
t  <- trans.bed.filt
#t      <- t[t$TDclass=='complete', ]
t   <- t[order(t$CDS.length, decreasing=T),]
t   <- t[!duplicated(t$TRA),]

# adding SsaCIG names and transcript number...
final.transcripts  <- data.table(t)[ , CIGnames:=paste(LOC, '.t', seq(1:nrow(.SD)), sep=''), by=LOC]
final.transcripts$CIGnames  <- gsub('XLOC', gene.prefix, final.transcripts$CIGnames)


# adding ref min and max coordinates for LOCUS and TRANSCRIPTS
gtf.filt  <- gtf[gtf$LOC %in% final.transcripts$LOC, ]
gtf.loc.region  <- as.data.table(gtf.filt)[ , list(ref.start=min(.SD$V4), ref.stop=max(.SD$V5)), by=LOC]
final.transcripts$LOCstart.ref  <- gtf.loc.region$ref.start[match(final.transcripts$LOC, gtf.loc.region$LOC)]
final.transcripts$LOCstop.ref   <- gtf.loc.region$ref.stop[match(final.transcripts$LOC, gtf.loc.region$LOC)]

gtf.trans.region  <- as.data.table(gtf.filt)[ , list(ref.start=min(.SD$V4), ref.stop=max(.SD$V5)), by=TRA]
final.transcripts$TRAstart.ref  <- gtf.trans.region$ref.start[match(final.transcripts$TRA, gtf.trans.region$TRA)]
final.transcripts$TRAstop.ref   <- gtf.trans.region$ref.stop[ match(final.transcripts$TRA, gtf.trans.region$TRA)]
final.transcripts  <- as.data.frame(final.transcripts)

#final.transcripts  <- final.transcripts[, -c(3,4)]

##################################################################
# Her skal vi    fitrere                                         #
##################################################################

all.dat  <- filter.longest(data=final.transcripts, support=FALSE, complete=FALSE, only.longest=FALSE)
all.support.dat  <- filter.longest(data=final.transcripts, support=TRUE, complete=FALSE, only.longest=FALSE)
longest.support.dat  <- filter.longest(data=final.transcripts, support=TRUE, complete=FALSE, only.longest=TRUE)
longest.all.dat  <- filter.longest(data=final.transcripts, support=FALSE, complete=FALSE, only.longest=TRUE)
longest.complete.support  <- filter.longest(data=final.transcripts, support=TRUE, complete=TRUE, only.longest=TRUE)
longest.complete.all  <- filter.longest(data=final.transcripts, support=FALSE, complete=TRUE, only.longest=TRUE)

cat(
  paste('Total number of transcripts:', nrow(all.dat), '\n'),
  paste('Total number of transcripts with homology support:',nrow(all.support.dat), '\n'),
  paste('Locus number all:', nrow(longest.all.dat), '\n'),
  paste('Locus number with homology support:', nrow(longest.support.dat), '\n'),
  paste('Locus number complete genes:', nrow(longest.complete.all), '\n'),
  paste('Locus number complete genes with homology support:', nrow(longest.complete.support), '\n'), 
  file=paste(out.dir, paste(gene.prefix, 'Transcript_and_gene_numbers.txt', sep='_')))

output.table(all.dat, 'All.annot.tab')
output.table(longest.all.dat, 'Longest_All.annot.tab')
output.table(longest.support.dat, 'Longest_Support.annot.tab')
output.table(longest.complete.support, 'Longest_Complete_Support.annot.tab')
output.table(longest.complete.all, 'Longest_All_Complete.annot.tab')

# print out evidence plot
png(paste(out.dir, 'Evidence_distribution_Transcripts.png', sep=''), units='cm', height=10, width=10, res=300, bg='transparent')
barplot(table(final.transcripts$evidence), cex.axis=0.5, cex.names=0.5) # per transcript
dev.off()

png(paste(out.dir, 'Evidence_distribution_Longest.png', sep=''), units='cm', height=10, width=10, res=300, bg='transparent')
barplot(table(longest.all.dat$evidence), cex.axis=0.5, cex.names=0.5) # per transcript
dev.off()

cat('\nReading in fasta files...might take some time....go get coffee!\n')
# extract fasta...
require(seqinr)
td.pep  <- read.fasta(paste(out.dir, 'TD_concatenated.pep', sep=''), as.string=T)
td.cds  <- read.fasta(paste(out.dir, 'TD_concatenated.cds', sep=''), as.string=T)

output.fasta(filter.table=longest.support.dat, file.name.tag='Longest_HomologySupport')

###### checking...



#################### NOT IN USE BELOW ###################################

#t.longestmodel[grep("ElCIG021149", t.longestmodel$CIGnames),]
#t.longestmodel[grep("ElCIG021988", t.longestmodel$CIGnames),]
#longest.support.dat[grep("ElCIG021988", longest.support.dat$CIGnames),]

### extra

# trans.bed.filt.sort  <- trans.bed.filt[order(trans.bed.filt$TRA.length, decreasing=T),]
# head(trans.bed.filt.sort)
# t.longestmodel  <- trans.bed.filt.sort[!duplicated(trans.bed.filt.sort$TRA),]
# dim(t.longestmodel)
# 
# # tst
# g.longestmodel  <- t.longestmodel[!duplicated(t.longestmodel$CIG.LOC),]
# dim(g.longestmodel)




# ##################################################
# #  identifying multi-prot hit ORFs in TD output  #
# ##################################################
# 
# require(data.table)
# test.dt  <- data.table(t.cds)
# #str(test.dt)
# 
# # identifying multi-protein hit transcripts
# cat('Identifying multi-protein hit transcripts...\n')
# test.dt$LOC  <- as.factor(test.dt$LOC) 
# test.dt  <- test.dt[!is.na(test.dt$target), ]
# test  <- test.dt[, list(Prot.hits=sum(!duplicated(target)), LOC=LOC[1]),by=V1]
# multihit.LOC  <- test[, max(Prot.hits),by=LOC]
# multihit.LOC.c  <- multihit.LOC$LOC[multihit.LOC$V1>1]
# test.dt.multihit  <- test.dt[test.dt$LOC %in% multihit.LOC.c ,]
# test.dt.multihit  <- test.dt.multihit[order(test.dt.multihit$evalue, decreasing=F) , ]
# 
# multi.frame.per.TR  <- test.dt.multihit[ , { mydata  <- copy(.SD)
#                                           frame  <- mydata$frame[!duplicated(target)]
#                                           list(frames=sum(!duplicated(frame)))},
#                                   by=V1]
# 
# ###################################################################
# #   gj??r frame shift effect i alle IKKE-multihit transcripter.... #
# ###################################################################
# 
# cat('Identifying Frame-Shift transcripts in single-protein hit transcript...\n')
# single.hit  <- test.dt[!test.dt$LOC %in% test.dt.multihit$LOC,]
# single.frame.per.TR  <- single.hit[ , list(frames=sum(!duplicated(frame)), LOC=LOC[1]), by=V1]
# single.frame.per.LOC  <- single.frame.per.TR[, list(frames=max(frames)), by=LOC]
# 
# 
# ###############################################
# ## Filter TD-table...and giving output files ##
# ###############################################
# 
# ## klassifisere og filtrere TRANS
# # 1) FS+med multi proteinhit (sannsynligvis tull) --> flagge
# # 2) FS innen samme protein --> flagge
# # 3) TD class (complete etc) --> flagge
# # 4) Velge lengste orf-model per transcript
# 
# 
# ##--> flagg FS+multiple proteiner
# remove.FS.diff.prot  <- multi.frame.per.TR$V1[multi.frame.per.TR$frames>1]
# t.cds$multiprotFS   <- rep(FALSE, nrow(t.cds))
# t.cds$multiprotFS[t.cds$V1 %in% remove.FS.diff.prot]  <- TRUE
# 
# 
# ##--> FLAGGING TRANS with FS+single proteiner
# flag.FS.same.prot  <- single.frame.per.TR$V1[single.frame.per.TR$frames>1]
# t.cds$singleprotFS  <- rep(FALSE,nrow(t.cds))
# t.cds$singleprotFS[t.cds$V1 %in% flag.FS.same.prot]  <- TRUE


# ## --> class TD ORF
# splt.bed  <- strsplit(bed[-1], '\t|;')
# bed.model  <-  gsub('ID=', 'cds.', sapply(splt.bed, '[', 4))
# bed.class  <-  sapply(strsplit(sapply(splt.bed, '[', 6), 'type:|_len'), '[', 2)
# t.cds$TDclass  <- rep(NA, nrow(t.cds))
# t.cds$TDclass  <- bed.class[match(t.cds$TRmodel, bed.model)]
# 
# # remove strand eror
# t.cds$ref.strand <-  ref.strand.info$V7[match(as.character(t.cds$V1), ref.strand.info$TRA)]
# head(t.cds)
# 
# dim(t.cds)
# t.cds.strandfilt  <- t.cds[!(t.cds$ref.strand=='+' &  t.cds$V7 == '-' | t.cds$ref.strand=='-' &  t.cds$V7 == '-'), ]
# dim(t.cds.strandfilt)

