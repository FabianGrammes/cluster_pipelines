# Rscript to collect transdecoder results into a SQLite DB
      
# INPUT variables:

# 0: Annotation data
# 1: merged.gtf 
# 2: TD_concatenated.bed (Transdecoder .bed file)
# 3: TD_concatenated.pfam (Transdecoder .pfam file)
# 4: BLASTP_concatenated.tab (BLAST P hit file; keep column order in mind)
# 5: topBlast.txt (from Blast2GO exported)
# 6: GO annotation (from Blast2GO exported)
# 7: SNP
# 8: Transpososns
# 9: KEGG annotation Gene id 2 K(from KAAS)
# 10: KEGG K 2 path ko
# 11: KEGG path ko 2 name
# 12: DB name (f.eks SsaCIG_annot.001.db )


# Example:
# $Rscript makeDB.R merged.gtf TD_concatenated.bed TD_concatenated.pfam BLASTP_concatenated.tab Test.db 
      
#---------------------------------------------------------------------
# read arguments
args <- commandArgs(TRUE)
ann.file <- args[1]
gtf.file <- args[2]
bed.file <- args[3]
pfam.file <- args[4]
blast.go.file <- args[5]
go.file <- args[6]
snp.file <- args[7]
transpos <- args[8]
knr.file <- args[9]
kpath.file <- args[10]
kname.file <- args[11]
db.name <- args[12]

# Check if the files exist, STOP if not:
for( i in 1:11){
    if(file.exists(args[i])){
        cat("FILE:", i, args[i], "\n")
    }else{
        stop(paste(args[i], "Does not exist!\n"))
    }
}

library(data.table)
cat("Starting....\n")
#=====================================================================
# 1) The ann data
#=====================================================================
ann <- fread( ann.file )
annot <- data.frame(gene_id = ann$CIG.LOC,
                    old_id = ann$LOC,
                    transcript = ann$TRA,
                    cds = gsub("\\|", ":", ann$TDmodel),
                    scaffold = ann$scaff,
                    stringsAsFactors = FALSE)

rm(ann)

# REPORT:
cat("-----------------------------------------------------------------\n")
cat("STEP 1: parsing .annot\n")
cat("dimensions: ", dim(annot), "\n")
cat("STEP 1: parsing .annot --> FINISHED\n")

#=====================================================================
# 2) Parse the merged.gtf from cufflinks
#=====================================================================
# Use data.table to read in the gtf file since it's so big..
    
gtf <-  fread(gtf.file)

split.col9  <-  strsplit(gtf$V9, '\\;')
gtf$LOC <- gsub('gene_id ', '', sapply(split.col9 , '[', 1))
gtf$LOC <- gsub('"', "", gtf$LOC)
gtf$TRA <- gsub(' transcript_id ', '', sapply(split.col9 , '[', 2))
gtf$TRA <- gsub('"', "", gtf$TRA)
gtf$EX <-  gsub('exon_number ', '', sapply(split.col9 , '[', 3))
gtf$EX <- gsub('"', "", gtf$EX)
gtf$TSS <- gsub('tss_id ', '', sapply(split.col9 , '[', 5))
gtf$TSS <- gsub('"', "", gtf$TSS)

setnames(gtf, c("seqname", "source", "feature", "start", 
                "end", "score","strand", "frame", "attributes", 
                "gene", "transcript", "exon_number", "tss_id"))

gtf[,source:=NULL]
gtf[,feature:=NULL]
gtf[,frame:=NULL]
gtf[,score:=NULL]
gtf[,seqname:=NULL]
gtf[,gene:=NULL]  # new
gtf[,attributes:=NULL]
gtf[,tss_id:=NULL]

gtf <- as.data.frame(gtf)
gtf <- gtf[gtf$transcript %in% annot$transcript,] # new 

# REPORT:
cat("-----------------------------------------------------------------\n")
cat("STEP 2: parsing .gtf\n")
cat("dimensions: ", dim(gtf), "\n")
cat("STEP 2: parsing .gtf --> FINISHED\n")
#=====================================================================
# 3) parse the .bed file
#=====================================================================
bed <-  read.table(bed.file, header = FALSE, stringsAsFactors=FALSE, 
                   fill = TRUE)
bed <- bed[!grepl('track', bed[,1]),]

split.col4 <- strsplit(bed$V4, '\\;')
bed$FRAG <- gsub('ID=', '', sapply(split.col4 , '[', 1))
bed$STATUS <- gsub( '_len', '', sapply(
    strsplit(
        sapply(split.col4, "[", 3),
        ":"), 
    "[", 2))

bed$frame <- NA
bed$frame[bed$V7 %in% seq(1, 1000000, by=3)]  <- '1'
bed$frame[bed$V7 %in% seq(2, 1000000, by=3)]  <- '2'
bed$frame[bed$V7 %in% seq(3, 1000000, by=3)]  <- '3'
                                        #table(bed$frame)

cds <- data.frame('cds' =  gsub("\\|", ":",bed$FRAG), 
                  'cds_status' = bed$STATUS, 
                  'cds_start' = bed$V7,
                  'cds_end' = bed$V8, 
                  'cds_strand' = bed$V6, 
                  'cds_frame' = bed$frame,
                  stringsAsFactors = FALSE )    

cds <- cds[cds$cds %in% annot$cds,] # new

rm(bed)
# REPORT:
cat("-----------------------------------------------------------------\n")
cat("STEP 3: parsing .bed\n")
cat("dimensions: ", dim(cds), "\n")
cat("STEP 3: parsing .bed --> FINISHED\n")
#=====================================================================
# 4) parse the pfam file
#=====================================================================
pfam <-  read.table(pfam.file, stringsAsFactors=FALSE, fill =TRUE, 
                    header=FALSE,  comment.char = "/", colClasses=c(rep("character", 4), 
                                                           rep("NULL",24)))

pfam$V3 <- NULL
pfam$V1 <- NULL
                                        # split.id <- sapply(strsplit(pfam$V4, "\\|"), "[", 1)
                                        # 'transcript' = split.id,
pfam <- data.frame(
    'cds' = gsub("\\|", ":",pfam$V4),
    'pfam_id' = pfam$V2,
    stringsAsFactors = FALSE
    )

pfam <- pfam[pfam$cds %in% annot$cds,] # new

# REPORT:
cat("-----------------------------------------------------------------\n")
cat("STEP 4: parsing .pfam\n")
cat("dimensions: ", dim(pfam), "\n")
cat("STEP 4: parsing .pfam --> FINISHED\n")

#=====================================================================
# 6 parse the BLAST2GO topBlast hits and TRANSPOSONS
#===================================================================== 
blast.go <- read.csv(blast.go.file, header=TRUE, row.names = 1,
                     stringsAsFactors=FALSE, fill=TRUE)

transposons <- read.table(transpos, header=FALSE, stringsAsFactors = FALSE)
transposons <- gsub(".t[0123456789]+","", as.vector(transposons[,1]))

hit.df <- data.frame(gene_id = blast.go[,1], 
                     gene_name = blast.go[,2], 
                     gene_symbol = blast.go[,3], 
                     hit_acc = blast.go[,4], 
                     hit_eval = blast.go[,5],
                     hit_org = blast.go[,7]
                     )

hit.df$transposon <- as.numeric( hit.df$gene_id %in% transposons )


rm(blast.go, transposons)
# REPORT:
cat("-----------------------------------------------------------------\n")
cat("STEP 6: parsing .TopBlast\n")
cat("dimensions: ", dim(hit.df), "\n")
cat("STEP 6: parsing .TopBlast --> FINISHED\n")
#=====================================================================
# 7 parse the BLAST2GO GO hits
#=====================================================================
go <- fread(go.file, header=TRUE)
go[,'Term':=NULL]
go[,'Hit-Desc':=NULL]
setnames(go, c("gene_id","go_cat", "go_id"))
go <- as.data.frame(go, stringsAsFactors=FALSE)
go$gene_id <- gsub(".t[0123456789]", "", go$gene_id )
go$go_cat <- gsub("P", "BP",go$go_cat )
go$go_cat <- gsub("F", "MF",go$go_cat )
go$go_cat <- gsub("C", "CC",go$go_cat )
  
# REPORT:
cat("-----------------------------------------------------------------\n")
cat("STEP 7: parsing .TopBlast\n")
cat("dimensions: ", dim(go), "\n")
cat("STEP 7: parsing .TopBlast --> FINISHED\n")

#=====================================================================
# 8 parse the SNP data
#=====================================================================
snp <- read.table( snp.file, sep="\t", header = FALSE, stringsAsFactors=FALSE )
colnames(snp) = c("snp_id", "scaffold", "pos")

# REPORT:
cat("-----------------------------------------------------------------\n")
cat("STEP 8: parsing .snp\n")
cat("dimensions: ", dim(snp), "\n")
cat("STEP 8: parsing .snp --> FINISHED\n")

#=====================================================================
# 9 parse the KEGG data
#=====================================================================
kegg.knr <- read.table( file = knr.file, sep="\t", header=TRUE,
                       stringsAsFactors = FALSE)
kegg.kpath <- read.table( file = kpath.file, sep="\t", header=TRUE,
                         stringsAsFactors = FALSE)
kegg.kname <- read.table( file = kname.file, sep="\t", header=TRUE,
                         stringsAsFactors = FALSE)


# REPORT:
cat("-----------------------------------------------------------------\n")
cat("STEP 10: parsing .kegg\n")
cat("dimensions kegg K numbers: ", dim(kegg.knr), "\n")
cat("dimensions kegg K pathways: ", dim(kegg.kpath), "\n")
cat("dimensions kegg K pathway names: ", dim(kegg.kname), "\n")
cat("STEP 10: parsing .keggg --> FINISHED\n")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# MAKE the SQLite DB
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(RSQLite)

                                        #setup the DB
db <- dbConnect(SQLite(), dbname= db.name)

                                        # insert Tables into DB
db.finished <- vector()
db.finished[1] <- dbWriteTable(conn = db, name = "GTF", value = gtf, row.names=FALSE) 
db.finished[2] <- dbWriteTable(conn = db, name = "CDS", value = cds, row.names=FALSE) 
db.finished[3] <- dbWriteTable(conn = db, name = "PFAM", value = pfam, row.names=FALSE) 
db.finished[4] <- dbWriteTable(conn = db, name = "GENES", value = annot, row.names=FALSE) 
db.finished[5] <- dbWriteTable(conn = db, name = "HITS", value = hit.df, row.names=FALSE) 
db.finished[6] <- dbWriteTable(conn = db, name = "GO", value = go, row.names=FALSE) 
db.finished[7] <- dbWriteTable(conn = db, name = "SNP", value = snp, row.names=FALSE) 
# Kegg stuff
db.finished[8] <- dbWriteTable(conn = db, name = "KEGGK", value = kegg.knr, row.names=FALSE)
db.finished[9] <- dbWriteTable(conn = db, name = "KEGGP", value =  kegg.kpath, row.names=FALSE)
db.finished[10] <- dbWriteTable(conn = db, name = "KEGGN", value = kegg.kname, row.names=FALSE) 

# FINAL CHECK
cat("-----------------------------------------------------------------\n")
if(all(db.finished) == TRUE){
    cat("STEP 9: making .DB --> FINISHED\n")
}else{
    cat("ERROR: DB is not complete! Table(s) are missing\n")
}
