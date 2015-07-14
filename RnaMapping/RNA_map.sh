#!/bin/bash
module load anaconda
module list 
date

# sh script_test.sh -d -s -m

#-------------------------------------------------------------------------------
# Read script arguments and check if files/dirs exist

while [[ $# > 0 ]]
do
key="$1"

case $key in
    -d|--dirin)  # .fastq path
    DIRIN="$2"
    if [ ! -d "$DIRIN" ]; then
	echo 'ERROR: Directory' $DIRIN 'Does not exist!'
	exit 1
    fi
    shift # past argument
    ;;
    -s|--samplesheet) # OPTIONAL: illumina sample sheet
    SHEET="$2"
    shift # past argument
    ;;
    -g|--genome)      # OPTIONAL: path to genome
    GENOME="$2"
    shift # past argument
    ;;
    -m|--mastersheet) # file name of the mastersheet
    MASTER="$2"
    shift # past argument
    ;;
    --gtf)            # OPTIONAL: path to .gtf
    GTF="$2" 
    shift # past argument
    ;;
    -r|--read)        # long / short 
    READ="$2"
    shift # past argument
    ;;
    --execute)        # Only used for testing!
    EXECUTE="$2"
    shift # past argument
    ;;
esac
shift # past argument or value
done

#-------------------------------------------------------------------------------
echo '--------------------------------------------------------------------------'
# Default set to Genome: 
if [ -z "$GENOME" ] 
then 
    GENOME=/mnt/users/fabig/Ssa_genome/CIG_3.6v2_chrom-NCBI/STAR_index
fi
# check if Genome file exists
if [ ! -d "$GENOME" ]; then
   echo 'ERROR: File' $GENOME 'Does not exist!'
   exit 1
else 
    echo 'FOUND: File' $GENOME
fi

# Default set to GTF: 
if [ -z "$GTF" ] 
then 
    GTF=/mnt/users/fabig/Ssa_genome/CIG_3.6v2_chrom-NCBI/GTF/Salmon_3p6_Chr_070715_All.filter.gtf
fi
# check if GTF file exists
if [ ! -f "$GTF" ]; then
    echo 'ERROR: File' $GTF 'Does not exist!'
    exit 1
else 
    echo 'FOUND: File' $GTF 
fi

if [ -f "$SHEET" ]; then
    echo 'FOUND: File' $SHEET
fi


#-------------------------------------------------------------------------------
# If an Illumina SampleSheet is provided parse it
if [ -n "$SHEET" ]
then 
    python /mnt/users/fabig/cluster_pipelines/RnaMapping/helper_scripts/SampleSheetParser.py -s $SHEET -o $MASTER
else
    if cat -v $MASTER | grep -q '\^M' 
    then
	 echo 'Converting line terminators'
	 sed 's/\r/\n/g' $MASTER > sheet.tmp
	 mv sheet.tmp $MASTER
    else
	echo 'Line terminators seem correct'
    fi
fi

#-------------------------------------------------------------------------------
# Arguments for STAR

# Number of samples
END=$(cat $MASTER | wc -l)
echo 'SAMPLES:' $END

# Calculation: how many cores to use for STAR
CORES=$(($END * 3))
if [ $CORES -gt 20 ]
then 
    CORES=20
fi

case $READ in 
    "short")
	STAR=STAR
     ;;
    "long")
	STAR=STARlong
     ;;
    "")
	echo "-r|read: You have to specify short/long"
        exit 1
     ;;
esac



# echo all input variables
echo '--------------------------------------------------------------------------'
echo 'Input arguments:'
echo $DIRIN
echo $GENOME
echo $GTF
echo $STAR
echo '--------------------------------------------------------------------------'

# Create the folder tree if it does not exist
mkdir -p {slurm,bash,fastq_trim,fastq_trim_pe,qc,star,count}



#===============================================================================



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# PART 1: Quality trim the reads (ARRAY JOB)

cat > bash/sbatch-trim.sh << EOF
#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --array=1-$END
#SBATCH --job-name=TRIMMER
#SBATCH --output=slurm/trim-%A_%a.out
  
module load cutadapt
module list
date
  
FILEBASE=\$(awk ' NR=='\$SLURM_ARRAY_TASK_ID' { print \$2 ; }' $MASTER)

R1=$DIRIN'/'\$FILEBASE'_R1_001.fastq.gz'
O1='fastq_trim/'\$FILEBASE'_R1_001.trim.fastq.gz'
R2=$DIRIN'/'\$FILEBASE'_R2_001.fastq.gz'
O2='fastq_trim/'\$FILEBASE'_R2_001.trim.fastq.gz'

#------READ1---------------------
adaptorR1=\$(awk ' NR=='\$SLURM_ARRAY_TASK_ID' { print \$4 ; }' $MASTER)
adaptorR1='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'\$adaptorR1'ATCTCGTATGCCGTCTTCTGCTTG'

cutadapt -a \$adaptorR1 -q 20 -O 8 -o \$O1 \$R1

#------READ2---------------------
adaptorR2='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'

cutadapt -a \$adaptorR2 -q 20 -O 8 -o \$O2 \$R2

EOF


#-------------------------------------------------------------------------------
# PART 2: Quality control (ARRAY JOB)

cat > bash/sbatch-qc.sh << EOF
#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --array=1-$END
#SBATCH --job-name=QC
#SBATCH --output=slurm/qc-%A_%a.out
  
module load fastqc
module list
date
  
FILEBASE=\$(awk ' NR=='\$SLURM_ARRAY_TASK_ID' { print \$2 ; }' $MASTER)

R1='fastq_trim/'\$FILEBASE'_R1_001.trim.fastq.gz'
R2='fastq_trim/'\$FILEBASE'_R2_001.trim.fastq.gz'

fastqc -o qc \$R1
fastqc -o qc \$R2

EOF

#-------------------------------------------------------------------------------
# PART 3: Remove reads where one of the pairs is too short

cat > bash/sbatch-pairs.sh << EOF
#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --job-name=paired
#SBATCH --array=1-$END
#SBATCH --output=slurm/trimPE-%A_%a.out
  
module load anaconda
module list
date
  
script=/mnt/users/fabig/cluster_pipelines/RnaMapping/helper_scripts/fastqPE.py
  
FILEBASE=\$(awk ' NR=='\$SLURM_ARRAY_TASK_ID' { print \$2 ; }' $MASTER)

R1='fastq_trim/'\$FILEBASE'_R1_001.trim.fastq.gz'
R2='fastq_trim/'\$FILEBASE'_R2_001.trim.fastq.gz'
  
python \$script --f1 \$R1 --f2 \$R2 --c 40 --p fastq_trim_pe

EOF

#-------------------------------------------------------------------------------
# PART 4: STAR (LOOP)

cat > bash/sbatch-star.sh << EOF
#!/bin/bash
#SBATCH --job-name=STAR
#SBATCH -n $CORES
#SBATCH -N 1
#SBATCH --output=slurm/star-%j.out
    

module load star
module list
date

for TASK in {1..$END}
do

FILEBASE=\$(awk ' NR=='\$TASK' { print \$2 ; }' $MASTER)

R1='fastq_trim_pe/'\$FILEBASE'_R1_001.trim.fastq.gz'
R2='fastq_trim_pe/'\$FILEBASE'_R2_001.trim.fastq.gz'

OUT=star/\$FILEBASE

$STAR --limitGenomeGenerateRAM 62000000000 \
--genomeDir $GENOME \
--readFilesCommand zcat \
--readFilesIn \$R1 \$R2 \
--outFileNamePrefix \$OUT \
--outSAMmode Full \
--outSAMtype BAM SortedByCoordinate \
--runThreadN $CORES \
--readMatesLengthsIn NotEqual 

echo "FILE --> " \$OUT " PROCESSED"
 
done

EOF

#-------------------------------------------------------------------------------
# PART 5: HTSeq (ARRAY JOB)

cat > bash/sbatch-htseq.sh << EOF
#!/bin/sh
#SBATCH -n 1     
#SBATCH --array=1-$END       
#SBATCH --job-name=HTseq  
#SBATCH --output=slurm/HTSeq-%A_%a.out

module load samtools anaconda
module list
date

FILEBASE=\$(awk ' NR=='\$SLURM_ARRAY_TASK_ID' { print \$2 ; }' $MASTER)
INC=star/\$FILEBASE'Aligned.sortedByCoord.out.bam'
INN=star/\$FILEBASE'Aligned.sortedByName.out.bam'
OUT=count/\$FILEBASE'.count'

samtools sort -n -o \$INN -T \$INN'.temp' -O bam \$INC
samtools view \$INN | htseq-count -q -s reverse - $GTF > \$OUT

echo \$OUT "FINISHED"
EOF

#-------------------------------------------------------------------------------
# PART 6: STAR quality control

cat > bash/sbatch-star-check.sh << EOF
#!/bin/bash
#SBATCH --job-name=STARstat
#SBATCH -n 1
#SBATCH --output=slurm/slurm-star-stats%j.out

module load R
module list
date

cd star
R CMD BATCH /mnt/users/fabig/cluster_pipelines/RnaMapping/helper_scripts/STAR_Log.R
cd ..

EOF

#-------------------------------------------------------------------------------
# PART 7: HTSeq quality control

cat > bash/sbatch-htseq-check.sh << EOF
#!/bin/bash
#SBATCH --job-name=HTstat
#SBATCH -n 1
#SBATCH --output=slurm/slurm-htseq-stats-%j.out

module load R
module list
date

cd count
R CMD BATCH /mnt/users/fabig/cluster_pipelines/RnaMapping/helper_scripts/HTseq_plot.R
cd ..

EOF

#===============================================================================
#===============================================================================
# SUBMIT JOBS TO SLURM


if [ "$EXECUTE" != "no" ] 
then

#-------------------------------------------------------------------------------
# run sbatch file
command="sbatch bash/sbatch-trim.sh"
TrimJob=$($command | awk ' { print $4 }')
for i in $(seq 1 $END)
do 
    job=$TrimJob'_'$i
    if [ -z $TrimJobArray ] 
    then 
	TrimJobArray=$job
    else 
	TrimJobArray=$TrimJobArray':'$job
    fi
done
echo '---------------'
echo '1) sequences submitted for trimming'
echo '   slurm ID:' $TrimJobArray

#-------------------------------------------------------------------------------

command="sbatch --dependency=afterok:$TrimJobArray bash/sbatch-qc.sh" # Double quotes are essential!
QcJob=$($command | awk ' { print $4 }')
echo '---------------'
echo '2) Trimmed sequneces submitted for quality control'
echo '   slurm ID' $QcJob

#-------------------------------------------------------------------------------

command="sbatch --dependency=afterok:$TrimJobArray bash/sbatch-pairs.sh" # Double quotes are essential!
PairJob=$($command | awk ' { print $4 }')
for i in $(seq 1 $END)
do 
    job=$PairJob'_'$i
    if [ -z $PairJobArray ] 
    then 
	PairJobArray=$job
    else 
	PairJobArray=$PairJobArray':'$job
    fi
done
echo '---------------'
echo '3) Trimmed sequneces submitted removing too short reads'
echo '   slurm ID' $PairJobArray

#-------------------------------------------------------------------------------

command="sbatch --dependency=afterok:$PairJobArray bash/sbatch-star.sh"
StarJob=$($command | awk ' { print $4 }')
echo '---------------'
echo '4) Trimmed sequneces submitted for mapping'
echo '   slurm ID' $StarJob

#-------------------------------------------------------------------------------

command="sbatch --dependency=afterok:$StarJob bash/sbatch-htseq.sh"
HtseqJob=$($command | awk ' { print $4 }')
for i in $(seq 1 $END)
do 
    job=$HtseqJob'_'$i
    if [ -z $HtseqJobArray ] 
    then 
	HtseqJobArray=$job
    else 
	HtseqJobArray=$HtseqJobArray':'$job
    fi
done
echo '---------------'
echo '5) Mapped sequences submitted for counting'
echo '   slurm ID' $HtseqJobArray

#-------------------------------------------------------------------------------

command="sbatch --dependency=afterok:$HtseqJobArray bash/sbatch-star-check.sh"
StarStatJob=$($command | awk ' { print $4 }')
echo '---------------'
echo '5) Checking Star Log stats'
echo '   slurm ID' $StarStatJob

#-------------------------------------------------------------------------------

command="sbatch --dependency=afterok:$HtseqJobArray bash/sbatch-htseq-check.sh"
HtseqStatJob=$($command | awk ' { print $4 }')
echo '---------------'
echo '6) Checking HTSeq counting stats'
echo '   slurm ID' $HtseqStatJob

else
    echo '=================='
    echo 'Nothing submitted!'
    echo '=================='
fi
