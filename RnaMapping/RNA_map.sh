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
    -d|--dirin)
    DIRIN="$2"
    if [ ! -d "$DIRIN" ]; then
	echo 'ERROR: Directory' $DIRIN 'Does not exist!'
	exit 1
    fi
    shift # past argument
    ;;
    -s|--samplesheet)
    SHEET="$2"
    if [ ! -f "$SHEET" ]; then
	echo 'ERROR: File' $SHEET 'Does not exist!'
	exit 1
    fi
    shift # past argument
    ;;
    -g|--genome)
    GENOME="$2"
    if [ ! -f "$SHEET" ]; then
	echo 'ERROR: File' $GENOME 'Does not exist!'
	exit 1
    fi
    shift # past argument
    ;;
    -m|--mastersheet)
    MASTER="$2"
    shift # past argument
    ;;
    --gtf)
    GTF="$2"
    shift # past argument
    ;;
    -r|--read)
    READ="$2"
    shift # past argument
    ;;
esac
shift # past argument or value
done
#-------------------------------------------------------------------------------

# Default set to Genome: 
if [ -z "$GENOME" ] 
then GENOME=/mnt/users/fabig/RNAseq/Ssa_genome/CIG_3.6v2_chrom-NCBI/STAR_index
fi
# Default set to GTF: 
if [ -z "$GTF" ] 
then GTF=/mnt/users/fabig/RNAseq/Ssa_genome/CIG_3.6v2_chrom-NCBI/GTF/Salmon_3p6_Chr_NCBI_230415.gtf
fi

# echo all input variables
echo '--------------------------------------------------------------------------'
echo 'Input arguments:'
echo $DIRIN
echo $SHEET
echo $GENOME
echo $GTF
echo '--------------------------------------------------------------------------'

# Create the folder tree if it does not exist
mkdir -p {slurm,bash,fastq_trim,fastq_trim_pe,qc,star,count}


# Parse the Illumina sample sheet
python /mnt/users/fabig/cluster_pipelines/RnaMapping/helper_scripts/SampleSheetParser.py -s $SHEET -o $MASTER


END=$(cat $MASTER | wc -l)

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
  
#------READ1---------------------
adaptorR1=\$(awk ' NR=='\$SLURM_ARRAY_TASK_ID' { print \$4 ; }' $MASTER)
adaptorR1='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'\$adaptorR1'ATCTCGTATGCCGTCTTCTGCTTG'

filebase=\$(awk ' NR=='\$SLURM_ARRAY_TASK_ID' { print \$5 ; }' $MASTER)

fileR1=$DIRIN'/'\$filebase'_L001_R1_001.fastq.gz'
fileO1='fastq_trim/'\$filebase'_L001_R1_001.trim.fastq.gz'

cutadapt -a \$adaptorR1 -q 20 -O 8 -o \$fileO1 \$fileR1

#------READ2---------------------
adaptorR2='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'

fileR2=$DIRIN'/'\$filebase'_L001_R2_001.fastq.gz'
fileO2='fastq_trim/'\$filebase'_L001_R2_001.trim.fastq.gz'

cutadapt -a \$adaptorR2 -q 20 -O 8 -o \$fileO2 \$fileR2

EOF

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
echo $TrimJobArray
echo '1) sequences submitted for trimming'

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
  
filebase=\$(awk ' NR=='\$SLURM_ARRAY_TASK_ID' { print \$5 ; }' $MASTER)

R1='fastq_trim/'\$filebase'_L001_R1_001.trim.fastq.gz'
R2='fastq_trim/'\$filebase'_L001_R2_001.trim.fastq.gz'

fastqc -o qc \$R1
fastqc -o qc \$R2

EOF

command="sbatch --dependency=afterok:$TrimJobArray bash/sbatch-qc.sh" # Double quotes are essential!
QcJob=$($command | awk ' { print $4 }')
echo '2) Trimmed sequneces submitted for quality control'

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
  
filebase=\$(awk ' NR=='\$SLURM_ARRAY_TASK_ID' { print \$5 ; }' $MASTER)

R1='fastq_trim/'\$filebase'_L001_R1_001.trim.fastq.gz'
R2='fastq_trim/'\$filebase'_L001_R2_001.trim.fastq.gz'
  
python $script --f1 \$R1 --f2 \$R2 --c 40 --p fastq_trim_pe

EOF
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
echo $PairJobArray

echo '3) Trimmed sequneces submitted removing too short reads'

#-------------------------------------------------------------------------------
# PART 4: STAR (LOOP)

# Calculation: how many cores to use for STAR
CORES=$(($END * 3))
if [ $CORES -gt 30 ]
then 
CORES=30
fi

# Which STAR command to use
if [ -z "$READ" ] 
then 
   STAR=STAR
else 
   STAR=STARlong
fi


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

filebase=\$(awk ' NR=='\$TASK' { print \$5 ; }' $MASTER)

R1='fastq_trim_pe/'\$filebase'_L001_R1_001.trim.fastq.gz'
R2='fastq_trim_pe/'\$filebase'_L001_R2_001.trim.fastq.gz'

OUT=star/\$filebase

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

command="sbatch --dependency=afterok:$PairJobArray bash/sbatch-star.sh"
StarJob=$($command | awk ' { print $4 }')
echo '4) Trimmed sequneces submitted for mapping'

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

filebase=\$(awk ' NR=='\$TASK' { print \$5 ; }' $MASTER)
IN=star/\$filebase'Aligned.sortedByCoord.out.bam'
OUT=count/\$filebase'.count'

samtools sort -n -O sam -T \$IN | htseq-count -q -s reverse - $GTF > \$OUT
echo "FINISHED"
EOF

command="sbatch --dependency=afterok:$StarJob bash/sbatch-htseq.sh"
HtseqJob=$($command | awk ' { print $4 }')
for i in $(seq 1 $END)
do 
    job=$TrimJob'_'$i
    if [ -z $HtseqJobArray ] 
    then 
	HtseqJobArray=$job
    else 
	HtseqJobArray=$HtseqJobArray':'$job
    fi
done
echo $HtseqJobArray
echo '5) Mapped sequences submitted for counting'


#-------------------------------------------------------------------------------
# PART 6: STAR quality control

cat > bash/sbatch-star-check.sh << EOF
#!/bin/bash
#SBATCH --job-name=STARstat
#SBATCH -n $CORES
#SBATCH -N 1
#SBATCH --output=SLURM/slurm-star-stats%j.out

module load R
module list
date

cd star
R CMD BATCH /mnt/users/fabig/cluster_pipelines/RnaMapping/helper_scripts/STAR_Log.R
cd ..

EOF

command="sbatch --dependency=afterok:$HtseqJobArray bash/sbatch-star-check.sh"
StarStatJob=$($command | awk ' { print $4 }')
echo '5) Checking Star Log stats'


#-------------------------------------------------------------------------------
# PART 7: HTSeq quality control

cat > bash/sbatch-htseq-check.sh << EOF
#!/bin/bash
#SBATCH --job-name=HTstat
#SBATCH -n $CORES
#SBATCH -N 1
#SBATCH --output=SLURM/slurm-htseq-stats-%j.out

module load R
module list
date

cd count
R CMD BATCH /mnt/users/fabig/cluster_pipelines/RnaMapping/helper_scripts/HTseq_plot.R
cd ..

EOF

command="sbatch --dependency=afterok:$HtseqJobArray bash/sbatch-htseq-check.sh"
StarStatJob=$($command | awk ' { print $4 }')
echo '6) Checking HTSeq counting stats'
