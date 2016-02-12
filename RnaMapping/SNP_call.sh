
#!/bin/bash

echo ''
echo $(date)
echo 'SNP_call.sh Version 1.0.0'

#-------------------------------------------------------------------------------
# Read script arguments and check if files/dirs exist

while [[ $# > 0 ]]
do
key="$1"

case $key in
    --sj)
	SJ="$2" # File path to all folders containing splice junction files (STAR folder)
	shift   # if more than 1 folder seperate by comma.
	;;
    --GenFa)      # OPTIONAL: path to genome
	GENFA="$2"
	shift # past argument
	;;
    --GenDir)      # OPTIONAL: path to the newgenome
	GENDIR="$2"
	shift # past argument
	;;
    --FaDir)      # .fastq path
	FADIR="$2"
	shift # past argument
	;;
    --read)      # Read length
	RLEN="$2"
	shift # past argument
	;;
    --gtf)            # OPTIONAL: path to .gtf
	GTF="$2" 
	shift # past argument
	;;
    -m|--mastersheet) # file name of the mastersheet
	MASTER="$2"
	shift # past argument
	;;
    --execute)        # Only used for testing; use --execute no
	EXECUTE="$2"
	shift # past argument
	;;
esac
shift # past argument or value
done

# Create the folder tree if it does not exist
mkdir -p {slurm,bash,star_2pass,star_genome,mapp_summary}

#-------------------------------------------------------------------------------
# Check line terminators in MASTER
if cat -v $MASTER | grep -q '\^M' 
then
    echo 'Converting line terminators'
    sed 's/\r/\n/g' $MASTER > sheet.tmp
    mv sheet.tmp $MASTER
else
    echo 'Line terminators seem correct'
fi


#-------------------------------------------------------------------------------
# Number of samples
END=$(sed '1d' $MASTER | wc -l) # skip hearder line

if [ "$RLEN" -ge 200 ]; then
    STAR=STARlong
elif [ "$RLEN" -le 200 ]; then
    STAR=STAR
fi

#------------------------------------------------------------------------------

# echo all input variables
echo ''
echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
echo 'Input arguments:'
echo '-----------------------'
echo 'Location of .fastq files:'
echo $FADIR
echo 'Genome .fasta:'
echo $GENFA
echo '.gtf file'
echo $GTF
echo 'STAR command:'
echo $STAR
echo '-----------------------'

echo 'Number of samples= ' $END
echo 'FIRST sample:' $(awk ' NR=='2' {OFS="\t"; print; }' $MASTER)
echo 'LAST sample:' $(awk ' NR=='$END+1' {OFS="\t"; print; }' $MASTER)
if [ "$CPFOLDER" != "no" ]; then
    echo '-----------------------'
    echo 'Copy to common foder:'
    echo $COMMON/$CPFOLDER
fi
echo ''
echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'


# ==============================================================================
# 1) Junctions
cat > bash/snp_call-SJ.sh << EOF
#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --job-name=JUNCTIONS
#SBATCH --output=slurm/snp_call-SJ-%A.out
  
module list
date
 
# convert 2 array
OIFS=\$IFS; IFS=",";
SJa=($SJ) 

# Join all SJ files
for ((i=0; i<\${#SJa[@]}; ++i))
    do
	SJdir=\${SJa[\$i]}
	awk '\$7>1' \$SJdir/*SJ.out.tab >> star_2pass/SJ_all.tab
done

# Filter the joined file
awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} {if(\$5>0){print \$1,\$2,\$3,strChar[\$4]}}' star_2pass/SJ_all.tab | sort | uniq > star_2pass/SJ_in.tab

echo '==>>FINISHED'
EOF

# ==============================================================================
# 2) Remake the STAR genome index considering the new splice junctions

cat > bash/snp_call-StarIdx.sh << EOF
#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --job-name=StarIdx
#SBATCH -n 10
#SBATCH -N 1
#SBATCH --mem=60G
#SBATCH --output=slurm/snp_call-StarIdx-%A.out

module load star
module list
date

mkdir -p $GENDIR

STAR --runMode genomeGenerate \
--runThreadN 10 \
--limitGenomeGenerateRAM 62000000000 \
--genomeChrBinNbits 12 \
--genomeDir $GENDIR \
--genomeFastaFiles $GENFA \
--sjdbOverhang \$(($RLEN-1)) \
--sjdbFileChrStartEnd star_2pass/SJ_in.tab \
--sjdbGTFfile $GTF

echo '==>>FINISHED'

EOF

# eventually kill the old STAR directory after this

# ==============================================================================
# 3) ReRun STAR

cat > bash/snp_call-Star.sh << EOF
#!/bin/bash
#SBATCH --job-name=STAR2nd
#SBATCH -n 20
#SBATCH -N 1
#SBATCH --mem=60G
#SBATCH --array=1-$END%5
#SBATCH --output=slurm/snp_call-Star-%A_%a.out
    
module load star
module list
date

SAMPLE=\$(awk ' NR=='\$SLURM_ARRAY_TASK_ID+1' { print \$1 ; }' $MASTER)
FILEBASE=\$(awk ' NR=='\$SLURM_ARRAY_TASK_ID+1' { print \$2 ; }' $MASTER)

R1A=( \$( ls $FADIR\$FILEBASE*_R1_*.fastq.gz ))
R2A=( \$( ls $FADIR\$FILEBASE*_R2_*.fastq.gz ))

R1=\$(printf ",%s" "\${R1A[@]}")
R2=\$(printf ",%s" "\${R2A[@]}")
R1=\${R1:1}
R2=\${R2:1}

#-------------------------------------------------------------------------------
# BAM read groups:
RGa=() # Array to hold the read group string
for ((i=0; i<\${#R1A[@]}; i++))
do
    printf "%s\t%s\n" \$( printf "ID:L%03d" \$((\$i+1)) ) \$( basename \${R1A[\$i]} .trim.fastq.gz ) >> ReadGroup_summary.txt
    RGa+=\$(printf "ID:L%03d PL:illumina LB:\$SAMPLE SM:\$SAMPLE , " \$(($i+1)))
done

RG=\$( printf "%s" "\${RGa[@]}" ) 			
RG=\$(echo \$RG | sed 's/ ,\$//g' ) 
#-------------------------------------------------------------------------------

OUT=star_2pass/\$FILEBASE'-2pass-'

echo "Running  --> " \$R1 \$R2

# Run STAR
$STAR --limitGenomeGenerateRAM 62000000000 \
--genomeDir $GENDIR \
--readFilesCommand zcat \
--readFilesIn \$R1 \$R2 \
--outFileNamePrefix \$OUT \
--outSAMmode Full \
--outSAMtype BAM Unsorted \
--runThreadN 20 \
--readMatesLengthsIn NotEqual \
--outSAMattrRGline \$RG

echo "FILE --> " \$OUT " PROCESSED"

EOF


# ==============================================================================
# 4) Add read groups, sort, mark duplicates, and create index

cat > bash/snp_call-mDupl.sh << EOF
#!/bin/bash
#SBATCH --job-name=picard
#SBATCH -n 2
#SBATCH --array=1-$END
#SBATCH --output=slurm/snp_call-picard-%A_%a.out

module load samtools picard
module list
date

#-------------------------------------------------------------------------------
FILEBASE=\$(awk ' NR=='\$SLURM_ARRAY_TASK_ID+1' { print \$2 ; }' $MASTER)
BAM=star_2pass/\$FILEBASE'-2pass-Aligned.out.bam'
BAMC=star_2pass/\$FILEBASE'-2pass-Aligned.out.sortedByCoord.bam'
PIC=star_2pass/\$FILEBASE'-2pass-Aligned.out.sortedByCoord.dedupped.bam'
PIM=star_2pass/\$FILEBASE'.dedupped.metrics'
#-------------------------------------------------------------------------------

samtools sort -n -o \$BAMC -T \$BAM'.temp' -O bam \$BAM
rm \$BAM

picard MarkDuplicates I=\$BAMC O=\$PIC CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=\$PIM 
rm \$BAMC

EOF

# ==============================================================================


# Unfinished stuff from TIM
cat > bash/snp_call-TIM.sh << EOF
#SBATCH -N 1
#SBATCH -n 5
module load gatk/3.5
set -o nounset   # Prevent unset variables from been used.
set -o errexit   # Exit if error occurs

## RNA-seq variant calling as described in GATK best practises https://www.broadinstitute.org/gatk/guide/article?id=3891
## Best practies scripted below where Last updated on 2015-12-07 11:08:30 in link above

## Usage test case: sbatch SNP_call.sh ref BAM_in

## Set variables
threads=5 
BAM_input=$2
VCF_temp=${RANDOM}.vcf
reference=$1
VCF_out=test/$(basename $BAM_input .bam).vcf

## Test variables
#BAM_input=test.bam
#VCF_temp=${RANDOM}.vcf
#reference=test.fasta
#VCF_out=$(basename $BAM_input .bam).vcf


## Do the variant calling on 2-pass recalibrated BAM.
## -recoverDanglingHeads is default in gatk 3.5, but was not in previos versions. Use 3.5! 

echo "Running commnds:"
echo
echo gatk -T HaplotypeCaller \
	-R $reference \
	-I $BAM_input \
	-dontUseSoftClippedBases \
	-stand_call_conf 20.0 \
	-stand_emit_conf 20.0 \
	--num_cpu_threads_per_data_thread $threads \
	-o $VCF_temp

## Do RNA-seq specific filtering
echo gatk -T VariantFiltration \
	-R $reference \
	-V $VCF_temp \
	-window 35 -cluster 3 \
	-filterName FS -filter "FS > 30.0" \
	-filterName QD -filter "QD < 2.0" \
	-o $VCF_out

gatk -T HaplotypeCaller \
	-R $reference \
	-I $BAM_input \
	-dontUseSoftClippedBases \
	-stand_call_conf 20.0 \
	-stand_emit_conf 20.0 \
	--num_cpu_threads_per_data_thread $threads \
	-o $VCF_temp

## Do RNA-seq specific filtering
gatk -T VariantFiltration \
	-R $reference \
	-V $VCF_temp \
	-window 35 -cluster 3 \
	-filterName FS -filter "FS > 30.0" \
	-filterName QD -filter "QD < 2.0" \
	-o $VCF_out

EOF
