
#!/bin/bash

echo ''
echo $(date)
echo 'SNP_call.sh Version 1.0.2'

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
    --RunType)
        RTYP="$2"  # Accepts makeIdx/mapp
	shift
	;;
esac
shift # past argument or value
done


if [ "$RTYP" != "makeIdx" ] && [ "$RTYP" != "mapp" ]
then
    echo "ERROR --RunType has to be set to makeIdx/mapp"
    exit 1
fi
	
#-------------------------------------------------------------------------------
# CHECK IF THE MAIN VARIABLES HAVE BEEN SET
    
# check if .fasta file exists
if [ ! -f "$GENFA" ]; then
    echo 'ERROR: File' $GENFA 'Does not exist!'
    exit 1
else 
    echo 'FOUND: File' $GENFA
fi

# check if .gtf file exists
if [ ! -f "$GTF" ]; then
    echo 'ERROR: File' $GTF 'Does not exist!'
    exit 1
else 
    echo 'FOUND: File' $GTF 
fi


   
# Dependent on the --makeIdx option (IDX) check if the main input exists
if [ "$RTYP" == "makeIdx" ]
then
    # check if Read length is set
    if [ -z "$RLEN" ]; then
	echo 'ERROR: You have to specify read length at --read'
	exit 1
    fi
    
    mkdir -p {slurm,bash,$GENDIR}
    # echo all input variables
    echo ''
    echo 'Executing' $EXECUTE
    echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    echo 'Input arguments:'
    echo '-----------------------'
    echo 'Genome .fasta:'
    echo $GENFA
    echo '.gtf file'
    echo $GTF
    echo 'Read length'
    echo $RLEN
    echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

else  

    # check if MASTER file exists
    if [ ! -f "$MASTER" ]
    then
	echo 'ERROR: File' $MASTER 'Does not exist!'
	exit 1
    else 
	echo 'FOUND: File' $MASTER 
    fi

    # Create the folder tree if it does not exist
    mkdir -p {slurm,bash,star_2pass,mapp_summary,gatk}

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

    # Determine which STAR command to use (e.g STAR/STARlong)
    if [ $RLEN -ge 300 ]
    then
	STAR="STARlong"
    else
	STAR="STAR"
    fi

    #-------------------------------------------------------------------------------
    # Number of samples
    END=$(sed '1d' $MASTER | wc -l) # skip hearder line

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
    echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    
fi

#-------------------------------------------------------------------------------






if [ "$RTYP" == "makeIdx"  ]  ## Only use the commands to make the STAR index
then
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
SJa=($SJ) 
SJa=\$(echo \$SJa | tr "," "\n")

# Join all SJ files
for i in \$SJa
do
    echo '++>>' $i
    for ii in \$(ls \$i/*SJ.out.tab)
    do
	echo \$ii
	awk '\$7>1' \$ii >> $GENDIR/SJ_all.tab
    done
done

# Filter the joined file
awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} {if(\$5>0){print \$1,\$2,\$3,strChar[\$4]}}' $GENDIR/SJ_all.tab | sort | uniq > $GENDIR/SJ_in.tab

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


STAR --runMode genomeGenerate \
--runThreadN 10 \
--limitGenomeGenerateRAM 62000000000 \
--genomeChrBinNbits 12 \
--genomeDir $GENDIR \
--sjdbGTFtagExonParentTranscript Parent \
--genomeFastaFiles $GENFA \
--sjdbOverhang \$(($RLEN-1)) \
--sjdbFileChrStartEnd $GENDIR/SJ_in.tab \
--sjdbGTFfile $GTF

echo '==>>FINISHED'

EOF

else 

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



cat > bash/snp_call-StarCheck.sh << EOF
#!/bin/bash
#SBATCH --job-name=STARstat
#SBATCH -n 1
#SBATCH --output=slurm/snp_call-StarStats%j.out

module load R
module list
date

cd star_2pass
R CMD BATCH /mnt/users/fabig/cluster_pipelines/RnaMapping/helper_scripts/STAR_Log.R
cd ..

EOF


# ==============================================================================
# 4) Add read groups, sort, mark duplicates, and create index

cat > bash/snp_call-mDupl.sh << EOF
#!/bin/bash
#SBATCH --job-name=picard
#SBATCH -n 5
#SBATCH --mem=10G
#SBATCH --array=1-$END%20
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

samtools sort -@5 -m 2G -o \$BAMC -T \$BAM'.temp' -O bam \$BAM
#rm \$BAM
echo '==> Done sorting'


picard MarkDuplicates I=\$BAMC O=\$PIC CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=\$PIM 
#rm \$BAMC
echo '==> Done marking duplicates'
EOF

# ==============================================================================
# 5) Split'N'Trim and reassign mapping qualities

cat > bash/snp_call-splitNtrim.sh << EOF
#!/bin/bash
#SBATCH --job-name=splitNtrim
#SBATCH -n 1
#SBATCH --array=1-$END
#SBATCH --output=slurm/snp_call-splitNtrim-%A_%a.out

module load gatk/3.5
module list 
date

FILEBASE=\$(awk ' NR=='\$SLURM_ARRAY_TASK_ID+1' { print \$2 ; }' $MASTER)
BAM_in=star_2pass/\$FILEBASE'-2pass-Aligned.out.sortedByCoord.dedupped.bam'
BAM_cig=gatk/\$FILEBASE'-2pass-Aligned.out.sortedByCoord.dedupped.splitCig.bam'

gatk SplitNCigarReads -R $GENFA \
-I \$BAM_in \
-o \$BAM_cig \
-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 \
-U ALLOW_N_CIGAR_READS

echo '==> Done splitting cigars'
#rm \$BAM_in
EOF



cat > bash/snp_call-Variant.sh << EOF
#!/bin/bash
#SBATCH --job-name=variant
#SBATCH -n 10
#SBATCH --nodes=1
#SBATCH --array=1-$END%10
#SBATCH --output=slurm/snp_call-variant-%A_%a.out

module load gatk/3.5
module list
date

# set -o nounset   # Prevent unset variables from been used.
set -o errexit   # Exit if error occurs

## RNA-seq variant calling as described in GATK best practises https://www.broadinstitute.org/gatk/guide/article?id=3891
## Best practies scripted below where Last updated on 2015-12-07 11:08:30 in link above
## Usage test case: sbatch SNP_call.sh ref BAM_in

## Set variables
FILEBASE=\$(awk ' NR=='\$SLURM_ARRAY_TASK_ID+1' { print \$2 ; }' $MASTER)
BAM_cig=gatk/\$FILEBASE'-2pass-Aligned.out.sortedByCoord.dedupped.splitCig.bam'
VCF_out=gatk/\$(basename \$BAM_cig .bam).vcf
VCF_temp=\$VCF_out'.tmp'

threads=10


## Do the variant calling on 2-pass recalibrated BAM.
## -recoverDanglingHeads is default in gatk 3.5, but was not in previos versions. Use 3.5! 

echo "==> HaplotypeCaller"

gatk -T HaplotypeCaller \
-R $GENFA -I \$BAM_cig \
-dontUseSoftClippedBases \
-stand_call_conf 20.0 \
-stand_emit_conf 20.0 \
--num_cpu_threads_per_data_thread \$threads \
-o \$VCF_temp
echo '==> Done HaplotypeCaller'

## Do RNA-seq specific filtering
echo "==> VariantFiltration"
gatk -T VariantFiltration -R $GENFA -V \$VCF_temp \
-window 35 -cluster 3 \
-filterName FS -filter "FS>30.0" \
-filterName QD -filter "QD<2.0" \
-o \$VCF_out
echo '==> VariantFiltration'

echo "==> Finished"
EOF

fi
# ==============================================================================

# SCRIPT SUBMISSION

if [ "$EXECUTE" != "no" ] 
then
    if [ "$RTYP" == "makeIdx" ] 
    then
	#-------------------------------------------------------------------------------
	# 1) splice junctions
	command='sbatch bash/snp_call-SJ.sh'
	SJjob=$($command | awk ' { print $4 }')
	echo '---------------'
	echo ' Splice junctions'
	echo ' slurm ID:' $SJjob
	
	#-------------------------------------------------------------------------------
	# 2) STAR index
	command="sbatch --dependency=afterok:$SJjob bash/snp_call-StarIdx.sh"
	IDXjob=$($command | awk ' { print $4 }')
	echo '---------------'
	echo ' Making STARindex'
	echo ' slurm ID' $IDXjob
    else
	#-------------------------------------------------------------------------------
	# 1) 2nd Round STAR
	command="sbatch bash/snp_call-Star.sh"
	STARjob=$($command | awk ' { print $4 }')
	echo '---------------'
	echo ' 2nd round mapping'
	echo ' slurm ID' $STARjob

	#-------------------------------------------------------------------------------
	command="sbatch --dependency=afterok:$STARjob bash/snp_call-StarCheck.sh"
	CHECKjob=$($command | awk ' { print $4 }')
	echo '---------------'
	echo ' star check processing'
	echo ' slurm ID' $CHECKjob

	#-------------------------------------------------------------------------------
	# 1) Picard
	command="sbatch --dependency=afterok:$STARjob bash/snp_call-mDupl.sh"
	PICjob=$($command | awk ' { print $4 }')
	echo '---------------'
	echo ' picard processing'
	echo ' slurm ID' $PICjob
    
    
	#-------------------------------------------------------------------------------
	command="sbatch --dependency=afterok:$PICjob bash/snp_call-splitNtrim.sh"
	SPLITjob=$($command | awk ' { print $4 }')
	echo '---------------'
	echo ' splitNtrim processing'
	echo ' slurm ID' $SPLITjob

        #-------------------------------------------------------------------------------
	command="sbatch --dependency=afterok:$SPLITjob bash/snp_call-Variant.sh"
	VARjob=$($command | awk ' { print $4 }')
	echo '---------------'
	echo ' Variant calling+filtering'
	echo ' slurm ID' $VARjob
    fi   
else
    echo '=================='
    echo 'Nothing submitted!'
    echo '=================='
fi
