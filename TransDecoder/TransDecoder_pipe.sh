#!/bin/bash
echo ''
date


#-------------------------------------------------------------------------------
# Read script arguments and check if files/dirs exist
while [[ $# > 0 ]]
do
key="$1"
case $key in
-f|--fasta) # Input transcriptome fasta, full path (inputFASTA)
FASTA="$2"
if [ ! -f "$FASTA"  ]; then
    echo 'ERROR: File' $FASTA 'Does not exist!'
    exit 1
fi
shift # past argument
;;
-s|--split) # Number of chunks to split the transcriptome into.
SPLIT="$2"
if [ -z "$SPLIT" ]; then
    echo 'ERROR: You have to specify the number of chunks'
    exit 1
fi
shift # past argument
;;
-b|--blastdb) 
BLASTDB="$2"  # path to BLASTP DB
shift # past argument
;;
-p|--pfamdb) # PFAM DB
PFAMDB="$2"   
shift # past argument
;;
 --execute) # Only used for testing!
EXECUTE="$2"
shift # past argument
;;
esac
shift # past argument or value
done

#-------------------------------------------------------------------------------
echo '--------------------------------------------------------------------------'
# Create folder structure
mkdir -p {slurm,bash,chunks,blast,pfam,td-final}

# assign variable for the longest ORF chunk path
CHUNKDIR=chunks/longest_orfs.pep_chunk_

# assign relative .fasta name
RELFA=$(basename $FASTA)

#-------------------------------------------------------------------------------
echo 'Input arguments:'
echo $FASTA
echo $BLASTDB
echo $PFAMDB
#-------------------------------------------------------------------------------
echo '--------------------------------------------------------------------------'
# part 1: Run TransDecoder.

cat > bash/TD1-run.sh << EOF 
#!/bin/bash
#SBATCH --job-name=TD.run:1of6
#SBATCH -n 3
#SBATCH --output=slurm/TD_predictORF-%A.out
        
module load transdecoder/2.0.1

TD=/local/genome/packages/transdecoder/2.0.1/TransDecoder.LongOrfs

\$TD -t ${FASTA} -m 30

EOF

#-------------------------------------------------------------------------------
# part 2: split .fasta in X chunks
 
cat > bash/TD2-split.sh << EOF 
#!/bin/bash
#SBATCH --job-name=TD.split:2of6
#SBATCH -n 1
#SBATCH --output=slurm/TD_splitORF-%A.out

module load cigene-tools/1.0

fastasplit -f $RELFA.transdecoder_dir/longest_orfs.pep -o chunks -c $SPLIT

EOF

#-------------------------------------------------------------------------------
# part 3: Run BlastP on each chunk 
cat > bash/TD3-blastp.sh << EOF 
#!/bin/sh
#SBATCH -n 2 # NB CPUS
#SBATCH --job-name=TD.blastp:3of6
#SBATCH --array=1-$SPLIT # array number
#SBATCH --output=slurm/TD_blastp-%A-%a.out

module load blast+/2.2.28
  
SGE_TASK_ID=\$(expr \$SLURM_ARRAY_TASK_ID - 1)
printf -v ID '%07g' \$SGE_TASK_ID              # will produce feks 0000001 (matches extension produced by fastasplit)

blastp -query $CHUNKDIR\$ID \
  -db $BLASTDB \
  -outfmt 6 -num_threads 2 \
  -evalue 1e-10  -max_target_seqs 1 \
  -out blast/TD_blastp-\$ID
EOF

#-------------------------------------------------------------------------------
# part 4: Run PFAM
cat > bash/TD4-pfam.sh << EOF 
#!/bin/sh
#SBATCH -n 2 
#SBATCH --job-name=TD.pfam:4of6
#SBATCH --array=1-$SPLIT
#SBATCH --output=slurm/TD_pfam-%A-%a.out

module load hmmer/3.1b1  
  
SGE_TASK_ID=\$(expr \$SLURM_ARRAY_TASK_ID - 1)
printf -v ID '%07g' \$SGE_TASK_ID              # will produce feks 0000001 (matches extension produced by fastasplit)

hmmscan --cpu 2 --domtblout pfam/TD_pfam-\$ID \
 $PFAMDB \
 $CHUNKDIR\$ID
EOF

#-------------------------------------------------------------------------------
# part 5: Join chunks  
cat > bash/TD5-join.sh << EOF 
#!/bin/sh
#SBATCH -n 1 # NB CPUS
#SBATCH --job-name=TD.join:5of6
#SBATCH --output=slurm/TD_join-%A.out

module load R anaconda

cat blast/TD_blastp-* > blast/TD_blastp_conc.tab
cat pfam/TD_pfam-* > pfam/TD_pfam_conc.tab 

# Fix some column issues
pyscr=/mnt/users/fabig/cluster_pipelines/TransDecoder/helper_scripts/pfam-space-to-tab.py
python \$pyscr pfam/TD_pfam_conc.tab > pfam/TD_pfam_conc.fixcolsep.tab

EOF

#-------------------------------------------------------------------------------
# part 6: TransDecoder 
cat > bash/TD6-run2.sh << EOF 
#!/bin/bash
#SBATCH --job-name=TD.run2:6of6
#SBATCH -n 3
#SBATCH --output=slurm/TD_run2-%A.out

module load transdecoder/2.0.1

TD=/local/genome/packages/transdecoder/2.0.1/TransDecoder.Predict

\$TD -t $FASTA --retain_pfam_hits pfam/TD_pfam_conc.tab \
  --retain_blastp_hits blast/TD_blastp_conc.tab \
  --retain_long_orfs 120

rm -r $RELFA'.transdecoder_dir'

# rename the terrible long dir
mkdir -p td-models
mv \$(ls | grep 'transdecoder' | grep -v 'transdecoder_dir') td-models

EOF

#-------------------------------------------------------------------------------
# part 7: Filter results
cat > bash/TD7-filt.sh << EOF 
#!/bin/bash
#SBATCH --job-name=TD.filt:7of7
#SBATCH -n 1
#SBATCH --output=slurm/TD_filt-%A.out

module load R

FSCR=/mnt/users/fabig/cluster_pipelines/TransDecoder/helper_scripts/TransDecoder_filtering.R

Rscript \$FSCR $FASTA td-models/$RELFA'.transdecoder.bed' blast/TD_blastp_conc.tab pfam/TD_pfam_conc.fixcolsep.tab 

EOF

#===============================================================================
# Submit the jobs to SLURM

if [ "$EXECUTE" != "no" ]
then


command="sbatch bash/TD1-run.sh"
td1=$($command | awk ' { print $4 }')
echo '1) TransDecoder part 1 submitted' $td1

#---
command="sbatch --dependency=$td1 bash/TD2-split.sh"  
td2=$($command | awk ' { print $4 }')
echo '2) FASTA split submitted' $td2

#----
command="sbatch --dependency=afterok:$td2 bash/TD3-blastp.sh"  
td3=$($command | awk ' { print $4 }')
echo '3) BLASTP submitted' $td3

#---
command="sbatch --dependency=afterok:$td2 bash/TD4-pfam.sh"  
td4=$($command | awk ' { print $4 }')
echo '4) PFAM submitted' $td4

#---
command="sbatch --dependency=afterok:$td3:$td4 bash/TD5-join.sh"  
td5=$($command | awk ' { print $4 }')
echo '5) Files joined' $td5

#---
command="sbatch --dependency=afterok:$td5 bash/TD6-run2.sh"  
td6=$($command | awk ' { print $4 }')
echo '6) TransDecoder part 2' $td6

#---
command="sbatch --dependency=afterok:$td6 bash/TD7-filt.sh"  
td7=$($command | awk ' { print $4 }')
echo '7) Filter results' $td7

echo '=================='
echo 'ALL SUBMITTED   :)'
echo '=================='

else

echo '==================='
echo 'Nothing submitted:('
echo '==================='

fi
