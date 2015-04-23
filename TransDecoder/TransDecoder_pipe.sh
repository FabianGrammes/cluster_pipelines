#!/bin/bash

# Input variables for the script: 
# 1 = Trinity transcriptome .fasta file
# 2 = How many chunks to run
# 3 = BLASTP DB
# 4 = PFAM DB (/local/genome/packages/transdecoder/rel16JAN2014/pfam)

inputFASTA=$1   # .fasta full path
JOBsplit=$2
blastDB=$3
pfamDB=$4
projectWD=$5

TransDwd=$(basename $inputFASTA)

# create project folder and use as working directory
mkdir $projectWD
cd $projectWD
mkdir slurm

#-------------------------------------------------------------------------------
# part 1: Run TransDecoder.

cat > sbatchTD1.sh << EOF 
#!/bin/bash
#SBATCH --job-name=TransDecoder
#SBATCH -n 3
#SBATCH --output=slurm/TransDecoder-%j.out
        
module load transdecoder/2.0.1

TD=/local/genome/packages/transdecoder/2.0.1/TransDecoder.LongOrfs

\$TD -t ${inputFASTA} -m 30
EOF

# run sbatch file
command="sbatch sbatchTD1.sh"
latest_id=$($command | awk ' { print $4 }')

echo '1) TransDecoder part 1 submitted'

#-------------------------------------------------------------------------------
# part 2: split .fasta in X chunks
 
cat > sbatchSPLIT.sh << EOF 
#!/bin/bash
#SBATCH --job-name=Split
#SBATCH -n 2
#SBATCH --output=slurm/FastaSplit-%j.out

module load cigene-tools/1.0

mkdir transcriptchunks
fastasplit -f ${TransDwd}.transdecoder_dir/longest_orfs.pep -o transcriptchunks -c ${JOBsplit}
EOF

command="sbatch --dependency=$latest_id sbatchSPLIT.sh"  # run sbatch script
latest_id=$($command | awk ' { print $4 }')

echo '2) FASTA split submitted'

#-------------------------------------------------------------------------------
# part 3: Run BlastP on each chunk 

# make blastp bash
mkdir BLAST/

cat > sbatchBLASTP.sh << EOF 
#!/bin/sh
#SBATCH -n 2 # NB CPUS
#SBATCH -J BLASTP # JOB NAME
#SBATCH --array=1-${JOBsplit} # array number
#SBATCH --output=slurm/BlastP-%j-%a.out

module load blast+/2.2.28
  
SGE_TASK_ID=\$(expr \$SLURM_ARRAY_TASK_ID - 1)
printf -v ID '%07g' \$SGE_TASK_ID              # will produce feks 0000001 (matches extension produced by fastasplit)

blastp -query transcriptchunks/longest_orfs.pep_chunk_\$ID \
  -db ${blastDB} \
  -outfmt 6 -num_threads 2 \
  -evalue 1e-10  -max_target_seqs 1 -out BLAST/BLAST_out_\$ID
EOF

command="sbatch --dependency=afterok:$latest_id sbatchBLASTP.sh"  # run sbatch script
blast_job=$($command | awk ' { print $4 }')

echo '3) BLASTP submitted'

#-------------------------------------------------------------------------------
# part 4: Run PFAM

mkdir PFAM/

cat > sbatchPFAM.sh << EOF 
#!/bin/sh
#SBATCH -n 2 # NB CPUS
#SBATCH -J PFAM # JOB NAME
#SBATCH --array=1-${JOBsplit}# array number
#SBATCH --output=slurm/Pfam-%j-%a.out

module load hmmer/3.1b1  
  
SGE_TASK_ID=\$(expr \$SLURM_ARRAY_TASK_ID - 1)
printf -v ID '%07g' \$SGE_TASK_ID              # will produce feks 0000001 (matches extension produced by fastasplit)

hmmscan --cpu 2 --domtblout PFAM/PFAM_out_\$ID \
 $pfamDB \
 transcriptchunks/longest_orfs.pep_chunk_\$ID
EOF

command="sbatch --dependency=afterok:$latest_id sbatchPFAM.sh"  # run sbatch script
pfam_job=$($command | awk ' { print $4 }')

echo '4) PFAM submitted'

#-------------------------------------------------------------------------------
# part 5: Join chunks  
cat > sbatchJOIN.sh << EOF 
#!/bin/sh
#SBATCH -n 1 # NB CPUS
#SBATCH -J Join # JOB NAME
#SBATCH --output=slurm/Join-%j.out

module load R

cat BLAST/BLAST_out_* > BLASTP_concatenated.tab # this must be run in the wdir
cat PFAM/PFAM_out_* > PFAM_concatenated.tab # this must be run in the wdir
EOF

command="sbatch --dependency=afterok:$pfam_job:$blast_job sbatchJOIN.sh"  # run sbatch script
latest_id=$($command | awk ' { print $4 }')

echo '5) Files joined'

#-------------------------------------------------------------------------------
# part 6: Integrating the Blast and Pfam search results into coding region selection

cat > sbatchTD2.sh << EOF 
#!/bin/bash
#SBATCH --job-name=TransDecoder2
#SBATCH -n 3
#SBATCH --output=slurm/TransDecoder2-%j.out
        
module load transdecoder/2.0.1

TD=/local/genome/packages/transdecoder/2.0.1/TransDecoder.Predict

\$TD -t $inputFASTA --retain_pfam_hits PFAM_concatenated.tab \
 --retain_blastp_hits BLASTP_concatenated.tab \
 --retain_long_orfs 120
EOF

command="sbatch --dependency=afterok:$latest_id sbatchTD2.sh"  # run sbatch script
latest_id=$($command | awk ' { print $4 }')

echo 'ALL GOOD :)'
