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
    -d|--db) 
	DBFILE="$2"  # path to BLASTP DB
	if [ ! -f "$DBFILE"  ]; then
	    echo 'ERROR: File' $DBFILE 'Does not exist!'
	    exit 1
	fi
	shift # past argument
	;;
    --mart) 
	MART="$2"  # path to BLASTP DB
	if [ -z "$MART"  ]; then
	    MART=plants_mart_28,description,external_gene_id
	fi
	shift # past argument
	;;
    --swiss) 
	SWISS="$2"  # path to BLASTP DB
	if [ -z "$SWISS"  ]; then
	    echo 'ERROR: Need a SWISSPROTDB'
	    exit 1
	fi
	shift # past argument
	;;
    --execute) # Only used for testing!
	EXECUTE="$2"
	shift # past argument
	;;
esac
shift # past argument or value
done

#------------------------------------------------------------------------
# Check the DB file 

if cat -v $DBFILE | grep -q '\^M'
then
    echo 'Converting line terminators'
    sed 's/\r/\n/g' $DBFILE > sheet.tmp
    mv sheet.tmp $DBFILE
else
    echo 'Line terminators seem correct'
fi


#--------------------------------------------------------------------------
# echo all input variables

# Number of samples
END=$(sed '1d' $DBFILE | wc -l) # skip header line

# FIX the R input arguments
Marray=(${MART//,/ })

# assign relative .fasta name
RELFA=$(basename $FASTA)

echo '--------------------------------------------------------------------------'
echo 'Input arguments:'
echo $FASTA
echo 'Chunks:' $SPLIT
echo $SWISS
echo '--------'
echo 'MART arguments:'
for i in "${!Marray[@]}"
do
    echo "${Marray[i]}"
done
echo '--------'
echo 'Number of proteoms' $END
echo 'FIRST proteome:' $(awk ' NR=='2' { print $1 ; }' $DBFILE)
echo 'LAST proteome:' $(awk ' NR=='$END+1' { print $1 ; }' $DBFILE)
echo '--------------------------------------------------------------------------'

# CREATE FOLDER STRUCTURE
mkdir -p {ann-db,ann-chunks,ann-hits,bash,slurm}

#------------------------------------------------------------------------
# DOWNLOAD ensembl information

# importnat
cat > bash/ANN1-wget.sh << EOF 
#!/bin/bash
#SBATCH --job-name=ANN.wget:XofX
#SBATCH -n 1
#SBATCH --output=slurm/ANN_wget-%A.out

# wget
for N in {2..$END}
do 
  FP=\$(awk 'NR=='\$N' {print \$3\$2}' $DBFILE)
  OP='ann-db/'\$(awk 'NR=='\$N' {print \$2}' $DBFILE)
  echo \$FP
  wget -O \$OP \$FP
done

# unzip
cd ann-db

for i in \$(ls *.fa.gz)
do 
  gunzip \$i
done 

for fa in \$(ls *.fa)
do
  size=\$(cat \$fa | grep '^>' | wc -l)
  echo \$fa \$size >> DB.size
done

cd ..

EOF

#------------------------------------------------------------------------
# BLASTDB

cat > bash/ANN2-bdb.sh << EOF 
#!/bin/sh
#SBATCH --job-name=ANN.bdb:XofX
#SBATCH -n 1
#SBATCH --output=slurm/ANN_blastdb-%A.out

module load blast+/2.2.28   
    
cd ann-db

cat \$(ls *.fa) >> collection_all.fa

makeblastdb -in collection_all.fa -dbtype prot

cd ..
EOF

#------------------------------------------------------------------------
# SPLIT input fasta

cat > bash/ANN3-split.sh << EOF
#!/bin/sh
#SBATCH --job-name=ANN.split:XofX
#SBATCH -n 1
#SBATCH --output=slurm/ANN_split-%A.out
  
module load cigene-tools/1.0   

fastasplit -f $FASTA -o ann-chunks -c $SPLIT

echo 'All Good'
EOF

#------------------------------------------------------------------------
# ERROR blast does not take the -db argument from variable
# Set DBS variable outside EOF
DBS="' "'"'$SWISS'" ''"ann-db/collection_all.fa"'" '"

# BLASTP
cat > bash/ANN4-blastp.sh << EOF 
#!/bin/sh
#SBATCH -n 2 # NB CPUS
#SBATCH --job-name=ANN.blastp:XofX
#SBATCH --array=1-$SPLIT # array number
#SBATCH --output=slurm/ANN_blastp-%A-%a.out

module load blast+/2.2.28
  
SGE_TASK_ID=\$(expr \$SLURM_ARRAY_TASK_ID - 1)
printf -v ID '%07g' \$SGE_TASK_ID              # will produce feks 0000001 (matches extension produced by fastasplit)

CHUNK=ann-chunks/$RELFA'_chunk_'\$ID

blastp -query \$CHUNK \
-db $DBS \
-outfmt "6 std qlen slen nident" \
-num_threads 2 \
-evalue 1e-5 \
-max_hsps_per_subject 1 \
-max_target_seqs 20 \
-out ann-hits/hits_chunk_\$ID

echo 'All Good'
EOF

#------------------------------------------------------------------------
# GENE NAMES

cat > bash/ANN5-mart.sh << EOF 
#!/bin/sh
#SBATCH -n 1 # NB CPUS
#SBATCH --job-name=ANN.mart:XofX
#SBATCH --array=1-$SPLIT # array number
#SBATCH --output=slurm/ANN_mart-%A-%a.out

module load R  

RSCR=/mnt/users/fabig/cluster_pipelines/TransDecoder/helper_scripts/Annotation_biomaRt.R

SGE_TASK_ID=\$(expr \$SLURM_ARRAY_TASK_ID - 1)
printf -v ID '%07g' \$SGE_TASK_ID 

Rscript \$RSCR ann-hits/hits_chunk_\$ID ann-hits/hits_chunk_\$ID.names $DBFILE ${Marray[0]} ${Marray[1]} ${Marray[2]}

echo 'All Good'
EOF

#-----------------------------------------------------------------------
# SUMMARIZE

cat > bash/ANN6-summary.sh << EOF 
#!/bin/sh
#SBATCH -n 1 # NB CPUS
#SBATCH --job-name=ANN.sum
#SBATCH --output=slurm/ANN_summary-%A.out

module load R

# join the files
awk 'FNR==1 && NR!=1{next;}{print}' \$(ls ann-hits/*.names) > ann-hits/All_names.txt

RSCR=/mnt/users/fabig/cluster_pipelines/TransDecoder/helper_scripts/sum-annot.R
Rscript \$RSCR ann-hits/All_names.txt ann-hits/All_names.pdf
EOF



#===============================================================================
# Submit the jobs to SLURM

if [ "$EXECUTE" != "no" ]
then

    command="sbatch bash/ANN1-wget.sh"
    td1=$($command | awk ' { print $4 }')
    echo '1) Annotation: fetching proteoms' $td1

    #---
    command="sbatch --dependency=$td1 bash/ANN2-bdb.sh"  
    td2=$($command | awk ' { print $4 }')
    echo '2) Annotation: stting up DBs' $td2
    
    #----
    command="sbatch --dependency=afterok:$td2 bash/ANN3-split.sh"  
    td3=$($command | awk ' { print $4 }')
    echo '3) splitting files' $td3

    #---
    command="sbatch --dependency=afterok:$td3 bash/ANN4-blastp.sh"  
    td4=$($command | awk ' { print $4 }')
    echo '4) BLASTP submitted' $td4

    #---
    command="sbatch --dependency=afterok:$td3:$td4 bash/ANN5-mart.sh"  
    td5=$($command | awk ' { print $4 }')
    echo '5) Running biomaRt' $td5

    #---
    command="sbatch --dependency=afterok:$td5 bash/ANN6-summary.sh"  
    td6=$($command | awk ' { print $4 }')
    echo '5) Summarizing' $td6

    echo '=================='
    echo 'ALL SUBMITTED   :)'
    echo '=================='

else

    echo '==================='
    echo 'Nothing submitted:('
    echo '==================='

fi
