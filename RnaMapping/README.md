# RNAmap

`RNAmap` is a bash script that automates the processes of:
1. Trimming
2. Quality Control
3. Mapping
4. Read counting per gene

The script currently works only if the following conditions are met
(possible to extend the script for other settings):
- Libraries were prepared using a *TruSeq stranded Kit*
- fastq data was generated on a MiSeq

## Example


```bash
script=/mnt/users/fabig/cluster_pipelines/RnaMapping/RNAmap.sh
fdir=/mnt/SeqData2/MiSeq_data/150612_M02210_0008_000000000-ACGHJ/Data/Intensities/BaseCalls
sdir=/mnt/SeqData2/MiSeq_data/150612_M02210_0008_000000000-ACGHJ/Data/Intensities/BaseCalls/SampleSheet.csv

bash $script -d $fdir -s $sdir -m SampleMaster.txt
```

## RNAmap options

- `-d|--dirin`: Full path of the dierectory that holds the .fasq.gz files
- `-s|--samplesheet`:  Full path to the Illumina sample sheet (.`.csv` format)
- `-g|--genome`: *Optional*, defaults to the latest version of the _Salmon salar_ genome which is
located at: `/mnt/backup2/users/aquagenome/Ssa/Synteny/CIGENE-3.6-unmasked/STAR/reference/cigene3.6v2-unmasked`
- `-m|--mastersheet`: Name of the master-sheet, this will be created in the current
dierectory and contains all relevant sample information. Information is extracted from
the SampleSheet and arranged in a more user friendly way.
- `--gtf`: *Optional*, defaults to the latest _Salmon salar_ gtf which islocated at:
`/mnt/users/fabig/Ssa_transcriptome/cigene3.6_chrom/gtf/Salmon_3p6_Chr_051214.gtf`


`RNAmap` creates a folder tree in the current dierectory:

- `slurm` folder to collect the slurm reports. *Note* each slurm report lists:
   - date of job submission
   - version of the used module 
- `bash` folder to collect the bash scripts, genesrated and submitted through `RNAmap`.
- `fastq_trim`folder for the adaptor and quality trimmed .fastq files.
- `qc` folder for the quality control reports of the quality trimmed .fastq files.
- `star` folder contains read alignments to the genome in `.BAM` format. Reads were aligned using the program `STAR`
- `count` folder contains read counts per gene genrated by the `HTSeq-count` script. 
