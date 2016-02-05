# Notes for SNP_call.sh development

## FASTA preparation
GATK needs a seqence dictionay and an indexed fasta to work. Create with
```
module load picard 
module load samtools 
sbatch --wrap "picard CreateSequenceDictionary R=umd3.1_rna-seq.fa O=umd3.1_rna-seq.dict"
sbatch --wrap "samtools faidx umd3.1_rna-seq.fa"
```
## Problems
At present time, fails with `INFO  15:58:12,385 MicroScheduler -   -> 10211039 reads (96.28% of total) failing MappingQualityUnavailableFilter `
