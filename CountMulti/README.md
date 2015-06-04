# count_multi.py

Script to count multimapping reads in a name sorted sam/bam file. The scripts currently
only consideres reads that mapp to 2 positions in the genome. These reads are normally
labled as anbigous (in the HTSeq script), however here they are counted.



### Run Example
In progress..
```bash
$Rscript makeDB3.R merged.gtf TD_concatenated.bed TD_concatenated.pfam BLASTP_concatenated.tab Test.db 
```


### HELP

_Input:_

1. `.gtf file`
2. `.sam file (name sorted!)` 
3. `output file`
