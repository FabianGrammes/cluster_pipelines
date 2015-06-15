# count_multi.py

Python scripts modified from the HTSeq count script. The script is
intended to work **only** on name sorted `.sam ` files. It maps and
counts reads that map to exactly 2 different genomic locations. 



### Run Example
In progress..
```bash
python count_multi.py -g Salmon_3p6_Chr_051214.gtf -i example.sam -o example_out.count
```


The counts actually add up; however it's imoprtant to keep in mind
that :
* `__Single_hit:Not_aligned` does NOT count to the to `__Total_reads`
* `__Double_hit:Not_aligned` does count to the to `__Total_reads`



### HELP

_Input:_

1. `.gtf file`
2. `.sam file (name sorted!)` 
3. `output file`
