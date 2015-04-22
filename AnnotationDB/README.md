# make Annotation DB 

This is a pieline to produce a SQLite DB from a bunch of Gene related
results


### Run Example

```bash
$Rscript makeDB3.R merged.gtf TD_concatenated.bed TD_concatenated.pfam BLASTP_concatenated.tab Test.db 
```


### HELP

_Input:_

0. Annotation data
1. `merged.gtf` 
2. `TD_concatenated.bed` (Transdecoder .bed file)
3. `TD_concatenated.pfam` (Transdecoder .pfam file)
4. `BLASTP_concatenated.tab` (BLAST P hit file; keep column order in mind)
5. `topBlast.txt` (from Blast2GO exported)
6. `GO annotation` (from Blast2GO exported)
7. `SNP`
8. `Transpososns`
9.  KEGG annotation Gene id 2 K(from KAAS)
10. KEGG K 2 path ko
11. KEGG path ko 2 name
12. DB name (f.eks `SsaCIG_annot.001.db` )

