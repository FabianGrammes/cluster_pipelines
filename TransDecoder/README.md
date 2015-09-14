# TransDecoder_pipe.sh

This is a pieline to run `TransDecoder` from the `Trinity` suite. The pipeline splits the BLAST and PFAM jobs automatically.

## input
Script expects:

1. A transcriptome `.fasta` file (full path)
2. Number of chunks to split into
3. A path to the swissprot blast DB ( currently
   `/mnt/NOBACKUP/db/swissport/swissprot` )
4. A path to the pfam DB (currently
   `/local/genome/packages/transdecoder/rel16JAN2014/pfam/Pfam-AB.hmm.bin`)
   


**Example**

First we generate a small subset of a transcriptome, which we run
subsequently through the pipeline. 

```bash
fasta=/mnt/users/fabig/Side_Projects/strawberry_2014/trinity/CT3Trinity.expr.fasta
cat $fasta | head -20000 > test_transcriptome.fa

scr=/mnt/users/fabig/cluster_pipelines/TransDecoder/TransDecoder_pipe.sh

bdb=/mnt/NOBACKUP/db/swissport/swissprot
pdb=/local/genome/packages/transdecoder/rel16JAN2014/pfam/Pfam-AB.hmm.bin

bash $scr -f test_transcriptome.fa -s 8 -b $bdb -p $pdb
```

## Short pipeline summary:

1. Run `TransDecoder.LongOrfs` to identify ORFs. running on option `-m 30`.
2. Split the TransDecoder output file `longest_orfs.pep` into X chnks.
3. Run `blastp` on each chunk.
4. Run `hmmer` (PFAM) on each chunk.
5. Merge the chunks.
6. Run `TransDecoder.Predict` using the results from step 1, 3 and 4.

## Output

2 Folders will be generated in your current path:

1. `td-out/` Contains the final TransDeoder models
2. `td-final/` Contains the (probably) most relevant files:
   1. `LOC_longest_all.annot.tab`: Table containing the best TD models
      per gene
   2. `LOC_longest_support.annot.tab`Table containing the best TD models
      per gene, that have some support (e.g BLASTP or PFAM hit or both).
   3. `TRA_longest_all.annot.tab`
   4. `TRA_longest_support.annot.tab`
   5. `<yor transcritome>.transdec.filt.cds` filtered for the best
      models (per gene), can be used for annotation ...
   6. `<yor transcritome>..fa.transdec.filt.fa` same as above
   7. `<yor transcritome>..fa.transdec.filt.pep`


# Annotation_pipe.sh

Pipeline to find *gene names*  for genes.

## Input
Script expects:

1. A transcriptome `.fasta` file containing (predicted) protein
   sequences.
2. Number of chunks to split into
3. DB file see below



### Input file
A `.tab` delimited file with a header (see below)

Species | Proteome |  FTP  | Mart |  Prefix 
------|---------|----|----|------
Arabidopsis.thaliana   | Arabidopsis_thaliana.TAIR10.28.pep.all.fa |
ensembl FTP | athaliana_eg_gene | ATMG

- *Species*  `Arabidopsis.thaliana` (Species name, generic and
  specific name should be seperated by a `.` NO space! )
- *Proteome* `Arabidopsis_thaliana.TAIR10.28.pep.all.fa` ()
- *FTP*
	`ftp://ftp.ensemblgenomes.org/pub/release-28/plants/fasta/arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.28.pep.all.fa.gz`
	(For ensembl proteoms the list can be found at `http://www.ensembl.org/info/data/ftp/index.html`)
- *Mart* `athaliana_eg_gene`
- *Prefix* `ATMG` 
