# TransDecoder_pipe.sh

This is a pieline to run `TransDecoder` from the `Trinity` suite. The
pipeline splits the BLAST and PFAM jobs automatically.

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
Script arguments:

- `-/--fasta`: A transcriptome `.fasta` file containing (predicted) protein
   sequences.
- `-s/--split`: Number of chunks to split into
- `-d/-db`: Input DB file (see below)
- `--mart`: Specific arguments to BioMart. 3 Arguments seperated
  by. Defaults to `plants_mart_28,description,external_gene_id`
- `--swiss`: Path to the swissprot db.
- `--execute`: If set to <no> all folders/scripts will be genearted,
  but none of the jobs will be executed.


## Requirments

In order for all subscripts to work the following `R` packages have to
be installed in your R environment:

- `bioconductor`
- `biomaRt`
- `plyr`


### Input DB file
A `.tab` delimited file with a header (see below), *important* keep
a file header with lower case column names as in the example below. 

species | proteome |  ftp  | mart |  prefix 
------|---------|----|----|------
Arabidopsis.thaliana   | Arabidopsis_thaliana.TAIR10.28.pep.all.fa.gz | ensembl FTP | athaliana_eg_gene | ATH

- *Species*  `Arabidopsis.thaliana` (Species name, generic and
  specific name should be seperated by a `.` NO space! )
- *Proteome* `Arabidopsis_thaliana.TAIR10.28.pep.all.fa.gz` ()
- *FTP*
	`ftp://ftp.ensemblgenomes.org/pub/release-28/plants/fasta/arabidopsis_thaliana/pep/`
	(For ensembl proteoms the list can be found at `http://www.ensembl.org/info/data/ftp/index.html`)
- *Mart* `athaliana_eg_gene`
- *Prefix* choose a *three* capital letter, unique identifyer for the species. Do
  NEVER use *gi* since this is reserved for swissprot. 

### Example

```sh
fasta=Test_Protein.fasta
db=Ann-DBs.txt
swiss=/mnt/NOBACKUP/db/swissport/swissprot

script=~/cluster_pipelines/TransDecoder/Annotation_pipe.sh

bash $script -f $fasta -s 5 -d $db --mart plants_mart_28,description,external_gene_id --swiss $swiss --execute no
```
