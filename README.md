# PIPELINES

Collection of analysis pipelines for the cigene cluster.

## TransDecoder_pipe.sh

This is a pieline to run `TransDecoder` from the `Trinity` suite. The pipeline splits the BLAST and PFAM jobs automatically.

### Run exsample
Script expects:

1. A transcriptome `.fasta` file (full path)
2. Number of chunks to split into
3. A path to the blast DB ( currently `/mnt/NOBACKUP/db/swissport/swissprot` )
4. A path to the pfam DB (currently `/local/genome/packages/transdecoder/rel16JAN2014/pfam/Pfam-AB.hmm.bin`)
5. A name for the project

**Example**

```bash
fasta=/mnt/users/fabig/Side_Projects/strawberry_2014/trinity/CT3Trinity.expr.fasta
bdb=/mnt/NOBACKUP/db/swissport/swissprot
pdb=/local/genome/packages/transdecoder/rel16JAN2014/pfam/Pfam-AB.hmm.bin

sh ~/cluster_pipelines/TransDecoder_pipe.sh $fasta 8 $bdb $pdb TD_pipeline

```

### Short pipeline summary:
1. Run `TransDecoder.LongOrfs` to identify ORFs. running on option `-m 30`.
2. Split the TransDecoder output file `longest_orfs.pep` into X chnks.
3. Run `blastp` on each chunk.
4. Run `hmmer` (PFAM) on each chunk.
5. Merge the chunks.
6. Run `TransDecoder.Predict` using the results from step 1, 3 and 4.




