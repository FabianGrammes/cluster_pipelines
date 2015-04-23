# Ortholog prediction

This is a pieline to identify ortholog genes in A.salmon. Only one
gene candidate at a time

### Short pipeline summary:
1. Run `ortho_pipeline.R` this will identify potential orthologos in
   the A.salmon genome.
2. Manually kick out A.salmon gene IDs that are clear outliers (as
   seen in the phylogenetic tree).
3. Re-compute the phylogentic trees by running
   `ortho_pipeline-ReRun.R`
4. Collect final information in a table by running `collect_sets.R`

### Run Example

**Step 1**
```bash
module load blast+ anaconda mafft
Rscript ortho_pipeline.R ACO1 hgnc_symbol hsapiens_gene_ensembl Iron/trees
```
**Step 3**
```bash
module load blast+ anaconda mafft
Rscript ortho_pipeline-ReRun.R ACO1 Iron/trees
```
**Step 4**
```bash
Rscript collect_sets.R Iron/trees iron Iron/iron.set.txt
```

### HELP

#### ortho_pipeline.R

_Input:_

- Gene identifier; any unique gene identifier used within ensembl is
  possible (f.eks gene Symbol:`ACO1` or ensembl gene
  id:`ENSGALG00000008900`).
- Type of gene identifier used in Argument 1 (above). In case of a
  ensembl gene id this has to be `ensembl_gene_id`
- Organism to which the genes should be ortholog (autimatically gets
  all fish orthologs). F.eks chicken `ggallus_gene_ensemb`. Check the
  RBioMart documentation.
- Folder name to collect the results

_Output:_
 The script will produce the following files (* as specified by argument 4 e.g. test/vitello )
- *.blast.tab  = blastp Hits ortholog sequences VS A.salmon translated
  protein sequences (from TransDecoder)
- *.spike.fa   = .fasta file containing all ortholog sequences,
  including Ssa sequences idientified by blast
- *.spike.algn = mafft amino acid alignment file
- *.tree1.pdf  = phylogenetic tree figure (rooted, if a Lamprey
  sequence was found among the orthologous)
- *.tree1.pdf = tree file
- *. ssa.ids = Ssa gene ids of the ortolog ssa sequences. If there are
  clear outliers (identified in the tree), delete the gene id of the
  outlier from the file and run `ortho_pipeline-ReRun.R`.

#### ortho_pipeline-ReRun.R

_Input:_

- Gene identifier; for re-computing the phylogeny

- Folder name were the results from `ortho_pipeline.R` are located
