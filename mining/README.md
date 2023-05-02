# TPS mining

Mining of the TPS sequences from databases using Pfam HMMs with Snakemake pipeline.

The source code of the pipeline can be found in `Snakefile`.

## Used databases

* UniParc (UniProt)
* Phytozome
* 1KP
* TSA

## Used Pfam HMMs

* [PF01397](https://pfam.xfam.org/family/PF01397) (Terpene synthase, N-terminal domain)
* [PF03936](https://pfam.xfam.org/family/PF03936) (Terpene synthase family, metal binding domain)
* [PF19086](https://pfam.xfam.org/family/PF19086) (Terpene synthase family 2, C-terminal metal binding)
* [PF13243](https://pfam.xfam.org/family/PF13243) (Squalene-hopene cyclase C-terminal domain)
* [PF13249](https://pfam.xfam.org/family/PF13249) (Squalene-hopene cyclase N-terminal domain)
* [PF00494](https://pfam.xfam.org/family/PF00494) (Squalene/phytoene synthase)

## Mining description

Any sequence containing any of the Pfam HMMs was mined. 

## Output files description

Output files can be found in the folder `results/`.

There are `.tsv` and `.fasta` files for individual databases containing the mined sequences.

Then there are also `final.tsv` and `final.fasta` file containing all the mined sequences from all the databases.

To have only unique sequences, files `final_unique.tsv` and `final_unique.fasta` were created. These files were used for the phylogenetic analysis. In total, there were 191,476 sequences.