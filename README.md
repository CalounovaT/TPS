# TPS
Snakemake pipeline to mine terpene synthase sequences from databases

## Current databases
* UniParc (UniProt)
* Phytozome
* 1KP
* TSA

## Current Pfam HMMs
* [PF01397](https://pfam.xfam.org/family/PF01397) (Terpene synthase, N-terminal domain)
* [PF03936](https://pfam.xfam.org/family/PF03936) (Terpene synthase family, metal binding domain)
* [PF19086](https://pfam.xfam.org/family/PF19086) (Terpene synthase family 2, C-terminal metal binding)
* [PF13243](https://pfam.xfam.org/family/PF13243) (Squalene-hopene cyclase C-terminal domain)
* [PF13249](https://pfam.xfam.org/family/PF13249) (Squalene-hopene cyclase N-terminal domain)
* [PF00494](https://pfam.xfam.org/family/PF00494) (Squalene/phytoene synthase)

## Packages

* snakemake
* seqkit
* hmmer
* biopython
* transdecoder
* dplyr
* ggplot2

## How to run it
###
Clone the repository
```
git clone https://github.com/CalounovaT/TPS.git
```
Download *.protein.fa.gz files from Phytozome

Extract the .zip in the folder db/Phytozome, don't keep the direcory structure (use `unzip -j`)

Run the main script
```
./run_script.sh
```
## Output structure
```
TPS
├── README.md
├── Snakefile
├── config.yaml
├── map.py
├── run_script.sh
├── status.py
├── tsa_mapping3.tsv
├── db
│   ├── 1kp/                          - contains 1kp *translated-protein.fa.gz files
│   ├── Phytozome/                    - contains Phytozome *.protein.fa.gz files
│   ├── TSA/                          - contains TSA *.fsa_nt.gz files
│   ├── uniparc/                      - contains uniparc_active.fasta.gz file
│   └── uniparc_chunks/               - split uniparc_active.fasta.gz 
├── hmm_profiles
│   ├── PF01397.hmm (Terpene synthase, N-terminal domain)
│   ├── PF03936.hmm (Terpene synthase family, metal binding domain)
|   ├── PF19086.hmm (Terpene synthase family 2, C-terminal metal binding)
|   ├── PF13243.hmm (Squalene-hopene cyclase C-terminal domain)
|   ├── PF13249.hmm (Squalene-hopene cyclase N-terminal domain)
│   └── PF00494.hmm (Squalene/phytoene synthase)
├── hmm_results                       - contains results of hmmsearch on all files 
│   ├── 1kp/
│   ├── phytozome/
│   ├── tsa/
│   └── uniparc/
├── filtered_results                  - contains filtered files from db/ directory
│   ├── 1kp/
│   ├── phytozome/
│   ├── tsa/
│   └── uniprot/
└── results                           - contains TPS sequences filtered from databases
    ├── 1kp_output.fasta              - FASTA file with header on one line and sequence on the line below
    ├── 1kp_output_long.fasta         - classic FASTA file
    ├── 1kp_output.txt                - tsv format
    ├── phytozome_output.fasta
    ├── phytozome_output_long.fasta
    ├── phytozome_output.txt
    ├── TSA_output.fasta
    ├── TSA_output_long.fasta
    ├── uniprot_output.fasta
    |── uniprot_output.txt
    └── final.tsv                     - final tsv table with all results together

```
