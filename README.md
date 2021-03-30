# TPS
Snakemake pipeline to mine terpene synthase sequences from databases

## Current databases
* UniParc (UniProt)
* Phytozome
* 1KP
* TSA

## Packages

* snakemake
* seqkit
* hmmer
* biopython
* transdecoder
* dplyr
* ggplot2
* mechanize

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
## Output tructure
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
│   ├── PF01397.hmm
│   ├── PF03936.hmm
│   └── PF19086.hmm
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
    └── uniprot_output.txt

```
