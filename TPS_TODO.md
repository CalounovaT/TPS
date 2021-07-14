# TPS Project TODO

## Outline

- [x] use new Pfam HMMs for search
- [ ] delete duplicated sequences, keep primarly UniProt sequences
- [ ] use new models (based on Adelas db) to assign TPS class
- [ ] create filter to keep complete sequences; the filtering will be based on:
	- length of a sequence (too short sequences are not complete)
	- domain architecture (whether the sequence contains complete domains)
- [ ] count how many sequences there are after filtering; do some statistics
- [ ] annotate UniProt sequences so they can be added to Adelas db (with a script inspired by Martin)
- [ ] assign phylogenetic class based on a phylogenetic tree


## 1) New Pfam models

* to find new models I looked in Adelas db for tri-, sester- and sesquarterpenes
* using downloaded PfamA db I searched what domains these terpenes contain
* from the domains after discussion we selected 3 other domains
    * PF13243 &emsp; Squalene-hopene cyclase C-terminal domain &emsp; https://pfam.xfam.org/family/PF13243
    * PF13249 &emsp; Squalene-hopene cyclase N-terminal domain &emsp; https://pfam.xfam.org/family/PF13249
    * PF00494 &emsp; Squalene/phytoene synthase	&emsp; https://pfam.xfam.org/family/PF00494


## 2) Delete duplicated sequences

* can use seqkit rmdup to remove duplicated reads
* there is not stated in manual based on what criteria the sequences are removed... but based on what I tried it seems that it is based on order
* in that case putting UniProt results first should ensure they are not removed

## 3) Assign TPS class

* already implemented I just have to run it on new results (using HMMs created based on Adelas db, not the Terzyme ones)

## 4) Create filter to keep complete sequences

* remove sequences shorther than a certain threshold (probably something like 100 aa)
* look at the domain architecture (could be helpful - https://www.biostars.org/p/134579/, http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf, https://biopython.org/docs/1.75/api/Bio.SearchIO.HmmerIO.html)
    * do hmmscan against the set of the TPS Pfam HMMS
    * positions of where the  domain starts and ends corresponds to `envfrom, env to` in the output
    * how big part of the hmm was used corresponds to `hmmfrom  hmm to` in the output, that is the criteria based on which the filtering will be made (tells if the domain is complete)
        * one can find full length of the pfam model on the web of Pfam
            * PF01397
                * full length: 192
                * average length: 162.70
                * average coverage of the sequence by the domain: 31.19 %
            * PF03936
                * full length: 267
                * average length: 200.00
                * average coverage of the sequence by the domain: 43.13 %
            * PF19086
                * full length: 199
                * average length: 186.00
                * average coverage of the sequence by the domain: 46.25 %
            * PF13243
                * full length: 319
                * average length: 280.90
                * average coverage of the sequence by the domain: 44.12 %
            * PF13249
                * full length: 291
                * average length: 270.10
                * average coverage of the sequence by the domain: 40.48 %
            * PF00494
                * full length: 263
                * average length: 252.80
                * average coverage of the sequence by the domain: 74.27 %

        * if `hmmfrom and hmm to` is the same as the length, then the domain is complete
        * otherwise have to decide...
    * IF sequence contains at least one complete domain -> complete (however sequence containing just a single domain may not be functional)
    * ELSE -> incomplete
    * can do some statistics of domain architectures
    * can compare the architectures to some architecture of TPS from well annotated genome
    * <sup><sub>can look at secondary structures of the domains (CDD db), predict secondary structure of the sequences and compare if they contain the parts of the domains (if small helix on one end is missing it probably is not that important)</sub></sup>

## 5) Count how many sequences survived filtering + statistics

* how many sequences survived
* length distribution
* how many sequences are in the groups

## 6) annotate the UniProt sequences

* look if there is all necessary information so the sequence can be added to Adelas db
* contact UniProt people and ask what is the best way to retrieve such data
* needed parts:
    * UniProt ID
    * Name
    * Sequence
    * Species
    * Kingdom (will look it up based on the specie)
    * Type (determined based on the formula)
    * Substrate
    * Substrate ChEBI ID
    * Cofactors
    * Product
    * Product formula
    * Product SMILES
    * Product ChEBI ID
    * Cyclic (determined using rdkit)
* https://www.uniprot.org/uniprot/?query=id:%22Q9SAK2%22&format=tab&columns=id,entry_name,organism,protein_names,sequence,chebi,chebi(Catalytic%20activity),catalytic_activity,chebi-id,chebi(Cofactor)

## 7) assign phylogenetic group

* use some phylogenetic tree to assign them




