# TPS scoring

Scoring of the TPS sequences using phylogenetic analysis and knowledge about characteristic TPS sequences to prioretize the sequences for further analysis.

## Used sequences
* mined tps (`tps_mining.fasta` - corresponds to file `final_unique.fasta` from the mining part)
* characterized sesqTPS and diTPS from our TPS database (`tps_characterized.fasta`)

## Phylogenetic tree construction:
* MSA using MAFFT
* MSA trimming using trimAl
* phylogenetic tree construction using Fasttree2

## Distance calculation in the tree and scoring

### Distance calculation
For all uncharacterized sequences, the distance to closest characterized sequence in the tree was calculated (sum of the lengths of the branches, i.e. the closest characterized sequence has the shortest distance). 

### Scoring

The scoring should reflect the probability that the sequence is sesqTPS or diTPS and also prioritize the sequences that are more distant from the characterized sequences in the tree.

The sequences were classified into types (mono/sesq/di/tri) using custom class HMMs.

There were following scoring criteria:
* `completness`: presence of methionin as a first amino acid
* `type`: predicted type
* `length`: length of the sequence similar to characterized sesqTPS and diTPS (using information about the predicted type)
* `distance`: normalized distance to closest characterized sequence in the tree
* `architecture`: Pfam domain architecture (using information about the predicted type)

All the criteria had value between 0 and 1. The criteria are described in more detail below and exact code can be found in `Snakefile`.

The final score was calculated followingly:
`total_score = completness + type + length + 2*distance + architecture`

The distance was multiplied by 2 to have more impact on the final score given the assumption that more distant sequences will give more interesting results.

#### Presence of methionin as a first amino acid
* 1 if methionin is present as a first amino acid
* 0 if methionin is not present as a first amino acid

#### Predicted type
* 1 if the predicted type is diTPS or sesqTPS
* 0 otherwise

#### Length of the sequence similar to characterized sesqTPS and diTPS
* score between 0 and 1 was determined using knowledge about the length distribution of characterized sesqTPS and diTPS
* if the length of predicted sesqTPS was similar to lengths of characterized sesqTPS, the score was closer to 1, otherwise it got lower values
* if the length of predicted diTPS was similar to lengths of characterized diTPS, the score was closer to 1, otherwise it got lower values
* if the sequence was not predicted as sesqTPS or diTPS but had length similar to what is observed for sesqTPS and diTPS, the score was closer to 1, otherwise it got lower values

#### Normalized distance to closest characterized sequence in the tree
* the distance was normalized using the maximum distance and minimum distance in the tree to get values between 0 and 1
* `distance_score(d) = (d-min_dist)/(max_dist-min_dist)`

#### Pfam domain architecture
* if the sequence was predicted as sesqTPS:
    * if the sequence had common architecture for sesqTPS: 1
    * if the sequence had less common but observed architecture for sesqTPS: 0.5
    * if the sequence had architecture not observed in sesqTPS: 0
* if the sequence was predicted as diTPS:
    * if the sequence had common architecture for diTPS: 1
    * if the sequence had less common but observed architecture for diTPS: 0.5
    * if the sequence had architecture not observed in diTPS: 0
* if the sequence was not predicted as sesqTPS or diTPS:
    * if the sequence had common architecture for sesqTPS or diTPS: 0.5
    * if the sequence had less common but observed architecture for sesqTPS or diTPS: 0.25
    * if the sequence had architecture not observed in sesqTPS or diTPS: 0

## Results

The final result is a table with all the sequences and their scores. The sequences are sorted by the total score in descending order. The table can be found in `annotated_scored_candidates.tsv`.