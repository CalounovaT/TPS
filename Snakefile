import os
import subprocess
import urllib.parse
import urllib.request
import Bio.SeqIO


configfile: "config.yaml"

#rule localrules: all, clean
# 1KP project   https://db.cngb.org/onekp/
rule all:
    input:
        "db/uniparc/uniparc_active.fasta.gz",
        "db/uniparc_chunks",
        "db/Phytozome",
        "db/1kp",
        "hmm_profiles/PF01397.hmm",
        "hmm_profiles/PF03936.hmm",
        "hmm_profiles/PF19086.hmm",
        "results/uniprot_output.txt",
        "results/1kp_output.txt",
        "results/phytozome_output.txt",
        "results/TSA_output.txt",
        "results/final.tsv"

rule clean:
    shell:
        """
        rm -rf hmm_results logs hmm_combined.out output.txt uniparc_chunks uniprot UPIs upi.txt tmp.list 
        """

## DOWNLOADING THE DATABASES
## Phytozome database transcriptomes were downloaded manually from the website and the file is located in db/

rule download_uniparc:
    output:
        "db/uniparc/uniparc_active.fasta.gz"
    log:
        "logs/download_uniparc.log"
    params:
        memory="2"
    threads:
        1
    shell:
        """
        wget --directory-prefix "db/uniparc" "https://ftp.expasy.org/databases/uniprot/current_release/uniparc/uniparc_active.fasta.gz" &> {log}
        """

checkpoint download1KP:
    output:
        directory("db/1kp")
    log:
        "logs/download_1kp"
    params:
        memory="2"
    threads:
        1
    shell:
        """
        wget --directory-prefix "db/1kp" -r -l2 --no-parent -nd -A "translated-protein.fa.gz" ftp://parrot.genomics.cn/gigadb/pub/10.5524/100001_101000/100627/assemblies/
        """

checkpoint download_TSA:
    output:
        directory("db/TSA")
    log:
        "logs/download_TSA.log"
    params:
        memory="2"
    threads:
        2
    shell:
        """
        wget -r -np -nd --directory-prefix "db/TSA" ftp://ftp.ncbi.nlm.nih.gov/genbank/tsa/H/*fsa_nt.gz 2> {log}
        wget -r -np -nd --directory-prefix "db/TSA" ftp://ftp.ncbi.nlm.nih.gov/genbank/tsa/G/*fsa_nt.gz 2>> {log}
        wget -r -np -nd --directory-prefix "db/TSA" ftp://ftp.ncbi.nlm.nih.gov/genbank/tsa/I/*fsa_nt.gz 2>> {log}
        """

## DOWNLOADING THE HMM PROFILES
rule download_hmmprofiles:
    output:
        terpene_synth="hmm_profiles/PF01397.hmm", # Terpene_synth (PF01397)
        terpene_synth_C="hmm_profiles/PF03936.hmm", # Family: Terpene_synth_C (PF03936)
        terpene_synth_C2="hmm_profiles/PF19086.hmm" # Family: Terpene_syn_C_2 (PF19086) 
    log:
        "logs/download_hmmprofiles.log"
    params:
        memory="2"
    threads:
        1
    shell:
        """
        wget -O hmm_profiles/PF01397.hmm "http://pfam.xfam.org/family/PF01397/hmm" &> {log}
        wget -O hmm_profiles/PF03936.hmm "http://pfam.xfam.org/family/PF03936/hmm" &>> {log}
        wget -O hmm_profiles/PF19086.hmm "http://pfam.xfam.org/family/PF19086/hmm" &>> {log}
        """

## *****UNIPARC DATA*****

## SPLITTING THE UNIPARC DATABASE INTO SMALLER CHUNKS OF 10,000 SEQUENCES
# the db contains 374 059 950 sequences in 02/2021
# splitting by 10 000 should result in 37 406 files

# TODO: filter files with seqkit so they dont contain sequences longer than 100k 
# seqkit seq -M 100000
checkpoint split_uniparc:
    input:
        uniparc="db/uniparc/uniparc_active.fasta.gz"
    output:
        directory("db/uniparc_chunks"),
        filtered_uniparc="db/uniparc/uniparc_filtered.fasta"
    log:
        "logs/uniparc_split.log"
    params:
        memory="250" #TODO: change the value?
    threads:
        2
    shell:
        """
        seqkit seq -M 100000 {input.uniparc} -o {output.filtered_uniparc} 2> {log} # problems with hmm search occured when there were sequences longer than 100,000
        seqkit split {output.filtered_uniparc} -s 10000 --two-pass -O "db/uniparc_chunks" --threads {threads} 2>> {log}
        """
## SEARCHING SEQUENCES MATCHING THE HMM PROFILES
# running hmmsearch on the chunks - returns sequences IDs matching the provided hmm
rule hmmsearch_uniparc:
    input:
        terpene_synth="hmm_profiles/PF01397.hmm", # Terpene_synth (PF01397)
        terpene_synth_C="hmm_profiles/PF03936.hmm", # Family: Terpene_synth_C (PF03936)
        terpene_synth_C2="hmm_profiles/PF19086.hmm", # Family: Terpene_syn_C_2 (PF19086) 
        uniparc_chunk="db/uniparc_chunks/uniparc_filtered.part_{index}.fasta"
    output:
        hmm_tps="hmm_results/uniparc/hmm_tps_{index}.tsv",
        hmm_tpsc="hmm_results/uniparc/hmm_tpsc_{index}.tsv",
        hmm_tpsc2="hmm_results/uniparc/hmm_tpsc2_{index}.tsv"
    log:
        "logs/hmmsearch_{index}.log"
    params:
        memory="10" #TODO: change the value?
    threads:
        16
    shell:
        """
        hmmsearch --noali --tblout {output.hmm_tps} -E {config[e_value_threshold]} --cpu {threads} {input.terpene_synth} {input.uniparc_chunk} &> {log}
        hmmsearch --noali --tblout {output.hmm_tpsc} -E {config[e_value_threshold]} --cpu {threads} {input.terpene_synth_C} {input.uniparc_chunk} &>> {log}
        hmmsearch --noali --tblout {output.hmm_tpsc2} -E {config[e_value_threshold]} --cpu {threads} {input.terpene_synth_C2} {input.uniparc_chunk} &>> {log}
        """
def aggregate_uniparc_input(wildcards):
    checkpoints_output = checkpoints.split_uniparc.get(**wildcards).output[0] # there could be a problem with changing the dir structure - uniparc moved to db/uniparc
    indeces = glob_wildcards(os.path.join(checkpoints_output, "uniparc_filtered.part_{index}.fasta")).index
    completed = expand(os.path.join("hmm_results", "uniparc","hmm_tps_{index}.tsv"),index=indeces)
    completed2 = expand(os.path.join("hmm_results", "uniparc", "hmm_tpsc_{index}.tsv"), index=indeces)
    completed3 = expand(os.path.join("hmm_results", "uniparc", "hmm_tpsc2_{index}.tsv"), index=indeces)
    completed.extend(completed2)
    completed.extend(completed3)
    return completed    

# using cat {input} failed because input contained too many files so creating list with the files and then using xargs cat in rule combine_hmm
rule combine_uniparc_hmm_list:
    input:
        aggregate_uniparc_input
    output:
        "tmp.list"
    log:
        "logs/hmm_uniparc_list.log"
    params:
        memory="5"
    threads:
        1
    run:
        with open(output[0], "w") as output_handle:
             output_handle.write('\n'.join(input))  

## GETTING THE UPIs OF SEQUENCES MATCHING THE HMM PROFILES
rule combine_uniparc_hmm:
    input:
        "tmp.list"
    output:
        "hmm_combined.out"
    log:
        "logs/hmm_uniparc_combined.log"
    params:
        memory="5"
    threads:
        1
    shell:
        """
        cat {input} | xargs cat | grep -v '#' > {output} 2> {log}
        """

rule get_UPIs:
    input:
        upi_file="hmm_combined.out"
    output:
        upi_list="upi.txt"
    log:
        "logs/get_UPIs.out"
    params:
        memory="1"
    threads:
        1
    shell:
        """
        cut -d ' ' -f 1 {input.upi_file} > {output.upi_list} 2> {log}
        """
## SPLITIING THE UPIs FILE INTO SMALLER CHUNKS
checkpoint split_UPIs:
    input:
        upi_list="upi.txt"
    output:
        directory("UPIs")
    log:
        "logs/split_UPIs"
    params:
        memory="1"
    threads:
        1
    shell:
        """
        mkdir -p {output} 2> {log}
        split -d -l 10 {input} {output}/upi_ 2> {log}
        """

## GETTING THE UNIPROT SEQUENCES BASED ON THE UPIs        
rule get_uniprot_data:
    input:
        upi_file="UPIs/upi_{index}"
    output:
        uniprot_data="filtered_results/uniprot/seq_{index}"  
    log:
        "logs/get_uniprot_data_{index}.log" 
    params:
        memory="1"
    threads:
        1
    shell:
        """
        ./map.py {input} 2> {log}
        """

def aggregate_uniprot(wildcards):
    checkpoints_output = checkpoints.split_UPIs.get(**wildcards).output[0]
    indeces = glob_wildcards(os.path.join(checkpoints_output, "upi_{index}")).index
    completed = expand("filtered_results/uniprot/seq_{index}", index=indeces)
    return completed

## CREATING THE FINAL FILE WITH THE UNIPROT DATA  

rule uniprot_merge_list:
    input:
        aggregate_uniprot
    output:
        "uniprot.tmp.list"
    log:
        "logs/uniprot_list.log"
    params:
        memory="5"
    threads:
        1
    run:
        with open(output[0], "w") as output_handle:
             output_handle.write('\n'.join(input))  

rule uniprot_merge:
    input:
        "uniprot.tmp.list"
    output:
        tmp="results/uniprot_output_tmp.txt",
        tsv="results/uniprot_output.txt",
        short_fasta="results/uniprot_output.fasta"
    log:
        "logs/merge_uniprot_results.log"
    params:
        memory="1"
    threads:
        1
    shell:
        """
        mkdir -p results 2> {log}
        cat {input} | xargs cat > {output.tmp} 2>> {log}
        sed '/^<html>/,/^<\/html>/{{/^<html>/!{{/^<\/html>/!d}}}}' {output.tmp} | grep -v "<html>" | grep -v "</html>" | grep -v 'Entry' | sort | uniq > {output.tsv} 2>> {log}
        sed 's$^$>$' {output.tsv} | sed -E 's$\\t$\\n$2' > {output.short_fasta} 2>> {log}
        #cat {input} | xargs cat > {output.short_fasta} 2>> {log}
        #cat {input} | xargs cat | grep -v 'Entry' | sort | uniq > {output.tsv} 2>> {log}
        """

## *****1KP DATA*****

## SEARCHING SEQUENCES MATCHING THE HMM PROFILES
rule onekp_hmmsearch:
    input:
        terpene_synth="hmm_profiles/PF01397.hmm", # Terpene_synth (PF01397)
        terpene_synth_C="hmm_profiles/PF03936.hmm", # Family: Terpene_synth_C (PF03936)
        terpene_synth_C2="hmm_profiles/PF19086.hmm", # Family: Terpene_syn_C_2 (PF19086) 
        onekp_chunk="db/1kp/{index}-translated-protein.fa.gz"
    output:
        hmm_tps="hmm_results/1kp/hmm_tps_{index}.tsv",
        hmm_tpsc="hmm_results/1kp/hmm_tpsc_{index}.tsv",
        hmm_tpsc2="hmm_results/1kp/hmm_tpsc2_{index}.tsv",
        gunzipped_onekp_chunk="db/1kp/{index}-translated-protein.fa"
    log:
        "logs/1kp_hmmsearch_{index}.log"
    params:
        memory="10" #TODO: change the value?
    threads:
        16
    shell:
        """
        hmmsearch --noali --tblout {output.hmm_tps} -E {config[e_value_threshold]} --cpu {threads} {input.terpene_synth} {input.onekp_chunk} &> {log}
        hmmsearch --noali --tblout {output.hmm_tpsc} -E {config[e_value_threshold]} --cpu {threads} {input.terpene_synth_C} {input.onekp_chunk} &>> {log}
        hmmsearch --noali --tblout {output.hmm_tpsc2} -E {config[e_value_threshold]} --cpu {threads} {input.terpene_synth_C2} {input.onekp_chunk} &>> {log}
        gunzip -k {input.onekp_chunk}
        """
## GETTING THE SEQUENCES MATCHING THE HMM PROFILES
rule onekp_filter:
    input:
        tps_result="hmm_results/1kp/hmm_tps_{index}.tsv",
        tpsc_result="hmm_results/1kp/hmm_tpsc_{index}.tsv",
        tpsc2_result="hmm_results/1kp/hmm_tpsc2_{index}.tsv",
        seq_file="db/1kp/{index}-translated-protein.fa"
    output:
        "filtered_results/1kp/{index}_filtered.fa"
    log:
        "logs/1kp_filter_{index}.log"
    params:
        memory="1"
    threads:
        1
    run:
        # get the lines matching tps model
        tps_ids = []
        tpsc_ids = []
        tpsc2_ids = []
        with open(input[0], "r") as tps_input_handle:
            for line in tps_input_handle:
                if line.startswith('#') == False:
                    tps_ids.append(line.split()[0])
        with open(input[1], "r") as tpsc_input_handle:
            for line in tpsc_input_handle:
                if line.startswith('#') == False:
                    tpsc_ids.append(line.split()[0])
        with open(input[2], "r") as tpsc2_input_handle:
            for line in tpsc2_input_handle:
                if line.startswith('#') == False:
                    tpsc2_ids.append(line.split()[0])
        
        # find corresponding sequences
        with open(input[3], "r") as input_handle, open(output[0], "w") as output_handle:
            for record in Bio.SeqIO.parse(input_handle, "fasta"):
                if record.id in tps_ids or record.id in tpsc_ids or record.id in tpsc2_ids:
                    Bio.SeqIO.write(record, output_handle, "fasta")

## CREATING THE FINAL FILE WITH THE 1KP DATA  
def aggregate_1kp_input(wildcards):
    checkpoints_output = checkpoints.download1KP.get(**wildcards).output[0] # there could be a problem with changed dir structure - 1kp was moved to db/1kp
    indeces = glob_wildcards(os.path.join(checkpoints_output, "{index}-translated-protein.fa.gz")).index
    completed = expand(os.path.join("filtered_results", "1kp", "{index}_filtered.fa"),index=indeces)
    return completed

rule onekp_merge:
    input:
        aggregate_1kp_input
    output:
        tsv="results/1kp_output.txt", # file with the sequence on the same line
        short_fasta="results/1kp_output.fasta", # fasta file with sequence on one line
        long_fasta="results/1kp_output_long.fasta" # fasta file with sequence on multiple lines
    log:
        "logs/1kp_merge.log"
    params:
        memory="1"
    threads:
        1
    shell:
        """
        mkdir -p results
        cat {input} > {output.long_fasta} 2> {log}
        seqkit seq -w 0 {output.long_fasta} -o {output.short_fasta}
        awk '{{printf "%s%s",$0,NR%2?"\t":RS}}' {output.short_fasta} > {output.tsv} 2>> {log}
        """

## *****PHYTOZOME DATA*****

# Phytozome transcriptome sequences were manually downloaded from the website and are stored in the db/Phytozome directory
## GUNZIPPING THE FILES 
rule phytozome_gunzip_chunks:
    input:
        phytozome_chunk="db/Phytozome/{index}.protein.fa.gz"  
    output:
        gunzipped_chunk="db/Phytozome/{index}.protein.fa"     
    log:
        "logs/phytozome_gunzip_{index}.log"
    params:
        memory="5"
    shell:
        """
        gunzip -k {input.phytozome_chunk}
        """
## ADDING INFORMATION ABOUT THE SPECIE 
## the downloaded transcriptomes do not contain information about the specie in the header so I use the filename to add this information
rule phytozome_add_specie:
    input:
        phytozome_chunk="db/Phytozome/{index}.protein.fa"
    output:
        annotated_chunk="db/Phytozome/{index}.protein2.fa"
    log:
        "logs/phytozome_add_specie_{index}.log"
    params:
        memory="10"
    run:
        with open(input[0], "r") as input_handle, open(output[0], "w") as output_handle:
            for record in Bio.SeqIO.parse(input_handle, "fasta"):
                record.description += "\t" + wildcards.index
                Bio.SeqIO.write(record, output_handle, 'fasta')

## SEARCHING SEQUENCES MATCHING THE HMM PROFILES
checkpoint phytozome_hmmsearch:
    input:
        terpene_synth="hmm_profiles/PF01397.hmm", # Terpene_synth (PF01397)
        terpene_synth_C="hmm_profiles/PF03936.hmm", # Family: Terpene_synth_C (PF03936)
        terpene_synth_C2="hmm_profiles/PF19086.hmm", # Family: Terpene_syn_C_2 (PF19086) 
        phytozome_chunk="db/Phytozome/{index}.protein2.fa"
    output:
        hmm_tps="hmm_results/phytozome/hmm_tps_{index}.tsv",
        hmm_tpsc="hmm_results/phytozome/hmm_tpsc_{index}.tsv",
        hmm_tpsc2="hmm_results/phytozome/hmm_tpsc2_{index}.tsv"#,
    log:
        "logs/phytozome_hmmsearch_{index}.log"
    params:
        memory="10"
    shell:
        """
        hmmsearch --noali --tblout {output.hmm_tps} -E {config[e_value_threshold]} --cpu {threads} {input.terpene_synth} {input.phytozome_chunk} &> {log}
        hmmsearch --noali --tblout {output.hmm_tpsc} -E {config[e_value_threshold]} --cpu {threads} {input.terpene_synth_C} {input.phytozome_chunk} &>> {log}
        hmmsearch --noali --tblout {output.hmm_tpsc2} -E {config[e_value_threshold]} --cpu {threads} {input.terpene_synth_C2} {input.phytozome_chunk} &>> {log}
        """

## GETTING THE SEQUENCES MATCHING THE HMM PROFILES
rule phytozome_filter:
    input:
        tps_result="hmm_results/phytozome/hmm_tps_{index}.tsv",
        tpsc_result="hmm_results/phytozome/hmm_tpsc_{index}.tsv",
        tpsc2_result="hmm_results/phytozome/hmm_tpsc2_{index}.tsv",
        seq_file="db/Phytozome/{index}.protein2.fa"
    output:
        "filtered_results/phytozome/{index}_filtered.fa"
    log:
        "logs/phytozome_filter_{index}.log"
    params:
        memory="1"
    threads:
        1
    run:
        # get the lines matching tps model
        tps_ids = []
        tpsc_ids = []
        tpsc2_ids = []
        with open(input[0], "r") as tps_input_handle:
            for line in tps_input_handle:
                if line.startswith('#') == False:
                    tps_ids.append(line.split()[0])
        with open(input[1], "r") as tpsc_input_handle:
            for line in tpsc_input_handle:
                if line.startswith('#') == False:
                    tpsc_ids.append(line.split()[0])
        with open(input[2], "r") as tpsc2_input_handle:
            for line in tpsc2_input_handle:
                if line.startswith('#') == False:
                    tpsc2_ids.append(line.split()[0])

        # find corresponding sequences
        with open(input[3], "r") as input_handle, open(output[0], "w") as output_handle:
            for record in Bio.SeqIO.parse(input_handle, "fasta"):
                if record.id in tps_ids or record.id in tpsc_ids or record.id in tpsc2_ids:
                    Bio.SeqIO.write(record, output_handle, "fasta")

## CREATING THE FINAL FILE WITH THE PHYTOZOME DATA  
def aggregate_phytozome_input(wildcards):
    indeces = glob_wildcards("db/Phytozome/{index}.protein.fa.gz").index
    completed = expand(os.path.join("filtered_results", "phytozome","{index}_filtered.fa"),index=indeces)
    return completed

rule phytozome_merge:
    input:
        aggregate_phytozome_input
    output:
        tsv="results/phytozome_output.txt",
        short_fasta="results/phytozome_output.fasta",
        long_fasta="results/phytozome_output_long.fasta"
    log:
        "logs/phytozome_merge.log"
    params:
        memory="1"
    threads:
        1
    shell:
        """
        mkdir -p results
        cat {input} > {output.long_fasta} 2> {log}
        seqkit seq -w 0 {output.long_fasta} -o {output.short_fasta}
        awk '{{printf "%s%s",$0,NR%2?"\t":RS}}' {output.short_fasta} > {output.tsv} 2>> {log}
        """

## *****TSA DATA*****

## PREDICTING PROTEIN SEQUENCES WITH TRANSDECODER 
## downloaded files from TSA are nucleotide sequences
rule transdecoder_TSA:
    input:
        transcripts_file="db/TSA/tsa.{ID}.{version}.fsa_nt.gz"
    output:
        nt_file="db/TSA/tsa.{ID}.{version}.fsa_nt",
        pep_file="db/TSA/{ID}.{version}.pep"
    log:
        "logs/transdecoder_{ID}.{version}.log"
    params:
        memory="16"
    threads:
        2
    shell:
        """
        TransDecoder.LongOrfs -t {input} &> {log}
        TransDecoder.Predict --no_refine_starts -t {input} &>> {log} # some files failed when running them without --no_refine_starts, also some files failed because they contain too few sequences so no protein sequence is predicted
        mv tsa.{wildcards.ID}.{wildcards.version}.fsa_nt {output.nt_file} 2>> {log}
        mv tsa.{wildcards.ID}.{wildcards.version}.fsa_nt.transdecoder.pep {output.pep_file} 2>> {log}

        # TransDecoder creates many files in the workinf dir however removing the files during the run of this rule can lead to some errors
        rm -r tsa.{wildcards.ID}.{wildcards.version}.fsa_nt.transdecoder*
        #rm pipeliner*
        """
## ADDING INFORMATION ABOUT THE SPECIE 
## the downloaded files do not contain information about the specie in the header so I use mapping table downloaded from the website which I adjusted for easier searching
rule TSA_add_specie:
    input:
        seq_file="db/TSA/{ID}.{version}.pep",
        map_file="tsa_mapping3.tsv"
    output:
        "db/TSA/{ID}.{version}_2.pep" 
    log:
        "logs/TSA_add_specie_{ID}.{version}.log"
    params:
        memory="2"
    threads:
        2
    run:
        with open(input[0], "r") as input_handle, open(output[0], "w") as output_handle, open(input[1], "r") as mapping_handle:
            ID = wildcards.ID
            if ID == 'GGYI': # this ID is not on the website where mapping could be found for some reason so there is just tab
                specie = ''
            else:
                specie = [" ".join(line.split()[1:]) for line in mapping_handle if line.split()[0] == ID][0]
            for record in Bio.SeqIO.parse(input_handle, "fasta"):
                record.description += "\t" + specie
                Bio.SeqIO.write(record, output_handle, 'fasta')  

## SEARCHING SEQUENCES MATCHING THE HMM PROFILES                
rule TSA_hmmsearch:
    input:
        terpene_synth="hmm_profiles/PF01397.hmm", # Terpene_synth (PF01397)
        terpene_synth_C="hmm_profiles/PF03936.hmm", # Family: Terpene_synth_C (PF03936)
        terpene_synth_C2="hmm_profiles/PF19086.hmm",
        tsa_chunk="db/TSA/{ID}.{version}_2.pep" 
    output:
        hmm_tps="hmm_results/tsa/hmm_tps_{ID}.{version}.tsv",
        hmm_tpsc="hmm_results/tsa/hmm_tpsc_{ID}.{version}.tsv",
        hmm_tpsc2="hmm_results/tsa/hmm_tpsc2_{ID}.{version}.tsv"
    log:
        "logs/tsa_hmmsearch_{ID}.{version}.log"
    params:
        memory="10"
    threads:
        16
    shell:
        """
        hmmsearch --noali --tblout {output.hmm_tps} -E {config[e_value_threshold]} --cpu {threads} {input.terpene_synth} {input.tsa_chunk} &> {log}
        hmmsearch --noali --tblout {output.hmm_tpsc} -E {config[e_value_threshold]} --cpu {threads} {input.terpene_synth_C} {input.tsa_chunk} &>> {log}
        hmmsearch --noali --tblout {output.hmm_tpsc2} -E {config[e_value_threshold]} --cpu {threads} {input.terpene_synth_C2} {input.tsa_chunk} &>> {log}
        """     

rule TSA_filter:
    input:
        tps_result="hmm_results/tsa/hmm_tps_{ID}.{version}.tsv",
        tpsc_result="hmm_results/tsa/hmm_tpsc_{ID}.{version}.tsv",
        tpsc2_result="hmm_results/tsa/hmm_tpsc2_{ID}.{version}.tsv",
        seq_file="db/TSA/{ID}.{version}_2.pep" 
    output:
        "filtered_results/tsa/{ID}.{version}_filtered.fa"
    log:
        "logs/tsa_filter_{ID}.{version}.log"
    params:
        memory="1"
    threads:
        1
    run:
        # get the lines matching tps model
        tps_ids = []
        tpsc_ids = []
        tpsc2_ids = []
        with open(input[0], "r") as tps_input_handle:
            for line in tps_input_handle:
                if line.startswith('#') == False:
                    tps_ids.append(line.split()[0])
        with open(input[1], "r") as tpsc_input_handle:
            for line in tpsc_input_handle:
                if line.startswith('#') == False:
                    tpsc_ids.append(line.split()[0])
        with open(input[2], "r") as tpsc2_input_handle:
            for line in tpsc2_input_handle:
                if line.startswith('#') == False:
                    tpsc2_ids.append(line.split()[0])

        # find corresponding sequences
        with open(input[3], "r") as input_handle, open(output[0], "w") as output_handle:
            for record in Bio.SeqIO.parse(input_handle, "fasta"):
                if record.id in tps_ids or record.id in tpsc_ids or record.id in tpsc2_ids:
                    Bio.SeqIO.write(record, output_handle, "fasta")



def aggregate_TSA_input(wildcards):
    ids = glob_wildcards("db/TSA/tsa.{ID}.{version}.fsa_nt.gz").ID #"TSA/tsa.{ID}.{version}.fsa_nt.gz" 
    versions=glob_wildcards("db/TSA/tsa.{ID}.{version}.fsa_nt.gz").version
    completed = expand(os.path.join("filtered_results", "tsa", "{ID}.{version}_filtered.fa"),zip,ID=ids, version=versions)
    return completed

rule make_TSA_list:
    input:
        aggregate_TSA_input
    output:
        "tsa.tmp.list"
    log:
        "logs/tsa_results_list.log"
    params:
        memory="5"
    threads:
        1
    run:
        with open(output[0], "w") as output_handle:
             output_handle.write('\n'.join(input))  

rule TSA_merge:
    input:
        "tsa.tmp.list"
    output:
        tsv="results/TSA_output.txt",
        short_fasta="results/TSA_output.fasta",
        long_fasta="results/TSA_output_long.fasta"
    log:
        "logs/tsa_combined.log"
    params:
        memory="5"
    threads:
        1
    shell:
        """
        mkdir -p results 2> {log}
        cat {input} | xargs cat > {output.long_fasta} 2>> {log}
        seqkit seq -w 0 {output.long_fasta} -o {output.short_fasta} 2>> {log}
        awk '{{printf "%s%s",$0,NR%2?"\t":RS}}' {output.short_fasta} > {output.tsv} 2>> {log}
        """

## MERGE ALL RESULTS INTO A SINGLE FILE
rule merge_results:
    input:
        onekp="results/1kp_output.txt",
        phytozome="results/phytozome_output.txt",
        uniprot="results/uniprot_output.txt",
        tsa="results/TSA_output.txt"
    output:
        onekp="results/1kp_output_f.tsv",
        phytozome="results/phytozome_output_f.tsv",
        uniprot="results/uniprot_output_f.tsv",
        tsa="results/TSA_output_f.tsv",
        fin="results/final.tsv"
    log:
        "logs/merge_results.log"
    params:
        memory="1"
    threads:
        1
    shell:
        """
        sed -r "s#(>scaffold-([A-Z]{{4}}-[0-9]{{7}})-([A-Za-z_]+).+)\t(.+)#\\2\t\\1\t\\3\t1kp\t\\4#" {input.onekp} > {output.onekp} 2> {log}
        sed -r 's#>(.+)\t(.+)\t(.+)#\\1\t\\1\t\\2\tphytozome\t\\3#' {input.phytozome} > {output.phytozome} 2>> {log}
        sed -r 's#>((.+) .+ .+ .+ .+ .+ .+)\t(.+)\t(.+)#\\2\t\\1\t\\3\tTSA\t\\4#' {input.tsa} > {output.tsa} 2>> {log}
        sed -r 's#(.+)\t(.+)\t(.+)#\\1\t\\1\t\\2\tuniprot\t\\3#' {input.uniprot} > {output.uniprot} 2>> {log}
        cat {output.onekp} {output.phytozome} {output.uniprot} {output.tsa} > {output.fin} 2>> {log}
        """
