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
        "db/uniparc_active.fasta.gz",
        "hmm_profiles/PF01397.hmm",
        "hmm_profiles/PF03936.hmm",
        "uniparc_chunks",
        "hmm_combined.out",
        "upi.txt",
        "output.txt",
        "1kp",
        "1kp_output.txt"
        #"1kp_hmm_combined.out"

rule clean:
    shell:
        """
        rm -rf hmm_results logs hmm_combined.out output.txt uniparc_chunks uniprot UPIs upi.txt tmp.list 
        """
    
rule download_uniparc:
    output:
        "db/uniparc_active.fasta.gz"
    log:
        "logs/download_uniparc.log"
    params:
        memory="2"
    threads:
        1
    shell:
        """
        wget --directory-prefix db "https://ftp.expasy.org/databases/uniprot/current_release/uniparc/uniparc_active.fasta.gz" &> {log}
        """

checkpoint download1KP:
    output:
        directory("1kp")
    log:
        "logs/download_1kp"
    params:
        memory="2"
    threads:
        1
    shell:
        """
        wget --directory-prefix 1kp -r -l2 --no-parent -nd -A "translated-protein.fa.gz" ftp://parrot.genomics.cn/gigadb/pub/10.5524/100001_101000/100627/assemblies/
        """

rule download_hmmprofiles:
    output:
        terpene_synth="hmm_profiles/PF01397.hmm", # Terpene_synth (PF01397)
        terpene_synth_C="hmm_profiles/PF03936.hmm" # Family: Terpene_synth_C (PF03936)
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
        """

# split the database into smaller chunks so it an be searched with hmmsearch - max is 100K?
# the db contains 374 059 950 sequences in 02/2021
# splitting by 10 000 should result in 37 406 files
# splitting by 50 000 should result in 7 482 files

# TODO: filter files with seqkit so they dont contain sequences longer than 100k 
# seqkit seq -M 100000
checkpoint split_uniparc:
    input:
        uniparc="db/uniparc_active.fasta.gz"
    output:
        directory("uniparc_chunks")
    log:
        "logs/uniparc_split.log"
    params:
        memory="15" #TODO: change the value?
    threads:
        2
    shell:
        """
        seqkit split {input.uniparc} -s 10000 --two-pass -O "uniparc_chunks" --threads 2 &> {log}
        """

# run hmmsearch on the chunks - returns sequences IDs matching the provided hmm
rule hmmsearch:
    input:
        terpene_synth="hmm_profiles/PF01397.hmm", # Terpene_synth (PF01397)
        terpene_synth_C="hmm_profiles/PF03936.hmm", # Family: Terpene_synth_C (PF03936)
        uniparc_chunk="uniparc_chunks/uniparc_active.part_{index}.fasta"
    output:
        hmm_tps="hmm_results/hmm_tps_{index}.tsv",
        hmm_tpsc="hmm_results/hmm_tpsc_{index}.tsv"
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
        """

rule onekp_hmmsearch:
    input:
        terpene_synth="hmm_profiles/PF01397.hmm", # Terpene_synth (PF01397)
        terpene_synth_C="hmm_profiles/PF03936.hmm", # Family: Terpene_synth_C (PF03936)
        onekp_chunk="1kp/{index}-translated-protein.fa.gz"
    output:
        hmm_tps="1kp_hmm_results/hmm_tps_{index}.tsv",
        hmm_tpsc="1kp_hmm_results/hmm_tpsc_{index}.tsv",
        gunzipped_onekp_chunk="1kp/{index}-translated-protein.fa"
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
        gunzip -k {input.onekp_chunk}
        """

def aggregate_input(wildcards):
    checkpoints_output = checkpoints.split_uniparc.get(**wildcards).output[0]
    indeces = glob_wildcards(os.path.join(checkpoints_output, "uniparc_active.part_{index}.fasta")).index
    completed = expand(os.path.join("hmm_results", "hmm_tps_{index}.tsv"),index=indeces)
    completed2 = expand(os.path.join("hmm_results", "hmm_tpsc_{index}.tsv"), index=indeces)
    completed.extend(completed2)
    return completed    

# using cat {input} failed because input contained too many files so creating list with the files and then using xargs cat in rule combine_hmm
rule combine_hmm_list:
    input:
        aggregate_input
    output:
        "tmp.list"
    log:
        "logs/hmm_list.log"
    params:
        memory="5"
    threads:
        1
    run:
        with open(output[0], "w") as output_handle:
             output_handle.write('\n'.join(input))  

rule combine_hmm:
    input:
        "tmp.list"
    output:
        "hmm_combined.out"
    log:
        "logs/hmm_combined.log"
    params:
        memory="5"
    threads:
        1
    shell:
        """
        cat {input} | xargs cat | grep -v '#' > {output} 2> {log}
        """

#def aggregate_1kp_input(wildcards):
#    checkpoints_output = checkpoints.download1KP.get(**wildcards).output[0]
#    indeces = glob_wildcards(os.path.join(checkpoints_output, "{index}-translated-protein.fa.gz")).index
#    completed = expand(os.path.join("1kp_hmm_results", "hmm_tps_{index}.tsv"),index=indeces)
#    completed2 = expand(os.path.join("1kp_hmm_results", "hmm_tpsc_{index}.tsv"), index=indeces)
#    completed.extend(completed2)
#    return completed
#
## using cat {input} failed because input contained too many files so creating list with the files and then using xargs cat in rule combine_hmm
#rule onekp_combine_hmm_list:
#    input:
#        aggregate_1kp_input
#    output:
#        "1kp_tmp.list"
#    log:
#        "logs/1kp_hmm_list.log"
#    params:
#        memory="5"
#    threads:
#        1
#    run:
#        with open(output[0], "w") as output_handle:
#             output_handle.write('\n'.join(input))
#
#rule onekp_combine_hmm:
#    input:
#        "1kp_tmp.list"
#    output:
#        "1kp_hmm_combined.out"
#    log:
#        "logs/1kp_hmm_combined.log"
#    params:
#        memory="5"
#    threads:
#        1
#    shell:
#        """
#        cat {input} | xargs cat | grep -v '#' > {output} 2> {log}
#        """

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
        
rule get_uniprot_data:
    input:
        upi_file="UPIs/upi_{index}"
    output:
        uniprot_data="uniprot/seq_{index}"  
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
    completed = expand("uniprot/seq_{index}", index=indeces)
    return completed

rule merge_results:
    input:
        aggregate_uniprot
    output:
        "output.txt"
    log:
        "logs/merge_results.log"
    params:
        memory="1"
    threads:
        1
    shell:
        """
        cat {input} | grep -v 'Entry' | sort | uniq > {output} 2> {log}
        """  

rule onekp_filter:
    input:
        tps_result="1kp_hmm_results/hmm_tps_{index}.tsv",
        tpsc_result="1kp_hmm_results/hmm_tpsc_{index}.tsv",
        seq_file="1kp/{index}-translated-protein.fa"
    output:
        "1kp_filtered_results/{index}_filtered.fa"
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
        with open(input[0], "r") as tps_input_handle:
            for line in tps_input_handle:
                if line.startswith('#') == False:
                    tps_ids.append(line.split()[0])
        with open(input[1], "r") as tpsc_input_handle:
            for line in tpsc_input_handle:
                if line.startswith('#') == False:
                    tpsc_ids.append(line.split()[0])
        
        # find corresponding sequences
        with open(input[2], "r") as input_handle, open(output[0], "w") as output_handle:
            for record in Bio.SeqIO.parse(input_handle, "fasta"):
                if record.id in tps_ids or record.id in tpsc_ids:
                    Bio.SeqIO.write(record, output_handle, "fasta")


def aggregate_1kp_input(wildcards):
    checkpoints_output = checkpoints.download1KP.get(**wildcards).output[0]
    indeces = glob_wildcards(os.path.join(checkpoints_output, "{index}-translated-protein.fa.gz")).index
    completed = expand(os.path.join("1kp_filtered_results", "{index}_filtered.fa"),index=indeces)
    #completed2 = expand(os.path.join("1kp_hmm_results", "hmm_tpsc_{index}.tsv"), index=indeces)
    #completed.extend(completed2)
    return completed

rule onekp_merge:
    input:
        aggregate_1kp_input
    output:
        "1kp_output.txt"
    log:
        "logs/1kp_merge.log"
    params:
        memory="1"
    threads:
        1
    shell:
        """
        cat {input} > 1kp_output_tmp.txt 2> {log}
        seqkit seq -w 0 1kp_output_tmp.txt -o 1kp_output_tmp2.txt
        awk '{{printf "%s%s",$0,NR%2?"\t":RS}}' 1kp_output_tmp2.txt > {output} 2>> {log}
        rm -f 1kp_output_tmp.txt
        rm -f 1kp_output_tmp2.txt
        """
