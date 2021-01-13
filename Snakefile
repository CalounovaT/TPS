

configfile: "config.yaml"

rule localrules: all, clean

rule all:
  input:
        "db/uniparc_active.fasta",
        "hmm_profiles/PF01397.hmm",
        "hmm_profiles/PF03936.hmm",
        "hmm_results"

rule clean:
    shell:
        """
        rm -r hmm_results log
        """
    
rule download_uniparc:
    output:
        "db/uniparc_active.fasta"
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
        wget --directory-prefix hmm_profiles -O PF01397.hmm "http://pfam.xfam.org/family/PF01397/hmm" &> {log}
        wget --directory-prefix hmm_profiles -O PF03936.hmm "http://pfam.xfam.org/family/PF03936/hmm" &>> {log}
        """

# split the database into smaller chunks so it an be searched with hmmsearch - max is 100K?
rule split_uniparc:
    input:
        uniparc="db/uniparc_active.fasta"
    output:
        directory("uniparc_chunks")
    log:
        "logs/uniparc_split.log"
    params:
        memory="2" #TODO: change the value?
    threads:
        2
    shell:
        """
        seqkit split {input.uniparc} -s 10000 --two-pass -O "uniparc_chunks" --threads 2 2> {log}
        """

# run hmmsearch on the chunks - returns sequences IDs matcching the provided hmm
rule hmmsearch:
    input:
        terpene_synth="hmm_profiles/PF01397.hmm", # Terpene_synth (PF01397)
        terpene_synth_C="hmm_profiles/PF03936.hmm", # Family: Terpene_synth_C (PF03936)
        uniparc_chunk="uniparc_chunks/{chunk}"
    output:
        "hmm_results/hmm_{chunk}"
    log:
        "logs/hmmsearch_{chunk}"
    params:
        memory="10" #TODO: change the value?
    threads:
        4
    shell:
        """
        hmmsearch --noali --tblout {output} -E {config[e_value_threshold]} {input.terpene_synth} {input.uniparc_chunk}
        """

# create desired output: sequence ID; AA sequence; species name; TPS class
#rule:
#input:
#output:
#log:
#params:

