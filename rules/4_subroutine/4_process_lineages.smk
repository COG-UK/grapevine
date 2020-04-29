configfile: workflow.current_basedir + "/config.yaml"

import os
import pandas as pd

##### Configuration #####

config["publish_path"] = os.path.abspath(config["publish_path"])

LINEAGES = config["lineages"].split()
OUTGROUPS = config["lineage_specific_outgroups"].split()

lineage_to_outgroup_map = {}
for i,lin in enumerate(LINEAGES):
    lineage_to_outgroup_map[lin] = OUTGROUPS[i].replace("/","_")

print("lineages", LINEAGES)
print("outgroups", OUTGROUPS)

##### Target rules #####

rule all:
    input:
        config["output_path"] + "/4/all_traits.csv",
        config["output_path"] + "/4/cog_gisaid_full.tree"

rule iq_tree:
    input:
        lineage_fasta = config["output_path"] + "/4/lineage_{lineage}.fasta"
    params:
        lineage = "{lineage}",
        outgroup = lambda wildcards: lineage_to_outgroup_map[wildcards.lineage]
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.tree"
    log:
        config["output_path"] + "/logs/4_iq_tree_{lineage}.log"
    resources: mem_per_cpu=10000
    shell:
        """
        echo "{params.outgroup} {input.lineage_fasta} {params.lineage}"
        iqtree -m HKY -bb 1000 -czb \
        -o \"{params.outgroup}\" \
        -cptime 300 \
        -nt AUTO \
        -s {input.lineage_fasta} &> {log}

        RESULT=$?
        if [ $RESULT -eq 0 ]
        then
          datafunk repair_names \
            --fasta {input.lineage_fasta} \
            --tree {input.lineage_fasta}.treefile \
            --out {output.tree} &>> {log}
        fi
        """

rule phylotype_tree:
    input:
        tree = rules.iq_tree.output.tree
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.phylotyped.tree"
    params:
        lineage="{lineage}",
        collapse=5E-6,
        threshold=2E-5,
    log:
        config["output_path"] + "/logs/4_phylotype_{lineage}.log"
    shell:
        """
        clusterfunk phylotype \
        --format newick \
        --collapse_to_polytomies {params.collapse} \
        --threshold {params.threshold} \
        --prefix {params.lineage}_1 \
        --input {input.tree} \
        --output {output.tree} &> {log}
        """


rule annotate_tree:
    input:
        tree = rules.phylotype_tree.output.tree,
        metadata = config["metadata"]
    params:
        lineage = "{lineage}",
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.phylotyped.annotated.tree"
    log:
        config["output_path"] + "/logs/4_annotate_{lineage}.log"
    shell:
        """
        clusterfunk annotate_tips \
          --in-metadata {input.metadata} \
          --trait-columns lineage \
          country uk_lineage \
          --index-column sequence_name \
          --boolean-for-trait country='UK' country='UK' country=UK \
          --boolean-trait-names country_uk country_uk_acctran country_uk_deltran\
          --input {input.tree} \
          --output {output.tree} &> {log}
        """

rule acctran_ancestral_reconstruction:
    input:
        tree = rules.annotate_tree.output.tree
    params:
        lineage = "{lineage}",
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.phylotyped.annotated.acc.tree"
    log:
        config["output_path"] + "/logs/4_acctran_ancestral_reconstruction_{lineage}.log"
    shell:
        """
        clusterfunk ancestral_reconstruction \
        --traits country_uk_acctran \
        --acctran \
        --ancestral_state False \
        --input {input.tree} \
        --output {output.tree} &> {log}
        """

rule deltran_ancestral_reconstruction:
    input:
        tree = rules.acctran_ancestral_reconstruction.output.tree
    params:
        lineage = "{lineage}",
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.phylotyped.annotated.acc.del.tree"
    log:
        config["output_path"] + "/logs/4_deltran_ancestral_reconstruction_{lineage}.log"
    shell:
        """
        clusterfunk ancestral_reconstruction \
        --traits country_uk_deltran \
        --deltran \
        --ancestral_state False \
        --input {input.tree} \
        --output {output.tree} &> {log}
        """

rule push_lineage_to_tips:
    input:
        tree = rules.deltran_ancestral_reconstruction.output.tree
    params:
        lineage = "{lineage}",
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.phylotyped.annotated.acc.del.uk_lineages.tree"
    log:
        config["output_path"] + "/logs/4_push_lineage_to_tips_{lineage}.log"
    shell:
        """
        clusterfunk push_annotations_to_tips \
          --traits uk_lineage \
          --stop-where-trait country_uk_acctran=False \
          --input {input.tree} \
          --output {output.tree} &> {log}
        """

rule label_acctran_introductions:
    input:
        tree = rules.push_lineage_to_tips.output.tree
    params:
        lineage = "{lineage}",
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.phylotyped.annotated.acc.del.uk_lineages.acc_labelled.tree"
    log:
        config["output_path"] + "/logs/4_label_acctran_introductions_{lineage}.log"
    shell:
        """
        clusterfunk label_transitions \
          --trait country_uk_acctran \
          --to True \
          --transition-name acc_lineage \
          --transition-prefix {params.lineage}_ \
          --include_root \
          --input {input.tree} \
          --output {output.tree} &> {log}
        """

rule label_deltran_introductions:
    input:
        tree = rules.label_acctran_introductions.output.tree
    params:
        lineage = "{lineage}",
        outdir = config["publish_path"] + "/COG_GISAID",
        prefix = config["publish_path"] + "/COG_GISAID/cog_gisaid_"
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.phylotyped.annotated.acc.del.uk_lineages.acc_labelled.del_labelled.tree"
    log:
        config["output_path"] + "/logs/4_label_deltran_introductions_{lineage}.log"
    shell:
        """
        clusterfunk label_transitions \
          --trait country_uk_deltran \
          --to True \
          --transition-name del_lineage \
          --transition-prefix {params.lineage}_ \
          --input {input.tree} \
          --output {output.tree} &> {log}

        mkdir -p {params.outdir}
        cp {output.tree} {params.prefix}_lineage_{params.lineage}.tree
        """

rule graft:
    input:
         # not sure how to pass this as a space separated list below. Also assuming the order here matches lineages
        trees = expand(config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.phylotyped.annotated.acc.del.uk_lineages.acc_labelled.del_labelled.tree", lineage=LINEAGES),
        guide_tree="/path/to/guide.tree" # path to guide tree here
    params:
        outdir = config["publish_path"] + "/COG_GISAID",
        prefix = config["publish_path"] + "/COG_GISAID/"
    output:
        tree = config["output_path"] + "/4/cog_gisaid_full.tree"
    log:
        config["output_path"] + "/logs/4_graft.log"
    shell:
        """
        clusterfunk graft \
        --full-graft \
        --scions {input.trees} \
        --scion_annotation_name scion_lineage \
        --annotate_scions {LINEAGES} \
        --input {input.guide_tree} \
        --output {output.tree} &> {log}
        
        mkdir -p {params.outdir}
        cp {output.tree} {params.prefix}
        """

rule output_annotations:
    input:
        tree = rules.label_deltran_introductions.output.tree,
    params:
        lineage = "{lineage}",
    output:
        traits = config["output_path"] + "/4/{lineage}/traits.csv"
    log:
        config["output_path"] + "/logs/4_output_annotations_{lineage}.log"
    shell:
        """
        clusterfunk extract_tip_annotations \
          --traits country lineage uk_lineage acc_lineage del_lineage phylotype \
          --input {input.tree} \
          --output {output.traits} &> {log}
        """

rule combine_traits_files:
    input:
        list_df = expand(config["output_path"] + "/4/{lineage}/traits.csv", lineage=LINEAGES)
    output:
        traits = config["output_path"] + "/4/all_traits.csv"
    log:
        config["output_path"] + "/logs/4_combine_traits_files.log"
    run:
        """
        result = pd.concat({input.list_df})
        result.to_csv({output}, index=False)
        """

