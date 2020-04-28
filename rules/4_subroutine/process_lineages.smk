configfile: workflow.current_basedir + "/config.yaml"

import datetime
import os

date = datetime.date.today()

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
        expand(config["output_path"] + "/4/{lineage}/traits.csv", lineage=LINEAGES),
        expand(config["output_path"] + "/4/{lineage}/cut_out_trees_done", lineage=LINEAGES)

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

rule annotate_tree:
    input:
        tree = rules.iq_tree.output.tree,
        metadata = config["metadata"]
    params:
        lineage = "{lineage}",
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.annotated.tree"
    log:
        config["output_path"] + "/logs/4_annotate_{lineage}.log"
    shell:
        """
        clusterfunk annotate_tips \
          --in-metadata {input.metadata} \
          --trait-columns lineage \
          country uk_lineage \
          --index-column sequence_name \
          --boolean-for-trait country='UK' \
          --boolean-trait-names country_uk \
          --input {input.tree} \
          --format newick \
          --output {output.tree} &> {log}
        """

rule ancestral_reconstruction:
    input:
        tree = rules.annotate_tree.output.tree
    params:
        lineage = "{lineage}",
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.annotated.acc.tree"
    log:
        config["output_path"] + "/logs/4_ancestral_reconstruction_{lineage}.log"
    shell:
        """
        clusterfunk ancestral_reconstruction \
        --traits country_uk \
        --acctran \
        --ancestral_state False \
        --input {input.tree} \
        --output {output.tree} &> {log}
        """

rule push_lineage_to_tips:
    input:
        tree = rules.ancestral_reconstruction.output.tree
    params:
        lineage = "{lineage}",
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.annotated.acc.uk_lineages.tree"
    log:
        config["output_path"] + "/logs/4_push_lineage_to_tips_{lineage}.log"
    shell:
        """
        clusterfunk push_annotations_to_tips \
          --traits uk_lineage \
          --stop-where-trait country_uk=False \
          --input {input.tree} \
          --output {output.tree} &> {log}
        """

rule label_introductions:
    input:
        tree = rules.push_lineage_to_tips.output.tree
    params:
        lineage = "{lineage}",
        outdir = config["publish_path"] + "/COG_GISAID",
        prefix = config["publish_path"] + "/COG_GISAID/cog_gisaid_%s" %date
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.annotated.acc.uk_lineages.labelled.tree"
    log:
        config["output_path"] + "/logs/4_label_introductions_{lineage}.log"
    shell:
        """
        clusterfunk label_transitions \
          --trait country_uk \
          --to True \
          --transition-name acc_lineage \
          --transition-prefix {params.lineage}_ \
          --input {input.tree} \
          --output {output.tree} &> {log}

        mkdir -p {params.outdir}
        cp {output.tree} {params.prefix}_lineage_{params.lineage}.tree
        """

rule cut_out_trees:
    input:
        tree = rules.label_introductions.output.tree
    params:
        lineage = "{lineage}",
        outdir = config["output_path"] + "/4/{lineage}/trees",
        pubdir = config["publish_path"] + "/COG_GISAID/acc_lineages/{lineage}",
    output:
        config["output_path"] + "/4/{lineage}/cut_out_trees_done"
    log:
        config["output_path"] + "/logs/4_cut_out_trees_{lineage}.log"
    threads: 8
    shell:
        """
        clusterfunk prune \
          --extract \
          --trait acc_lineage \
          --input {input.tree} \
          --threads {threads} \
          --output {params.outdir} &> {log}

        if [ ! -z "$(ls -A {params.outdir})" ]
        then
          mkdir -p {params.pubdir}
          cp {params.outdir}/*.tree {params.pubdir}/
        fi

        touch {output}
        """

rule output_annotations:
    input:
        tree = rules.label_introductions.output.tree,
    params:
        lineage = "{lineage}",
    output:
        traits = config["output_path"] + "/4/{lineage}/traits.csv"
    log:
        config["output_path"] + "/logs/4_output_annotations_{lineage}.log"
    shell:
        """
        clusterfunk extract_tip_annotations \
          --traits country lineage uk_lineage acc_lineage \
          --input {input.tree} \
          --output {output.traits} &> {log}
        """

