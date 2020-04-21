configfile: workflow.current_basedir + "/config.yaml"

import datetime

date = datetime.date.today()

##### Configuration #####

if config.get("output_path"):
    config["output_path"] = config["output_path"].rstrip("/")
else:
    config["output_path"] = "analysis"

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
        expand(config["output_path"] + "/4/{lineage}/trees/traits.csv", lineage=LINEAGES)

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
    threads: 16
    shell:
        """
        echo "{params.outgroup} {input.lineage_fasta} {params.lineage}"
        iqtree -m GTR+G -bb 1000 -czb \
        -o \"{params.outgroup}\" \
        -cptime 300 \
        -ntmax {threads} \
        -s {input.lineage_fasta} &> {log}
        mv {input.lineage_fasta}.treefile {output.tree}
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
          --boolean-for-trait country=UK-.* \
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
        pubdir = config["publish_path"] + "/COG_GISAID/acc_lineages",
    output:
        config["output_path"] + "/4/{lineage}/trees/cut_out_trees_done"
    log:
        config["output_path"] + "/logs/4_cut_out_trees_{lineage}.log"
    shell:
        """
        clusterfunk prune \
          --extract \
          --trait acc_lineage \
          --input {input.tree} \
          --output {params.outdir} &> {log}
        touch {output}

        mkdir -p {params.pubdir}
        cp {params.outdir}/*.tree {params.pubdir}/
        """

rule output_annotations:
    input:
        tree = rules.label_introductions.output.tree,
        trees_cut = rules.cut_out_trees.output
    params:
        lineage = "{lineage}",
    output:
        traits = config["output_path"] + "/4/{lineage}/trees/traits.csv"
    log:
        config["output_path"] + "/logs/4_output_annotations_{lineage}.log"
    shell:
        """
        clusterfunk extract_tip_annotations \
          --traits country lineage uk_lineage acc_lineage \
          --input {input.tree} \
          --output {output.traits} &> {log}
        """

