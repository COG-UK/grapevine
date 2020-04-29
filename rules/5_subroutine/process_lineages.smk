configfile: workflow.current_basedir + "/config.yaml"

import os

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
        expand(config["output_path"] + "/5/{lineage}/cut_out_trees_done", lineage=LINEAGES)


rule annotate_tree:
    input:
        tree = rules.phylotype_tree.output.tree, #need to change
        metadata = config["metadata"]
    params:
        lineage = "{lineage}",
    output:
        tree = config["output_path"] + "/5/{lineage}/cog_gisaid_{lineage}.annotated.tree"
    log:
        config["output_path"] + "/logs/5_annotate_{lineage}.log"
    shell:
        """
        clusterfunk annotate_tips \
          --in-metadata {input.metadata} \
          --trait-columns uk_lineage \
          --index-column sequence_name \
          --input {input.tree} \
          --output {output.tree} &> {log}
        """


rule cut_out_trees:
    input:
        tree = rules.annotate_tree.output.tree
    params:
        lineage = "{lineage}",
        outdir = config["output_path"] + "/5/{lineage}/trees",
        pubdir = config["publish_path"] + "/COG_GISAID/acc_lineages/{lineage}",
    output:
        config["output_path"] + "/5/{lineage}/cut_out_trees_done"
    log:
        config["output_path"] + "/logs/5_cut_out_trees_{lineage}.log"
    threads: 8
    shell:
        """
        clusterfunk prune \
          --extract \
          --trait uk_lineage \
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

