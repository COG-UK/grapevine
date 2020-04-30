configfile: workflow.current_basedir + "/config.yaml"

import os

##### Configuration #####

config["publish_path"] = os.path.abspath(config["publish_path"])

LINEAGES = config["lineages"].split()

print("lineages", LINEAGES)

##### Target rules #####

rule all:
    input:
         expand(config["output_path"] + "/5/{lineage}/cut_out_trees_done", lineage=LINEAGES)

rule annotate_tree:
    input:
         tree=config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.phylotyped.tree",
         metadata=config["metadata"]
    params:
          lineage="{lineage}",
    output:
          tree=config["output_path"] + "/5/{lineage}/cog_gisaid_{lineage}.annotated.tree"
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




"""
This is a checkpoint so the dag will get remade after this run. At that point aggregate_input_csv (below) will
ask for trait files with wildcards that match the number of uk lineage (i). These were made in this step and everything
should flow from there
"""
checkpoint cut_out_trees:
    input:
         tree=rules.annotate_tree.output.tree
    params:
          lineage="{lineage}",
          outdir=config["output_path"] + "/5/{lineage}/trees",
          pubdir=config["publish_path"] + "/COG_GISAID/trees",
    output:
          directory(config["output_path"] + "/5/{lineage}/trees")
    log:
       config["output_path"] + "/logs/5_cut_out_trees_{lineage}.log"
    threads: 40
    shell:
         """
         clusterfunk prune \
           --extract \
           --trait uk_lineage \
           --input {input.tree} \
           --threads {threads} \
           --output {params.outdir} &> {log}
         """

rule phylotype_cut_trees:
    input:
         tree=config["output_path"] + "/5/{lineage}/trees/uk_lineage_UK{i}.tree"
    output:
          tree=config["output_path"] + "/5/{lineage}/phylotyped_trees/uk_lineage_UK{i}.tree"
    params:
          lineage="{lineage}",
          collapse=5E-6,
          threshold=2E-5,
    log:
       config["output_path"] + "/logs/5_phylotype_{lineage}_UK{i}.log"
    shell:
         """
         clusterfunk phylotype \
         --format newick \
         --threshold {params.threshold} \
         --prefix UK{i}_1 \
         --input {input.tree} \
         --output {output.tree} &> {log}
         """

rule get_uk_phylotypes_csv:
    input:
         tree=config["output_path"] + "/5/{lineage}/phylotyped_trees/uk_lineage_UK{i}.tree"
    params:
          lineage="{lineage}",
    output:
          traits=config["output_path"] + "/5/{lineage}/phylotyped_trees/uk_lineage_UK{i}.csv"
    log:
       config["output_path"] + "/logs/5_traits_{lineage}_UK{i}.log"
    shell:
         """
         clusterfunk extract_tip_annotations \
           --traits country phylotype \
           --input {input.tree} \
           --output {output.traits} &> {log}
         """


def aggregate_input_csv(wildcards):
    checkpoint_output_directory = checkpoints.cut_out_trees.get(**wildcards).output[0]
    lineage = wildcards.lineage
    required_files = expand(config["output_path"] + "/5/%s/phylotyped_trees/uk_lineage_UK{i}.csv" % (lineage),
                            i=glob_wildcards(os.path.join(checkpoint_output_directory, "uk_lineage_UK{i}.tree")).i)

    return (required_files)


rule combine_phylotypes_csv:
    input:
         aggregate_input_csv
    output:
          phylotype_csv=config["output_path"] + "/5/{lineage}/UK_phylotypes.csv"
    log:
       config["output_path"] + "/logs/5_traits_{lineage}_combine_phylotype_csv.log"
    run:
        """
        result = pd.concat({input})
        result.to_csv({output}, index=False)
        """


rule merge_with_metadata:
    input:
        rules.combine_phylotypes_csv.ouput.traits.phylotype_csv
    shell:
         """
         echo "We still need to implement this to merge the UK phylotypes into the metatdata csv"
         """