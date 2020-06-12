configfile: workflow.current_basedir + "/config.yaml"

import os
import pandas as pd


##### Configuration #####

config["publish_path"] = os.path.abspath(config["publish_path"])

LINEAGES = config["lineages"].split()[1:]

print("lineages", LINEAGES)

##### Target rules #####

rule all:
    input:
        config["output_path"] + "/5/cog_gisaid.lineages.with_all_traits.with_phylotype_traits.csv",
        config["output_path"] + "/5/cog_gisaid_full.tree.nexus",

rule annotate_tree:
    input:
        tree=config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.annotated.acc.del.acc_labelled.del_labelled.acc_merged.del_merged.tree",
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
        pubdir=config["export_path"] + "/trees/uk_lineages",
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

          echo "Look at all those trees:" >> {log}
          ls {output}/* >> {log}
          if [ ! -z "$(ls -A {params.outdir})" ]
          then
            mkdir -p {params.pubdir}
            cp {params.outdir}/*.tree {params.pubdir}/
          fi
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
        i="{i}"
    log:
        config["output_path"] + "/logs/5_phylotype_{lineage}_UK{i}.log"
    shell:
        """
        clusterfunk phylotype \
        --threshold {params.threshold} \
        --prefix UK{params.i}_1 \
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
    print(checkpoints.cut_out_trees.get(**wildcards).output[0])
    lineage = wildcards.lineage
    required_files = expand( "%s/5/%s/phylotyped_trees/uk_lineage_UK{i}.csv" %(config["output_path"],lineage),
                            i=glob_wildcards(os.path.join(checkpoint_output_directory, "uk_lineage_UK{i}.tree")).i)
    return (required_files)

def aggregate_input_trees(wildcards):
    checkpoint_output_directory = checkpoints.cut_out_trees.get(**wildcards).output[0]
    print(checkpoints.cut_out_trees.get(**wildcards).output[0])
    lineage = wildcards.lineage
    required_files = expand( "%s/5/%s/phylotyped_trees/uk_lineage_UK{i}.tree" %(config["output_path"],lineage),
                            i=glob_wildcards(os.path.join(checkpoint_output_directory, "uk_lineage_UK{i}.tree")).i)
    return (sorted(required_files))

def aggregate_input_labels(wildcards):
    checkpoint_output_directory = checkpoints.cut_out_trees.get(**wildcards).output[0]
    print(checkpoints.cut_out_trees.get(**wildcards).output[0])
    labels = expand( "UK{i}",i=glob_wildcards(os.path.join(checkpoint_output_directory, "uk_lineage_UK{i}.tree")).i)
    return (sorted(labels))

rule combine_phylotypes_csv:
    input:
        files=aggregate_input_csv
    output:
        phylotype_csv=config["output_path"] + "/5/{lineage}/UK_phylotypes.csv"
    log:
        config["output_path"] + "/logs/5_traits_{lineage}_combine_phylotype_csv.log"
    run:
        dfs = [pd.read_csv(x) for x in input.files]
        result = pd.concat(dfs)
        result.to_csv(output[0], index=False)

rule combine_lineage_csv:
    input:
        expand(config["output_path"] + "/5/{lineage}/UK_phylotypes.csv",lineage=LINEAGES)
    output:
        phylotype_csv=config["output_path"] + "/5/UK_phylotypes.csv"
    log:
        config["output_path"] + "/logs/5_combine_lineage_csv.log"
    run:
        dfs = [pd.read_csv(x) for x in input]
        result = pd.concat(dfs)
        result.to_csv(output[0], index=False)



# This should create the full tree.
rule graft_lineages:
    input:
        scions = expand(config["output_path"] + "/5/{lineage}/cog_gisaid_{lineage}.annotated.tree", lineage=sorted(LINEAGES)),
        guide_tree = config["guide_tree"]
    params:
        lineages = sorted(LINEAGES),
    output:
        tree = config["output_path"] + "/5/cog_gisaid_full.no_phylotypes.tree.nexus",
    log:
        config["output_path"] + "/logs/5_graft_lineages.log"
    shell:
        """
        clusterfunk graft \
        --scions {input.scions} \
        --scion-annotation-name scion_lineage \
        --annotate-scions {params.lineages} \
        --input {input.guide_tree} \
        --output {output.tree} &> {log}
        """


# rule get_private_nexus_tree:
#     input:
#         tree = rules.graft_lineages.output.tree
#     output:
#         tree = '/5/cog_gisaid_full.tree.private.nexus'
#     log:
#         config["output_path"] + "/logs/5_publish_private_nexus_tree.log"
#     shell:
#         """
#         clusterfunk annotate_lineages \
#           --input {input.tree} \
#           --output {output.tree} \
#           --trait lineage
#         """

rule merge_with_metadata:
    input:
        metadata = config["metadata"],
        traits = rules.combine_lineage_csv.output.phylotype_csv
    output:
        metadata = config["output_path"] + "/5/cog_gisaid.lineages.with_all_traits.with_phylotype_traits.csv"
    log:
        config["output_path"] + "/logs/5_merge_with_metadata.log"
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.traits} \
          --index-column sequence_name \
          --join-on taxon \
          --new-columns phylotype \
          --out-metadata {output.metadata} &> {log}
        """

rule annotate_phylotypes:
    input:
        tree=rules.graft_lineages.output.tree,
        metadata = rules.merge_with_metadata.output.metadata
    output:
        annotated_tree = config["output_path"] + "/5/cog_gisaid_full.unordered.tree.nexus",
        ordered_tree = config["output_path"] + "/5/cog_gisaid_full.tree.nexus",
    log:
        config["output_path"] + "/logs/5_annotate_phylotypes.log"
    shell:
        """
        clusterfunk annotate_tips \
          --in-metadata {input.metadata} \
          --trait-columns phylotype \
          --index-column sequence_name \
          --input {input.tree} \
          --output {output.annotated_tree} &> {log}

        clusterfunk sort \
          --in-format nexus \
          -i {output.annotated_tree} \
          --out-format nexus \
          -o {output.ordered_tree} &>> {log}
        """
