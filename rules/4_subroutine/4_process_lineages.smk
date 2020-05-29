configfile: workflow.current_basedir + "/config.yaml"

import os
import pandas as pd

##### Configuration #####

config["publish_path"] = os.path.abspath(config["publish_path"])

LINEAGES = config["lineages"].split()[1:]
OUTGROUPS = config["lineage_specific_outgroups"].split()[1:]

lineage_to_outgroup_map = {}
for i,lin in enumerate(LINEAGES):
    lineage_to_outgroup_map[lin] = OUTGROUPS[i].replace("/","_")

print("lineages", LINEAGES)
print("outgroups", OUTGROUPS)

##### Target rules #####

rule all:
    input:
        traits=config["output_path"] + "/4/all_traits.csv",
        tree=config["output_path"] + "/4/cog_gisaid_full.tree.public.newick"


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
        iqtree -m HKY -czb \
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
        # iqtree -m HKY -bb 1000 -czb \

#
# rule subroutine_4_add_global_lineages_to_metadata:
#     input:
#         metadata = config["metadata"],
#         global_lineages = config["globallineages"],
#         new_global_lineages = config["output_path"] + "/2/normal_pangolin/lineage_report.csv"
#     output:
#         metadata_temp = temp(config["output_path"] + "/4/cog_gisaid.global.lineages.temp.csv"),
#         metadata = config["output_path"] + "/4/cog_gisaid.global.lineages.csv"
#     log:
#         config["output_path"] + "/logs/4_add_global_lineages_to_metadata.log"
#     shell:
#         """
#         fastafunk add_columns \
#           --in-metadata {input.metadata} \
#           --in-data {input.global_lineages} \
#           --index-column sequence_name \
#           --join-on taxon  \
#           --new-columns lineage lineage_support lineages_version \
#           --where-column lineage_support=UFbootstrap \
#           --out-metadata {output.metadata_temp} &>> {log}
#
#         fastafunk add_columns \
#           --in-metadata {output.metadata_temp} \
#           --in-data {input.new_global_lineages} \
#           --index-column sequence_name \
#           --join-on taxon \
#           --new-columns lineage lineage_support lineages_version \
#           --where-column lineage_support=UFbootstrap \
#           --out-metadata {output.metadata} &>> {log}
#         """


rule annotate_tree:
    input:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.tree",
        metadata = config["metadata"]
    params:
        lineage = "{lineage}",
        collapse=0.000005,

    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.annotated.tree"
    log:
        config["output_path"] + "/logs/4_annotate_{lineage}.log"
    shell:
        """
        clusterfunk annotate_tips \
          --in-metadata {input.metadata} \
          --trait-columns country \
          --index-column sequence_name \
          --boolean-for-trait country='UK' country='UK' country='UK' country='UK' \
          --boolean-trait-names country_uk country_uk_acctran country_uk_deltran\
          --in-format newick \
          --out-format nexus \
          --collapse_to_polytomies {params.collapse} \
          --input {input.tree} \
          --output {output.tree} &> {log}
        """

rule acctran_ancestral_reconstruction:
    input:
        tree = rules.annotate_tree.output.tree
    params:
        lineage = "{lineage}",
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.annotated.acc.tree"
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
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.annotated.acc.del.tree"
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

# rule push_lineage_to_tips:
#     input:
#         tree = rules.deltran_ancestral_reconstruction.output.tree
#     params:
#         lineage = "{lineage}",
#     output:
#         tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.annotated.acc.del.uk_lineages.tree"
#     log:
#         config["output_path"] + "/logs/4_push_lineage_to_tips_{lineage}.log"
#     shell:
#         """
#         clusterfunk push_annotations_to_tips \
#           --traits uk_lineage \
#           --stop-where-trait country_uk_acctran=False \
#           --input {input.tree} \
#           --output {output.tree} &> {log}
#         """

rule label_acctran_introductions:
    input:
        tree = rules.deltran_ancestral_reconstruction.output.tree
    params:
        lineage = "{lineage}",
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.annotated.acc.del.acc_labelled.tree"
    log:
        config["output_path"] + "/logs/4_label_acctran_introductions_{lineage}.log"
    shell:
        """
        clusterfunk label_transitions \
          --trait country_uk_acctran \
          --to True \
          --transition-name acc_introduction \
          --transition-prefix {params.lineage}_ \
          --include_root \
          --stubborn \
          --input {input.tree} \
          --output {output.tree} &> {log}
        """

rule label_deltran_introductions:
    input:
        tree = rules.label_acctran_introductions.output.tree
    params:
        lineage = "{lineage}",
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.annotated.acc.del.acc_labelled.del_labelled.tree"
    log:
        config["output_path"] + "/logs/4_label_deltran_introductions_{lineage}.log"
    shell:
        """
        clusterfunk label_transitions \
          --trait country_uk_deltran \
          --to True \
          --transition-name del_introduction \
          --transition-prefix {params.lineage}_ \
          --stubborn \
          --input {input.tree} \
          --output {output.tree} &> {log}

        """

rule merge_sibling_acc_introduction:
    input:
        tree = rules.label_deltran_introductions.output.tree
    params:
        lineage = "{lineage}",

    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.annotated.acc.del.acc_labelled.del_labelled.acc_merged.tree"
    log:
        config["output_path"] + "/logs/4_merge_acctran_introductions_{lineage}.log"
    shell:
        """
        clusterfunk merge_transitions \
          --trait-to-merge acc_introduction \
          --merged-trait-name acc_lineage \
          --max-merge 1 \
          --prefix {params.lineage}_ \
          --input {input.tree} \
          --output {output.tree} &> {log}
        """

rule merge_sibling_del_introduction:
    input:
        tree = rules.merge_sibling_acc_introduction.output.tree
    params:
        lineage = "{lineage}",
        outdir = config["publish_path"] + "/COG_GISAID",
        prefix = config["publish_path"] + "/COG_GISAID/cog_gisaid"
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.annotated.acc.del.acc_labelled.del_labelled.acc_merged.del_merged.tree"
    log:
        config["output_path"] + "/logs/4_merge_deltran_introductions_{lineage}.log"
    shell:
        """
        clusterfunk merge_transitions \
          --trait-to-merge del_introduction \
          --merged-trait-name del_lineage \
          --max-merge 1 \
          --prefix {params.lineage}_ \
          --input {input.tree} \
          --output {output.tree} &> {log}

        mkdir -p {params.outdir}
        cp {output.tree} {params.prefix}_lineage_{params.lineage}.tree
        """

# now outputs a newick tree for publication. The private nexus tree will be made in stage 5.
rule graft:
    input:
         # not sure how to pass this as a space separated list below. Also assuming the order here matches lineages
        scions = expand(config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.annotated.acc.del.acc_labelled.del_labelled.acc_merged.del_merged.tree", lineage=sorted(LINEAGES)),
        guide_tree = config["guide_tree"]
    params:
        lineages = sorted(LINEAGES),
    output:
        tree = config["output_path"] + '/4/cog_gisaid_full.tree.public.newick',
    log:
        config["output_path"] + "/logs/4_graft.log"
    shell:
        """
        clusterfunk graft \
        --scions {input.scions} \
        --scion_annotation_name scion_lineage \
        --annotate_scions {params.lineages} \
        --input {input.guide_tree} \
        --out-format newick \
        --output {output.tree} &> {log}
        """
# not needed since graft outputs a newick tree now.
# rule get_public_newick_tree:
#     input:
#         tree = rules.graft.output.tree
#     output:
#         tree = '/4/cog_gisaid_full.tree.public.newick'
#     log:
#         config["output_path"] + "/logs/4_publish_public_newick_tree.log"
#     shell:
#         """
#         clusterfunk reformat \
#           --input {input.tree} \
#           --output {output.tree} \
#           --out-format newick \
#           --in-format nexus &> {log}
#         """



rule output_annotations:
    input:
        tree = rules.merge_sibling_del_introduction.output.tree,
    params:
        lineage = "{lineage}",
    output:
        traits = config["output_path"] + "/4/{lineage}/traits.csv"
    log:
        config["output_path"] + "/logs/4_output_annotations_{lineage}.log"
    shell:
        """
        clusterfunk extract_tip_annotations \
          --traits country acc_introduction acc_lineage del_introduction del_lineage \
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
        dfs = [pd.read_csv(x) for x in input.list_df]
        result = pd.concat(dfs)
        result.to_csv(output[0], index=False)
