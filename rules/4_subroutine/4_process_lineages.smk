configfile: workflow.current_basedir + "/config.yaml"

import os
import pandas as pd

##### Configuration #####

config["publish_path"] = os.path.abspath(config["publish_path"])

LINEAGES = config["lineages"].split()[1:]
OUTGROUPS = config["lineage_specific_outgroups"].split()[1:]

lineage_to_outgroup_map = {}
for i,lin in enumerate(LINEAGES):
    lineage_to_outgroup_map[lin] = OUTGROUPS[i]

print("lineages", LINEAGES)
print("outgroups", OUTGROUPS)

##### Target rules #####

rule all:
    input:
        traits=config["output_path"] + "/4/all_traits.csv",
        tree=config["output_path"] + "/4/cog_gisaid_full.tree.public.newick"


rule fasttree:
    input:
        lineage_fasta = config["output_path"] + "/4/lineage_{lineage}.fasta"
    params:
        lineage = "{lineage}",
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.unrooted.tree",
    log:
        config["output_path"] + "/logs/4_fasttree_{lineage}.log"
    resources: mem_per_cpu=10000
    threads: 3
    shell:
        """
        echo "{input.lineage_fasta} {params.lineage}"

        export OMP_NUM_THREADS={threads}

        FastTreeMP -nosupport -nt {input.lineage_fasta} > {output.tree} 2> {log}
        """


rule root_tree:
    input:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.unrooted.tree",
    params:
        lineage = "{lineage}",
        outgroup = lambda wildcards: lineage_to_outgroup_map[wildcards.lineage]
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.tree",
    log:
        config["output_path"] + "/logs/4_root_tree_{lineage}.log",
    shell:
        """
        echo "{params.lineage} {params.outgroup}"

        clusterfunk root \
          --in-format newick \
          -i {input.tree} \
          --out-format newick \
          -o {output.tree} \
          --outgroup {params.outgroup} &>> {log}
        """


rule graft:
    input:
        # not sure how to pass this as a space separated list below. Also assuming the order here matches lineages
        scions = expand(config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.tree", lineage=sorted(LINEAGES)),
        guide_tree = config["guide_tree"]
    params:
        lineages = sorted(LINEAGES),
    output:
        grafted_tree = config["output_path"] + '/4/cog_gisaid_grafted.tree',
    log:
        config["output_path"] + "/logs/4_graft.log",
    run:
        if len(input.scions) > 1:
            shell("""
                clusterfunk graft \
                  --scions {input.scions} \
                  --scion-annotation-name scion_lineage \
                  --annotate-scions {params.lineages} \
                  --input {input.guide_tree} \
                  --in-format newick \
                  --out-format newick \
                  --output {output.grafted_tree} &> {log}
                """)

        else:
            shell("""
                clusterfunk reformat \
                  -i {input.scions} \
                  --in-format newick \
                  --out-format newick \
                  -o {output.grafted_tree} &> {log}
            """)

rule sort:
    input:
        grafted_tree = rules.graft.output.grafted_tree,
    output:
        sorted_tree = config["output_path"] + '/4/cog_gisaid_full.tree.public.newick',
    log:
        config["output_path"] + "/logs/4_sort.log",
    shell:
        """
        gotree rotate sort -i {input.grafted_tree} -o {output.sorted_tree} &> {log}
        """

rule annotate_tree:
    input:
        tree = rules.graft.output.grafted_tree,
        metadata = config["metadata"]
    params:
        collapse=0.000005,
    output:
        tree = config["output_path"] + "/4/cog_gisaid_grafted.annotated.tree"
    log:
        config["output_path"] + "/logs/4_annotate_tree.log"
    shell:
        """
        clusterfunk annotate_tips \
          --in-metadata {input.metadata} \
          --trait-columns country uk_lineage \
          --index-column sequence_name \
          --boolean-for-trait country='UK' country='UK' country='UK' country='UK' \
          --boolean-trait-names country_uk country_uk_acctran country_uk_deltran \
          --in-format newick \
          --out-format nexus \
          --collapse_to_polytomies {params.collapse} \
          --input {input.tree} \
          --output {output.tree} &> {log}
        """

rule acctran_ancestral_reconstruction:
    input:
        tree = rules.annotate_tree.output.tree
    output:
        tree = config["output_path"] + "/4/cog_gisaid_grafted.annotated.acc.tree"
    log:
        config["output_path"] + "/logs/4_acctran_ancestral_reconstruction.log"
    shell:
        """
        clusterfunk ancestral_reconstruction \
        --traits country_uk_acctran \
        --acctran \
        --ancestral-state False \
        --input {input.tree} \
        --output {output.tree} &> {log}
        """

rule deltran_ancestral_reconstruction:
    input:
        tree = rules.acctran_ancestral_reconstruction.output.tree
    output:
        tree = config["output_path"] + "/4/cog_gisaid_grafted.annotated.acc.del.tree"
    log:
        config["output_path"] + "/logs/4_deltran_ancestral_reconstruction.log"
    shell:
        """
        clusterfunk ancestral_reconstruction \
        --traits country_uk_deltran \
        --deltran \
        --ancestral-state False \
        --input {input.tree} \
        --output {output.tree} &> {log}
        """

rule label_acctran_introductions:
    input:
        tree = rules.deltran_ancestral_reconstruction.output.tree
    output:
        tree = config["output_path"] + "/4/cog_gisaid_grafted.annotated.acc.del.acc_labelled.tree"
    log:
        config["output_path"] + "/logs/4_label_acctran_introductions.log"
    shell:
        """
        clusterfunk label_transitions \
          --trait country_uk_acctran \
          --to True \
          --transition-name acc_introduction \
          --transition-prefix acc_trans_ \
          --include_root \
          --stubborn \
          --input {input.tree} \
          --output {output.tree} &> {log}
        """

rule label_deltran_introductions:
    input:
        tree = rules.label_acctran_introductions.output.tree
    output:
        tree = config["output_path"] + "/4/cog_gisaid_grafted.annotated.acc.del.acc_labelled.del_labelled.tree"
    log:
        config["output_path"] + "/logs/4_label_deltran_introductions.log"
    shell:
        """
        clusterfunk label_transitions \
          --trait country_uk_deltran \
          --to True \
          --transition-name del_introduction \
          --transition-prefix del_trans_ \
          --stubborn \
          --input {input.tree} \
          --output {output.tree} &> {log}
        """

rule merge_sibling_acc_introduction:
    input:
        tree = rules.label_deltran_introductions.output.tree
    output:
        tree = config["output_path"] + "/4/cog_gisaid_grafted.annotated.acc.del.acc_labelled.del_labelled.acc_merged.tree"
    log:
        config["output_path"] + "/logs/4_merge_acctran_introductions.log"
    shell:
        """
        clusterfunk merge_transitions \
          --trait-to-merge acc_introduction \
          --merged-trait-name acc_lineage \
          --merge-siblings \
          --input {input.tree} \
          --output {output.tree} &> {log}
        """

rule merge_sibling_del_introduction:
    input:
        tree = rules.merge_sibling_acc_introduction.output.tree
    params:
        outdir = config["publish_path"] + "/COG_GISAID",
    output:
        tree = config["output_path"] + "/4/cog_gisaid_grafted.annotated.acc.del.acc_labelled.del_labelled.acc_merged.del_merged.tree"
    log:
        config["output_path"] + "/logs/4_merge_deltran_introductions.log"
    shell:
        """
        clusterfunk merge_transitions \
          --trait-to-merge del_introduction \
          --merged-trait-name del_lineage \
          --merge-siblings \
          --input {input.tree} \
          --output {output.tree} &> {log}

        mkdir -p {params.outdir}
        """


    # shell:
    #     """
    #     clusterfunk graft \
    #       --scions {input.scions} \
    #       --scion-annotation-name scion_lineage \
    #       --annotate-scions {params.lineages} \
    #       --input {input.guide_tree} \
    #       --out-format newick \
    #       --output {output.grafted_tree} &> {log}
    #
    #     clusterfunk sort \
    #       --in-format newick \
    #       -i {output.grafted_tree} \
    #       --out-format newick \
    #       -o {output.ordered_tree} &>> {log}
    #     """


rule output_annotations:
    input:
        tree = rules.merge_sibling_del_introduction.output.tree,
    output:
        traits = config["output_path"] + "/4/all_traits.csv"
    log:
        config["output_path"] + "/logs/4_output_annotations.log"
    shell:
        """
        clusterfunk extract_tip_annotations \
          --traits country uk_lineage acc_introduction acc_lineage del_introduction del_lineage \
          --input {input.tree} \
          --output {output.traits} &> {log}
        """
