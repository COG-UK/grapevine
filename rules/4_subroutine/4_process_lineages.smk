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
        # iqtree -m HKY -bb 1000 -czb \



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
          --trait-columns lineage \
          country uk_lineage \
          --index-column sequence_name \
          --boolean-for-trait country='UK' country='UK' country='UK' country='UK' \
          --boolean-trait-names country_uk country_uk_acctran country_uk_deltran country_uk_maxtran\
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

rule maxtran_ancestral_reconstruction:
    input:
        tree = rules.acctran_ancestral_reconstruction.output.tree
    params:
        lineage = "{lineage}",
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.annotated.acc.max.tree"
    log:
        config["output_path"] + "/logs/4_maxtran_ancestral_reconstructionn_{lineage}.log"
    shell:
        """
        clusterfunk ancestral_reconstruction \
        --traits country_uk_maxtran \
        --maxtran-with-value True \
        --ancestral_state False \
        --input {input.tree} \
        --output {output.tree} &> {log}
        """

rule deltran_ancestral_reconstruction:
    input:
        tree = rules.maxtran_ancestral_reconstruction.output.tree
    params:
        lineage = "{lineage}",
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.annotated.acc.max.del.tree"
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
#         tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.annotated.acc.max.del.uk_lineages.tree"
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
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.annotated.acc.max.del.acc_labelled.tree"
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
        prefix = config["publish_path"] + "/COG_GISAID/cog_gisaid"
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.annotated.acc.max.del.acc_labelled.del_labelled.tree"
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

rule label_maxtran_introductions:
    input:
        tree = rules.label_deltran_introductions.output.tree
    params:
        lineage = "{lineage}",
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.annotated.acc.max.del.acc_labelled.del_labelled.max_labelled.tree"
    log:
        config["output_path"] + "/logs/4_label_acctran_introductions_{lineage}.log"
    shell:
        """
        clusterfunk label_transitions \
          --trait country_uk_maxtran \
          --to True \
          --transition-name max_lineage \
          --transition-prefix {params.lineage}_ \
          --include_root \
          --input {input.tree} \
          --output {output.tree} &> {log}
        """

# now outputs a newick tree for publication. The private nexus tree will be made in stage 5.
rule graft:
    input:
         # not sure how to pass this as a space separated list below. Also assuming the order here matches lineages
        scions = expand(config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.annotated.acc.max.del.acc_labelled.del_labelled.max_labelled.tree", lineage=sorted(LINEAGES)),
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
        tree = rules.label_maxtran_introductions.output.tree,
    params:
        lineage = "{lineage}",
    output:
        traits = config["output_path"] + "/4/{lineage}/traits.csv"
    log:
        config["output_path"] + "/logs/4_output_annotations_{lineage}.log"
    shell:
        """
        clusterfunk extract_tip_annotations \
          --traits country lineage uk_lineage acc_lineage del_lineage max_lineage\
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
