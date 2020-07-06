import os
import pandas as pd

LINEAGES = []
LINEAGES_df = pd.read_csv(config["lineage_splits"])
for i,row in LINEAGES_df.iterrows():
    LINEAGES.append(row["lineage"])


rule merge_and_create_new_uk_lineages:
    input:
        config["output_path"] + "/4/all_traits.csv"
    output:
        config["output_path"] + "/5/updated_traits.csv"
    log:
        config["output_path"] + "/logs/5_merge_and_create_new_uk_lineages.log"
    shell:
        """
        datafunk curate_lineages -i {input} -o {output} &> {log}
        """


rule step_5_generate_sankey_plot:
    input:
        old_traits = config["output_path"] + "/4/all_traits.csv",
        new_traits = config["output_path"] + "/5/updated_traits.csv",
    output:
        links = config["output_path"] + "/5/sankey_links.txt",
        plot = config["output_path"] + "/5/sankey.html"
    params:
        python_script = os.path.join(workflow.current_basedir, "../utilities/get_sankey_links.py"),
        R_script = os.path.join(workflow.current_basedir, "../utilities/plot_sankey.R")
    log:
        config["output_path"] + "/logs/5_generate_sankey_plot.log"
    shell:
        """
        python {params.python_script} {input.old_traits} {input.new_traits} {output.links} &> {log}
        Rscript {params.R_script} {output.links} {output.plot} &>> {log}
        """


rule five_update_global_lineage_metadata:
    input:
        metadata = config["output_path"] + "/3/cog_gisaid.lineages.csv",
        global_lineages = config["global_lineages"],
        new_global_lineages = config["output_path"] + "/2/normal_pangolin/lineage_report.csv"
    output:
        metadata_temp = temp(config["output_path"] + "/5/cog_gisaid.global.lineages.with_all_traits.temp.csv"),
        metadata = config["output_path"] + "/5/cog_gisaid.global.lineages.with_all_traits.csv"
    log:
        config["output_path"] + "/logs/5_five_update_global_lineage_metadata.log"
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.global_lineages} \
          --index-column sequence_name \
          --join-on taxon \
          --new-columns lineage lineage_support lineages_version \
          --where-column lineage_support=UFbootstrap \
          --out-metadata {output.metadata_temp} &> {log}

        fastafunk add_columns \
          --in-metadata {output.metadata_temp} \
          --in-data {input.new_global_lineages} \
          --index-column sequence_name \
          --join-on taxon \
          --new-columns lineage lineage_support lineages_version \
          --where-column lineage_support=UFbootstrap \
          --out-metadata {output.metadata} &>> {log}
        """


rule update_lineage_metadata:
    input:
        metadata = rules.five_update_global_lineage_metadata.output.metadata,
        traits = rules.output_annotations.output.traits,
        updated_lineages = rules.merge_and_create_new_uk_lineages.output
    output:
        traits_metadata = temp(config["output_path"] + "/5/cog_gisaid.lineages.with_traits.csv"),
        all_metadata = config["output_path"] + "/5/cog_gisaid.lineages.with_all_traits.csv"
    log:
        config["output_path"] + "/logs/5_update_metadata.log"
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.traits} \
          --index-column sequence_name \
          --join-on taxon \
          --new-columns acc_lineage del_lineage acc_introduction del_introduction \
          --out-metadata {output.traits_metadata} &> {log} ;

        fastafunk add_columns \
          --in-metadata {output.traits_metadata} \
          --in-data {input.updated_lineages} \
          --index-column sequence_name \
          --join-on taxon \
          --new-columns uk_lineage microreact_lineage \
          --out-metadata {output.all_metadata} &> {log}
        """

# rule run_5_subroutine_on_lineages:
#     input:
#         metadata = rules.update_lineage_metadata.output.all_metadata,
#     params:
#         path_to_script = workflow.current_basedir,
#         output_path = config["output_path"],
#         publish_path = config["publish_path"],
#         export_path = config["export_path"],
#     output:
#         metadata = config["output_path"] + "/5/cog_gisaid.lineages.with_all_traits.with_phylotype_traits.csv",
#         full_tree = config["output_path"] + "/5/cog_gisaid_full.tree.nexus",
#     log:
#         config["output_path"] + "/logs/5_run_5_subroutine_on_lineages.log"
#     threads: 40
#     shell:
#         """
#         snakemake --nolock \
#           --snakefile {params.path_to_script}/5_subroutine/5_process_lineages.smk \
#           --cores {threads} \
#           --configfile {params.path_to_script}/5_subroutine/config.yaml \
#           --config \
#           output_path={params.output_path} \
#           publish_path={params.publish_path} \
#           export_path={params.export_path} \
#           metadata={input.metadata} &> {log}
#         """

################################################################################
rule step_5_annotate_tree:
    input:
        tree = rules.merge_sibling_del_introduction.output.tree,
        metadata = rules.update_lineage_metadata.output.all_metadata,
    output:
        tree=config["output_path"] + "/5/cog_gisaid_grafted.annotated.tree",
    log:
        config["output_path"] + "/logs/5_annotate.log"
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
        tree=rules.step_5_annotate_tree.output.tree
    params:
        outdir=config["output_path"] + "/5/trees",
    output:
        directory(config["output_path"] + "/5/trees")
    log:
        config["output_path"] + "/logs/5_cut_out_trees.log"
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
        """


rule phylotype_cut_trees:
    input:
        tree=config["output_path"] + "/5/trees/uk_lineage_UK{i}.tree"
    output:
        tree=config["output_path"] + "/5/phylotyped_trees/uk_lineage_UK{i}.tree"
    params:
        collapse=5E-6,
        threshold=2E-5,
        i="{i}"
    log:
        config["output_path"] + "/logs/5_phylotype_UK{i}.log"
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
        tree=config["output_path"] + "/5/phylotyped_trees/uk_lineage_UK{i}.tree"
    output:
        traits=config["output_path"] + "/5/phylotyped_trees/uk_lineage_UK{i}.csv"
    log:
        config["output_path"] + "/logs/5_traits_UK{i}.log"
    shell:
        """
        clusterfunk extract_tip_annotations \
          --traits country phylotype \
          --input {input.tree} \
          --output {output.traits} &> {log}
        """


def aggregate_input_csv(wildcards):
    checkpoint_output_directory = checkpoints.cut_out_trees.get(**wildcards).output[0]
    # print(checkpoints.cut_out_trees.get(**wildcards).output[0])
    required_files = expand( "%s/5/phylotyped_trees/uk_lineage_UK{i}.csv" %(config["output_path"]),
                            i=glob_wildcards(os.path.join(checkpoint_output_directory, "uk_lineage_UK{i}.tree")).i)
    return (required_files)
#
# def aggregate_input_trees(wildcards):
#     checkpoint_output_directory = checkpoints.cut_out_trees.get(**wildcards).output[0]
#     print(checkpoints.cut_out_trees.get(**wildcards).output[0])
#     lineage = wildcards.lineage
#     required_files = expand( "%s/5/%s/phylotyped_trees/uk_lineage_UK{i}.tree" %(config["output_path"],lineage),
#                             i=glob_wildcards(os.path.join(checkpoint_output_directory, "uk_lineage_UK{i}.tree")).i)
#     return (sorted(required_files))
#
# def aggregate_input_labels(wildcards):
#     checkpoint_output_directory = checkpoints.cut_out_trees.get(**wildcards).output[0]
#     print(checkpoints.cut_out_trees.get(**wildcards).output[0])
#     labels = expand( "UK{i}",i=glob_wildcards(os.path.join(checkpoint_output_directory, "uk_lineage_UK{i}.tree")).i)
#     return (sorted(labels))


rule combine_phylotypes_csv:
    input:
        files=aggregate_input_csv
    output:
        phylotype_csv=config["output_path"] + "/5/UK_phylotypes.csv"
    log:
        config["output_path"] + "/logs/5_traits_combine_phylotype_csv.log"
    run:
        dfs = [pd.read_csv(x) for x in input.files]
        result = pd.concat(dfs)
        result.to_csv(output[0], index=False)


rule merge_with_metadata:
    input:
        metadata = rules.update_lineage_metadata.output.all_metadata,
        traits = rules.combine_phylotypes_csv.output.phylotype_csv
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
        tree=rules.step_5_annotate_tree.output.tree,
        metadata = rules.merge_with_metadata.output.metadata
    output:
        annotated_tree = config["output_path"] + "/5/cog_gisaid_full.tree.nexus",
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
        """

################################################################################

rule summarize_define_uk_lineages_and_cut_out_trees:
    input:
        annotated_tree = rules.annotate_phylotypes.output.annotated_tree,
        sankey_plot = rules.step_5_generate_sankey_plot.output.plot
    params:
        webhook = config["webhook"],
    log:
        config["output_path"] + "/logs/5_summarize_define_uk_lineages_and_cut_out_trees.log"
    shell:
        """
        echo "5_subroutine complete" &>> {log}

        echo '{{"text":"' > 5b_data.json
        echo "*Step 5: Generate UK lineage trees is complete*\\n" >> 5b_data.json
        echo "> _Look at this Sankey plot_: {input.sankey_plot}\\n" >> 5b_data.json
        echo '"}}' >> 5b_data.json
        echo "webhook {params.webhook}"
        curl -X POST -H "Content-type: application/json" -d @5b_data.json {params.webhook}
        """
