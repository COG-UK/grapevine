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


rule five_update_global_lineage_metadata:
    input:
        metadata = config["output_path"] + "/3/cog_gisaid.lineages.csv",
        global_lineages = config["global_lineages"]
    output:
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
          --out-metadata {output.metadata} &> {log}
        """


rule update_lineage_metadata:
    input:
        metadata = rules.five_update_global_lineage_metadata.output.metadata,
        traits = rules.run_4_subroutine_on_lineages.output.traits,
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
          --new-columns uk_lineage acc_lineage del_lineage acc_introduction del_introduction \
          --out-metadata {output.traits_metadata} &> {log} ;

        fastafunk add_columns \
          --in-metadata {output.traits_metadata} \
          --in-data {input.updated_lineages} \
          --index-column sequence_name \
          --join-on taxon \
          --new-columns uk_lineage microreact_lineage \
          --out-metadata {output.all_metadata} &> {log}
        """

rule run_5_subroutine_on_lineages:
    input:
        metadata = rules.update_lineage_metadata.output.all_metadata,
        lineage = config["lineage_splits"]
    params:
        path_to_script = workflow.current_basedir,
        output_path = config["output_path"],
        publish_path = config["publish_path"],
        prefix = config["output_path"] + "/5/lineage_",
        guide_tree = config["guide_tree"],
    output:
        metadata = config["output_path"] + "/5/cog_gisaid.lineages.with_all_traits.with_phylotype_traits.csv",
        full_tree = config["output_path"] + "/5/cog_gisaid_full.tree.nexus",
    log:
        config["output_path"] + "/logs/5_run_5_subroutine_on_lineages.log"
    threads: 40
    shell:
        """
        lineages=$(cat {input.lineage} | cut -f1 -d"," | tr '\\n' '  ')
        snakemake --nolock \
          --snakefile {params.path_to_script}/5_subroutine/5_process_lineages.smk \
          --cores {threads} \
          --configfile {params.path_to_script}/5_subroutine/config.yaml \
          --config \
          output_path={params.output_path} \
          publish_path={params.publish_path} \
          guide_tree="{params.guide_tree}" \
          lineages="$lineages" \
          metadata={input.metadata} &> {log}
        """


rule summarize_generate_report_and_cut_out_trees:
    input:
        #report = rules.generate_report.output
        metadata = rules.run_5_subroutine_on_lineages.output.metadata
    params:
        webhook = config["webhook"],
        outdir = config["publish_path"] + "/COG_GISAID",
        export_dir1 = config["export_path"] + "/trees/uk_lineages",
        export_dir2 = config["export_path"] + "/reports",
    log:
        config["output_path"] + "/logs/5_summarize_generate_report_and_cut_out_trees.log"
    shell:
        """
        mkdir -p {params.export_dir1}
        mkdir -p {params.export_dir2}

        if [ ! -z "$(ls -A {params.outdir}/trees)" ]
        then
          cp {params.outdir}/trees/* {params.export_dir1}/
        fi
        echo "> UK lineage trees have been published in _{params.outdir}/trees_ and _{params.export_dir1}_\\n" >> {log}
        echo ">\\n" >> {log}

        echo '{{"text":"' > 5b_data.json
        echo "*Step 5: Generate UK lineage trees is complete*\\n" >> 5b_data.json
        echo '"}}' >> 5b_data.json
        echo "webhook {params.webhook}"
        curl -X POST -H "Content-type: application/json" -d @5b_data.json {params.webhook}
        """

        # cat {log} >> 5b_data.json


#
# rule summarize_generate_report_and_cut_out_trees:
#     input:
#         #report = rules.generate_report.output
#         metadata = rules.run_5_subroutine_on_lineages.output.metadata
#     params:
#         webhook = config["webhook"],
#         outdir = config["publish_path"] + "/COG_GISAID",
#         export_dir1 = config["export_path"] + "/trees/uk_lineages",
#         export_dir2 = config["export_path"] + "/reports",
#     log:
#         config["output_path"] + "/logs/5_summarize_generate_report_and_cut_out_trees.log"
#     shell:
#         """
#         mkdir -p {params.export_dir1}
#         mkdir -p {params.export_dir2}
#
#         if [ ! -z "$(ls -A {params.outdir}/trees)" ]
#         then
#           cp {params.outdir}/trees/* {params.export_dir1}/
#         fi
#         echo "> UK lineage trees have been published in _{params.outdir}/trees_ and _{params.export_dir1}_\\n" >> {log}
#         echo ">\\n" >> {log}
#
#
#         echo '{{"text":"' > 5b_data.json
#         echo "*Step 5: Generate UK lineage trees is complete*\\n" >> 5_data.json
#         cat {log} >> 5b_data.json
#         echo '"}}' >> 5b_data.json
#         echo "webhook {params.webhook}"
#         curl -X POST -H "Content-type: application/json" -d @5b_data.json {params.webhook}
#         """

        #echo "*Step 5: Generate report and UK lineage trees is complete*\\n" >> 5_data.json
        #cp -r {input.report} {params.outdir}/
        #cp -r {input.report} {params.export_dir2}/
        #echo "> COG UK weekly report has been published in _{params.outdir}_ and _{params.export_dir2}_\\n" >> {log}
