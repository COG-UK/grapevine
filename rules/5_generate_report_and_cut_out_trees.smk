# rule update_metadata_with_uk_lineages:
#     pass
#
# rule run_subroutine_on_lineages:
#     input:
#         split_done = rules.split_based_on_lineages.output,
#         metadata = rules.combine_gisaid_and_cog.output.metadata,
#         lineage = config["lineage_splits"]
#     params:
#         path_to_script = workflow.current_basedir,
#         output_path = config["output_path"],
#         publish_path = config["publish_path"],
#         prefix = config["output_path"] + "/5/lineage_"
#     output:
#         config["output_path"] + "/5/trees_done"
#     log:
#         config["output_path"] + "/logs/5_run_subroutine_on_lineages.log"
#     threads: 16
#     shell:
#         """
#         lineages=$(cat {input.lineage} | cut -f1 -d"," | tr '\\n' '  ')
#         outgroups=$(cat {input.lineage} | cut -f2 -d"," | tr '\\n' '  ')
#         snakemake --nolock \
#           --snakefile {params.path_to_script}/5_subroutine/process_lineages.smk \
#           --cores {threads} \
#           --configfile {params.path_to_script}/5_subroutine/config.yaml \
#           --config \
#           output_path={params.output_path} \
#           publish_path={params.publish_path} \
#           lineages="$lineages" \
#           lineage_specific_outgroups="$outgroups" \
#           metadata={input.metadata} &> {log}
#
#         touch {output}
#         """
#
# rule summarize_generate_report_and_cut_out_trees:
#     input:
#         lineage = config["lineage_splits"],
#         trees_done = rules.run_subroutine_on_lineages.output,
#     params:
#         webhook = config["webhook"],
#         outdir = config["publish_path"] + "/COG_GISAID",
#     log:
#         config["output_path"] + "/logs/5_summarize_make_trees.log"
#     shell:
#         """
#         echo "> Trees have been published in {params.outdir}\\n" >> {log}
#
#         echo '{{"text":"' > 5b_data.json
#         echo "*Step 5: Construct and annotate trees completed*\\n" >> 5_data.json
#         cat {log} >> 5b_data.json
#         echo '"}}' >> 5b_data.json
#         echo "webhook {params.webhook}"
#         curl -X POST -H "Content-type: application/json" -d @5b_data.json {params.webhook}
#         #rm 5b_data.json
#         """