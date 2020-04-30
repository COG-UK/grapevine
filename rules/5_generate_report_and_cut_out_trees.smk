rule merge_and_create_new_uk_lineages:
    input:
        config["output_path"] + "/4/all_traits.csv"
    output:
        config["output_path"] + "/5/updated_traits.csv"
    log:
        config["output_path"] + "/logs/5_merge_and_create_new_uk_lineages.log"
    shell:
        """
        datafunk merge_lineages -i {input} -o {output} &> {log}
        """

rule update_metadata:
    input:
        metadata = rules.combine_gisaid_and_cog.output.metadata,
        traits = rules.run_4_subroutine_on_lineages.output,
        updated_lineages = rules.merge_and_create_new_uk_lineages.output
    output:
        traits_metadata = temp(config["output_path"] + "/4/cog_gisaid.with_traits.csv"),
        all_metadata = config["output_path"] + "/4/cog_gisaid.with_all_traits.csv"
    log:
        config["output_path"] + "/logs/4_update_metadata.log"
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.traits} \
          --index-column sequence_name \
          --join-on taxon \
          --new-columns uk_lineage acc_lineage del_lineage phylotype \
          --out-metadata {output.traits_metadata} &> {log}

        fastafunk add_columns \
          --in-metadata {output.traits_metadata} \
          --in-data {input.updated_lineages} \
          --index-column sequence_name \
          --join-on taxon \
          --new-columns uk_lineage \
          --out-metadata {output.all_metadata} &>> {log}
        """

rule run_5_subroutine_on_lineages:
    input:
        split_done = rules.split_based_on_lineages.output,
        metadata = rules.combine_gisaid_and_cog.output.metadata,
        lineage = config["lineage_splits"]
    params:
        path_to_script = workflow.current_basedir,
        output_path = config["output_path"],
        publish_path = config["publish_path"],
        prefix = config["output_path"] + "/5/lineage_"
    output:
        config["output_path"] + "/5/trees_done"
    log:
        config["output_path"] + "/logs/5_run_subroutine_on_lineages.log"
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
          lineages="$lineages" \
          lineage_specific_outgroups="$outgroups" \
          metadata={input.metadata} &> {log}

        touch {output}
        """

rule generate_report:
    pass

rule summarize_generate_report_and_cut_out_trees:
    input:
        lineage = config["lineage_splits"],
        trees_done = rules.run_subroutine_on_lineages.output,
    params:
        webhook = config["webhook"],
        outdir = config["publish_path"] + "/COG_GISAID",
    log:
        config["output_path"] + "/logs/5_summarize_make_trees.log"
    shell:
        """
        echo "> Trees have been published in {params.outdir}\\n" >> {log}

        echo '{{"text":"' > 5b_data.json
        echo "*Step 5: Construct and annotate trees completed*\\n" >> 5_data.json
        cat {log} >> 5b_data.json
        echo '"}}' >> 5b_data.json
        echo "webhook {params.webhook}"
        curl -X POST -H "Content-type: application/json" -d @5b_data.json {params.webhook}
        #rm 5b_data.json
        """