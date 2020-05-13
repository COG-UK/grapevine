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

rule update_metadata:
    input:
        metadata = config["output_path"] + "/3/cog_gisaid.csv",
        traits = rules.run_4_subroutine_on_lineages.output.traits,
        updated_lineages = rules.merge_and_create_new_uk_lineages.output
    params:
        export_dir = config["export_path"] + "/metadata",
        export_prefix = config["export_path"] + "/metadata/cog_global_" + config["date"],
        webhook = config["webhook"]
    output:
        traits_metadata = temp(config["output_path"] + "/5/cog_gisaid.with_traits.csv"),
        all_metadata = config["output_path"] + "/5/cog_gisaid.with_all_traits.csv"
    log:
        config["output_path"] + "/logs/5_update_metadata.log"
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.traits} \
          --index-column sequence_name \
          --join-on taxon \
          --new-columns uk_lineage acc_lineage del_lineage \
          --out-metadata {output.traits_metadata} &> {log} ;

        fastafunk add_columns \
          --in-metadata {output.traits_metadata} \
          --in-data {input.updated_lineages} \
          --index-column sequence_name \
          --join-on taxon \
          --new-columns uk_lineage \
          --out-metadata {output.all_metadata} &> {log}
        """

rule run_5_subroutine_on_lineages:
    input:
        metadata = rules.update_metadata.output.all_metadata,
        lineage = config["lineage_splits"]
    params:
        path_to_script = workflow.current_basedir,
        output_path = config["output_path"],
        publish_path = config["publish_path"],
        prefix = config["output_path"] + "/5/lineage_"
    output:
        metadata = config["output_path"] + "/5/cog_gisaid.with_all_traits.with_phylotype_traits.csv",
        full_tree = config["output_path"] + "/5/cog_gisaid_full.tree.nexus"
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
          lineages="$lineages" \
          metadata={input.metadata} &> {log}
        """


"""
the routine above outputs: /5/cog_gisaid.with_all_traits.with_phylotype_traits.csv

"""

rule uk_publish_cog:
    input:
        fasta = rules.uk_filter_low_coverage_sequences.output.fasta,
        metadata = rules.uk_update_metadata_lineages.output.metadata
    output:
        fasta = config["output_path"] + "/2/uk.regularized.fasta",
        metadata = config["output_path"] + "/2/uk.regularized.csv"
    log:
        config["output_path"] + "/logs/2_uk_output_cog.log"
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name sample_date epi_week \
                          country adm1 adm2 outer_postcode \
                          is_surveillance is_community is_hcw \
                          is_travel_history travel_history lineage special_lineage \
                          lineage_support uk_lineage \
          --where-column epi_week=edin_epi_week country=adm0 \
                         sample_date=received_date sample_date=collection_date \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --log-file {log} \
          --restrict
        """

rule uk_publish_cog_public:
    input:
        fasta = rules.uk_remove_duplicates.output.fasta,
        metadata = rules.uk_update_metadata_lineages.output.metadata,
        omit_list = rules.uk_filter_low_coverage_sequences.log
    output:
        fasta = config["output_path"] + "/2/uk.public.fasta",
        metadata = config["output_path"] + "/2/uk.public.csv"
    log:
        config["output_path"] + "/logs/2_uk_output_cog_public.log"
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name country adm1 \
                          sample_date epi_week lineage \
                          lineage_support \
          --where-column epi_week=edin_epi_week country=adm0 \
                         sample_date=received_date sample_date=collection_date \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --log-file {log} \
          --restrict

        fastafunk remove \
          --in-fasta {output.fasta} \
          --in-metadata {input.omit_list} \
          --out-fasta removed.fa
        mv removed.fa {output.fasta}
        """

rule publish_metadata:
    input:
        metadata = rules.run_5_subroutine_on_lineages.output.metadata,
    params:
        outdir = config["publish_path"] + "/COG_GISAID",
        prefix = config["publish_path"] + "/COG_GISAID/cog_gisaid",
        export_dir = config["export_path"] + "/metadata",
        export_prefix = config["export_path"] + "/metadata/cog_global_" + config["date"],
        webhook = config["webhook"]
    log:
        config["output_path"] + "/logs/5_publish_metadata.log"
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p {params.export_dir}
        cp {input.metadata} {params.prefix}_metadata.csv
        cp {input.metadata} {params.export_prefix}_metadata.csv
        echo "> Updated COG and GISAID metadata published to _{params.prefix}_metadata.csv_\\n" >> {log}
        echo "> and to _{params.export_prefix}_metadata.csv_\\n" >> {log}

        echo {params.webhook}

        echo '{{"text":"' > 5a_data.json
        echo "*Step 5: Updated complete metadata with UK lineages, acctrans and deltrans*\\n" >> 5a_data.json
        cat {log} >> 5a_data.json
        echo '"}}' >> 5a_data.json
        echo "webhook {params.webhook}"
        curl -X POST -H "Content-type: application/json" -d @5a_data.json {params.webhook}
        #rm 5a_data.json
        """


# rule generate_report:
#     input:
#         metadata = rules.run_5_subroutine_on_lineages.output.metadata
#     params:
#         name_stem = "UK_" + config["date"],
#         date = config["date"],
#         outdir = config["output_path"] + "/5/"
#     output:
#         report = config["output_path"] + "/5/UK_" + config["date"] + ".md",
#         figures = config["output_path"] + "/5/figures_" + config["date"],
#         summary = config["output_path"] + "/5/summary_files_" + config["date"]
#     log:
#         config["output_path"] + "/logs/5_generate_report.log"
#     shell:
#         """
#         generate_report --m {input.metadata} --w {params.date} --s {params.name_stem} --od {params.outdir} &> {log}
#         mv {params.name_stem}* {params.outdir}
#         """

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
        echo "*Step 5: Generate UK lineage trees is complete*\\n" >> 5_data.json
        cat {log} >> 5b_data.json
        echo '"}}' >> 5b_data.json
        echo "webhook {params.webhook}"
        curl -X POST -H "Content-type: application/json" -d @5b_data.json {params.webhook}
        #rm 5b_data.json
        """

        #echo "*Step 5: Generate report and UK lineage trees is complete*\\n" >> 5_data.json
        #cp -r {input.report} {params.outdir}/
        #cp -r {input.report} {params.export_dir2}/
        #echo "> COG UK weekly report has been published in _{params.outdir}_ and _{params.export_dir2}_\\n" >> {log}
