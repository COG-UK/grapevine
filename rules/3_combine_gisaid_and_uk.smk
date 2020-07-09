
rule gisaid_output_lineage_table:
    input:
        fasta =  config["GISAID_background_fasta"],
        metadata = config["GISAID_background_metadata"],
    output:
        fasta = config["output_path"] + "/3/gisaid.matched.fasta",
        metadata = config["output_path"] + "/3/gisaid.matched.lineages.csv",
    log:
        config["output_path"] + "/logs/3_gisaid_output_lineage_table.log"
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name country adm1 adm2 \
                          sample_date epi_week \
                          lineage uk_lineage \
          --where-column uk_omit=is_uk sample_date=covv_collection_date \
                                 epi_week=edin_epi_week country=edin_admin_0 \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --log-file {log} \
          --restrict
        """


rule combine_gisaid_and_cog:
    input:
        previous_stage = config["output_path"] + "/logs/2_summarize_pangolin_lineage_typing.log",
        gisaid_fasta = rules.gisaid_output_lineage_table.output.fasta,
        gisaid_metadata = rules.gisaid_output_lineage_table.output.metadata,
        uk_fasta = rules.uk_output_lineage_table.output.fasta,
        uk_metadata = rules.uk_output_lineage_table.output.metadata
    output:
        fasta = config["output_path"] + "/3/cog_gisaid.fasta",
        metadata = config["output_path"] + "/3/cog_gisaid.lineages.csv"
    log:
        config["output_path"] + "/logs/3_combine_gisaid_and_cog.log"
    shell:
        """
        fastafunk merge \
          --in-fasta {input.gisaid_fasta} {input.uk_fasta} \
          --in-metadata {input.gisaid_metadata} {input.uk_metadata} \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --index-column sequence_name \
          --log-file {log}
        """


rule summarize_combine_gisaid_and_cog:
    input:
        fasta = rules.combine_gisaid_and_cog.output.fasta,
        metadata = rules.combine_gisaid_and_cog.output.metadata,
    params:
        grapevine_webhook = config["grapevine_webhook"],
    log:
        config["output_path"] + "/logs/3_summarize_combine_gisaid_and_cog.log"
    shell:
        """
        echo "> Number of sequences in combined COG and GISAID matched files: $(cat {input.fasta} | grep ">" | wc -l)\\n" &>> {log}
        echo "> \\n" &>> {log}

        echo '{{"text":"' > 3_data.json
        echo "*Step 3: Combine COG-UK and GISAID data complete*\\n" >> 3_data.json
        cat {log} >> 3_data.json
        echo '"}}' >> 3_data.json
        echo "webhook {params.grapevine_webhook}"
        curl -X POST -H "Content-type: application/json" -d @3_data.json {params.grapevine_webhook}
        """
