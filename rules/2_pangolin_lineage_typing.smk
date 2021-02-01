

rule uk_normal_pangolin:
    input:
        previous_stage = config["output_path"] + "/logs/1_summarize_preprocess_uk.log",
        fasta = rules.uk_add_dups_to_lineageless.output.fasta,
    params:
        outdir = config["output_path"] + "/2/normal_pangolin",
        tmpdir = config["output_path"] + "/2/normal_pangolin/tmp"
    output:
        lineages = config["output_path"] + "/2/normal_pangolin/lineage_report.csv"
    log:
        config["output_path"] + "/logs/2_uk_normal_pangolin.log"
    shell:
        """
        pangolin {input.fasta} \
        -p \
        --outdir {params.outdir} \
        --tempdir {params.tmpdir}  >> {log} 2>&1
        """


rule uk_add_pangolin_lineages_to_metadata:
    input:
        metadata = rules.uk_add_previous_lineages_to_metadata.output.metadata,
        lineages = rules.uk_normal_pangolin.output.lineages
    output:
        metadata = config["output_path"] + "/2/uk.with_new_lineages.csv",
    log:
        config["output_path"] + "/logs/2_uk_add_normal_pangolin_lineages_to_metadata.log"
    resources: mem_per_cpu=20000
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.lineages} \
          --index-column sequence_name \
          --join-on taxon \
          --new-columns lineage \
          --out-metadata {output.metadata} &>> {log}
        """


rule uk_output_lineage_table:
    input:
        fasta = rules.uk_filter_omitted_sequences.output.fasta,
        metadata = rules.uk_add_pangolin_lineages_to_metadata.output.metadata
    output:
        fasta = config["output_path"] + "/2/uk.matched.fasta",
        metadata = config["output_path"] + "/2/uk.matched.lineages.csv"
    log:
        config["output_path"] + "/logs/2_uk_output_full_lineage_table.log"
    resources: mem_per_cpu=20000
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name country adm1 adm2 \
                          sample_date epi_week \
                          lineage uk_lineage \
          --where-column epi_week=edin_epi_week country=adm0 \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --log-file {log} \
          --restrict
        """


rule summarize_pangolin_lineage_typing:
    input:
        fasta = rules.uk_output_lineage_table.output.fasta,
        metadata = rules.uk_output_lineage_table.output.metadata,
    params:
        grapevine_webhook = config["grapevine_webhook"],
        json_path = config["json_path"],
        date = config["date"]
    log:
        config["output_path"] + "/logs/2_summarize_pangolin_lineage_typing.log"
    shell:
        """
        echo '{{"text":"' > {params.json_path}/2_data.json
        echo "*Step 2: {params.date} COG-UK pangolin typing complete*\\n" >> {params.json_path}/2_data.json
        echo '"}}' >> {params.json_path}/2_data.json
        echo "webhook {params.grapevine_webhook}"
        curl -X POST -H "Content-type: application/json" -d @{params.json_path}/2_data.json {params.grapevine_webhook}
        """
