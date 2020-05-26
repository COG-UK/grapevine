

rule uk_pangolin:
    input:
        previous_stage = config["output_path"] + "/logs/1_summarize_preprocess_uk.log",
        fasta = rules.uk_extract_lineageless.output.fasta,
    params:
        outdir = config["output_path"] + "/2/pangolin",
        tmpdir = config["output_path"] + "/2/pangolin/tmp"
    output:
        lineages = protected(config["output_path"] + "/2/pangolin/lineage_report.csv")
    log:
        config["output_path"] + "/logs/2_uk_pangolin.log"
    threads: 40
    shell:
        """
        pangolin {input.fasta} \
        --threads {threads} \
        --outdir {params.outdir} \
        --tempdir {params.tmpdir}  >> {log} 2>&1
        """


rule uk_add_pangolin_lineages_to_metadata:
    input:
        metadata = rules.uk_add_previous_uk_lineages_to_metadata.output.metadata,
        lineages = rules.uk_pangolin.output.lineages
    output:
        metadata = config["output_path"] + "/2/uk.with_new_lineages.csv",
    log:
        config["output_path"] + "/logs/2_uk_add_pangolin_lineages_to_metadata.log"
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.lineages} \
          --index-column sequence_name \
          --join-on taxon \
          --new-columns special_lineage \
          --where-column special_lineage=lineage\
          --out-metadata {output.metadata} &>> {log}
        """


rule uk_output_lineage_table:
    input:
        fasta = rules.uk_filter_low_coverage_sequences.output.fasta,
        metadata = rules.uk_add_pangolin_lineages_to_metadata.output.metadata
    output:
        fasta = config["output_path"] + "/2/uk.matched.fasta",
        metadata = config["output_path"] + "/2/uk.matched.lineages.csv"
    log:
        config["output_path"] + "/logs/2_uk_output_lineage_table.log"
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name country adm1 adm2 \
                          sample_date epi_week \
                          special_lineage uk_lineage \
          --where-column epi_week=edin_epi_week country=adm0 \
                         sample_date=received_date sample_date=collection_date \
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
        webhook = config["webhook"],
    log:
        config["output_path"] + "/logs/2_summarize_pangolin_lineage_typing.log"
    shell:
        """
        echo '{{"text":"' > 2_data.json
        echo "*Step 2: COG-UK pangolin typing complete*\\n" >> 2_data.json
        echo '"}}' >> 2_data.json
        echo "webhook {params.webhook}"
        curl -X POST -H "Content-type: application/json" -d @2_data.json {params.webhook}
        """
