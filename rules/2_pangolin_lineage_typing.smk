rule update_pangolin:
    output:
        temp("updated_pangolin")
    shell:
        """
        pip install --upgrade git+https://github.com/hCoV-2019/pangolin.git
        touch updated_pangolin
        """

rule uk_pangolin:
    input:
        fasta = rules.uk_filter_low_coverage_sequences.output,
        update = rules.update_pangolin.output
    params:
        outdir = config["output_path"] + "/2/pangolin"
    output:
        lineages = config["output_path"] + "/2/pangolin/lineage_report.csv"
    log:
        config["output_path"] + "/logs/2_uk_pangolin.log"
    threads: 32
    shell:
        """
        pangolin {input.fasta} \
        --threads {threads} \
        --outdir {params.outdir}
        """

rule uk_add_pangolin_lineages_to_metadata:
    input:
        metadata = rules.uk_remove_duplicates.output.metadata,
        lineages = rules.uk_pangolin.output.lineages
    output:
        metadata = config["output_path"] + "/2/uk_with_lineages.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_add_pangolin_lineages_to_metadata.log"
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.lineages} \
          --index-column sequence_name \
          --new-columns lineage bootstrap \
          --out-metadata {output.metadata} &> {log}
        """