import datetime

date = datetime.date.today()

rule combine_gisaid_and_cog:
    input:
        summary_gisaid = rules.gisaid_summarize_preprocess.log,
        summary_uk = rules.uk_summarize_preprocess.log,
        gisaid_fasta = rules.gisaid_combine_previous_and_new.output.fasta,
        gisaid_metadata = rules.gisaid_combine_previous_and_new.output.metadata,
        uk_fasta = rules.uk_output_cog.output.fasta,
        uk_metadata = rules.uk_output_cog.output.metadata
    params:
        outdir = config["publish_path"] + "/COG_GISAID",
        prefix = config["publish_path"] + "/COG_GISAID/cog_gisaid_%s" %date,
    output:
        fasta = config["output_path"] + "/3/cog_gisaid.fasta",
        metadata = config["output_path"] + "/3/cog_gisaid.csv"
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

        mkdir -p {params.outdir}
        cp {output.fasta} {params.prefix}_alignment.fasta
        cp {output.metadata} {params.prefix}_metadata.csv
        """