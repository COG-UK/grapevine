import datetime

date = datetime.date.today()

rule combine_gisaid_and_cog:
    input:
        summary_gisaid = rules.gisaid_summarize_preprocess.log,
        summary_uk = rules.uk_summarize_preprocess.log,
        summary_uk_pangolin = rules.uk_summarize_pangolin.log,
        gisaid_fasta = rules.gisaid_output_gisaid.output.fasta,
        gisaid_metadata = rules.gisaid_output_gisaid.output.metadata,
        uk_fasta = rules.uk_output_cog.output.fasta,
        uk_metadata = rules.uk_output_cog.output.metadata
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

rule summarize_combine:
    input:
        fasta = rules.combine_gisaid_and_cog.output.fasta,
        metadata = rules.combine_gisaid_and_cog.output.metadata,
    params:
        outdir = config["publish_path"] + "/COG_GISAID",
        prefix = config["publish_path"] + "/COG_GISAID/cog_gisaid_%s" %date,
    log:
        config["output_path"] + "/logs/3_summarize_combine.log"
    shell:
        """
        mkdir -p {params.outdir}
        cp {input.fasta} {params.prefix}_alignment.fasta
        cp {input.metadata} {params.prefix}_metadata.csv
        echo "> Combined COG and GISAID trimmed alignments published to {params.prefix}_alignment.fasta\n" >> {log}
        echo "> Combined COG and GISAID restricted metadata published to {params.prefix}_metadata.csv\n" >> {log}

        echo '{{"text":"' > 3_data.json
        echo "*Step 3: Combine COG-UK and GISAID data*\n" >> 3_data.json
        cat {log} >> 3_data.json
        echo '"}}' >> 3_data.json
        curl -X POST -H "Content-type: application/json" -d @1_data.json https://hooks.slack.com/services/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
        rm 3_data.json
        """