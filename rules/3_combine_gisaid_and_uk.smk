import datetime

date = datetime.date.today()

rule combine_gisaid_and_cog:
    input:
        previous_stage = config["output_path"] + "/logs/2_summarize_pangolin_lineage_typing.log",
        gisaid_fasta = config["output_path"] + "/0/gisaid.regularized.fasta",
        gisaid_metadata = config["output_path"] + "/0/gisaid.regularized.csv",
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
        """

rule summarize_combine_gisaid_and_cog:
    input:
        fasta = rules.combine_gisaid_and_cog.output.fasta,
        metadata = rules.combine_gisaid_and_cog.output.metadata,
    params:
        webhook = config["webhook"],
        outdir = config["publish_path"] + "/COG_GISAID",
        prefix = config["publish_path"] + "/COG_GISAID/cog_gisaid_%s" %date,
    log:
        config["output_path"] + "/logs/3_summarize_combine_gisaid_and_cog.log"
    shell:
        """
        mkdir -p {params.outdir}
        cp {input.fasta} {params.prefix}_alignment.fasta
        cp {input.metadata} {params.prefix}_metadata.csv
        echo "> Combined COG and GISAID trimmed alignments published to {params.prefix}_alignment.fasta\\n" >> {log}
        echo "> Combined COG and GISAID restricted metadata published to {params.prefix}_metadata.csv\\n" >> {log}
        echo "> Number of sequences in combined COG and GISAID matched files: $(cat {input.fasta} | grep ">" | wc -l)\\n" &>> {log}

        echo {params.webhook}
        echo '{{"text":"' > 3_data.json
        echo "*Step 3: Combine COG-UK and GISAID data complete*\\n" >> 3_data.json
        cat {log} >> 3_data.json
        echo '"}}' >> 3_data.json
        echo "webhook {params.webhook}"
        curl -X POST -H "Content-type: application/json" -d @3_data.json {params.webhook}
        #rm 3_data.json
        """