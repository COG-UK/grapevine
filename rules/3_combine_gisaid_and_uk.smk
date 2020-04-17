import datetime

date = datetime.date.today()

rule combine_gisaid_and_cog:
    input:
        gisaid_fasta = rules.gisaid_summarize_preprocess.output.fasta,
        gisaid_metadata = rules.gisaid_summarize_preprocess.output.fasta,
        uk_fasta = rules.uk_combine_previous_and_new.output.fasta,
        uk_metadata = rules.uk_combine_previous_and_new.output.metadata
    params:
        prefix = config["output_path"] + "/COG_GISAID/cog_gisaid_%s_" %date,
    output:
        fasta = config["output_path"] + "/3/cog_gisaid_%s.fasta" %date,
        metadata = config["output_path"] + "/3/cog_gisaid_%s.csv" %date
    log:
        config["output_path"] + "/logs/3_combine_gisaid_and_cog.log"
    shell:
        """
        fastafunk merge \
          --in-fasta {input.gisaid_} {input.uk_fasta} \
          --in-metadata {input.gisaid_metadata} {input.uk_metadata} \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --log-file {log}

        num_seqs=$(cat {output.fasta} | grep ">" | wc -l)
        cp {output.fasta} {params.prefix}$num_seqs.fasta
        cp {output.metadata} {params.prefix}$num_seqs.csv
        """