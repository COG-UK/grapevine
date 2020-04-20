import datetime

date = datetime.date.today()

rule combine_gisaid_and_cog:
    input:
        summary_gisaid = rules.gisaid_summarize_preprocess.log,
        summary_uk = rules.uk_summarize_preprocess.log,
        gisaid_fasta = rules.gisaid_combine_previous_and_new.output.fasta,
        gisaid_metadata = rules.gisaid_combine_previous_and_new.output.metadata,
        uk_fasta = rules.uk_combine_previous_and_new.output.fasta,
        uk_metadata = rules.uk_combine_previous_and_new.output.metadata
    params:
        prefix = "COG_GISAID/cog_gisaid_%s_" %date,
    output:
        fasta = config["output_path"] + "/3/cog_gisaid_%s.fasta" %date,
        metadata = config["output_path"] + "/3/cog_gisaid_%s.csv" %date
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

        mkdir -p COG_GISAID
        num_seqs=$(cat {output.fasta} | grep ">" | wc -l)
        cp {output.fasta} {params.prefix}$num_seqs.fasta
        cp {output.metadata} {params.prefix}$num_seqs.csv
        """