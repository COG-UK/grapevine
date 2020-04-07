import datetime

date = datetime.date.today()

rule process_uk_download:
    input:
        uk = config["uk_input"],
    output:
        config["output_path"] + "/uk_filtered.fasta"
    shell:
        """
        datafunk process_gisaid_sequence_data \
          -i \"{input.uk}\" \
          -o {output}
        """

rule uk_filter_short_sequences:
    input:
        fasta = rules.process_uk_download.output
    params:
        min_covg = config["min_covg"],
        min_length = config["min_length"]
    output:
        config["output_path"] + "/uk_filtered.covg_length_fitered.fasta"
    shell:
        """
        datafunk filter_fasta_by_covg_and_length \
          -i {input.fasta} \
          -o {output} \
          --min_length {params.min_length}
        """

rule uk_minimap2_to_reference:
    input:
        fasta = rules.uk_filter_short_sequences.output,
        reference = config["reference_fasta"]
    output:
        config["output_path"] + "/uk_filtered.covg_length_fitered.mapped.sam"
    shell:
        """
        minimap2 -a -x asm5 {input.reference} {input.fasta} > {output}
        """

rule uk_remove_insertions_and_trim:
    input:
        sam = rules.uk_minimap2_to_reference.output,
        reference = config["reference_fasta"]
    params:
        trim_start = config["trim_start"],
        trim_end = config["trim_end"],
    output:
        config["output_path"] + "/uk_filtered_fixed_trimmed.fasta"
    shell:
        """
        datafunk sam_2_fasta \
          -s {input.sam} \
          -r {input.reference} \
          -o {output} \
          -t [{params.trim_start}:{params.trim_end}] \
          --prefix_ref
        """

rule uk_filter_low_coverage_sequences:
    input:
        fasta = rules.uk_remove_insertions_and_trim.output
    params:
        min_covg = config["min_covg"]
    output:
        config["output_path"] + "/uk_%s_filtered.fasta" %date
    shell:
        """
        datafunk filter_fasta_by_covg_and_length \
          -i {input.fasta} \
          -o {output} \
          --min_covg {params.min_covg}
        """