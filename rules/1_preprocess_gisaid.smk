rule process_gisaid_download:
    input:
        gisaid = config["gisaid_input"],
        omitted = config["omitted_file"]
    output:
        temp(config["output_path"] + "/gisaid_filtered.fasta")
    shell:
        """
        datafunk process_gisaid_sequence_data \
          -i \"{input.gisaid}\" \
          -o {output} \
          -e \"{input.omitted}\" \
          --exclude_uk
        """

rule filter_short_sequences:
    input:
        fasta = rules.process_gisaid_download.output
    params:
        min_covg = config["min_covg"],
        min_length = config["min_length"]
    output:
        temp(config["output_path"] + "/gisaid_filtered.covg_length_fitered.fasta")
    shell:
        """
        datafunk filter_fasta_by_covg_and_length \
          -i {input.fasta} \
          -o {output} \
          --min_length {params.min_length}
        """

rule minimap2_to_reference:
    input:
        fasta = rules.filter_short_sequences.output,
        reference = config["reference_fasta"]
    output:
        temp(config["output_path"] + "/gisaid_filtered.covg_length_fitered.mapped.sam")
    shell:
        """
        minimap2 -a -x asm5 {input.reference} {input.fasta} > {output}
        """

rule remove_insertions_and_trim:
    input:
        sam = rules.minimap2_to_reference.output,
        reference = config["reference_fasta"]
    params:
        trim_start = config["trim_start"],
        trim_end = config["trim_end"],
    output:
        temp(config["output_path"] + "/filtered_fixed_trimmed_gisaid.fasta")
    shell:
        """
        datafunk sam_2_fasta \
          -s {input.sam} \
          -r {input.reference} \
          -o {output} \
          -t [{params.trim_start}:{params.trim_end}] \
          --prefix_ref
        """

rule filter_low_coverage_sequences:
    input:
        fasta = rules.remove_insertions_and_trim.output
    params:
        min_covg = config["min_covg"]
    output:
        config["output_path"] + "/filtered_fixed_gisaid.fasta"
    shell:
        """
        datafunk filter_fasta_by_covg_and_length \
          -i {input.fasta} \
          -o {output} \
          --min_covg {params.min_covg}
        """