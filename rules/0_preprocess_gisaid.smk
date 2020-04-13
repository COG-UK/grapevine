import datetime

date = datetime.date.today()

#rule gisaid_download_new_json:
#    input:
#    output:
#        json = temp(config["output_path"] + "/0/gisaid_dump.json")
#    log:
#        config["output_path"] + "/logs/0_gisaid_download_new_json.log"
#    shell:
#        """
#        """

rule gisaid_process_json:
    input:
        #json = rules.gisaid_download_new_json.output.json
        json = config["latest_gisaid_json"],
        gisaid_fasta = config["previous_gisaid_fasta"],
        gisaid_metadata = config["previous_gisaid_metadata"],
        gisaid_omitted = config["previous_omitted_file"],
    params:
        flags = config["gisaid_flags"]
    output:
        fasta = temp(config["output_path"] + "/0/gisaid_latest.fasta"),
        metadata = config["output_path"] + "/0/gisaid_%s_filtered.csv" %date
    log:
        config["output_path"] + "/logs/0_gisaid_process_json.log"
    shell:
        """
        datafunk process_gisaid_sequence_data \
          --input \"{input.json}\" \
          --output {output.fasta} \
          --exclude \"{input.gisaid_omitted}\" \
          {params.flags}

        datafunk gisaid_json_2_metadata \
          --new \"{input.json}\" \
          --csv \"{input.gisaid_metadata}\" \
          --output {output.metadata} \
          --exclude \"{input.gisaid_omitted}\" \
        """

rule gisaid_filter_short_sequences:
    input:
        fasta = rules.gisaid_process_json.output
    params:
        min_length = config["min_length"]
    output:
        temp(config["output_path"] + "/0/gisaid_latest.covg_length_fitered.fasta")
    log:
        config["output_path"] + "/logs/0_gisaid_filter_short_sequences.log"
    shell:
        """
        datafunk filter_fasta_by_covg_and_length \
          -i {input.fasta} \
          -o {output} \
          --min_length {params.min_length}
        """

rule gisaid_minimap2_to_reference:
    input:
        fasta = rules.gisaid_filter_short_sequences.output,
        reference = config["reference_fasta"]
    output:
        temp(config["output_path"] + "/0/gisaid_latest.covg_length_fitered.mapped.sam")
    log:
        config["output_path"] + "/logs/0_gisaid_minimap2_to_reference.log"
    shell:
        """
        minimap2 -a -x asm5 {input.reference} {input.fasta} > {output}
        """

rule gisaid_remove_insertions_and_trim:
    input:
        sam = rules.gisaid_minimap2_to_reference.output,
        reference = config["reference_fasta"]
    params:
        trim_start = config["trim_start"],
        trim_end = config["trim_end"],
    output:
        temp(config["output_path"] + "/0/gisaid_latest.covg_length_fitered.trimmed.fasta")
    log:
        config["output_path"] + "/logs/0_gisaid_remove_insertions_and_trim.log"
    shell:
        """
        datafunk sam_2_fasta \
          -s {input.sam} \
          -r {input.reference} \
          -o {output} \
          -t [{params.trim_start}:{params.trim_end}] \
          --prefix_ref
        """

rule gisaid_filter_low_coverage_sequences:
    input:
        fasta = rules.gisaid_remove_insertions_and_trim.output
    params:
        min_covg = config["min_covg"]
    output:
        config["output_path"] + "/0/gisaid_%s_filtered.fasta" %date
    log:
        config["output_path"] + "/logs/0_gisaid_filter_low_coverage_sequences.log"
    shell:
        """
        datafunk filter_fasta_by_covg_and_length \
          -i {input.fasta} \
          -o {output} \
          --min_covg {params.min_covg}
        """

rule gisaid_summarize_preprocess:
    input:
        previous_fasta = config["previous_gisaid_fasta"],
        latest_fasta = rules.gisaid_process_json.output.fasta,
        removed_short_fasta = rules.gisaid_filter_short_sequences.output,
        removed_low_covg_fasta = rules.gisaid_filter_low_coverage_sequences.output
    log:
        config["output_path"] + "/logs/0_summary_preprocess_gisaid.log"
    shell:
        """
        echo "Number of sequences in previous GISAID fasta: $(cat {input.previous_fasta} | grep ">" | wc -l)"
        echo "Number of sequences in latest GISAID download: $(cat {input.latest_fasta} | grep ">" | wc -l)"
        echo "Number of sequences after removing sequences <29000bps: $(cat {input.removed_short_fasta} | grep ">" | wc -l)"
        echo "Number of sequences after trimming and removing those with <95% coverage: $(cat {input.removed_low_covg_fasta} | grep ">" | wc -l)"
        """
