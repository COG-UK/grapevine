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
        fasta = config["previous_gisaid_fasta"],
        metadata = config["previous_gisaid_metadata"],
        omitted = config["previous_omitted_file"],
    params:
        flags = config["gisaid_flags"]
    output:
        fasta = temp(config["output_path"] + "/0/gisaid_latest.fasta"),
        metadata = temp(config["output_path"] + "/0/gisaid_latest.csv")
    log:
        config["output_path"] + "/logs/0_gisaid_process_json.log"
    shell:
        """
        datafunk process_gisaid_sequence_data \
          --input \"{input.json}\" \
          --output {output.fasta} \
          --exclude \"{input.omitted}\" \
          {params.flags}

        datafunk gisaid_json_2_metadata \
          --new \"{input.json}\" \
          --csv \"{input.metadata}\" \
          --output {output.metadata} \
          --exclude \"{input.omitted}\" \
        """

rule gisaid_unify_headers:
    input:
        fasta = rules.gisaid_process_json.output.fasta,
        metadata = rules.gisaid_process_json.output.metadata
    output:
        fasta = temp(config["output_path"] + "/0/gisaid_latest.unify_headers.fasta"),
        metadata = temp(config["output_path"] + "/0/gisaid_latest.unify_headers.csv")
    log:
        config["output_path"] + "/logs/0_gisaid_unify_headers.log"
    shell:
        """
        datafunk set_uniform_header \
          --input_fasta {input.fasta} \
          --input_metadata {input.metadata} \
          --output_fasta {output.fasta} \
          --output_metadata {output.metadata} \
          --log {log} \
          --gisaid
        """

rule gisaid_extract_new:
    input:
        fasta = rules.gisaid_unify_headers.output.fasta,
        metadata = rules.gisaid_unify_headers.output.metadata,
        previous_metadata = config["previous_gisaid_metadata"]
    output:
        fasta = temp(config["output_path"] + "/0/gisaid_latest.unify_headers.new.fasta"),
        metadata = temp(config["output_path"] + "/0/gisaid_latest.unify_headers.new.csv"),
    log:
        config["output_path"] + "/logs/0_gisaid_extract_new.log"
    shell:
        """
        fastafunk new \
          --input-fasta {input.fasta} \
          --input-metadata {input.previous_metadata} {input.metadata} \
          --output-fasta {output.fasta} \
          --output-metadata {output.metadata} \
          --date-column 'edin_date_stamp' \
          --log {log}
        """

rule gisaid_filter_short_sequences:
    input:
        fasta = rules.gisaid_extract_new.output.fasta
    params:
        min_length = config["min_length"]
    output:
        temp(config["output_path"] + "/0/gisaid_latest.unify_headers.new.length_fitered.fasta")
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
        temp(config["output_path"] + "/0/gisaid_latest.unify_headers.new.length_fitered.mapped.sam")
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
        temp(config["output_path"] + "/0/gisaid_latest.unify_headers.new.length_fitered.trimmed.fasta")
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
        temp(config["output_path"] + "/0/gisaid_latest.unify_headers.new.length_fitered.trimmed.low_covg_filtered.fasta")
    log:
        config["output_path"] + "/logs/0_gisaid_filter_low_coverage_sequences.log"
    shell:
        """
        datafunk filter_fasta_by_covg_and_length \
          -i {input.fasta} \
          -o {output} \
          --min_covg {params.min_covg}
        """

rule gisaid_pangolin:
    input:
        fasta = rules.gisaid_filter_low_coverage_sequences.output
    params:
        outdir = config["output_path"] + "/0/pangolin"
    output:
        config["output_path"] + "/0/pangolin/lineage_report.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_pangolin.log"
    threads: 8
    shell:
        """
        pangolin {input.fasta}
        --threads {threads} \
        --outdir {params.outdir} \
        """
# rule gisaid_add_pangolin_lineages_to_metadata:
#
# rule gisaid_combine_previous_and_new:
#     input:
#         previous_fasta =
#         previous_metadata =
#         new_fasta =
#         new_metadata =
#     output:
#         fasta =
#         metadata =
#     shell:
#         """
#         fastafunk merge
#         """
#
# rule gisaid_omit_sequences:
#     input:
#         fasta =
#         metadata =
#         omitted =
#     output:
#         fasta =
#         metadata =
#     shell:
#         """
#         fastafunk remove
#         """

rule gisaid_summarize_preprocess:
    input:
        previous_fasta = config["previous_gisaid_fasta"],
        latest_fasta = rules.gisaid_process_json.output.fasta,
        unify_headers_fasta = rules.gisaid_unify_headers.output.fasta,
        removed_short_fasta = rules.gisaid_filter_short_sequences.output,
        removed_low_covg_fasta = rules.gisaid_filter_low_coverage_sequences.output,
        final_fasta = rules.gisaid_filter_low_coverage_sequences.output,
        final_metadata = rules.gisaid_unify_headers.metadata
    params:
        prefix = config["output_path"] + "/gisaid_%s_" %date
    output:
        fasta = config["output_path"] + "/gisaid_%s_*.fasta" %date,
        metadata = config["output_path"] + "/gisaid_%s_*.csv" %date
    log:
        config["output_path"] + "/logs/0_summary_preprocess_gisaid.log"
    shell:
        """
        echo "Number of sequences in previous GISAID fasta: $(cat {input.previous_fasta} | grep ">" | wc -l)"
        echo "Number of sequences in latest GISAID download: $(cat {input.latest_fasta} | grep ">" | wc -l)"
        echo "Number of sequences after matching headers: $(cat {input.unify_headers_fasta} | grep ">" | wc -l)"
        echo "Number of sequences after removing sequences <29000bps: $(cat {input.removed_short_fasta} | grep ">" | wc -l)"
        echo "Number of sequences after trimming and removing those with <95% coverage: $(cat {input.removed_low_covg_fasta} | grep ">" | wc -l)"

        num_seqs=$(cat {input.final_fasta} | grep ">" | wc -l)
        cp {input.final_fasta} {params.prefix}$num_seqs.fasta
        cp {input.final_metadata} {params.prefix}$num_seqs.csv
        """
