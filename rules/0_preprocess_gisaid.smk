import datetime

date = datetime.date.today()

#rule gisaid_download_new_json:
#    input:
#    output:
#        json = config["output_path"] + "/0/gisaid_dump.json"
#    log:
#        config["output_path"] + "/logs/0_gisaid_download_new_json.log"
#    shell:
#        """
#        """

rule gisaid_process_json:
    input:
        #json = rules.gisaid_download_new_json.output.json
        json = config["latest_gisaid_json"],
        metadata = config["previous_gisaid_metadata"],
        omitted = config["previous_omitted_file"],
    params:
        flags = config["gisaid_flags"]
    output:
        fasta = config["output_path"] + "/0/gisaid_latest.fasta",
        metadata = config["output_path"] + "/0/gisaid_latest.csv"
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
        fasta = config["output_path"] + "/0/gisaid_latest.unify_headers.fasta",
        metadata = config["output_path"] + "/0/gisaid_latest.unify_headers.csv"
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
        fasta = config["output_path"] + "/0/gisaid_latest.unify_headers.new.fasta",
        metadata = config["output_path"] + "/0/gisaid_latest.unify_headers.new.csv",
    log:
        config["output_path"] + "/logs/0_gisaid_extract_new.log"
    shell:
        """
        fastafunk new \
          --in-fasta {input.fasta} \
          --in-metadata {input.previous_metadata} {input.metadata} \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --date-column 'edin_date_stamp' \
          --log {log}
        """

rule gisaid_filter_short_sequences:
    input:
        fasta = rules.gisaid_extract_new.output.fasta
    params:
        min_length = config["min_length"]
    output:
        fasta = config["output_path"] + "/0/gisaid_latest.unify_headers.new.length_fitered.fasta"
    log:
        config["output_path"] + "/logs/0_gisaid_filter_short_sequences.log"
    shell:
        """
        datafunk filter_fasta_by_covg_and_length \
          -i {input.fasta} \
          -o {output} \
          --min_length {params.min_length} &> {log}
        """

rule gisaid_minimap2_to_reference:
    input:
        fasta = rules.gisaid_filter_short_sequences.output.fasta,
        reference = config["reference_fasta"]
    output:
        sam = config["output_path"] + "/0/gisaid_latest.unify_headers.new.length_fitered.mapped.sam"
    log:
        config["output_path"] + "/logs/0_gisaid_minimap2_to_reference.log"
    shell:
        """
        minimap2 -a -x asm5 {input.reference} {input.fasta} -o {output} &> {log}
        """

rule gisaid_remove_insertions_and_trim:
    input:
        sam = rules.gisaid_minimap2_to_reference.output.sam,
        reference = config["reference_fasta"]
    params:
        trim_start = config["trim_start"],
        trim_end = config["trim_end"],
    output:
        fasta = config["output_path"] + "/0/gisaid_latest.unify_headers.new.length_fitered.trimmed.fasta"
    log:
        config["output_path"] + "/logs/0_gisaid_remove_insertions_and_trim.log"
    shell:
        """
        datafunk sam_2_fasta \
          -s {input.sam} \
          -r {input.reference} \
          -o {output} \
          -t [{params.trim_start}:{params.trim_end}] \
        """

rule gisaid_filter_low_coverage_sequences:
    input:
        fasta = rules.gisaid_remove_insertions_and_trim.output.fasta
    params:
        min_covg = config["min_covg"]
    output:
        fasta = config["output_path"] + "/0/gisaid_latest.unify_headers.new.length_fitered.trimmed.low_covg_filtered.fasta"
    log:
        config["output_path"] + "/logs/0_gisaid_filter_low_coverage_sequences.log"
    shell:
        """
        datafunk filter_fasta_by_covg_and_length \
          -i {input.fasta} \
          -o {output} \
          --min_covg {params.min_covg} &> {log}
        """

rule gisaid_pangolin:
    input:
        fasta = rules.gisaid_filter_low_coverage_sequences.output.fasta
    params:
        outdir = config["output_path"] + "/0/pangolin"
    output:
        lineages = config["output_path"] + "/0/pangolin/lineage_report.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_pangolin.log"
    threads: 8
    shell:
        """
        pangolin {input.fasta} \
        --threads {threads} \
        --outdir {params.outdir} > {log} 2>&1
        """

rule gisaid_add_pangolin_lineages_to_metadata:
    input:
        metadata = rules.gisaid_extract_new.output.metadata,
        lineages = rules.gisaid_pangolin.output.lineages
    output:
        metadata = config["output_path"] + "/0/gisaid_latest.unify_headers.new.lineages.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_add_pangolin_lineages_to_metadata.log"
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.lineages} \
          --index-column header \
          --new-columns lineage bootstrap \
          --out-metadata {output.metadata} &> {log}
        """

rule gisaid_add_epi_week:
    input:
        metadata = rules.gisaid_add_pangolin_lineages_to_metadata.output.metadata
    output:
        metadata = config["output_path"] + "/0/gisaid_latest.unify_headers.new.lineages.epi_week.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_add_epi_week.log"
    shell:
        """
        datafunk add_epi_week \
        --input_metadata {input.metadata} \
        --output_metadata {output.metadata} \
        --date_column covv_collection_date \
        --epi_column_name edin_epi_week &> {log}
        """

rule gisaid_combine_previous_and_new:
    input:
        previous_fasta = config["previous_gisaid_fasta"],
        previous_metadata = config["previous_gisaid_metadata"],
        new_fasta = rules.gisaid_filter_low_coverage_sequences.output.fasta,
        new_metadata = rules.gisaid_add_epi_week.output.metadata
    output:
        fasta = config["output_path"] + "/0/gisaid_latest.unify_headers.combined.fasta",
        metadata = config["output_path"] + "/0/gisaid_latest.unify_headers.combined.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_combine_previous_and_new.log"
    shell:
        """
        fastafunk merge \
          --in-fasta {input.previous_fasta} {input.new_fasta} \
          --in-metadata {input.previous_metadata} {input.new_metadata} \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --log-file {log}
        """
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

rule gisaid_counts_by_country:
    input:
        metadata = rules.gisaid_unify_headers.output.metadata
    output:
        config["output_path"] + "/0/gisaid_counts_by_country.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_counts_by_country.log"
    shell:
        """
        fastafunk count \
          --in-metadata {input.metadata} \
          --group-column edin_admin_0 \
          --log-file {output} &> {log}
        """

rule gisaid_summarize_preprocess:
    input:
        previous_fasta = config["previous_gisaid_fasta"],
        latest_fasta = rules.gisaid_process_json.output.fasta,
        unify_headers_fasta = rules.gisaid_unify_headers.output.fasta,
        new_fasta = rules.gisaid_extract_new.output.fasta,
        removed_short_fasta = rules.gisaid_filter_short_sequences.output,
        removed_low_covg_fasta = rules.gisaid_filter_low_coverage_sequences.output,
        final_fasta = rules.gisaid_filter_low_coverage_sequences.output,
        final_metadata = rules.gisaid_unify_headers.output.metadata
    params:
        prefix = config["output_path"] + "/GISAID/gisaid_%s_" %date,
        prefix_hack = config["output_path"] + "/0/gisaid_%s" %date
    output:
        fasta = config["output_path"] + "/0/gisaid_%s.fasta" %date,
        metadata = config["output_path"] + "/0/gisaid_%s.csv" %date
    log:
        config["output_path"] + "/logs/0_summary_preprocess_gisaid.log"
    shell:
        """
        echo "Number of sequences in previous GISAID fasta: $(cat {input.previous_fasta} | grep ">" | wc -l)" >> {log}
        echo "Number of sequences in latest GISAID download: $(cat {input.latest_fasta} | grep ">" | wc -l)" >> {log}
        echo "Number of sequences after matching headers: $(cat {input.unify_headers_fasta} | grep ">" | wc -l)" >> {log}
        echo "Number of new sequences: $(cat {input.new_fasta} | grep ">" | wc -l)" >> {log}
        echo "Number of sequences after removing sequences <29000bps: $(cat {input.removed_short_fasta} | grep ">" | wc -l)" >> {log}
        echo "Number of sequences after trimming and removing those with <95% coverage: $(cat {input.removed_low_covg_fasta} | grep ">" | wc -l)" >> {log}

        num_seqs=$(cat {input.final_fasta} | grep ">" | wc -l)
        cp {input.final_fasta} {params.prefix}$num_seqs.fasta
        cp {input.final_metadata} {params.prefix}$num_seqs.csv
        cp {input.final_fasta} {params.prefix_hack}.fasta
        cp {input.final_metadata} {params.prefix_hack}.csv
        """
