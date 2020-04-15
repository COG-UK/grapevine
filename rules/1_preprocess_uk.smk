configfile: workflow.current_basedir + "/config.yaml"

import datetime

date = datetime.date.today()

rule uk_unify_headers:
    input:
        fasta = config["latest_uk_fasta"],
        metadata = config["latest_uk_metadata"]
    output:
        fasta = temp(config["output_path"] + "/1/uk_latest.unify_headers.fasta"),
        metadata = temp(config["output_path"] + "/1/uk_latest.unify_headers.csv")
    log:
        config["output_path"] + "/logs/1_uk_unify_headers.log"
    shell:
        """
        datafunk set_uniform_header \
          --input_fasta {input.fasta} \
          --input_metadata {input.metadata} \
          --output_fasta {output.fasta} \
          --output_metadata {output.metadata} \
          --log {log} \
          --cog_uk
        """

rule uk_annotate_to_remove_duplicates:
    input:
        fasta = rules.uk_unify_headers.output.fasta
        metadata = rules.uk_unify_headers.output.metadata
    output:
        annotate_metadata = config["output_path"] + "/1/uk_%s_filtered.csv" %date,
    log:
        config["output_path"] + "/logs/1_uk_annotate_to_remove_duplicates.log"
    shell:
        """
        fastafunk annotate \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata}
          --out-metadata {output.annotate_metadata} \
          --log-file {log} \
          --add-cov-id
        """

rule uk_remove_duplicates:
    input:
        fasta = rules.uk_unify_headers.output.fasta,
        metadata = rules.uk_annotate_to_remove_duplicates.output.annotate_metadata
    output:
        temp(subsample_fasta = config["output_path"] + "/1/coguk.deduplicated.fasta")
    log:
        config["output_path"] + "/logs/1_uk_remove_duplicates.log"
    shell:
        """
        fastafunk subsample \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --group-column cov_id \
          --index-column header \
          --out-fasta {output.subsample_fasta} \
          --sample-size 1 \
          --select-by-min-column gaps &> {log}
        """

rule uk_filter_short_sequences:
    input:
        fasta = rules.uk_remove_duplicates.output
    params:
        min_covg = config["min_covg"],
        min_length = config["min_length"]
    output:
        temp(config["output_path"] + "/1/coguk.deduplicated.fixed_headers.covg_length_fitered.fasta")
    log:
        config["output_path"] + "/logs/1_uk_filter_short_sequences.log"
    shell:
        """
        datafunk filter_fasta_by_covg_and_length \
          -i {input.fasta} \
          -o {output} \
          --min_length {params.min_length} &> {log}
        """

rule uk_minimap2_to_reference:
    input:
        fasta = rules.uk_filter_short_sequences.output,
        reference = config["reference_fasta"]
    output:
        temp(config["output_path"] + "/1/coguk.deduplicated.fixed_headers.length_fitered.mapped.sam")
    log:
        config["output_path"] + "/logs/1_uk_minimap2_to_reference.log"
    shell:
        """
        minimap2 -a -x asm5 {input.reference} {input.fasta} > {output} 2> {log}
        """

rule uk_remove_insertions_and_trim:
    input:
        sam = rules.uk_minimap2_to_reference.output,
        reference = config["reference_fasta"]
    params:
        trim_start = config["trim_start"],
        trim_end = config["trim_end"],
    output:
        temp(config["output_path"] + "/1/coguk.deduplicated.fixed_headers.length_fitered.trimmed.fasta")
    log:
        config["output_path"] + "/logs/1_uk_remove_insertions_and_trim.log"
    shell:
        """
        datafunk sam_2_fasta \
          -s {input.sam} \
          -r {input.reference} \
          -o {output} \
          -t [{params.trim_start}:{params.trim_end}] \
          --log_inserts &> {log}
        """

rule uk_filter_low_coverage_sequences:
    input:
        fasta = rules.uk_remove_insertions_and_trim.output
    params:
        min_covg = config["min_covg"]
    output:
        config["output_path"] + "/1/uk_%s_filtered.fasta" %date
    log:
        config["output_path"] + "/logs/1_uk_filter_low_coverage_sequences.log"
    shell:
        """
        datafunk filter_fasta_by_covg_and_length \
          -i {input.fasta} \
          -o {output} \
          --min_covg {params.min_covg} &> {log}
        """

rule uk_summarize_preprocess:
    input:
        raw_fasta = config["latest_uk_fasta"]
        deduplicated_fasta = rules.uk_remove_duplicates.output.subsample_fasta,
        removed_short_fasta = rules.uk_filter_short_sequences.output,
        removed_low_covg_fasta = rules.uk_filter_low_coverage_sequences.output
    log:
        config["output_path"] + "/logs/1_summary_preprocess_uk.log"
    shell:
        """
        echo "Number of sequences in raw UK fasta: $(cat {input.raw_fasta} | grep ">" | wc -l)"
        echo "Number of sequences after deduplication: $(cat {input.deduplicated_fasta} | grep ">" | wc -l)"
        echo "Number of sequences after removing sequences <29000bps: $(cat {input.removed_short_fasta} | grep ">" | wc -l)"
        echo "Number of sequences after trimming and removing those with <95% coverage: $(cat {input.removed_low_covg_fasta} | grep ">" | wc -l)"
        """
