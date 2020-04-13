import datetime

date = datetime.date.today()

rule uk_annotate:
    input:
        fasta = config["uk_input"],
    output:
        annotate_metadata = config["output_path"] + "/1/coguk_annotated.csv",
        annotate_fasta = config["output_path"] + "/1/coguk_annotated.fasta",
        annotate_log = config["output_path"] + "/1/coguk_annotated.log",
    shell:
        """
        fastafunk annotate \
          --in-fasta {input.fasta} \
          --out-metadata {output.annotate_metadata} \
          --out-fasta {output.annotate_fasta} \
          --log-file {output.annotate_log} \
          --add-cov-id
        """

rule uk_remove_duplicates:
    input:
        fasta = config["uk_input"],
        metadata = rules.uk_annotate.output.annotate_metadata
    output:
        subsample_fasta = config["output_path"] + "/1/coguk_subsample.fasta",
        subsample_log = config["output_path"] + "/1/coguk_subsample.log"
    shell:
        """
        fastafunk subsample \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --group-column cov_id \
          --index-column header \
          --out-fasta {output.subsample_fasta} \
          --sample-size 1 \
          --select-by-min-column gaps &> {output.subsample_log}
        """

rule process_uk_download:
    input:
        uk = rules.uk_remove_duplicates.output.subsample_fasta
    output:
        config["output_path"] + "/1/uk_filtered.fasta"
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
        config["output_path"] + "/1/uk_filtered.covg_length_fitered.fasta"
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
        config["output_path"] + "/1/uk_filtered.covg_length_fitered.mapped.sam"
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
        config["output_path"] + "/1/uk_filtered_fixed_trimmed.fasta"
    shell:
        """
        datafunk sam_2_fasta \
          -s {input.sam} \
          -r {input.reference} \
          -o {output} \
          -t [{params.trim_start}:{params.trim_end}] \
          --prefix_ref \
          --log_inserts
        """

rule uk_filter_low_coverage_sequences:
    input:
        fasta = rules.uk_remove_insertions_and_trim.output
    params:
        min_covg = config["min_covg"]
    output:
        config["output_path"] + "/1/uk_%s_filtered.fasta" %date
    shell:
        """
        datafunk filter_fasta_by_covg_and_length \
          -i {input.fasta} \
          -o {output} \
          --min_covg {params.min_covg}
        """
