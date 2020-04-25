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
          --input-fasta {input.fasta} \
          --input-metadata {input.metadata} \
          --output-fasta {output.fasta} \
          --output-metadata {output.metadata} \
          --log {log} \
          --cog-uk

        sed --in-place=.tmp 's/United Kingdom/UK/g' {output.metadata}
        """

rule uk_add_epi_week:
    input:
        metadata = rules.uk_unify_headers.output.metadata
    output:
        metadata = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.csv"
    log:
        config["output_path"] + "/logs/1_uk_add_epi_week.log"
    shell:
        """
        datafunk add_epi_week \
        --input-metadata {input.metadata} \
        --output-metadata {output.metadata} \
        --date-column collection_date \
        --epi-column-name edin_epi_week &> {log}
        """

rule uk_annotate_to_remove_duplicates:
    input:
        fasta = rules.uk_unify_headers.output.fasta,
        metadata = rules.uk_add_epi_week.output.metadata
    output:
        metadata = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.annotated.csv"
    log:
        config["output_path"] + "/logs/1_uk_annotate_to_remove_duplicates.log"
    shell:
        """
        fastafunk annotate \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --out-metadata {output.metadata} \
          --log-file {log} \
          --add-cov-id
        """

rule uk_remove_duplicates:
    input:
        fasta = rules.uk_unify_headers.output.fasta,
        metadata = rules.uk_annotate_to_remove_duplicates.output.metadata
    output:
        fasta = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.fasta",
        metadata = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.csv"
    log:
        config["output_path"] + "/logs/1_uk_remove_duplicates.log"
    shell:
        """
        fastafunk subsample \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --group-column cov_id \
          --index-column sequence_name \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --sample-size 1 \
          --select-by-min-column gaps &> {log}
        """

rule uk_filter_short_sequences:
    input:
        fasta = rules.uk_remove_duplicates.output.fasta
    params:
        min_covg = config["min_covg"],
        min_length = config["min_length"]
    output:
        fasta = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.length_fitered.fasta"
    log:
        config["output_path"] + "/logs/1_uk_filter_short_sequences.log"
    shell:
        """
        datafunk filter_fasta_by_covg_and_length \
          -i {input.fasta} \
          -o {output.fasta} \
          --min-length {params.min_length} &> {log}
        """

rule uk_minimap2_to_reference:
    input:
        fasta = rules.uk_filter_short_sequences.output,
        reference = config["reference_fasta"]
    output:
        sam = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.length_fitered.mapped.sam"
    log:
        config["output_path"] + "/logs/1_uk_minimap2_to_reference.log"
    shell:
        """
        minimap2 -a -x asm5 {input.reference} {input.fasta} > {output.sam} 2> {log}
        """

rule uk_remove_insertions_and_trim_and_pad:
    input:
        sam = rules.uk_minimap2_to_reference.output.sam,
        reference = config["reference_fasta"]
    params:
        trim_start = config["trim_start"],
        trim_end = config["trim_end"],
        insertions = config["output_path"] + "/1/uk_insertions.txt"
    output:
        fasta = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.length_fitered.trimmed_alignment.fasta"
    log:
        config["output_path"] + "/logs/1_uk_remove_insertions_and_trim_and_pad.log"
    shell:
        """
        datafunk sam_2_fasta \
          -s {input.sam} \
          -r {input.reference} \
          -o {output.fasta} \
          -t [{params.trim_start}:{params.trim_end}] \
          --pad \
          --log-inserts &> {log}
        mv insertions.txt {params.insertions}
        """

rule uk_filter_low_coverage_sequences:
    input:
        fasta = rules.uk_remove_insertions_and_trim_and_pad.output.fasta
    params:
        min_covg = config["min_covg"]
    output:
        fasta = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.length_fitered.trimmed.low_covg_filtered.fasta"
    log:
        config["output_path"] + "/logs/1_uk_filter_low_coverage_sequences.log"
    shell:
        """
        datafunk filter_fasta_by_covg_and_length \
          -i {input.fasta} \
          -o {output.fasta} \
          --min-covg {params.min_covg} &> {log}
        """

rule uk_full_untrimmed_alignment:
    input:
        sam = rules.uk_minimap2_to_reference.output.sam,
        reference = config["reference_fasta"],
        omit_list = rules.uk_filter_low_coverage_sequences.log
    output:
        fasta = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.length_fitered.full_alignment.fasta"
    log:
        config["output_path"] + "/logs/1_uk_full_untrimmed_alignment.log"
    shell:
        """
        datafunk sam_2_fasta \
          -s {input.sam} \
          -r {input.reference} \
          -o {output.fasta} \
          &> {log}

        fastafunk remove \
          --in-fasta {output.fasta} \
          --in-metadata {input.omit_list} \
          --out-fasta removed.fa
        mv removed.fa {output.fasta}
        """

rule summarize_preprocess_uk:
    input:
        raw_fasta = config["latest_uk_fasta"],
        unify_headers_fasta = rules.uk_unify_headers.output.fasta,
        deduplicated_fasta = rules.uk_remove_duplicates.output.fasta,
        removed_short_fasta = rules.uk_filter_short_sequences.output.fasta,
        removed_low_covg_fasta = rules.uk_filter_low_coverage_sequences.output.fasta,
        full_alignment = rules.uk_full_untrimmed_alignment.output.fasta
    params:
        webhook = config["webhook"],
        outdir = config["publish_path"] + "/COG",
        prefix = config["publish_path"] + "/COG/cog_%s" %date
    log:
        config["output_path"] + "/logs/1_summarize_preprocess_uk.log"
    shell:
        """
        echo "> Number of sequences in raw UK fasta: $(cat {input.raw_fasta} | grep ">" | wc -l)\\n" &> {log}
        echo "> Number of sequences in raw UK fasta after unifying headers: $(cat {input.unify_headers_fasta} | grep ">" | wc -l)\\n" &>> {log}
        echo "> Number of sequences after deduplication: $(cat {input.deduplicated_fasta} | grep ">" | wc -l)\\n" &>> {log}
        echo "> Number of sequences after removing sequences <29000bps: $(cat {input.removed_short_fasta} | grep ">" | wc -l)\\n" &>> {log}
        echo "> Number of sequences after trimming and removing those with <95% coverage: $(cat {input.removed_low_covg_fasta} | grep ">" | wc -l)\\n" &>> {log}

        mkdir -p {params.outdir}
        cp {input.full_alignment} {params.prefix}_alignment.full.fasta
        echo "> Full untrimmed COG alignment published to {params.prefix}_alignment.full.fasta\\n" >> {log}
        cp {input.removed_low_covg_fasta} {params.prefix}_alignment.trimmed.fasta
        echo "> Trimmed COG alignment published to {params.prefix}_alignment.trimmed.fasta\\n" >> {log}

        echo '{{"text":"' > 1_data.json
        echo "*Step 1: COG-UK preprocessing complete*\\n" >> 1_data.json
        cat {log} >> 1_data.json
        echo '"}}' >> 1_data.json
        curl -X POST -H "Content-type: application/json" -d @1_data.json {params.webhook}
        rm 1_data.json
        """
