

rule uk_add_header_column:
    input:
        fasta = config["latest_uk_fasta"],
        metadata = config["latest_uk_metadata"]
    output:
        fasta = config["output_path"] + "/1/uk_latest.add_header.fasta",
        metadata = config["output_path"] + "/1/uk_latest.add_header.csv"
    log:
        config["output_path"] + "/logs/1_add_header_column.log"
    shell:
        """
        datafunk add_header_column \
        --input-fasta {input.fasta} \
        --input-metadata {input.metadata} \
        --output-metadata {output.metadata} \
        --output-fasta {output.fasta} \
        --cog-uk &> {log}
        """


rule uk_annotate_to_remove_duplicates:
    input:
        fasta = rules.uk_add_header_column.output.fasta,
        metadata = rules.uk_add_header_column.output.metadata
    output:
        metadata = config["output_path"] + "/1/uk_latest.add_header.annotated.csv"
    log:
        config["output_path"] + "/logs/1_uk_annotate_to_remove_duplicates.log"
    shell:
        """
        fastafunk annotate \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --out-metadata {output.metadata} \
          --log-file {log} \
          --add-cov-id \
          --index-column header &> {log}
        """


rule uk_remove_duplicates_covid_by_gaps:
    input:
        fasta = rules.uk_add_header_column.output.fasta,
        metadata = rules.uk_annotate_to_remove_duplicates.output.metadata
    output:
        fasta = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated_cov_id.fasta",
        metadata = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated_cov_id.csv"
    log:
        config["output_path"] + "/logs/1_uk_remove_duplicates.log"
    shell:
        """
        fastafunk subsample \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --group-column cov_id \
          --index-column header \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --sample-size 1 \
          --select-by-min-column gaps &> {log}
        """


rule uk_remove_duplicates_biosamplesourceid_by_date:
    input:
        fasta = rules.uk_remove_duplicates_covid_by_gaps.output.fasta,
        metadata = rules.uk_remove_duplicates_covid_by_gaps.output.metadata
    output:
        fasta = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated_cov_id_biosample_source_id.fasta",
        metadata = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated_cov_id_biosample_source_id.csv"
    log:
        config["output_path"] + "/logs/1_uk_remove_duplicates.log"
    shell:
        """
        fastafunk subsample \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --group-column biosample_source_id \
          --index-column header \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --sample-size 1 \
          --select-by-min-column collection_date &> {log}
        """


rule uk_unify_headers:
    input:
        fasta = rules.uk_remove_duplicates_biosamplesourceid_by_date.output.fasta,
        metadata = rules.uk_remove_duplicates_biosamplesourceid_by_date.output.metadata
    output:
        fasta = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated.unify_headers.fasta",
        metadata = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated.unify_headers.csv"
    log:
        config["output_path"] + "/logs/1_uk_unify_headers.log"
    shell:
        """
        datafunk set_uniform_header \
          --input-fasta {input.fasta} \
          --input-metadata {input.metadata} \
          --output-fasta {output.fasta} \
          --output-metadata {output.metadata} \
          --cog-uk  &> {log}

        sed --in-place=.tmp 's/United Kingdom/UK/g' {output.metadata}
        """


rule uk_add_epi_week:
    input:
        metadata = rules.uk_unify_headers.output.metadata
    output:
        metadata = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.csv",
        tmp_metadata = temp(config["output_path"] + "/1/uk_latest.unify_headers.epi_week.csv.tmp")
    log:
        config["output_path"] + "/logs/1_uk_add_epi_week.log"
    shell:
        """
        datafunk add_epi_week \
        --input-metadata {input.metadata} \
        --output-metadata {output.tmp_metadata} \
        --date-column received_date \
        --epi-column-name edin_epi_week &> {log}

        datafunk add_epi_week \
        --input-metadata {output.tmp_metadata} \
        --output-metadata {output.metadata} \
        --date-column collection_date \
        --epi-column-name edin_epi_week &>> {log}
        """


rule uk_minimap2_to_reference:
    input:
        fasta = rules.uk_unify_headers.output.fasta,
        reference = config["reference_fasta"]
    output:
        sam = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.mapped.sam"
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
        fasta = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.alignment.trimmed.fasta"
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
        fasta = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.trimmed.low_covg_filtered.fasta"
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
        fasta = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.alignment.full.fasta"
    log:
        config["output_path"] + "/logs/1_uk_full_untrimmed_alignment.log"
    shell:
        """
        datafunk sam_2_fasta \
          -s {input.sam} \
          -r {input.reference} \
          -o {output.fasta} \
          &> {log}
        """
        # fastafunk remove \
        #   --in-fasta {output.fasta} \
        #   --in-metadata {input.omit_list} \
        #   --out-fasta removed.fa
        # mv removed.fa {output.fasta}


rule run_snp_finder:
    input:
        fasta = rules.uk_full_untrimmed_alignment.output.fasta,
        snps = config["snps"]
    output:
        found = config["output_path"] + "/1/cog.snp_finder.csv",
    log:
        config["output_path"] + "/logs/1_run_snp_finder.log"
    shell:
        """
        datafunk snp_finder -a {input.fasta} -o {output.found} --snp-csv {input.snps} &> {log}
        """

rule add_snp_finder_result_to_metadata:
    input:
        snps = config["snps"],
        metadata = rules.uk_add_epi_week.output.metadata,
        new_data = rules.run_snp_finder.output.found
    output:
        metadata = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.with_snp_finder.csv"
    log:
        config["output_path"] + "/logs/1_add_snp_finder_result_to_metadata.log"
    shell:
        """
        columns=$(head -n1 {input.new_data} | cut -d',' -f2-)
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.new_data} \
          --index-column sequence_name \
          --join-on name \
          --new-columns "$columns" \
          --out-metadata {output.metadata} &>> {log}
        """

rule summarize_preprocess_uk:
    input:
        raw_fasta = config["latest_uk_fasta"],
        deduplicated_fasta = rules.uk_remove_duplicates.output.fasta,
        unify_headers_fasta = rules.uk_unify_headers.output.fasta,
        removed_low_covg_fasta = rules.uk_filter_low_coverage_sequences.output.fasta,
        full_alignment = rules.uk_full_untrimmed_alignment.output.fasta,
        full_metadata = rules.add_snp_finder_result_to_metadata.output.metadata
    params:
        webhook = config["webhook"],
        outdir = config["publish_path"] + "/COG",
        prefix = config["publish_path"] + "/COG/cog",
        export_dir = config["export_path"] + "/alignments",
        export_prefix = config["export_path"] + "/alignments/cog_" + config["date"]
    log:
        config["output_path"] + "/logs/1_summarize_preprocess_uk.log"
    shell:
        """
        echo "> Number of sequences in raw UK fasta: $(cat {input.raw_fasta} | grep ">" | wc -l)\\n" &> {log}
        echo "> Number of sequences after deduplication: $(cat {input.deduplicated_fasta} | grep ">" | wc -l)\\n" &>> {log}
        echo "> Number of sequences in raw UK fasta after unifying headers: $(cat {input.unify_headers_fasta} | grep ">" | wc -l)\\n" &>> {log}
        echo "> Number of sequences after trimming and removing those with <95% coverage: $(cat {input.removed_low_covg_fasta} | grep ">" | wc -l)\\n" &>> {log}
        echo ">\\n" >> {log}

        echo '{{"text":"' > 1_data.json
        echo "*Step 1: COG-UK preprocessing complete*\\n" >> 1_data.json
        cat {log} >> 1_data.json
        echo '"}}' >> 1_data.json
        echo "webhook {params.webhook}"
        curl -X POST -H "Content-type: application/json" -d @1_data.json {params.webhook}
        """



        # mkdir -p {params.outdir}
        # mkdir -p {params.export_dir}
        # cp {input.full_alignment} {params.prefix}_alignment.full.fasta
        # cp {input.full_alignment} {params.export_prefix}_alignment.full.fasta
        # echo "> Full untrimmed COG alignment published to _{params.prefix}_alignment.full.fasta_\\n" >> {log}
        # echo "> and to _{params.export_prefix}_alignment.full.fasta_\\n" >> {log}
        # echo ">\\n" >> {log}
        # cp {input.full_metadata} {params.prefix}_metadata.full.csv
        # cp {input.full_metadata} {params.export_prefix}_metadata.full.csv
        # echo "> Full COG only metadata published to _{params.prefix}_metadata.full.csv_\\n" >> {log}
        # echo "> and to _{params.export_prefix}_metadata.full.csv_\\n" >> {log}
        # echo ">\\n" >> {log}
        # cp {input.removed_low_covg_fasta} {params.prefix}_alignment.trimmed.fasta
        # echo "> Trimmed COG alignment published to _{params.prefix}_alignment.trimmed.fasta_\\n" >> {log}
