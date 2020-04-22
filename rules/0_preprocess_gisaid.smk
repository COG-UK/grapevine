import datetime

date = datetime.date.today()

rule gisaid_process_json:
    input:
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
        datafunk process_gisaid_data \
          --input-json \"{input.json}\" \
          --input-metadata {input.metadata} \
          --exclude-file \"{input.omitted}\" \
          --output-fasta {output.fasta} \
          --output-metadata {output.metadata} \
          {params.flags} &> {log}
        """

rule gisaid_remove_duplicates:
    input:
        fasta = rules.gisaid_process_json.output.fasta,
        metadata = rules.gisaid_process_json.output.metadata
    output:
        fasta = config["output_path"] + "/0/gisaid_latest.deduplicated.fasta",
        metadata = config["output_path"] + "/0/gisaid_latest.deduplicated.csv",
    log:
        config["output_path"] + "/logs/0_gisaid_remove_duplicates.log"
    shell:
        """
        fastafunk subsample \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --group-column covv_virus_name \
          --index-column edin_header \
          --out-fasta {output.fasta} \
          --sample-size 1 \
          --out-metadata {output.metadata} \
          --select-by-min-column edin_epi_week &> {log}
        """

rule gisaid_unify_headers:
    input:
        fasta = rules.gisaid_remove_duplicates.output.fasta,
        metadata = rules.gisaid_remove_duplicates.output.metadata
    output:
        fasta = config["output_path"] + "/0/gisaid_latest.unify_headers.fasta",
        metadata = config["output_path"] + "/0/gisaid_latest.unify_headers.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_unify_headers.log"
    shell:
        """
        datafunk set_uniform_header \
          --input-fasta {input.fasta} \
          --input-metadata {input.metadata} \
          --output-fasta {output.fasta} \
          --output-metadata {output.metadata} \
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
          --min-length {params.min_length} &> {log}
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
        insertions = config["output_path"] + "/0/gisaid_insertions.txt"
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
          --log-inserts &> {log}
        mv insertions.txt {params.insertions}
        """

rule gisaid_remove_insertions_and_pad:
    input:
        sam = rules.gisaid_minimap2_to_reference.output.sam,
        reference = config["reference_fasta"]
    params:
        trim_start = config["trim_start"],
        trim_end = config["trim_end"],
        insertions = config["output_path"] + "/0/gisaid_insertions.txt"
    output:
        fasta = config["output_path"] + "/0/gisaid_latest.unify_headers.new.length_fitered.padded.fasta"
    log:
        config["output_path"] + "/logs/0_gisaid_remove_insertions_and_pad.log"
    shell:
        """
        datafunk sam_2_fasta \
          -s {input.sam} \
          -r {input.reference} \
          -o {output} \
          -t [{params.trim_start}:{params.trim_end}] \
          --pad \
          --log-inserts &> {log}
        mv insertions.txt {params.insertions}
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
          --min-covg {params.min_covg} &> {log}
        """

rule gisaid_filter_low_coverage_sequences_padded:
    input:
        fasta = rules.gisaid_remove_insertions_and_pad.output.fasta
    params:
        min_covg = config["min_covg"]
    output:
        fasta = config["output_path"] + "/0/gisaid_latest.unify_headers.new.length_fitered.padded.low_covg_filtered.fasta"
    log:
        config["output_path"] + "/logs/0_gisaid_filter_low_coverage_sequences_padded.log"
    shell:
        """
        datafunk filter_fasta_by_covg_and_length \
          -i {input.fasta} \
          -o {output} \
          --min-covg {params.min_covg} &> {log}
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
    threads: 32
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
          --index-column sequence_name \
          --join-on taxon \
          --new-columns lineage UFbootstrap \
          --out-metadata {output.metadata} &> {log}
        """

rule gisaid_combine_previous_and_new:
    input:
        previous_fasta = config["previous_gisaid_fasta"],
        previous_metadata = config["previous_gisaid_metadata"],
        new_fasta = rules.gisaid_filter_low_coverage_sequences.output.fasta,
        new_metadata = rules.gisaid_add_pangolin_lineages_to_metadata.output.metadata
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
          --index-column sequence_name \
          --log-file {log}

        if [ -f info/to_omit_gisaid.txt ]; then
            fastafunk remove \
              --in-fasta {output.fasta} \
              --in-metadata info/to_omit_gisaid.txt \
              --out-fasta tmp.fa
            mv tmp.fa {output.fasta}
        fi
        """

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
        deduplicated_fasta = rules.gisaid_remove_duplicates.output.fasta,
        unify_headers_fasta = rules.gisaid_unify_headers.output.fasta,
        new_fasta = rules.gisaid_extract_new.output.fasta,
        removed_short_fasta = rules.gisaid_filter_short_sequences.output,
        removed_low_covg_fasta = rules.gisaid_filter_low_coverage_sequences.output,
        final_fasta = rules.gisaid_combine_previous_and_new.output.fasta,
        final_metadata = rules.gisaid_combine_previous_and_new.output.metadata
    log:
        config["output_path"] + "/logs/0_summary_preprocess_gisaid.log"
    shell:
        """
        echo "Number of sequences in previous GISAID fasta: $(cat {input.previous_fasta} | grep ">" | wc -l)" >> {log}
        echo "Number of sequences in latest GISAID download: $(cat {input.latest_fasta} | grep ">" | wc -l)" >> {log}
        echo "Number of deduplicated sequences: $(cat {input.latest_fasta} | grep ">" | wc -l)" >> {log}
        echo "Number of sequences after matching headers: $(cat {input.unify_headers_fasta} | grep ">" | wc -l)" >> {log}
        echo "Number of new sequences: $(cat {input.new_fasta} | grep ">" | wc -l)" >> {log}
        echo "Number of sequences after removing sequences <29000bps: $(cat {input.removed_short_fasta} | grep ">" | wc -l)" >> {log}
        echo "Number of sequences after trimming and removing those with <95% coverage: $(cat {input.removed_low_covg_fasta} | grep ">" | wc -l)" >> {log}

        curl -X POST -H ‘Content-type: application/json’ --data ‘{“text”:$(cat ${log})}’ https://hooks.slack.com/services/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
        """

rule gisaid_output_gisaid:
    input:
        fasta = rules.gisaid_combine_previous_and_new.output.fasta,
        metadata = rules.gisaid_combine_previous_and_new.output.metadata
    params:
        outdir = config["publish_path"] + "/GISAID",
        prefix = config["publish_path"] + "/GISAID/cog_%s" %date
    output:
        fasta = config["output_path"] + "/0/gisaid.regularized.fasta",
        metadata = config["output_path"] + "/0/gisaid.regularized.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_output_gisaid.log"
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name collection_date epi_week \
                          country adm1 adm2 outer_postcode \
                          is_surveillance is_community is_hcw \
                          is_travel_history travel_history lineage
                          lineage_support uk_lineage \
          --where-column collection_date=covv_collection_date epi_week=edin_epi_week \
                         country=edin_admin_0 travel_history=edin_travel lineage_support=ufbootstrap \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --log-file {log} \
          --restrict

        echo "Number of sequences in combined full GISAID fasta: $(cat {input.fasta} | grep ">" | wc -l)" >> {log}
        echo "Number of lines in combined full GISAID metadata: $(cat {input.metadata} | wc -l)" >> {log}
        echo "Published to {params.prefix}_alignment.full.fasta and {params.prefix}_metadata.full.fasta" >> {log}

        echo "Number of non-omitted sequences in GISAID fasta: $(cat {output.fasta} | grep ">" | wc -l)" >> {log}
        echo "Number of lines in regularized GISAID metadata: $(cat {output.metadata} | wc -l)" >> {log}
        echo "Published to {params.prefix}_alignment.fasta and {params.prefix}_metadata.fasta" >> {log}

        mkdir -p {params.outdir}
        cp {input.fasta} {params.prefix}_alignment.full.fasta
        cp {input.metadata} {params.prefix}_metadata.full.fasta
        cp {output.fasta} {params.prefix}_alignment.fasta
        cp {output.metadata} {params.prefix}_metadata.fasta

        curl -X POST -H ‘Content-type: application/json’ --data ‘{“text”:$(cat ${log})}’ https://hooks.slack.com/services/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
        """