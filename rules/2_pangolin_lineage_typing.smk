rule uk_extract_new:
    input:
        previous_stage = config["output_path"] + "/logs/1_summarize_preprocess_uk.log",
        fasta = rules.uk_filter_low_coverage_sequences.output,
        metadata = rules.uk_remove_duplicates.output.metadata,
        previous_metadata = config["previous_uk_metadata"]
    output:
        fasta = config["output_path"] + "/2/uk.new.fasta",
        metadata = config["output_path"] + "/2/uk.new.csv",
    log:
        config["output_path"] + "/logs/2_extract_new.log"
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

rule uk_pangolin:
    input:
        fasta = rules.uk_extract_new.output.fasta,
    params:
        outdir = config["output_path"] + "/2/pangolin"
    output:
        lineages = config["output_path"] + "/2/pangolin/lineage_report.csv"
    log:
        config["output_path"] + "/logs/2_uk_pangolin.log"
    threads: 40
    shell:
        """
        pangolin {input.fasta} \
        --threads {threads} \
        --outdir {params.outdir} > {log} 2>&1
        """

rule uk_add_previous_uk_lineages_to_metadata:
    input:
        previous_metadata = config["previous_uk_metadata"],
        metadata = rules.uk_remove_duplicates.output.metadata,
    output:
        metadata = config["output_path"] + "/2/uk.with_previous_lineages.csv",
    log:
        config["output_path"] + "/logs/2_uk_add_previous_uk_lineages_to_metadata.log"
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.previous_metadata} \
          --index-column sequence_name \
          --join-on sequence_name \
          --new-columns uk_lineage lineage lineage_support edin_date_stamp \
          --out-metadata {output.metadata} &>> {log}
        """

rule uk_add_pangolin_lineages_to_metadata:
    input:
        metadata = rules.uk_add_previous_uk_lineages_to_metadata.output.metadata,
        new_metadata = rules.uk_extract_new.output.metadata,
        lineages = rules.uk_pangolin.output.lineages
    output:
        metadata = config["output_path"] + "/2/uk.with_new_lineages.csv",
        tmp_metadata = temp(config["output_path"] + "/2/uk.with_new_lineages.csv.tmp")
    log:
        config["output_path"] + "/logs/2_uk_add_pangolin_lineages_to_metadata.log"
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.new_metadata} \
          --index-column sequence_name \
          --join-on sequence_name \
          --new-columns edin_date_stamp \
          --out-metadata {output.tmp_metadata} &> {log}

        fastafunk add_columns \
          --in-metadata {output.tmp_metadata} \
          --in-data {input.lineages} \
          --index-column sequence_name \
          --join-on taxon \
          --new-columns lineage lineage_support \
          --where-column lineage_support=UFbootstrap \
          --out-metadata {output.metadata} &>> {log}
        """

rule uk_output_cog:
    input:
        fasta = rules.uk_filter_low_coverage_sequences.output.fasta,
        metadata = rules.uk_add_pangolin_lineages_to_metadata.output.metadata
    output:
        fasta = config["output_path"] + "/2/uk.regularized.fasta",
        metadata = config["output_path"] + "/2/uk.regularized.csv"
    log:
        config["output_path"] + "/logs/2_uk_output_cog.log"
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name sample_date epi_week \
                          country adm1 adm2 outer_postcode \
                          is_surveillance is_community is_hcw \
                          is_travel_history travel_history lineage \
                          lineage_support uk_lineage \
          --where-column epi_week=edin_epi_week country=adm0 \
                         sample_date=reception_date sample_date=collection_date \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --log-file {log} \
          --restrict
        """

rule uk_output_cog_public:
    input:
        fasta = rules.uk_remove_duplicates.output.fasta,
        metadata = rules.uk_add_pangolin_lineages_to_metadata.output.metadata
    output:
        fasta = config["output_path"] + "/2/uk.public.fasta",
        metadata = config["output_path"] + "/2/uk.public.csv"
    log:
        config["output_path"] + "/logs/2_uk_output_cog_public.log"
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name country adm1 \
                          sample_date epi_week lineage \
                          lineage_support \
          --where-column epi_week=edin_epi_week country=adm0 \
                         sample_date=reception_date sample_date=collection_date \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --log-file {log} \
          --restrict
        """

rule summarize_pangolin_lineage_typing:
    input:
        fasta = rules.uk_output_cog.output.fasta,
        metadata = rules.uk_output_cog.output.metadata,
        public_fasta = rules.uk_output_cog_public.output.fasta,
        public_metadata = rules.uk_output_cog_public.output.metadata,
        full_metadata = rules.uk_add_pangolin_lineages_to_metadata.output.metadata
    params:
        webhook = config["webhook"],
        outdir = config["publish_path"] + "/COG",
        prefix = config["publish_path"] + "/COG/cog_"
    log:
        config["output_path"] + "/logs/2_summarize_pangolin_lineage_typing.log"
    shell:
        """
        mkdir -p {params.outdir}
        cp {input.full_metadata} {params.prefix}_metadata.full.csv
        echo "> Full COG metadata published to {params.prefix}_metadata.full.csv\\n" >> {log}

        cp {input.fasta} {params.prefix}_alignment.matched.fasta
        cp {input.metadata} {params.prefix}_metadata.matched.csv
        echo "> Matched COG fasta and restricted metadata published to {params.prefix}_alignment.matched.fasta and {params.prefix}_metadata.matched.csv\\n" >> {log}
        echo "> Number of sequences in matched COG files: $(cat {input.fasta} | grep ">" | wc -l)\\n" &>> {log}

        cp {input.public_fasta} {params.prefix}_sequences.public.fasta
        cp {input.public_metadata} {params.prefix}_metadata.public.csv
        echo "> Public unaligned COG fasta and restricted metadata published to {params.prefix}_sequences.public.fasta and {params.prefix}_metadata.public.csv\\n" >> {log}
        echo "> Number of sequences in public COG files: $(cat {input.public_fasta} | grep ">" | wc -l)\\n" &>> {log}

        echo '{{"text":"' > 2_data.json
        echo "*Step 2: COG-UK pangolin typing complete*\\n" >> 2_data.json
        cat {log} >> 2_data.json
        echo '"}}' >> 2_data.json
        echo "webhook {params.webhook}"
        curl -X POST -H "Content-type: application/json" -d @2_data.json {params.webhook}
        #rm 2_data.json
        """