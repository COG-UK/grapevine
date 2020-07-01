

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
        config["output_path"] + "/logs/1_uk_filter_duplicates_bycovid.log"
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


rule uk_add_epi_week:
    input:
        metadata = rules.uk_remove_duplicates_covid_by_gaps.output.metadata
    output:
        metadata = config["output_path"] + "/1/uk_latest.epi_week.csv",
        tmp_metadata = temp(config["output_path"] + "/1/uk_latest.epi_week.csv.tmp")
    log:
        config["output_path"] + "/logs/1_uk_add_epi_week.log"
    shell:
        """
        datafunk add_epi_week \
        --input-metadata {input.metadata} \
        --output-metadata {output.tmp_metadata} \
        --date-column received_date \
        --epi-week-column-name edin_epi_week \
        --epi-day-column-name edin_epi_day &> {log}

        datafunk add_epi_week \
        --input-metadata {output.tmp_metadata} \
        --output-metadata {output.metadata} \
        --date-column collection_date \
        --epi-week-column-name edin_epi_week \
        --epi-day-column-name edin_epi_day &>> {log}
        """


rule uk_annotate_to_remove_duplicates_by_biosample:
    input:
        metadata = rules.uk_add_epi_week.output.metadata,
    output:
        metadata = config["output_path"] + "/1/uk_latest.epi_week.annotated2.csv",
    log:
        config["output_path"] + "/logs/1_uk_annotate_to_remove_duplicates_by_biosample.log",
    run:
        import pandas as pd

        df = pd.read_csv(input.metadata)

        edin_biosample = []

        for i,row in df.iterrows():
            if pd.isnull(row['biosample_source_id']):
                edin_biosample.append(row['cov_id'])
            else:
                edin_biosample.append(row['biosample_source_id'])

        df['edin_biosample'] = edin_biosample
        df.to_csv(output.metadata, index=False)



rule uk_remove_duplicates_biosamplesourceid_by_date:
    input:
        fasta = rules.uk_remove_duplicates_covid_by_gaps.output.fasta,
        metadata = rules.uk_annotate_to_remove_duplicates_by_biosample.output.metadata
    output:
        fasta = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated_cov_id_biosample_source_id.fasta",
        metadata = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated_cov_id_biosample_source_id.csv"
    log:
        config["output_path"] + "/logs/1_uk_filter_duplicates_by_biosample.log"
    shell:
        """
        fastafunk subsample \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --group-column edin_biosample \
          --index-column header \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --sample-size 1 \
          --select-by-min-column edin_epi_day &> {log}
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


rule uk_mask_1:
    input:
        fasta = rules.uk_remove_insertions_and_trim_and_pad.output.fasta,
        mask = config["gisaid_mask_file"]
    output:
        fasta = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.alignment.trimmed.masked.fasta",
    log:
        config["output_path"] + "/logs/1_uk_mask_1.log"
    shell:
        """
        datafunk mask \
          --input-fasta {input.fasta} \
          --output-fasta {output.fasta} \
          --mask-file \"{input.mask}\" 2> {log}
        """


rule uk_filter_low_coverage_sequences:
    input:
        fasta = rules.uk_mask_1.output.fasta
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


rule uk_mask_2:
    input:
        fasta = rules.uk_full_untrimmed_alignment.output.fasta,
        mask = config["gisaid_mask_file"]
    output:
        fasta = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.alignment.full.masked.fasta"
    log:
        config["output_path"] + "/logs/1_uk_mask_2.log"
    shell:
        """
        datafunk mask \
          --input-fasta {input.fasta} \
          --output-fasta {output.fasta} \
          --mask-file \"{input.mask}\" 2> {log}
        """


rule run_snp_finder:
    input:
        fasta = rules.uk_mask_2.output.fasta,
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
        metadata = rules.uk_unify_headers.output.metadata,
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


rule uk_del_finder:
    input:
        fasta = rules.uk_mask_2.output.fasta,
        dels = config["dels"]
    output:
        metadata = config["output_path"] + "/1/cog.del_finder.csv",
    log:
        config["output_path"] + "/logs/1_uk_del_finder.log"
    shell:
        """
        datafunk del_finder \
            -i {input.fasta} \
            --deletions-file {input.dels} \
            --genotypes-table {output.metadata} &> {log}
        """


rule uk_add_del_finder_result_to_metadata:
    input:
        dels = config["dels"],
        metadata = rules.add_snp_finder_result_to_metadata.output.metadata,
        new_data = rules.uk_del_finder.output.metadata
    output:
        metadata = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.with_snp_finder.with_del_finder.csv"
    log:
        config["output_path"] + "/logs/1_add_del_finder_result_to_metadata.log"
    shell:
        """
        columns=$(head -n1 {input.new_data} | cut -d',' -f2-)
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.new_data} \
          --index-column sequence_name \
          --join-on sequence_name \
          --new-columns "$columns" \
          --out-metadata {output.metadata} &>> {log}
        """


"""
Instead of new sequences (as determined by a date stamp), it might be more robust
to extract sequences for lineage typing that don't currently have an associated
lineage designation in the metadata file.
"""
rule uk_add_previous_lineages_to_metadata:
    input:
        metadata = rules.uk_add_del_finder_result_to_metadata.output.metadata,
        previous_metadata = config["previous_uk_metadata"],
        global_lineages = config["global_lineages"]
    output:
        metadata_temp = temp(config["output_path"] + "/1/uk.with_previous_lineages.temp.csv"),
        metadata = config["output_path"] + "/1/uk.with_previous_lineages.csv",
    log:
        config["output_path"] + "/logs/1_uk_add_previous_lineages_to_metadata.log"
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.previous_metadata} \
          --index-column sequence_name \
          --join-on sequence_name \
          --new-columns uk_lineage lineage lineage_support lineages_version edin_date_stamp \
          --out-metadata {output.metadata_temp} &>> {log}

        fastafunk add_columns \
          --in-metadata {output.metadata_temp} \
          --in-data {input.global_lineages} \
          --index-column sequence_name \
          --join-on taxon \
          --new-columns lineage lineage_support lineages_version \
          --where-column lineage_support=UFbootstrap \
          --out-metadata {output.metadata} &>> {log}
        """


rule uk_extract_lineageless:
    input:
        fasta = rules.uk_filter_low_coverage_sequences.output,
        metadata = rules.uk_add_previous_lineages_to_metadata.output.metadata,
    output:
        fasta = config["output_path"] + "/1/uk.new.pangolin_lineages.fasta",
    log:
        config["output_path"] + "/logs/1_extract_lineageless.log"
    run:
        from Bio import SeqIO
        import pandas as pd

        fasta_in = SeqIO.index(str(input.fasta), "fasta")
        df = pd.read_csv(input.metadata)

        sequence_record = []

        with open(str(output.fasta), 'w') as fasta_out:
            for i,row in df.iterrows():
                if pd.isnull(row['lineage']):
                    sequence_name = row['sequence_name']
                    if sequence_name in fasta_in:
                        if sequence_name not in sequence_record:
                            record = fasta_in[sequence_name]
                            fasta_out.write('>' + record.id + '\n')
                            fasta_out.write(str(record.seq) + '\n')
                            sequence_record.append(sequence_name)


rule summarize_preprocess_uk:
    input:
        raw_fasta = config["latest_uk_fasta"],
        deduplicated_fasta_by_covid = rules.uk_remove_duplicates_covid_by_gaps.output.fasta,
        deduplicated_fasta_by_biosampleid = rules.uk_remove_duplicates_biosamplesourceid_by_date.output.fasta,
        unify_headers_fasta = rules.uk_unify_headers.output.fasta,
        removed_low_covg_fasta = rules.uk_filter_low_coverage_sequences.output.fasta,
        full_alignment = rules.uk_full_untrimmed_alignment.output.fasta,
        full_metadata = rules.add_snp_finder_result_to_metadata.output.metadata,
        lineageless_fasta = rules.uk_extract_lineageless.output.fasta
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
        echo "> Number of sequences after deduplication by cov_id: $(cat {input.deduplicated_fasta_by_covid} | grep ">" | wc -l)\\n" &>> {log}
        echo "> Number of sequences after deduplication by bio_sample_id: $(cat {input.deduplicated_fasta_by_biosampleid} | grep ">" | wc -l)\\n" &>> {log}
        echo "> Number of sequences after unifying headers: $(cat {input.unify_headers_fasta} | grep ">" | wc -l)\\n" &>> {log}
        echo "> Number of sequences after trimming and removing those with <95% coverage: $(cat {input.removed_low_covg_fasta} | grep ">" | wc -l)\\n" &>> {log}
        echo "> Number of new sequences passed to Pangolin for typing: $(cat {input.lineageless_fasta} | grep ">" | wc -l)\\n" &>> {log}
        echo ">\\n" >> {log}

        echo '{{"text":"' > 1_data.json
        echo "*Step 1: COG-UK preprocessing complete*\\n" >> 1_data.json
        cat {log} >> 1_data.json
        echo '"}}' >> 1_data.json
        echo "webhook {params.webhook}"
        curl -X POST -H "Content-type: application/json" -d @1_data.json {params.webhook}
        """
