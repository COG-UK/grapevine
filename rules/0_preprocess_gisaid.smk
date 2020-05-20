import pandas as pd

rule gisaid_process_json:
    input:
        json = config["latest_gisaid_json"],
        omitted = config["previous_omitted_file"],
        metadata = config["previous_gisaid_metadata"]
    output:
        fasta = config["output_path"] + "/0/gisaid.fasta",
        metadata = config["output_path"] + "/0/gisaid.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_process_json.log"
    shell:
        """
        datafunk process_gisaid_data \
          --input-json \"{input.json}\" \
          --input-metadata \"{input.metadata}\" \
          --exclude-file \"{input.omitted}\" \
          --output-fasta {output.fasta} \
          --output-metadata {output.metadata} \
          --exclude-undated &> {log}
        """

rule gisaid_remove_duplicates:
    input:
        fasta = rules.gisaid_process_json.output.fasta,
        metadata = rules.gisaid_process_json.output.metadata
    output:
        fasta = config["output_path"] + "/0/gisaid.RD.fasta",
        metadata = config["output_path"] + "/0/gisaid.RD.csv"
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
# ^^^ gives you a subsample_omit column in your metadata - True if omit it


rule gisaid_unify_headers:
    input:
        fasta = rules.gisaid_remove_duplicates.output.fasta,
        metadata = rules.gisaid_remove_duplicates.output.metadata
    output:
        fasta = config["output_path"] + "/0/gisaid.RD.UH.fasta",
        metadata = config["output_path"] + "/0/gisaid.RD.UH.csv"
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


rule gisaid_counts_by_country:
    input:
        metadata = rules.gisaid_unify_headers.output.metadata
    output:
        counts = config["output_path"] + "/0/gisaid_counts_by_country.csv",
        published_counts = config["publish_path"] + "/GISAID/gisaid_counts_by_country.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_counts_by_country.log"
    shell:
        """
        fastafunk count \
          --in-metadata {input.metadata} \
          --group-column edin_admin_0 \
          --log-file {output.counts} &> {log}

        cp {output.counts} {output.published_counts}
        """


rule gisaid_extract_new:
    input:
        fasta = rules.gisaid_unify_headers.output.fasta,
        metadata = rules.gisaid_unify_headers.output.metadata,
        previous_metadata = config["previous_gisaid_metadata"]
    output:
        fasta = config["output_path"] + "/0/gisaid.RD.UH.new.fasta",
        metadata = config["output_path"] + "/0/gisaid.RD.UH.new.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_extract_new.log"
    shell:
        """
        fastafunk new \
          --in-fasta {input.fasta} \
          --in-metadata {input.previous_metadata} {input.metadata} \
          --index-column sequence_name \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --date-column 'edin_date_stamp' \
          --log {log}
        """
# ^^^ get's new things based on sequence_name if you don't specific an --index-column


rule gisaid_filter_1:
    input:
        fasta = rules.gisaid_extract_new.output.fasta
    params:
        min_covg = config["min_covg"],
        min_length = config["min_length"]
    output:
        fasta = config["output_path"] + "/0/gisaid.RD.new.UH.filt1.fasta"
    log:
        config["output_path"] + "/logs/0_gisaid_filter_1.log"
    shell:
        """
        datafunk filter_fasta_by_covg_and_length \
          -i {input.fasta} \
          -o {output} \
          --min-length {params.min_length} \
          --min-covg {params.min_covg} &> {log}
        """


rule gisaid_minimap2_to_reference:
    input:
        fasta = rules.gisaid_filter_1.output.fasta,
        reference = config["reference_fasta"]
    output:
        sam = config["output_path"] + "/0/gisaid.new.mapped.sam"
    log:
        config["output_path"] + "/logs/0_gisaid_minimap2_to_reference.log"
    shell:
        """
        minimap2 -a -x asm5 {input.reference} {input.fasta} -o {output} &> {log}
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
        fasta = config["output_path"] + "/0/gisaid.RD.new.UH.filt1.mapped.fasta"
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
          --log-all-inserts &> {log}

        mv insertions.txt {params.insertions}
        """


rule gisaid_filter_2:
    input:
        fasta = rules.gisaid_remove_insertions_and_pad.output.fasta
    output:
        fasta = config["output_path"] + "/0/gisaid.RD.new.UH.filt.mapped.filt2.fasta"
    params:
        min_covg = config["min_covg"]
    log:
        config["output_path"] + "/logs/0_gisaid_filter_2.log"
    shell:
        """
        datafunk filter_fasta_by_covg_and_length \
          -i {input.fasta} \
          -o {output.fasta} \
          --min-covg {params.min_covg} &> {log}
        """


rule gisaid_mask_1:
    input:
        fasta = rules.gisaid_filter_2.output.fasta,
        mask = config["gisaid_mask_file"]
    output:
        fasta = config["output_path"] + "/0/gisaid.RD.new.UH.filt.mapped.filt2.masked.fasta",
        published_fasta = config["publish_path"] + "/GISAID/gisaid.trimmed_alignment.new.fasta"
    shell:
        """
        datafunk mask \
          --input-fasta {input.fasta} \
          --output-fasta {output.fasta} \
          --mask-file \"{input.mask}\"

        cp {output.fasta} {output.published_fasta}
        """


rule gisaid_pangolin:
    input:
        fasta = rules.gisaid_mask_1.output.fasta
    params:
        outdir = config["output_path"] + "/0/pangolin",
        tmpdir = config["output_path"] + "/0/pangolin/tmp"
    output:
        lineages = protected(lineages = config["output_path"] + "/0/pangolin/lineage_report.csv")
    log:
        config["output_path"] + "/logs/0_gisaid_pangolin.log"
    threads: 40
    shell:
        """
        pangolin {input.fasta} \
        --threads {threads} \
        --tempdir {params.tmp} \
        --outdir {params.outdir} > {log} 2>&1
        """


# rule update_pangolin_lineages:
#     log:
#         config["output_path"] + "/logs/2_update_pangolin_lineages.log"
#     shell:
#         """
#         pip install --upgrade git+https://github.com/hCoV-2019/spangolin.git
#         """


# rule gisaid_special_pangolin:
#     input:
#         fasta = rules.gisaid_mask_1.output.fasta
#     params:
#         outdir = config["output_path"] + "/0/pangolin"
#         tmp = config["output_path"] + "/0/pangolin/tmp"
#     output:
#         lineages = config["output_path"] + "/0/pangolin/lineage_report.csv"
#     log:
#         config["output_path"] + "/logs/0_gisaid_special_pangolin.log"
#     threads: 40
#     shell:
#         """
#         pangolin {input.fasta} \
#         --threads {threads} \
#         --outdir {params.outdir} \
#         --tempdir {params.tmp} > {log} 2>&1
#         """


rule gisaid_add_pangolin_lineages_to_metadata:
    input:
        metadata = rules.gisaid_unify_headers.output.metadata,
        lineages = rules.gisaid_pangolin.output.lineages
    output:
        metadata = config["output_path"] + "/0/gisaid.RD.new.UH.filt.mapped.filt2.lineages.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_add_pangolin_lineages_to_metadata.log"
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.lineages} \
          --index-column sequence_name \
          --join-on taxon \
          --new-columns special_lineage lineage_support \
          --where-column lineage_support=UFbootstrap special_lineage=lineage \
          --out-metadata {output.metadata} &> {log}
        """


rule gisaid_combine_previous_and_new:
    input:
        previous_fasta = config["previous_gisaid_fasta"],
        new_fasta = rules.gisaid_mask_1.output.fasta,
        new_metadata = rules.gisaid_add_pangolin_lineages_to_metadata.output.metadata
    output:
        fasta = config["output_path"] + "/0/gisaid.combined.fasta",
        metadata = config["output_path"] + "/0/gisaid.combined.csv"
    params:
        outdir = config["publish_path"] + "/GISAID"
    log:
        config["output_path"] + "/logs/0_gisaid_combine_previous_and_new.log"
    shell:
        """
        fastafunk merge \
          --in-fasta {input.previous_fasta} {input.new_fasta} \
          --in-metadata {input.new_metadata} \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --index-column sequence_name \
          --log-file {log}
        """


rule gisaid_update_metadata_lineages:
    input:
        metadata = rules.gisaid_combine_previous_and_new.output.metadata
    output:
        metadata = config["output_path"] + "/0/gisaid.combined.updated.csv",
    log:
        config["output_path"] + "/logs/0_gisaid_update_metadata_lineages.log"
    run:
        df = pd.read_csv(input.metadata)
        lineages = []
        for i,row in df.iterrows():
            if row['special_lineage']:
                lineages.append(str(row['special_lineage']).replace(".X","").replace(".Y",""))
            else:
                lineages.append(row['lineage'])
        df['lineage'] = lineages
        df.to_csv(output.metadata, index=False)


# rule gisaid_publish_full_metadata:
#     input:
#         metadata = rules.gisaid_update_metadata_lineages.output.metadata
#     output:
#         published_metadata = config["publish_path"] + "/GISAID/gisaid.metadata.full.csv"
#     log:
#         config["output_path"] + "/logs/0_gisaid_publish_full_metadata.log"
#     shell:
#         """
#         cp {input.metadata} {output.published_metadata} 2> {log}
#         """


rule gisaid_mask_2:
    input:
        fasta = rules.gisaid_combine_previous_and_new.output.fasta,
        mask = config["gisaid_mask_file"]
    output:
        fasta = config["output_path"] + "/0/gisaid.full.masked.fasta",
    log:
        config["output_path"] + "/logs/0_gisaid_mask_2.log"
    shell:
        """
        datafunk mask \
          --input-fasta {input.fasta} \
          --output-fasta {output.fasta} \
          --mask-file \"{input.mask}\" 2> {log}
        """


# rule gisaid_publish_full_alignment:
#     input:
#         fasta = rules.gisaid_mask_2.output.fasta
#     output:
#         published_fasta = config["publish_path"] + "/GISAID/gisaid.trimmed_alignment.full.fasta"
#     log:
#         config["output_path"] + "/logs/0_gisaid_publish_full_alignment.log"
#     shell:
#         """
#         cp {input.fasta} {output.published_fasta} 2> {log}
#         """


rule gisaid_distance_QC:
    input:
        fasta = rules.gisaid_mask_2.output.fasta,
        metadata = rules.gisaid_update_metadata_lineages.output.metadata
    log:
        config["output_path"] + "/logs/0_gisaid_distance_QC.log"
    output:
        table = config["output_path"] + "/0/QC_distances.tsv",
    shell:
        """
        datafunk distance_to_root \
          --input-fasta {input.fasta} \
          --input-metadata {input.metadata} &> {log}

        mv distances.tsv {output.table}
        mv distances.png {output.plot}
        """


rule gisaid_output_lineage_table:
    input:
        fasta = rules.gisaid_mask_2.output.fasta,
        metadata = rules.gisaid_update_metadata_lineages.output.metadata
    output:
        fasta = config["output_path"] + "/0/gisaid.matched.fasta",
        metadata = config["output_path"] + "/0/gisaid.matched.lineages.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_output_lineage_table.log"
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name country adm1 adm2 \
                          sample_date epi_week \
                          lineage special_lineage uk_lineage \
          --where-column uk_omit=is_uk sample_date=covv_collection_date \
                                 epi_week=edin_epi_week country=edin_admin_0 \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --log-file {log} \
          --restrict
        """


# rule gisaid_output_gisaid:
#     input:
#         fasta = rules.gisaid_mask_2.output.fasta,
#         metadata = rules.gisaid_update_metadata_lineages.output.metadata
#     output:
#         fasta = config["output_path"] + "/0/gisaid.matched.fasta",
#         metadata = config["output_path"] + "/0/gisaid.matched.csv",
#         published_fasta = config["publish_path"] + "/GISAID/gisaid.trimmed_alignment.matched.fasta",
#         published_metadata = config["publish_path"] + "/GISAID/gisaid.metadata.matched.csv"
#     log:
#         config["output_path"] + "/logs/0_gisaid_output_gisaid.log"
#     shell:
#         """
#         fastafunk fetch \
#           --in-fasta {input.fasta} \
#           --in-metadata {input.metadata} \
#           --index-column sequence_name \
#           --filter-column sequence_name sample_date epi_week \
#                           country adm1 adm2 outer_postcode \
#                           is_surveillance is_community is_hcw \
#                           is_travel_history travel_history special_lineage \
#                           lineage_support uk_lineage lineage \
#           --where-column uk_omit=is_uk sample_date=covv_collection_date epi_week=edin_epi_week \
#                          country=edin_admin_0 travel_history=edin_travel \
#           --out-fasta {output.fasta} \
#           --out-metadata {output.metadata} \
#           --log-file {log} \
#           --restrict
#
#         cp {output.fasta} {output.published_fasta}
#         cp {output.metadata} {output.published_metadata}
#         """

rule summarize_preprocess_gisaid:
    input:
        latest_fasta = rules.gisaid_process_json.output.fasta,
        deduplicated_fasta = rules.gisaid_remove_duplicates.output.fasta,
        unify_headers_fasta = rules.gisaid_unify_headers.output.fasta,
        new_fasta = rules.gisaid_extract_new.output.fasta,
        removed_short_fasta = rules.gisaid_filter_1.output,
        removed_low_covg_fasta = rules.gisaid_filter_2.output.fasta,
        new_fasta_masked = rules.gisaid_mask_1.output.fasta,
        # full_fasta = rules.gisaid_publish_full_alignment.output.published_fasta,
        # full_metadata = rules.gisaid_publish_full_metadata.output.published_metadata,
        # matched_fasta = rules.gisaid_output_gisaid.output.fasta,
        # matched_metadata = rules.gisaid_output_gisaid.output.metadata,
        matched_fasta = rules.gisaid_output_lineage_table.output.fasta,
        matched_lineage_table = rules.gisaid_output_lineage_table.output.metadata,
    params:
        prefix = config["publish_path"] + "/GISAID/gisaid",
        webhook = config["webhook"]
    log:
        config["output_path"] + "/logs/0_summarize_preprocess_gisaid.log"
    shell:
        """
        echo "Number of sequences in latest GISAID download: $(cat {input.latest_fasta} | grep '>' | wc -l)\\n" >> {log}
        echo "Number of sequence after deduplicating: $(cat {input.deduplicated_fasta} | grep '>' | wc -l)\\n" >> {log}
        echo "Number of sequences after matching headers: $(cat {input.unify_headers_fasta} | grep '>' | wc -l)\\n" >> {log}
        echo "Number of new sequences: $(cat {input.new_fasta} | grep '>' | wc -l)\\n" >> {log}
        echo "Number of sequences after removing sequences <29000bps and with <95%% coverage: $(cat {input.removed_short_fasta} | grep '>' | wc -l)\\n" >> {log}
        echo "Number of sequences after mapping and removing those with <95%% coverage remaining: $(cat {input.removed_low_covg_fasta} | grep '>' | wc -l)\\n" >> {log}
        echo "> \\n" >> {log}
        echo "> Number of sequences in matched GISAID files: $(cat {input.matched_fasta} | grep '>' | wc -l)\\n" >> {log}
        echo ">\\n" >> {log}
        echo "> Counts by country published to {params.prefix}_counts_by_country.csv\\n" >> {log}
        echo '{{"text":"' > 0_data.json
        echo "*Step 0: GISAID preprocessing complete*\\n" >> 0_data.json
        cat {log} >> 0_data.json
        echo '"}}' >> 0_data.json
        echo 'webhook {params.webhook}'
        curl -X POST -H "Content-type: application/json" -d @0_data.json {params.webhook}
        # rm 0_data.json
        """

        # echo "> Full masked and trimmed GISAID alignment published to {params.prefix}.trimmed_alignment.full.fasta\\n" >> {log}
        # echo "> Full GISAID metadata published to {params.prefix}.metadata.full.fasta\\n" >> {log}
        # echo "> Matched GISAID fasta and restricted metadata published to {params.prefix}.trimmed_alignment.matched.fasta and {params.prefix}.metadata.matched.csv\\n" >> {log}
