

rule uk_add_lineage_information_back_to_master_metadata:
    input:
        metadata = rules.uk_update_metadata_lineages.output.metadata,
        lineage_data = rules.run_5_subroutine_on_lineages.output.metadata
    output:
        metadata = config["output_path"] + "/6/uk.master.csv"
    log:
        config["output_path"] + "/logs/6_uk_update_master_metadata.log"
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.lineage_data} \
          --index-column sequence_name \
          --join-on sequence_name \
          --new-columns uk_lineage acc_lineage del_lineage phylotype \
          --out-metadata {output.metadata} &> {log}
        """


rule publish_COG_master_metadata:
    input:
        metadata = rules.uk_add_lineage_information_back_to_master_metadata.output.metadata
    output:
        metadata = config["publish_path"] + "/COG/master.csv"
    log:
        config["output_path"] + "/logs/6_publish_COG_master_metadata.log"
    shell:
        """
        cp {input.metadata} {output.metadata} &> {log}
        """


rule publish_gisaid_master_metadata:
    input:
        metadata = config["output_path"] + "/0/gisaid.combined.updated.csv"
    output:
        metadata = config["publish_path"] + "/GISAID/master.csv"
    log:
        config["output_path"] + "/logs/6_publish_gisaid_master_metadata.log"
    shell:
        """
        cp {input.metadata} {output.metadata} &> {log}
        """


rule publish_unaligned_cog_sequences:
    input:
        fasta = rules.uk_unify_headers.output.fasta
    output:
        fasta = config["export_path"] + "/alignments/cog_" + config["date"] + '_all.fasta'
    log:
        config["output_path"] + "/logs/6_publish_unaligned_cog_sequences.log"
    shell:
        """
        cp {input.fasta} {output.fasta} 2> {log}
        """


rule publish_full_aligned_cog_data:
    input:
        fasta = rules.uk_full_untrimmed_alignment.output.fasta,
        metadata = rules.uk_add_lineage_information_back_to_master_metadata.output.metadata
    output:
        fasta = config["export_path"] + "/alignments/cog_" + config["date"] + '_all_alignment.fasta',
        metadata = config["export_path"] + "/alignments/cog_" + config["date"] + '_all_metadata.csv'
    log:
        config["output_path"] + "/logs/6_publish_full_aligned_cog_data.log"
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column adm0 adm1 adm2 biosample_source_id central_sample_id \
                          collection_date flowcell_id flowcell_type instrument_make instrument_model \
                          layout_insert_length layout_read_length library_layout_config library_name \
                          library_primers library_protocol library_selection library_seq_kit \
                          library_seq_protocol library_source library_strategy meta.artic.primers \
                          meta.artic.protocol metadata published_as received_date root_sample_id \
                          run_group run_name sample_type_collected sample_type_received sampling_strategy \
                          secondary_accession secondary_identifier sequencing_org sequencing_org_code \
                          sequencing_submission_date sequencing_uuid source_age source_sex start_time \
                          submission_org submission_org_code submission_user swab_site header sequence_name \
                          length missing gaps cov_id subsample_omit edin_epi_week d614g \
          --where-column country=adm0 \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --restrict &> {log}
        """


rule publish_filtered_aligned_cog_data:
    input:
        fasta = rules.uk_filter_low_coverage_sequences.output.fasta,
        metadata = rules.uk_add_lineage_information_back_to_master_metadata.output.metadata
    output:
        fasta = config["export_path"] + "/alignments/cog_" + config["date"] + '_alignment.fasta',
        metadata = config["export_path"] + "/alignments/cog_" + config["date"] + '_metadata.csv'
    log:
        config["output_path"] + "/logs/6_publish_filtered_aligned_cog_data.log"
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
                          lineage_support uk_lineage acc_lineage del_lineage phylotype \
          --where-column epi_week=edin_epi_week country=adm0 \
                         sample_date=received_date sample_date=collection_date \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --restrict &> {log}
        """


rule publish_full_annotated_tree_and_metadata:
    input:
        tree = rules.run_5_subroutine_on_lineages.output.full_tree,
        cog_fasta = rules.uk_filter_low_coverage_sequences.output.fasta,
        cog_metadata = rules.uk_add_lineage_information_back_to_master_metadata.output.metadata,
        gisaid_fasta = config["output_path"] + "/0/gisaid.full.masked.fasta",
        gisaid_metadata = config["output_path"] + "/0/gisaid.combined.updated.csv"
    params:
        intermediate_cog_fasta = config["output_path"] + "/6/cog.publish_full_annotated_tree_and_metadata.temp.cog.fasta",
        intermediate_cog_metadata = config["output_path"] + "/6/cog.publish_full_annotated_tree_and_metadata.temp.cog.csv",
        intermediate_gisaid_fasta = config["output_path"] + "/6/gisaid.publish_full_annotated_tree_and_metadata.temp.gisaid.fasta",
        intermediate_gisaid_metadata = config["output_path"] + "/6/gisaid.publish_full_annotated_tree_and_metadata.temp.gisaid.csv",
        intermediate_combined_fasta = config["output_path"] + "/6/gisaid.publish_full_annotated_tree_and_metadata.temp.combined.fasta",
        intermediate_combined_metadata = config["output_path"] + "/6/gisaid.publish_full_annotated_tree_and_metadata.temp.combined.csv",
    output:
        tree = config["export_path"] + "/trees/cog_global_" + config["date"] + '_tree.nexus',
        metadata = config["export_path"] + "/trees/cog_global_" + config["date"] + '_metadata.csv',
        fasta = config["output_path"] + "/6/cog_global.fasta"
    log:
        config["output_path"] + "/logs/6_publish_full_annotated_tree_and_metadata.log"
    shell:
        """
        cp {input.tree} {output.tree} &> {log}


        fastafunk fetch \
          --in-fasta {input.cog_fasta} \
          --in-metadata {input.cog_metadata} \
          --index-column sequence_name \
          --filter-column sequence_name sample_date epi_week \
                          country adm1 adm2 outer_postcode \
                          is_surveillance is_community is_hcw \
                          is_travel_history travel_history lineage \
                          lineage_support uk_lineage acc_lineage del_lineage phylotype \
          --where-column epi_week=edin_epi_week country=adm0 \
                         sample_date=received_date sample_date=collection_date \
          --out-fasta {params.intermediate_cog_fasta} \
          --out-metadata {params.intermediate_cog_metadata} \
          --restrict &>> {log}


        fastafunk fetch \
          --in-fasta {input.gisaid_fasta} \
          --in-metadata {input.gisaid_metadata} \
          --index-column sequence_name \
          --filter-column sequence_name sample_date epi_week \
                          country adm1 adm2 outer_postcode \
                          is_surveillance is_community is_hcw \
                          is_travel_history travel_history lineage \
                          lineage_support uk_lineage acc_lineage del_lineage phylotype \
          --where-column uk_omit=is_uk sample_date=covv_collection_date epi_week=edin_epi_week \
                         country=edin_admin_0 travel_history=edin_travel \
          --out-fasta {params.intermediate_gisaid_fasta} \
          --out-metadata {params.intermediate_gisaid_metadata} \
          --restrict &>> {log}


        fastafunk merge \
          --in-fasta {params.intermediate_gisaid_fasta} {params.intermediate_cog_fasta} \
          --in-metadata {params.intermediate_gisaid_metadata} {params.intermediate_cog_metadata} \
          --out-fasta {params.intermediate_combined_fasta} \
          --out-metadata {params.intermediate_combined_metadata} \
          --index-column sequence_name &>> {log}

        fastafunk fetch \
          --in-fasta {params.intermediate_combined_fasta} \
          --in-metadata {params.intermediate_combined_metadata} \
          --index-column sequence_name \
          --filter-column sequence_name sample_date epi_week \
                          country adm1 adm2 outer_postcode \
                          is_surveillance is_community is_hcw \
                          is_travel_history travel_history lineage \
                          lineage_support uk_lineage acc_lineage del_lineage phylotype \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --restrict &>> {log}
        """


rule publish_public_cog_data:
    input:
        public_tree = rules.run_4_subroutine_on_lineages.output.public_tree,
        fasta = rules.uk_unify_headers.output.fasta,
        metadata = rules.uk_add_lineage_information_back_to_master_metadata.output.metadata
    output:
        public_tree = config["export_path"] + "/public/cog_global_" + config["date"] + "_tree.newick",
        fasta = config["export_path"] + "/public/cog_" + config["date"] + ".fasta",
        metadata = config["export_path"] + "/public/cog_" + config["date"] + "_metadata.csv"
    log:
        config["output_path"] + "/logs/6_publish_public_cog_data.log"
    shell:
        """
        cp {input.public_tree} {output.public_tree} &> {log}

        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name country adm1 \
                          sample_date epi_week lineage \
                          lineage_support \
          --where-column epi_week=edin_epi_week country=adm0 \
                         sample_date=received_date sample_date=collection_date \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --restrict &>> {log}
        """


rule summarize_publish:
    input:
        GISAID_meta_master = rules.publish_gisaid_master_metadata.output.metadata,
        COG_meta_master = rules.publish_COG_master_metadata.output.metadata,

        COG_seq_all = rules.publish_unaligned_cog_sequences.output.fasta,
        COG_seq_all_aligned = rules.publish_full_aligned_cog_data.output.fasta,
        COG_meta_all_aligned = rules.publish_full_aligned_cog_data.output.metadata,

        COG_seq_all_aligned_filtered = rules.publish_filtered_aligned_cog_data.output.metadata,
        COG_meta_all_aligned_filtered = rules.publish_filtered_aligned_cog_data.output.metadata,

        COG_GISAID_nexus_tree = rules.publish_full_annotated_tree_and_metadata.output.tree,
        COG_GISAID_meta = rules.publish_full_annotated_tree_and_metadata.output.metadata,

        public_COG_GISAID_newick_tree = rules.publish_public_cog_data.output.public_tree,
        public_COG_GISAID_seq_all = rules.publish_unaligned_cog_sequences.output.fasta,
        public_COG_meta = rules.publish_public_cog_data.output.metadata
    params:
        webhook = config["webhook"],
        uk_trees_path = config["export_path"] + "/trees/uk_lineages/",
    log:
        config["output_path"] + "/logs/6_summarize_publish.log"
    shell:
        """
        echo "> Gisaid master metadata published to {input.GISAID_meta_master}\\n" >> {log}
        echo "> COG master metadata published to {input.COG_meta_master}\\n" >> {log}
        echo "> \\n" >> {log}
        echo "> Unaligned (deduplicated, clean headers) COG sequences published to {input.COG_seq_all}\\n" >> {log}
        echo "> \\n" >> {log}
        echo "> Aligned (deduplicated, clean headers) COG sequences published to {input.COG_seq_all_aligned}\\n" >> {log}
        echo "> Matching metadata published to {input.COG_meta_all_aligned}\\n" >> {log}
        echo "> \\n" >> {log}
        echo "> Filtered, aligned COG sequences published to {input.COG_seq_all_aligned_filtered}\\n" >> {log}
        echo "> Matching metadata with lineage information published to {input.COG_meta_all_aligned_filtered}\\n" >> {log}
        echo "> \\n" >> {log}
        echo "> Full, annotated tree published to {input.COG_GISAID_nexus_tree}\\n" >> {log}
        echo "> Matching metadata published to {input.COG_GISAID_meta}\\n" >> {log}
        echo "> UK lineage subtrees published in {params.uk_trees_path}\\n" >> {log}
        echo "> \\n" >> {log}
        echo "> Public tree published to {input.public_COG_GISAID_newick_tree}\\n" >> {log}
        echo "> Associated unaligned sequences published to {input.public_COG_GISAID_seq_all}\\n" >> {log}
        echo "> Matching metadata with public fields only published to {input.public_COG_meta}\\n" >> {log}
        echo "> \\n" >> {log}

        echo '{{"text":"' > 6_data.json
        echo "*Step 6: publish data complete*\\n" >> 6_data.json
        cat {log} >> 6_data.json
        echo '"}}' >> 6_data.json
        echo 'webhook {params.webhook}'
        curl -X POST -H "Content-type: application/json" -d @6_data.json {params.webhook}
        """










#
