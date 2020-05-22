


rule uk_add_lineage_information_back_to_master_metadata:
    input:
        metadata = config["output_path"] + "/2/uk.with_new_lineages.csv",
        uk_lineage_data = config["output_path"] + "/5/cog_gisaid.lineages.with_all_traits.with_phylotype_traits.csv",
        global_lineage_data = config["global_lineages"]
    output:
        metadata_temp = temp(config["output_path"] + "/6/temp.uk.master.csv")
        metadata = config["output_path"] + "/6/uk.master.csv"
    log:
        config["output_path"] + "/logs/6_uk_add_lineage_information_back_to_master_metadata.log"
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.uk_lineage_data} \
          --index-column sequence_name \
          --join-on sequence_name \
          --new-columns uk_lineage microreact_lineages acc_lineage del_lineage acc_introduction del_introduction phylotype \
          --out-metadata {output.metadata_temp} &> {log}

        fastafunk add_columns \
          --in-metadata {output.metadata_temp} \
          --in-data {input.global_lineage_data} \
          --index-column sequence_name \
          --join-on taxon  \
          --new-columns lineage lineage_support lineages_version \
          --where-column lineage_support=UFbootstrap \
          --out-metadata {output.metadata} &>> {log}
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


rule gisaid_add_lineage_information_back_to_master_metadata:
    input:
        metadata = config["output_path"] + "/0/gisaid.combined.csv",
        global_lineage_data = config["global_lineages"]
    output:
        metadata = config["output_path"] + "/6/gisaid.master.csv"
    log:
        config["output_path"] + "/logs/6_gisaid_add_lineage_information_back_to_master_metadata.log"
    shell:
        """
      fastafunk add_columns \
        --in-metadata {input.metadata} \
        --in-data {input.global_lineage_data} \
        --index-column sequence_name \
        --join-on taxon  \
        --new-columns lineage lineage_support lineages_version \
        --where-column lineage_support=UFbootstrap \
        --out-metadata {output.metadata} &>> {log}
        """


rule publish_gisaid_master_metadata:
    input:
        metadata = rules.gisaid_add_lineage_information_back_to_master_metadata.output.metadata,
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
        fasta = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated.unify_headers.fasta",
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
        fasta = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.alignment.full.masked.fasta",
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
        fasta = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.trimmed.low_covg_filtered.fasta",
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


rule combine_cog_gisaid:
    input:
        cog_fasta = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.trimmed.low_covg_filtered.fasta",
        cog_metadata = rules.uk_add_lineage_information_back_to_master_metadata.output.metadata,
        gisaid_fasta = config["output_path"] + "/0/gisaid.full.masked.filtered.fasta",
        gisaid_metadata = config["output_path"] + "/0/gisaid.combined.csv"
    params:
        intermediate_cog_fasta = config["output_path"] + "/6/cog.combine_cog_gisaid.temp.cog.fasta",
        intermediate_cog_metadata = config["output_path"] + "/6/cog.combine_cog_gisaid.temp.cog.csv",
        intermediate_gisaid_fasta = config["output_path"] + "/6/gisaid.combine_cog_gisaid.temp.gisaid.fasta",
        intermediate_gisaid_metadata = config["output_path"] + "/6/gisaid.combine_cog_gisaid.temp.gisaid.csv",
    output:
        fasta = config["output_path"] + "/6/gisaid.combine_cog_gisaid.combined.fasta",
        metadata = config["output_path"] + "/6/gisaid.combine_cog_gisaid.combined.csv",
    log:
        config["output_path"] + "/logs/6_combine_cog_gisaid.log"
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.cog_fasta} \
          --in-metadata {input.cog_metadata} \
          --index-column sequence_name \
          --filter-column covv_accession_id \
                          sequence_name sample_date epi_week \
                          country adm1 adm2 submission_org_code \
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
          --filter-column covv_accession_id \
                          sequence_name sample_date epi_week \
                          country adm1 adm2 submission_org_code \
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
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --index-column sequence_name &>> {log}
        """


rule publish_full_annotated_tree_and_metadata:
    input:
        newick_tree = config["output_path"] + "/4/cog_gisaid_full.tree.public.newick",
        annotated_tree = config["output_path"] + "/5/cog_gisaid_full.tree.nexus",
        combined_fasta = rules.combine_cog_gisaid.output.fasta,
        combined_metadata = rules.combine_cog_gisaid.output.metadata,
    output:
        newick_tree = config["export_path"] + "/trees/cog_global_" + config["date"] + '_tree.newick',
        annotated_tree = config["export_path"] + "/trees/cog_global_" + config["date"] + '_tree.nexus',
        metadata = config["export_path"] + "/trees/cog_global_" + config["date"] + '_metadata.csv',
        fasta = config["output_path"] + "/6/cog_global.fasta",
    log:
        config["output_path"] + "/logs/6_publish_full_annotated_tree_and_metadata.log"
    shell:
        """
        cp {input.annotated_tree} {output.annotated_tree} &> {log}
        cp {input.newick_tree} {output.newick_tree} &>> {log}


        fastafunk fetch \
          --in-fasta {input.combined_fasta} \
          --in-metadata {input.combined_metadata} \
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


rule publish_cog_gisaid_data_for_lineage_release_work:
    input:
        combined_fasta = rules.combine_cog_gisaid.output.fasta,
        combined_metadata = rules.combine_cog_gisaid.output.metadata,
    output:
        fasta = config["export_path"] + "/lineage_release/cog_gisaid.fasta",
        metadata = config["export_path"] + "/lineage_release/cog_gisaid.csv",
    log:
        config["output_path"] + "/logs/6_publish_cog_gisaid_data_for_lineage_release_work.log"
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.combined_fasta} \
          --in-metadata {input.combined_metadata} \
          --index-column sequence_name \
          --filter-column covv_accession_id sequence_name country travel_history sample_date epi_week \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --restrict &>> {log}
        """


rule publish_public_cog_data:
    input:
        public_tree = config["output_path"] + "/4/cog_gisaid_full.tree.public.newick",
        fasta = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated.unify_headers.fasta",
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


rule publish_microreact_specific_output:
    input:
        newick_tree = config["output_path"] + "/4/cog_gisaid_full.tree.public.newick",
        metadata = rules.combine_cog_gisaid.output.metadata,
        fasta = rules.combine_cog_gisaid.output.fasta,
    output:
        newick_tree = config["export_path"] + "/microreact/cog_global_" + config["date"] + '_tree.newick',
        public_metadata = config["export_path"] + "/microreact/cog_global_" + config["date"] + '_metadata_public.csv',
        private_metadata = config["export_path"] + "/microreact/cog_global_" + config["date"] + '_metadata_private.csv',
        fasta1 = temp(config["output_path"] + "/6/cog_global_microreact1.fasta"),
        fasta2 = temp(config["output_path"] + "/6/cog_global_microreact2.fasta")
    log:
        config["output_path"] + "/logs/6_publish_microreact_specific_output.log"
    shell:
        """
        cp {input.newick_tree} {output.newick_tree}

        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name sample_date epi_week \
                          country adm1 adm2 submission_org_code lineage \
                          lineage_support uk_lineage \
          --out-fasta {output.fasta1} \
          --out-metadata {output.public_metadata} \
          --restrict &>> {log}

        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name sample_date epi_week \
                          country adm1 adm2 submission_org_code \
                          is_hcw travel_history \
                          lineage lineage_support uk_lineage d614g \
          --out-fasta {output.fasta2} \
          --out-metadata {output.private_metadata} \
          --restrict &>> {log}
        """


rule summarize_publish:
    input:
        GISAID_meta_master = rules.publish_gisaid_master_metadata.output.metadata,
        COG_meta_master = rules.publish_COG_master_metadata.output.metadata,

        COG_seq_all = rules.publish_unaligned_cog_sequences.output.fasta,
        COG_seq_all_aligned = rules.publish_full_aligned_cog_data.output.fasta,
        COG_meta_all_aligned = rules.publish_full_aligned_cog_data.output.metadata,

        COG_seq_all_aligned_filtered = rules.publish_filtered_aligned_cog_data.output.fasta,
        COG_meta_all_aligned_filtered = rules.publish_filtered_aligned_cog_data.output.metadata,

        COG_GISAID_nexus_tree = rules.publish_full_annotated_tree_and_metadata.output.annotated_tree,
        COG_GISAID_meta = rules.publish_full_annotated_tree_and_metadata.output.metadata,

        public_COG_GISAID_newick_tree = rules.publish_public_cog_data.output.public_tree,
        public_COG_GISAID_seq_all = rules.publish_unaligned_cog_sequences.output.fasta,
        public_COG_meta = rules.publish_public_cog_data.output.metadata,

        microreact_tree = rules.publish_microreact_specific_output.output.newick_tree,
        microreact_public_metadata = rules.publish_microreact_specific_output.output.public_metadata,
        microreact_private_metadata = rules.publish_microreact_specific_output.output.private_metadata,

        lineage_report_fasta = rules.publish_cog_gisaid_data_for_lineage_release_work.output.fasta,
        lineage_report_metadata = rules.publish_cog_gisaid_data_for_lineage_release_work.output.csv,
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
        echo "> Newick tree for microreact published to {input.microreact_tree}\\n" >> {log}
        echo "> Public metadata for microreact published to {input.microreact_public_metadata}\\n" >> {log}
        echo "> Private metadata for microreact published to {input.microreact_private_metadata}\\n" >> {log}
        echo "> \\n" >> {log}

        echo '{{"text":"' > 6_data.json
        echo "*Step 6: publish data complete*\\n" >> 6_data.json
        cat {log} >> 6_data.json
        echo '"}}' >> 6_data.json
        echo 'webhook {params.webhook}'
        curl -X POST -H "Content-type: application/json" -d @6_data.json {params.webhook}
        """










#
