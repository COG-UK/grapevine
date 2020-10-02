


rule uk_add_lineage_information_back_to_master_metadata:
    input:
        metadata = config["output_path"] + "/2/uk.with_new_lineages.csv",
        uk_lineage_data = config["output_path"] + "/5/cog_gisaid.lineages.with_all_traits.with_phylotype_traits.csv",
        global_lineage_data = config["global_lineages"],
        new_global_lineages = config["output_path"] + "/2/normal_pangolin/lineage_report.csv"
    output:
        metadata_temp1 = temp(config["output_path"] + "/7/temp1.uk.master.csv"),
        metadata_temp2 = temp(config["output_path"] + "/7/temp2.uk.master.csv"),
        metadata = config["output_path"] + "/7/uk.master.csv"
    log:
        config["output_path"] + "/logs/7_uk_add_lineage_information_back_to_master_metadata.log"
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.uk_lineage_data} \
          --index-column sequence_name \
          --join-on sequence_name \
          --new-columns uk_lineage microreact_lineage acc_lineage del_lineage acc_introduction del_introduction phylotype \
          --out-metadata {output.metadata_temp1} &> {log}

        fastafunk add_columns \
          --in-metadata {output.metadata_temp1} \
          --in-data {input.global_lineage_data} \
          --index-column sequence_name \
          --join-on taxon  \
          --new-columns lineage lineage_support lineages_version \
          --where-column lineage_support=probability lineages_version=pangoLEARN_version \
          --out-metadata {output.metadata_temp2} &>> {log}

        fastafunk add_columns \
          --in-metadata {output.metadata_temp2} \
          --in-data {input.new_global_lineages} \
          --index-column sequence_name \
          --join-on taxon  \
          --new-columns lineage lineage_support lineages_version \
          --where-column lineage_support=probability lineages_version=pangoLEARN_version \
          --out-metadata {output.metadata} &>> {log}
        """


rule publish_COG_master_metadata:
    input:
        fasta = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.trimmed.low_covg_filtered.omissions_filtered.fasta",
        metadata = rules.uk_add_lineage_information_back_to_master_metadata.output.metadata
    output:
        metadata_master = config["publish_path"] + "/COG/master.csv",
        metadata_report = config["publish_path"] + "/COG/report_metadata.csv",
        fasta = temp(config["output_path"] + "/7/7_publish_COG_master_metadata.temp.fasta")
    log:
        config["output_path"] + "/logs/7_publish_COG_master_metadata.log"
    shell:
        """
        cp {input.metadata} {output.metadata_master} &> {log}

        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name sample_date epi_week \
                          country adm1 adm2 \
                          lineage lineage_support lineages_version uk_lineage acc_lineage del_lineage phylotype acc_introduction del_introduction \
                          sequencing_org_code sequencing_submission_date \
          --where-column epi_week=edin_epi_week country=adm0 \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata_report} \
          --restrict &>> {log}
        """


# rule gisaid_add_lineage_information_back_to_master_metadata:
#     input:
#         metadata = config["GISAID_background_metadata"],
#         global_lineage_data = config["global_lineages"]
#     output:
#         metadata = config["output_path"] + "/7/gisaid.master.csv"
#     log:
#         config["output_path"] + "/logs/7_gisaid_add_lineage_information_back_to_master_metadata.log"
#     shell:
#         """
#       fastafunk add_columns \
#         --in-metadata {input.metadata} \
#         --in-data {input.global_lineage_data} \
#         --index-column sequence_name \
#         --join-on taxon  \
#         --new-columns lineage lineage_support lineages_version \
#         --where-column lineage_support=UFbootstrap \
#         --out-metadata {output.metadata} &>> {log}
#         """


# rule publish_gisaid_master_metadata:
#     input:
#         metadata = rules.gisaid_add_lineage_information_back_to_master_metadata.output.metadata,
#     output:
#         metadata = config["publish_path"] + "/GISAID/master.csv"
#     log:
#         config["output_path"] + "/logs/7_publish_gisaid_master_metadata.log"
#     shell:
#         """
#         cp {input.metadata} {output.metadata} &> {log}
#         """


rule publish_unaligned_cog_sequences:
    input:
        fasta = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated.unify_headers.fasta",
    output:
        fasta = config["export_path"] + "/alignments/cog_" + config["date"] + '_all.fasta'
    log:
        config["output_path"] + "/logs/7_publish_unaligned_cog_sequences.log"
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
        config["output_path"] + "/logs/7_publish_full_aligned_cog_data.log"
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column \
              adm0 adm1 adm2 adm2_private biosample_source_id \
              central_sample_id collected_by collection_date end_time \
              flowcell_id flowcell_type instrument_make instrument_model \
              layout_insert_length layout_read_length \
              library_adaptor_barcode library_layout_config library_name library_primers library_protocol \
              library_selection library_seq_kit library_seq_protocol library_source library_strategy \
              meta.artic.primers meta.artic.protocol meta.epi.cluster meta.investigation.cluster \
              meta.investigation.name meta.investigation.site metric.ct.1.ct_value metric.ct.1.test_kit \
              metric.ct.1.test_platform metric.ct.1.test_target metric.ct.2.ct_value metric.ct.2.test_kit         \
              metric.ct.2.test_platform metric.ct.2.test_target metric.ct.max_ct \
              metric.ct.min_ct metric.ct.num_tests \
              published_as received_date root_sample_id run_group run_name \
              sample_type_collected sample_type_received secondary_accession secondary_identifier \
              sequencing_org sequencing_org_code sequencing_submission_date sequencing_uuid \
              source_age source_sex start_time \
              submission_org submission_org_code submission_user swab_site \
              header sequence_name length missing gaps cov_id sample_date subsample_omit edin_epi_week d614g n439k p323l del_1605_3 \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --restrict &> {log}
        """


rule publish_filtered_aligned_cog_data:
    input:
        fasta = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.trimmed.low_covg_filtered.omissions_filtered.fasta",
        metadata = rules.uk_add_lineage_information_back_to_master_metadata.output.metadata
    output:
        fasta = config["export_path"] + "/alignments/cog_" + config["date"] + '_alignment.fasta',
        metadata = config["export_path"] + "/alignments/cog_" + config["date"] + '_metadata.csv'
    log:
        config["output_path"] + "/logs/7_publish_filtered_aligned_cog_data.log"
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name secondary_identifier sample_date epi_week \
                          country adm1 adm2 outer_postcode \
                          is_surveillance is_community is_hcw \
                          is_travel_history travel_history lineage \
                          lineage_support uk_lineage acc_lineage del_lineage phylotype \
          --where-column epi_week=edin_epi_week country=adm0 \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --restrict &> {log}
        """


rule combine_cog_gisaid:
    input:
        cog_fasta = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.trimmed.low_covg_filtered.omissions_filtered.fasta",
        cog_metadata = rules.uk_add_lineage_information_back_to_master_metadata.output.metadata,
        gisaid_fasta = config["GISAID_background_fasta"],
        gisaid_metadata = config["GISAID_background_metadata"]
    params:
        intermediate_cog_fasta = config["output_path"] + "/7/cog.combine_cog_gisaid.temp.cog.fasta",
        intermediate_cog_metadata = config["output_path"] + "/7/cog.combine_cog_gisaid.temp.cog.csv",
        intermediate_gisaid_fasta = config["output_path"] + "/7/gisaid.combine_cog_gisaid.temp.gisaid.fasta",
        intermediate_gisaid_metadata = config["output_path"] + "/7/gisaid.combine_cog_gisaid.temp.gisaid.csv",
    output:
        fasta = config["output_path"] + "/7/gisaid.combine_cog_gisaid.combined.fasta",
        metadata = config["output_path"] + "/7/gisaid.combine_cog_gisaid.combined.csv",
    log:
        config["output_path"] + "/logs/7_combine_cog_gisaid.log"
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.cog_fasta} \
          --in-metadata {input.cog_metadata} \
          --index-column sequence_name \
          --filter-column covv_accession_id central_sample_id biosample_source_id secondary_identifier \
                          sequence_name sample_date epi_week \
                          country adm1 adm2 submission_org_code \
                          is_surveillance is_community is_hcw \
                          is_travel_history travel_history lineage lineage_support lineages_version \
                          uk_lineage microreact_lineage acc_lineage del_lineage acc_introduction del_introduction phylotype d614g n439k p323l del_1605_3 \
          --where-column epi_week=edin_epi_week country=adm0 \
          --out-fasta {params.intermediate_cog_fasta} \
          --out-metadata {params.intermediate_cog_metadata} \
          --restrict &>> {log}


        fastafunk fetch \
          --in-fasta {input.gisaid_fasta} \
          --in-metadata {input.gisaid_metadata} \
          --index-column sequence_name \
          --filter-column covv_accession_id central_sample_id biosample_source_id secondary_identifier \
                          sequence_name sample_date epi_week \
                          country adm1 adm2 submission_org_code \
                          is_surveillance is_community is_hcw \
                          is_travel_history travel_history lineage lineage_support lineages_version \
                          uk_lineage microreact_lineage acc_lineage del_lineage acc_introduction del_introduction phylotype d614g n439k p323l del_1605_3 \
          --where-column adm1=edin_admin_1 travel_history=edin_travel \
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


rule publish_updated_global_lineages:
    input:
        combined_fasta = rules.combine_cog_gisaid.output.fasta,
        combined_metadata = rules.combine_cog_gisaid.output.metadata,
    output:
        fasta = config["publish_path"] + "/COG_GISAID/cog_gisaid.fasta",
        temp_metadata_1 = temp(config["publish_path"] + "/COG_GISAID/global_lineages_temp1.csv"),
        temp_metadata_2 = temp(config["publish_path"] + "/COG_GISAID/global_lineages_temp2.csv"),
        temp_metadata_3 = temp(config["publish_path"] + "/COG_GISAID/global_lineages_temp3.csv"),
        metadata = config["publish_path"] + "/COG_GISAID/global_lineages.csv",
    log:
        config["output_path"] + "/logs/7_publish_updated_global_lineages.log"
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.combined_fasta} \
          --in-metadata {input.combined_metadata} \
          --index-column sequence_name \
          --filter-column sequence_name lineage lineage_support lineages_version \
          --out-fasta {output.fasta} \
          --out-metadata {output.temp_metadata_1} \
          --restrict &>> {log}

        sed '1s/sequence_name/taxon/' {output.temp_metadata_1} > {output.temp_metadata_2}
        sed '1s/lineage_support/probability/' {output.temp_metadata_2} > {output.temp_metadata_3}
        sed '1s/lineages_version/pangoLEARN_version/' {output.temp_metadata_3} > {output.metadata}
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
        fasta = config["output_path"] + "/7/cog_global.fasta",
    log:
        config["output_path"] + "/logs/7_publish_full_annotated_tree_and_metadata.log"
    shell:
        """
        cp {input.annotated_tree} {output.annotated_tree} &> {log}
        cp {input.newick_tree} {output.newick_tree} &>> {log}


        fastafunk fetch \
          --in-fasta {input.combined_fasta} \
          --in-metadata {input.combined_metadata} \
          --index-column sequence_name \
          --filter-column sequence_name secondary_identifier sample_date epi_week \
                          country adm1 adm2 outer_postcode \
                          is_surveillance is_community is_hcw \
                          is_travel_history travel_history lineage \
                          lineage_support uk_lineage acc_lineage del_lineage acc_introduction del_introduction phylotype \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --restrict &>> {log}
        """


rule publish_civet_data:
    input:
        cog_all_fasta = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.alignment.full.masked.fasta",
        cog_fasta = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.trimmed.low_covg_filtered.omissions_filtered.fasta",
        cog_metadata = rules.uk_add_lineage_information_back_to_master_metadata.output.metadata,
        combined_metadata = rules.combine_cog_gisaid.output.metadata,
        combined_fasta = rules.combine_cog_gisaid.output.fasta,
        tree = config["output_path"] + "/5/cog_gisaid_full.tree.nexus",
    output:
        cog_all_fasta = config["export_path"] + "/civet/cog_" + config["date"] + '_alignment_all.fasta',
        cog_all_metadata = config["export_path"] + "/civet/cog_" + config["date"] + '_metadata_all.csv',
        cog_fasta = config["export_path"] + "/civet/cog_" + config["date"] + '_alignment.fasta',
        cog_metadata = config["export_path"] + "/civet/cog_" + config["date"] + '_metadata.csv',
        combined_metadata = config["export_path"] + "/civet/cog_global_" + config["date"] + '_metadata.csv',
        combined_fasta = config["export_path"] + "/civet/cog_global_" + config["date"] + '_alignment.fasta',
        tree = config["export_path"] + "/civet/cog_global_"  + config["date"] +  "_tree.nexus",
    log:
        config["output_path"] + "/logs/7_publish_civet_data.log"
    shell:
        """
        cp {input.tree} {output.tree} &>> {log}

        fastafunk fetch \
          --in-fasta {input.cog_all_fasta} \
          --in-metadata {input.cog_metadata} \
          --index-column sequence_name \
          --filter-column central_sample_id biosample_source_id sequence_name secondary_identifier \
                          sample_date epi_week \
                          country adm1 adm2 outer_postcode \
                          is_surveillance is_community is_hcw \
                          is_travel_history travel_history lineage \
                          lineage_support uk_lineage acc_lineage del_lineage phylotype \
          --where-column epi_week=edin_epi_week country=adm0 \
          --out-fasta {output.cog_all_fasta} \
          --out-metadata {output.cog_all_metadata} \
          --restrict &>> {log}

        fastafunk fetch \
          --in-fasta {input.cog_fasta} \
          --in-metadata {input.cog_metadata} \
          --index-column sequence_name \
          --filter-column central_sample_id biosample_source_id sequence_name secondary_identifier \
                          sample_date epi_week \
                          country adm1 adm2 outer_postcode \
                          is_surveillance is_community is_hcw \
                          is_travel_history travel_history lineage \
                          lineage_support uk_lineage acc_lineage del_lineage phylotype \
          --where-column epi_week=edin_epi_week country=adm0 \
          --out-fasta {output.cog_fasta} \
          --out-metadata {output.cog_metadata} \
          --restrict &>> {log}

        fastafunk fetch \
          --in-fasta {input.combined_fasta} \
          --in-metadata {input.combined_metadata} \
          --index-column sequence_name \
          --filter-column central_sample_id biosample_source_id sequence_name secondary_identifier \
                          sample_date epi_week \
                          country adm1 adm2 outer_postcode \
                          is_surveillance is_community is_hcw \
                          is_travel_history travel_history lineage \
                          lineage_support uk_lineage acc_lineage del_lineage phylotype \
          --where-column epi_week=edin_epi_week country=adm0 \
          --out-fasta {output.combined_fasta} \
          --out-metadata {output.combined_metadata} \
          --restrict &>> {log}
        """


rule publish_uk_lineage_specific_fasta_and_metadata_files:
    input:
        step_5_done = rules.summarize_define_uk_lineages_and_cut_out_trees.log,
        tree = config["output_path"] + "/5/trees/uk_lineage_UK{i}.tree",
        metadata = rules.publish_full_annotated_tree_and_metadata.output.metadata,
        fasta = rules.publish_full_annotated_tree_and_metadata.output.fasta
    params:
        temp_fasta = temp(config["output_path"] + "/7/uk_lineage_UK{i}.fasta")
    output:
        fasta = config["export_path"] + "/trees/uk_lineages/uk_lineage_UK{i}.fasta",
        metadata = config["export_path"] + "/trees/uk_lineages/uk_lineage_UK{i}.csv",
        tree = config["export_path"] + "/trees/uk_lineages/uk_lineage_UK{i}.newick"
    log:
        config["output_path"] + "/logs/7_publish_uk_lineage_specific_fasta_and_metadata_files_uk{i}.log"
    shell:
        """
        cp {input.tree} {output.tree}

        fastafunk extract \
          --in-fasta {input.fasta} \
          --in-tree {input.tree} \
          --out-fasta {output.fasta} &>> {log}

          fastafunk fetch \
          --in-fasta {output.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name secondary_identifier sample_date epi_week \
                          country adm1 adm2 outer_postcode \
                          is_surveillance is_community is_hcw \
                          is_travel_history travel_history lineage \
                          lineage_support uk_lineage acc_lineage del_lineage acc_introduction del_introduction phylotype \
          --out-fasta {params.temp_fasta} \
          --out-metadata {output.metadata} \
          --restrict &>> {log}
        """



def aggregate_input_publish_trees_logs(wildcards):
    checkpoint_output_directory = checkpoints.get_uk_lineage_samples.get(**wildcards).output[0]
    print(checkpoints.get_uk_lineage_samples.get(**wildcards).output[0])
    required_files = expand( "%s/logs/7_publish_uk_lineage_specific_fasta_and_metadata_files_uk{i}.log" %(config["output_path"]),
                            i=glob_wildcards(os.path.join(checkpoint_output_directory, "UK{i}.samples.txt")).i)
    return (sorted(required_files))

# X = glob_wildcards(config["output_path"] + "/5/trees/uk_lineage_UK{i}.tree").i
# csvs = expand(config["export_path"] + "/trees/uk_lineages/uk_lineage_UK{i}.csv", i = X)

rule summarize_publish_uk_lineage_specific_fasta_and_metadata_files:
    input:
        logs=aggregate_input_publish_trees_logs,
        step_5_done = rules.summarize_define_uk_lineages_and_cut_out_trees.log,
    log:
        config["output_path"] + "/logs/7_summarize_publish_uk_lineage_specific_fasta_and_metadata_files.log"
    shell:
        """
        echo "done" > {log}
        """


rule publish_cog_gisaid_data_for_lineage_release_work:
    input:
        combined_fasta = rules.combine_cog_gisaid.output.fasta,
        combined_metadata = rules.combine_cog_gisaid.output.metadata,
    output:
        fasta = config["publish_path"] + "/lineage_release/cog_gisaid.fasta",
        metadata = config["publish_path"] + "/lineage_release/cog_gisaid.csv",
    log:
        config["output_path"] + "/logs/7_publish_cog_gisaid_data_for_lineage_release_work.log"
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
        metadata = rules.uk_add_lineage_information_back_to_master_metadata.output.metadata,
        alignment = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.trimmed.low_covg_filtered.omissions_filtered.fasta",
        unmasked_alignment = rules.uk_get_unmasked_alignment.output.fasta,
    output:
        public_tree = config["export_path"] + "/public/cog_global_" + config["date"] + "_tree.newick",
        fasta = config["export_path"] + "/public/cog_" + config["date"] + "_all.fasta",
        metadata = config["export_path"] + "/public/cog_" + config["date"] + "_metadata.csv",
        alignment = config["export_path"] + "/public/cog_" + config["date"] + "_alignment.fasta",
        unmasked_alignment = config["export_path"] + "/public/cog_" + config["date"] + "_unmasked_alignment.fasta",
    log:
        config["output_path"] + "/logs/7_publish_public_cog_data.log"
    shell:
        """
        cp {input.public_tree} {output.public_tree} &> {log}
        cp {input.alignment} {output.alignment} &>> {log}
        cp {input.unmasked_alignment} {output.unmasked_alignment} &>> {log}

        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name country adm1 \
                          sample_date epi_week lineage \
                          lineage_support \
          --where-column epi_week=edin_epi_week country=adm0 \
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
        public_tree = config["export_path"] + "/microreact/cog_global_" + config["date"] + '_tree_public.newick',
        private_tree = config["export_path"] + "/microreact/cog_global_" + config["date"] + '_tree_private.newick',
        public_metadata = config["export_path"] + "/microreact/cog_global_" + config["date"] + '_metadata_public.csv',
        private_metadata = config["export_path"] + "/microreact/cog_global_" + config["date"] + '_metadata_private.csv',
        temp_public_metadata = temp(config["output_path"] + "/7/cog_global_microreact1.csv"),
        fasta1 = temp(config["output_path"] + "/7/cog_global_microreact1.fasta"),
        fasta2 = temp(config["output_path"] + "/7/cog_global_microreact2.fasta")
    log:
        config["output_path"] + "/logs/7_publish_microreact_specific_output.log"
    shell:
        """
        cp {input.newick_tree} {output.private_tree} &>> {log}

        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name sample_date epi_week \
                          country adm1 adm2 submission_org_code lineage \
                          lineage_support uk_lineage primary_uk_lineage d614g n439k p323l \
          --where-column primary_uk_lineage=microreact_lineage \
          --out-fasta {output.fasta1} \
          --out-metadata {output.temp_public_metadata} \
          --restrict &>> {log}

        datafunk anonymize_microreact \
          --input-tree {input.newick_tree} \
          --input-metadata {output.temp_public_metadata} \
          --output-tree {output.public_tree} \
          --output-metadata {output.public_metadata} &>> {log}

        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name sample_date epi_week \
                          country adm1 adm2 submission_org_code \
                          is_hcw travel_history \
                          lineage lineage_support uk_lineage primary_uk_lineage d614g n439k p323l \
          --where-column primary_uk_lineage=microreact_lineage \
          --out-fasta {output.fasta2} \
          --out-metadata {output.private_metadata} \
          --restrict &>> {log}
        """


# rule publish_time_trees:
#     input:
#         config["output_path"] + "/logs/6_summarize_treetime.log"
#     params:
#         treedir = config["output_path"] + "/6/",
#         outdir = config["export_path"] + "/trees/uk_lineages/"
#     log:
#         config["output_path"] + "/logs/7_publish_time_trees.log"
#     shell:
#         """
#         for DIR in {params.treedir}trees/*timetree/
#         do
#             FILECOUNT=$(ls $DIR | wc -l)
#             if [ $FILECOUNT -gt 0 ]
#             then
#                 cp -r $DIR {params.outdir}
#             fi
#         done &>> {log}
#         """


rule publish_full_report:
    input:
        metadata = rules.publish_COG_master_metadata.output.metadata_report,
    params:
        date = config["date"],
        template = config["latex_template"],
        publish_dir = config["export_path"] + "/reports/"
    output:
        outdir = directory(config["output_path"] + "/7/full_report/"),
    log:
        config["output_path"] + "/logs/7_publish_report_full.log"
    shell:
        """
        generate_report --m {input.metadata} \
            --w {params.date} \
            --s UK_report \
            --od {output.outdir} &>>{log}

        pandoc {output.outdir}UK_report.md \
            -V linkcolor:blue \
            -V geometry:a4paper \
            -V geometry:margin=2cm \
            -V mainfont="Helvetica Neue" \
            -V monofont="Helvetica Neue" \
            -V fontsize=9pt \
            --template={params.template} \
            --latex-engine=pdflatex \
            -o {output.outdir}UK_report.pdf &>> {log}

        sed -i.bak "s/(.*\/figures/(figures/g" {output.outdir}UK_report.md 2>> {log}

        mkdir -p {params.publish_dir} &>> {log}
        cp {output.outdir}UK_report.pdf {output.outdir}UK_report.md {params.publish_dir} &>> {log}
        cp -r {output.outdir}figures/ {params.publish_dir} &>> {log}
        cp -r {output.outdir}summary_files/ {params.publish_dir} &>> {log}
        """


rule publish_full_pillar_2_report:
    input:
        metadata = rules.publish_COG_master_metadata.output.metadata_report,
    params:
        date = config["date"],
        template = config["latex_template"],
        publish_dir = config["export_path"] + "/reports/full_report_pillar_2/"
    output:
        outdir = directory(config["output_path"] + "/7/full_report_pillar_2/"),
    log:
        config["output_path"] + "/logs/7_publish_report_pillar_2_full.log"
    shell:
        """
        generate_report --m {input.metadata} \
            --w {params.date} \
            --s UK_report \
            --od {output.outdir} \
            --p2 True &>> {log}

        pandoc {output.outdir}UK_report.md \
            -V linkcolor:blue \
            -V geometry:a4paper \
            -V geometry:margin=2cm \
            -V mainfont="Helvetica Neue" \
            -V monofont="Helvetica Neue" \
            -V fontsize=9pt \
            --template={params.template} \
            --latex-engine=pdflatex \
            -o {output.outdir}UK_report_pillar_2.pdf &>> {log}

        sed -i.bak "s/(.*\/figures/(figures/g" {output.outdir}UK_report.md 2>> {log}

        mkdir -p {params.publish_dir} &>> {log}
        cp {output.outdir}UK_report_pillar_2.pdf {params.publish_dir} &>> {log}
        cp {output.outdir}UK_report.md {params.publish_dir}UK_report_pillar_2.md &>> {log}
        cp -r {output.outdir}figures/ {params.publish_dir} &>> {log}
        cp -r {output.outdir}summary_files/ {params.publish_dir} &>> {log}
        """


rule publish_adm1_reports:
    input:
        metadata = rules.publish_COG_master_metadata.output.metadata_report,
    params:
        date = config["date"],
        template = config["latex_template"],
        publish_dir = config["export_path"] + "/reports/adm1_reports/"
    wildcard_constraints:
        adm1 = "[A-Z][a-z].*"
    output:
        outdir = directory(config["output_path"] + "/7/adm1_reports/{adm1}/"),
    log:
        config["output_path"] + "/logs/7_publish_report_{adm1}.log",
    shell:
        """
        generate_report --m {input.metadata} \
            --w {params.date} \
            --s {wildcards.adm1} \
            --adm {wildcards.adm1} \
            --od {output.outdir} &>> {log}

        pandoc {output.outdir}{wildcards.adm1}.md \
            -V linkcolor:blue \
            -V geometry:a4paper \
            -V geometry:margin=2cm \
            -V mainfont="Helvetica Neue" \
            -V monofont="Helvetica Neue" \
            -V fontsize=9pt \
            --template={params.template} \
            --latex-engine=pdflatex \
            -o {output.outdir}{wildcards.adm1}.pdf &>> {log}

        sed -i.bak "s/(.*\/figures/(figures/g" {output.outdir}{wildcards.adm1}.md 2>> {log}

        mkdir -p {params.publish_dir}{wildcards.adm1} &>> {log}
        cp {output.outdir}{wildcards.adm1}.pdf {output.outdir}{wildcards.adm1}.md {params.publish_dir}{wildcards.adm1}/ &>> {log}
        cp -r {output.outdir}figures/ {params.publish_dir}{wildcards.adm1}/ &>> {log}
        cp -r {output.outdir}summary_files/ {params.publish_dir}{wildcards.adm1}/ &>> {log}
        """


rule publish_sc_reports:
    input:
        metadata = rules.publish_COG_master_metadata.output.metadata_report,
    params:
        date = config["date"],
        template = config["latex_template"],
        publish_dir = config["export_path"] + "/reports/regional_reports/"
    wildcard_constraints:
        sc = "[A-Z]{4}|[A-Z]{4}_[A-Z]{4}"
    output:
        outdir = directory(config["output_path"] + "/7/regional_reports/{sc}/"),
    log:
        config["output_path"] + "/logs/7_publish_report_{sc}.log",
    shell:
        """
        generate_report --m {input.metadata} \
            --w {params.date} \
            --s report_{wildcards.sc} \
            --sc {wildcards.sc} \
            --od {output.outdir} &>> {log}

        pandoc {output.outdir}report_{wildcards.sc}.md \
            -V linkcolor:blue \
            -V geometry:a4paper \
            -V geometry:margin=2cm \
            -V mainfont="Helvetica Neue" \
            -V monofont="Helvetica Neue" \
            -V fontsize=9pt \
            --template={params.template} \
            --latex-engine=pdflatex \
            -o {output.outdir}report_{wildcards.sc}.pdf &>> {log}

        sed -i.bak "s/(.*\/figures/(figures/g" {output.outdir}report_{wildcards.sc}.md 2>> {log}

        mkdir -p {params.publish_dir}results/results_{wildcards.sc} &>> {log}
        cp {output.outdir}report_{wildcards.sc}.pdf {output.outdir}report_{wildcards.sc}.md {params.publish_dir} &>> {log}
        cp -r {output.outdir}/figures/ {params.publish_dir}results/results_{wildcards.sc}/ &>> {log}
        cp -r {output.outdir}summary_files/ {params.publish_dir}results/results_{wildcards.sc}/ &>> {log}
        """


ADM1 = ['England', 'Scotland', 'Wales', 'Northern_Ireland']
SC = ['BIRM','CAMB','EDIN','EXET','GLAS','GSTT','LIVE','LOND','LOND_BART','LOND_UCLH','NIRE','NORT','NORW','NOTT','OXON','PHEC','PHWC','PORT','SANG','SHEF']
REPORTS = ['full', 'pillar_2_full'] + ADM1 + SC

rule summarize_publish_reports:
    input:
        logs = expand(config["output_path"] + "/logs/7_publish_report_{X}.log", X = REPORTS)
    log:
        config["output_path"] + "/logs/7_summarize_publish_reports.log"
    shell:
        """
        echo "Done" > {log}
        """


rule summarize_publish:
    input:
        reports_log = rules.summarize_publish_reports.log,

        # GISAID_meta_master = rules.publish_gisaid_master_metadata.output.metadata,
        COG_meta_master = rules.publish_COG_master_metadata.output.metadata_master,

        COG_seq_all = rules.publish_unaligned_cog_sequences.output.fasta,
        COG_seq_all_aligned = rules.publish_full_aligned_cog_data.output.fasta,
        COG_meta_all_aligned = rules.publish_full_aligned_cog_data.output.metadata,

        COG_seq_all_aligned_filtered = rules.publish_filtered_aligned_cog_data.output.fasta,
        COG_meta_all_aligned_filtered = rules.publish_filtered_aligned_cog_data.output.metadata,

        COG_GISAID_nexus_tree = rules.publish_full_annotated_tree_and_metadata.output.annotated_tree,
        COG_GISAID_meta = rules.publish_full_annotated_tree_and_metadata.output.metadata,
        updated_global_lineages = rules.publish_updated_global_lineages.output.metadata,

        public_COG_GISAID_newick_tree = rules.publish_public_cog_data.output.public_tree,
        public_COG_GISAID_seq_all = rules.publish_unaligned_cog_sequences.output.fasta,
        public_COG_meta = rules.publish_public_cog_data.output.metadata,

        microreact_public_tree = rules.publish_microreact_specific_output.output.public_tree,
        microreact_private_tree = rules.publish_microreact_specific_output.output.private_tree,
        microreact_public_metadata = rules.publish_microreact_specific_output.output.public_metadata,
        microreact_private_metadata = rules.publish_microreact_specific_output.output.private_metadata,

        lineage_report_fasta = rules.publish_cog_gisaid_data_for_lineage_release_work.output.fasta,
        lineage_report_metadata = rules.publish_cog_gisaid_data_for_lineage_release_work.output.metadata,

        uk_lineage_fasta_csv_summary = rules.summarize_publish_uk_lineage_specific_fasta_and_metadata_files.log,

        # log_uk_lineage_timetrees = rules.publish_time_trees.log,

        log_civet = rules.publish_civet_data.log,
    params:
        date = config["date"],
        grapevine_webhook = config["grapevine_webhook"],
        export_path = config["export_path"],
        json_path = config["json_path"],
        uk_trees_path = config["export_path"] + "/trees/uk_lineages/",
        local_civet_path = config["export_path"] + "/civet/",
        reports_path = config["export_path"] + "/reports/",
        pillar_2_reports_path = config["export_path"] + "/reports/full_report_pillar_2/"
    log:
        config["output_path"] + "/logs/7_summarize_publish.log"
    shell:
        """
        echo "> Reports published to {params.reports_path}\\n" >> {log}
        echo "> Pillar_2 full report published to {params.pillar_2_reports_path}\\n" >> {log}
        echo "> \\n" >> {log}
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
        echo "> Public tree for microreact published to {input.microreact_public_tree}\\n" >> {log}
        echo "> Public metadata for microreact published to {input.microreact_public_metadata}\\n" >> {log}
        echo "> Private metadata for microreact published to {input.microreact_private_metadata}\\n" >> {log}
        echo "> \\n" >> {log}
        echo "> Data (there is an alignment too) for Pangolin lineage releases published to {input.lineage_report_metadata}\\n" >> {log}
        echo "> \\n" >> {log}
        echo "> Data for Civet published to {params.local_civet_path}\\n" >> {log}
        echo '{{"text":"' > {params.json_path}/7_data.json
        echo "*Step 7: publish data to {params.export_path} complete*\\n" >> {params.json_path}/7_data.json
        cat {log} >> {params.json_path}/7_data.json
        echo '"}}' >> {params.json_path}/7_data.json
        echo 'webhook {params.grapevine_webhook}'
        curl -X POST -H "Content-type: application/json" -d @{params.json_path}/7_data.json {params.grapevine_webhook}
        """
        # echo "> Gisaid master metadata published to {input.GISAID_meta_master}\\n" >> {log}


rule postpublish_cp_civet_data:
    input:
        cog_all_fasta = rules.publish_civet_data.output.cog_all_fasta,
        cog_all_metadata = rules.publish_civet_data.output.cog_all_metadata,
        cog_fasta = rules.publish_civet_data.output.cog_fasta,
        cog_metadata = rules.publish_civet_data.output.cog_metadata,
        combined_metadata = rules.publish_civet_data.output.combined_metadata,
        combined_fasta = rules.publish_civet_data.output.combined_fasta,
        tree = rules.publish_civet_data.output.tree,
    output:
        cog_all_fasta = "/cephfs/covid/bham/civet-cat/cog_alignment_all.fasta",
        cog_all_metadata = "/cephfs/covid/bham/civet-cat/cog_metadata_all.csv",
        cog_fasta = "/cephfs/covid/bham/civet-cat/cog_alignment.fasta",
        cog_metadata = "/cephfs/covid/bham/civet-cat/cog_metadata.csv",
        combined_metadata = "/cephfs/covid/bham/civet-cat/cog_global_metadata.csv",
        combined_fasta = "/cephfs/covid/bham/civet-cat/cog_global_alignment.fasta",
        tree = "/cephfs/covid/bham/civet-cat/cog_global_tree.nexus",
    log:
        config["output_path"] + "/logs/7_postpublish_cp_civet_data.log"
    shell:
        """
        cp {input.cog_all_fasta} {output.cog_all_fasta} &> {log}
        cp {input.cog_all_metadata} {output.cog_all_metadata} &>> {log}
        cp {input.cog_fasta} {output.cog_fasta} &>> {log}
        cp {input.cog_metadata} {output.cog_metadata} &>> {log}
        cp {input.combined_metadata} {output.combined_metadata} &>> {log}
        cp {input.combined_fasta} {output.combined_fasta} &>> {log}
        cp {input.tree} {output.tree} &>> {log}
        """

rule postpublish_rsync_phylogenetics_data:
    input:
        publishdone = rules.summarize_publish.log
    params:
        date = config["date"],
        parsed_date = config["date"].replace('-', ''),
        export_path = config["export_path"],
    log:
        config["output_path"] + "/logs/7_postpublish_rsync_phylogenetics_data.log"
    shell:
        """
        rsync -r {params.export_path}/ /cephfs/covid/bham/results/phylogenetics/{params.parsed_date}/
        ln -sfn /cephfs/covid/bham/results/phylogenetics/{params.parsed_date} /cephfs/covid/bham/results/phylogenetics/latest
        """

rule summarize_postpublish:
    input:
        civet_log = config["output_path"] + "/logs/7_postpublish_cp_civet_data.log",
        rsync_log = config["output_path"] + "/logs/7_postpublish_rsync_phylogenetics_data.log",
    params:
        date = config["date"],
        phylopipe_webhook = config["phylopipe_webhook"],
        phylopipe_token = config["phylopipe_token"],
        json_path = config["json_path"],
        reports_path = config["export_path"] + "/reports/",
        report_upload_log = config["output_path"] + "/logs/7_upload_report.log",
    log:
        config["output_path"] + "/logs/7_summarize_postpublish.log"
    shell:
        """
        echo "> Phylogenetics pipeline output published to \`/cephfs/covid/bham/results/phylogenetics/latest/\`\\n" >> {log}
        echo "> \\n" >> {log}
        echo "> Reports published to \`reports/\`\\n" >> {log}
        echo "> \\n" >> {log}
        echo "> Unaligned (deduplicated, clean headers) COG sequences published to \`alignments/cog_{params.date}_all.fasta\`\\n" >> {log}
        echo "> \\n" >> {log}
        echo "> Aligned (deduplicated, clean headers) COG sequences published to \`alignments/cog_{params.date}_all_alignment.fasta\`\\n" >> {log}
        echo "> Matching metadata published to \`alignments/cog_{params.date}_all_metadata.csv\`\\n" >> {log}
        echo "> \\n" >> {log}
        echo "> Filtered, aligned COG sequences published to \`alignments/cog_{params.date}_alignment.fasta\`\\n" >> {log}
        echo "> Matching metadata with lineage information published to \`alignments/cog_{params.date}_metadata.csv\`\\n" >> {log}
        echo "> \\n" >> {log}
        echo "> Full, annotated tree published to \`trees/cog_{params.date}_tree.nexus\`\\n" >> {log}
        echo "> Matching metadata published to \`trees/cog_global_{params.date}_metadata.csv\`\\n" >> {log}
        echo "> UK lineage subtrees published in \`trees/uk_lineages/\`\\n" >> {log}
        echo "> \\n" >> {log}
        echo "> Public tree published to \`public/cog_global_{params.date}_tree.newick\`\\n" >> {log}
        echo "> Associated unaligned sequences published to \`alignments/cog_{params.date}_all.fasta\`\\n" >> {log}
        echo "> Matching metadata with public fields only published to \`public/cog_{params.date}_metadata.csv\`\\n" >> {log}
        echo "> \\n" >> {log}
        echo "> Public tree for microreact published to \`microreact/cog_global_{params.date}_tree_public.newick\`\\n" >> {log}
        echo "> Public metadata for microreact published to \`microreact/cog_global_{params.date}_metadata_public.csv\`\\n" >> {log}
        echo "> Private metadata for microreact published to \`microreact/cog_global_{params.date}_metadata_private.csv\`\\n" >> {log}
        echo "> \\n" >> {log}
        echo "> Data for Civet published to \`/cephfs/covid/bham/civet-cat/\`\\n" >> {log}

        echo '{{"text":"' > {params.json_path}/publish_data_to_consortium.json
        echo "*Phylogenetic pipeline complete*\\n" >> {params.json_path}/publish_data_to_consortium.json
        cat {log} >> {params.json_path}/publish_data_to_consortium.json
        echo '"}}' >> {params.json_path}/publish_data_to_consortium.json
        echo 'webhook {params.phylopipe_webhook}'
        curl -X POST -H "Content-type: application/json" -d @{params.json_path}/publish_data_to_consortium.json {params.phylopipe_webhook}

        ln -sfn /cephfs/covid/bham/raccoon-dog/{params.date} /cephfs/covid/bham/raccoon-dog/previous

        curl -F file=@{params.reports_path}UK_report.pdf -F "initial_comment=UK report for the latest phylogenetics pipeline run" -F channels=CVACDLCKA -H "Authorization: Bearer {params.phylopipe_token}" https://slack.com/api/files.upload &> {params.report_upload_log}
        """
        # THREAD_TS=`grep -oh -E '\"ts\":\"[0-9]+\.[0-9]+\"' {params.report_upload_log} | sed 's/"ts":"//' | sed 's/"//'`
        #
        # curl -F file=@{params.reports_path}adm1_reports/England/England.pdf -F channels=CVACDLCKA thread_ts=${{THREAD_TS}} -H "Authorization: Bearer {params.phylopipe_token}" https://slack.com/api/files.upload &>> {params.report_upload_log}
        # curl -F file=@{params.reports_path}adm1_reports/Northern_Ireland/Northern_Ireland.pdf -F channels=CVACDLCKA thread_ts=${{THREAD_TS}} -H "Authorization: Bearer {params.phylopipe_token}" https://slack.com/api/files.upload &>> {params.report_upload_log}
        # curl -F file=@{params.reports_path}adm1_reports/Scotland/Scotland.pdf -F channels=CVACDLCKA thread_ts=${{THREAD_TS}} -H "Authorization: Bearer {params.phylopipe_token}" https://slack.com/api/files.upload &>> {params.report_upload_log}
        # curl -F file=@{params.reports_path}adm1_reports/Wales/Wales.pdf -F channels=CVACDLCKA thread_ts=${{THREAD_TS}} -H "Authorization: Bearer {params.phylopipe_token}" https://slack.com/api/files.upload &>> {params.report_upload_log}



#
