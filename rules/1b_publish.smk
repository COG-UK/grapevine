rule clean_and_publish_cog_geography:
    input:
        metadata = rules.uk_type_dels.output.metadata,
        fasta = rules.uk_trim_alignment.output.fasta,
    params:
        geog_script = os.path.join(workflow.current_basedir, "../utilities/geography_cleaning.py"),
        geog_util_dir = os.path.join(workflow.current_basedir, "../utilities/geography_utils"),
        outdir = config["output_path"] + "/1b/geography/",
    output:
        geography_metadata = config["output_path"] + "/1b/geography/geography.csv",
        geography_log = config["output_path"] + "/1b/geography/log_file.txt",
        locations = config["output_path"] + "/1b/geography/new_unclean_locations.csv",
        postcodes = config["output_path"] + "/1b/geography/new_unclean_postcodes.csv",
        bad_seqs = config["output_path"] + "/1b/geography/sequences_with_incompatible_locs.csv",
        junkfasta = temp(config["output_path"] + "/1b/geography/junk.fasta"),
        metadata_temp = config["output_path"] + "/1b/geography/metadata.junkcsv",
    resources: mem_per_cpu=20000
    log:
        config["output_path"] + "/logs/1b_clean_and_publish_cog_geography.log",
    shell:
        """
        mkdir -p {params.outdir}

        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column central_sample_id sequence_name sample_date epi_week \
                          adm0 adm1 adm2 adm2_private \
          --out-fasta {output.junkfasta} \
          --out-metadata {output.metadata_temp} \
          --restrict &>> {log}

        python {params.geog_script} \
          --metadata {output.metadata_temp} \
          --country-col adm0 \
          --adm1-col adm1 \
          --adm2-col adm2 \
          --outer-postcode-col adm2_private \
          --mapping-utils-dir {params.geog_util_dir} \
          --outdir {params.outdir} 2>> {log}
        """

rule publish_COG_master_metadata:
    input:
        fasta = rules.uk_trim_alignment.output.fasta,
        metadata = rules.uk_type_dels.output.metadata,
        geography_metadata = rules.clean_and_publish_cog_geography.output.geography_metadata,
    output:
        metadata_master = config["publish_path"] + "/COG/master.csv",
        metadata_report = config["publish_path"] + "/COG/report_metadata.csv",
        metadata_report_temp = config["output_path"] + "/1b/report_metadata_temp.csv",
        fasta = temp(config["output_path"] + "/1b/1b_publish_COG_master_metadata.temp.fasta")
    log:
        config["output_path"] + "/logs/1b_publish_COG_master_metadata.log"
    resources: mem_per_cpu=20000
    shell:
        """
        cp {input.metadata} {output.metadata_master} &> {log}

        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column central_sample_id sequence_name sample_date epi_week \
                          country adm1 adm2 \
                          lineage lineage_support lineages_version \
                          sequencing_org_code sequencing_submission_date \
          --where-column epi_week=edin_epi_week country=adm0 \
          --out-fasta {output.fasta} \
         --out-metadata {output.metadata_report_temp} \
          --restrict &>> {log}

        fastafunk add_columns \
          --in-metadata {output.metadata_report_temp} \
          --in-data {input.geography_metadata} \
          --index-column sequence_name \
          --join-on sequence_name \
          --new-columns adm2 \
          --out-metadata {output.metadata_report} &>> {log}
        """


rule publish_unaligned_cog_sequences:
    input:
        fasta = rules.uk_remove_duplicates_root_biosample_by_gaps.output.fasta,
    output:
        fasta = config["export_path"] + "/alignments/cog_" + config["date"] + '_all.fasta'
    log:
        config["output_path"] + "/logs/1b_publish_unaligned_cog_sequences.log"
    shell:
        """
        cp {input.fasta} {output.fasta} 2> {log}
        """


rule publish_full_aligned_cog_data:
    input:
        fasta = rules.uk_trim_alignment.output.fasta,
        metadata = rules.uk_type_dels.output.metadata
    output:
        fasta = config["export_path"] + "/alignments/cog_" + config["date"] + '_all_alignment.fasta',
        metadata = config["export_path"] + "/alignments/cog_" + config["date"] + '_all_metadata.csv'
    log:
        config["output_path"] + "/logs/1b_publish_full_aligned_cog_data.log"
    resources: mem_per_cpu=20000
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column \
              country adm1 adm2 outer_postcode biosample_source_id \
              central_sample_id collected_by collection_date end_time \
              flowcell_id flowcell_type instrument_make instrument_model is_surveillance \
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
              header sequence_name length missing gaps cov_id sample_date subsample_omit epi_week \
              d614g n439k p323l a222v y453f n501y t1001i p681h q27stop del_21765_6 del_1605_3 del_1605_3 \
          --where-column epi_week=edin_epi_week country=adm0 outer_postcode=adm2_private \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --restrict &> {log}
        """

rule publish_filtered_aligned_cog_data:
    input:
        fasta = rules.uk_trim_alignment.output.fasta,
        metadata = rules.uk_type_dels.output.metadata,
    output:
        fasta = config["export_path"] + "/alignments/cog_" + config["date"] + '_alignment.fasta',
        metadata = config["export_path"] + "/alignments/cog_" + config["date"] + '_metadata.csv'
    log:
        config["output_path"] + "/logs/1b_publish_filtered_aligned_cog_data.log"
    resources: mem_per_cpu=20000
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
                          lineage_support\
          --where-column epi_week=edin_epi_week country=adm0 outer_postcode=adm2_private \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --restrict &> {log}
        """


# TODO: split variants into AAs and synonymous SNPs (and indels...), etc.
rule combine_variants:
    input:
        COG_variants = rules.uk_get_variants.output.variants,
        COG_fasta = rules.uk_trim_alignment.output.fasta,
        GISAID_variants = config["GISAID_variants"],
        GISAID_fasta = config["GISAID_background_fasta"],
    output:
        variants = config["output_path"] + "/1b/combined_variants.csv",
    log:
        config["output_path"] + "/logs/1b_combine_variants.log"
    run:
        from Bio import SeqIO

        COG_fasta = SeqIO.index(str(input.COG_fasta), "fasta")
        GISAID_fasta = SeqIO.index(str(input.GISAID_fasta), "fasta")

        first = True
        with open(str(input.COG_variants), "r") as COG_variants_in, \
         open(str(input.GISAID_variants), "r") as GISAID_variants_in, \
         open(str(output.variants), "w") as variants_out:
            for line in GISAID_variants_in:
                if first:
                    variants_out.write(line)
                    first = False
                    continue

                sample = line.strip().split(",")[0]

                if sample in GISAID_fasta:
                    variants_out.write(line)

            first = True
            for line in COG_variants_in:
                if first:
                    first = False
                    continue

                sample = line.strip().split(",")[0]

                if sample in COG_fasta:
                    variants_out.write(line)


rule combine_cog_gisaid:
    input:
        cog_fasta = rules.uk_trim_alignment.output.fasta,
        cog_metadata = rules.uk_type_dels.output.metadata,
        gisaid_fasta = config["GISAID_background_fasta"],
        gisaid_metadata = config["GISAID_background_metadata"]
    params:
        intermediate_cog_fasta = config["output_path"] + "/1b/cog.combine_cog_gisaid.temp.cog.fasta",
        intermediate_cog_metadata = config["output_path"] + "/1b/cog.combine_cog_gisaid.temp.cog.csv",
        intermediate_gisaid_fasta = config["output_path"] + "/1b/gisaid.combine_cog_gisaid.temp.gisaid.fasta",
        intermediate_gisaid_metadata = config["output_path"] + "/1b/gisaid.combine_cog_gisaid.temp.gisaid.csv",
    output:
        fasta = config["output_path"] + "/1b/gisaid.combine_cog_gisaid.combined.fasta",
        metadata = config["output_path"] + "/1b/gisaid.combine_cog_gisaid.combined.csv",
    log:
        config["output_path"] + "/logs/1b_combine_cog_gisaid.log"
    resources: mem_per_cpu=20000
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.cog_fasta} \
          --in-metadata {input.cog_metadata} \
          --index-column sequence_name \
          --filter-column covv_accession_id central_sample_id biosample_source_id secondary_identifier \
                          sequence_name sample_date epi_week \
                          country adm1 adm2 outer_postcode adm2_raw adm2_source nuts1 region latitude longitude location \
                          submission_org_code is_surveillance is_community is_hcw \
                          is_travel_history travel_history \
                          lineage lineage_support lineages_version \
                          d614g n439k p323l a222v y453f n501y t1001i p681h q27stop del_21765_6 del_1605_3 del_1605_3 \
                          source_age source_sex sample_type_collected sample_type_received swab_site \
                          ct_n_ct_value ct_n_test_kit ct_n_test_platform ct_n_test_target \
          --where-column epi_week=edin_epi_week country=adm0 outer_postcode=adm2_private \
          --out-fasta {params.intermediate_cog_fasta} \
          --out-metadata {params.intermediate_cog_metadata} \
          --restrict &>> {log}

        fastafunk fetch \
          --in-fasta {input.gisaid_fasta} \
          --in-metadata {input.gisaid_metadata} \
          --index-column sequence_name \
          --filter-column covv_accession_id central_sample_id biosample_source_id secondary_identifier \
                          sequence_name sample_date epi_week \
                          country adm1 adm2 outer_postcode adm2_raw adm2_source nuts1 region latitude longitude location \
                          submission_org_code is_surveillance is_community is_hcw \
                          is_travel_history travel_history \
                          lineage lineage_support lineages_version \
                          d614g n439k p323l a222v y453f n501y t1001i p681h q27stop del_21765_6 del_1605_3 del_1605_3 \
                          source_age source_sex sample_type_collected sample_type_received swab_site \
                          ct_n_ct_value ct_n_test_kit ct_n_test_platform ct_n_test_target \
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


rule make_metadata_dir_outputs:
    input:
        combined_fasta = rules.combine_cog_gisaid.output.fasta,
        combined_metadata = rules.combine_cog_gisaid.output.metadata,
        clean_COG_geography = rules.clean_and_publish_cog_geography.output.geography_metadata,
        combined_variants = rules.combine_variants.output.variants,
    output:
        geography_metadata = config["export_path"] + "/metadata/cog_global_" + config["date"] + "_geography.csv",
        geography_metadata_temp = temp(config["output_path"] + "/1b/make_metadata_dir_geog.junkcsv1"),
        geography_metadata_temp2 = temp(config["output_path"] + "/1b/make_metadata_dir_geog.junkcsv2"),
        junkfasta3 = temp(config["output_path"] + "/1b/make_metadata_dir_consort.junkfasta3"),

        public_metadata = config["export_path"] + "/metadata/cog_global_" + config["date"] + "_public.csv",

        consortium_metadata = config["export_path"] + "/metadata/cog_global_" + config["date"] + "_consortium.csv",
        consortium_metadata_temp = temp(config["output_path"] + "/1b/make_metadata_dir_consort.junkcsv1"),
        consortium_metadata_temp2 = temp(config["output_path"] + "/1b/make_metadata_dir_consort.junkcsv2"),
        junkfasta1 = temp(config["output_path"] + "/1b/make_metadata_dir_consort.junkfasta1"),
        junkfasta2 = temp(config["output_path"] + "/1b/make_metadata_dir_consort.junkfasta2"),

        variants_metadata = config["export_path"] + "/metadata/cog_global_" + config["date"] + '_mutations.csv',
        variants_metadata_temp = temp(config["output_path"] + "/1b/make_metadata_dir_variants.junkcsv"),
        junkfasta4 = temp(config["output_path"] + "/1b/make_metadata_dir_variants.junkfasta4"),
    log:
        config["output_path"] + "/logs/1b_make_metadata_dir_outputs.log"
    resources: mem_per_cpu=20000
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.combined_fasta} \
          --in-metadata {input.combined_metadata} \
          --index-column sequence_name \
          --filter-column sequence_name sample_date \
          --out-fasta {output.junkfasta4} \
          --out-metadata {output.variants_metadata_temp} \
          --restrict &>> {log}

        fastafunk add_columns \
          --in-metadata {output.variants_metadata_temp} \
          --in-data {input.combined_variants} \
          --index-column sequence_name \
          --join-on query \
          --new-columns variants \
          --out-metadata {output.variants_metadata} 2>> {log}

        fastafunk fetch \
          --in-fasta {input.combined_fasta} \
          --in-metadata {input.combined_metadata} \
          --index-column sequence_name \
          --filter-column central_sample_id sequence_name sample_date epi_week \
                          country adm1 adm2 \
          --out-fasta {output.junkfasta3} \
          --out-metadata {output.geography_metadata_temp} \
          --restrict &>> {log}

        fastafunk add_columns \
          --in-metadata {output.geography_metadata_temp} \
          --in-data {input.clean_COG_geography} \
          --index-column sequence_name \
          --join-on sequence_name \
          --new-columns adm1 adm2 outer_postcode adm2_raw adm2_source nuts1 region latitude longitude location \
          --out-metadata {output.geography_metadata_temp2} 2>> {log}

        sed '1s/nuts1/NUTS1/' {output.geography_metadata_temp2} > {output.geography_metadata} 2>> {log}

        fastafunk fetch \
          --in-fasta {input.combined_fasta} \
          --in-metadata {input.combined_metadata} \
          --index-column sequence_name \
          --filter-column sequence_name sample_date epi_week \
                          country adm1 \
                          is_surveillance \
                          is_travel_history travel_history lineage \
                          lineage_support \
          --out-fasta {output.junkfasta1} \
          --out-metadata {output.public_metadata} \
          --restrict &>> {log}

        fastafunk fetch \
          --in-fasta {input.combined_fasta} \
          --in-metadata {input.combined_metadata} \
          --index-column sequence_name \
          --filter-column sequence_name gisaid_id sample_date epi_week submission_org_code \
                          country adm1 adm2 outer_postcode adm2_raw adm2_source nuts1 region latitude longitude location \
                          source_age source_sex sample_type_collected sample_type_received swab_site \
                          ct_n_ct_value ct_n_test_kit ct_n_test_platform ct_n_test_target \
                          is_surveillance \
                          is_travel_history travel_history \
                          lineage lineage_support \
                          d614g n439k p323l a222v y453f n501y t1001i p681h q27stop del_21765_6 del_1605_3 \
          --where-column gisaid_id=covv_accession_id \
          --out-fasta {output.junkfasta2} \
          --out-metadata {output.consortium_metadata_temp} \
          --restrict &>> {log}

        fastafunk add_columns \
          --in-metadata {output.consortium_metadata_temp} \
          --in-data {input.clean_COG_geography} \
          --index-column sequence_name \
          --join-on sequence_name \
          --new-columns adm1 adm2 adm2_raw adm2_source nuts1 outer_postcode region latitude longitude location \
          --out-metadata {output.consortium_metadata_temp2} &>> {log}

        sed '1s/nuts1/NUTS1/' {output.consortium_metadata_temp2} > {output.consortium_metadata} 2>> {log}
        """


rule MINK:
    input:
        combined_variants = rules.make_metadata_dir_outputs.output.variants_metadata,
        report_metadata = rules.publish_COG_master_metadata.output.metadata_report,
    params:
        date = config["date"],
        variant_list = config["MINK_variants"],
        AA_filter_script = os.path.join(workflow.current_basedir, "../utilities/get_AA_mutations_only.py"),
    output:
        COG_AAs_only = config["output_path"] + "/1b/COG_AA_variants.csv",
        MINK_results_dir = directory(config["output_path"] + "/1b/MINK_report/"),
        MINK_report = config["output_path"] + "/1b/MINK_report/MINK_variant_report.pdf",
    log:
        config["output_path"] + "/logs/1b_MINK.log"
    shell:
        """
        python {params.AA_filter_script} {input.combined_variants} {output.COG_AAs_only} 2> {log}

        mink \
          --metadata-file {input.report_metadata} \
          --snp-file {output.COG_AAs_only} \
          --date-data {params.date} \
          --snp-list $(cat {params.variant_list} | tr "\n" "," | sed "s/,$//") \
          --flag-fastest \
          --outdir {output.MINK_results_dir} 2>> {log}

        sed "s|figures/|{output.MINK_results_dir}figures/|g" {output.MINK_results_dir}/mutation_report.md > {output.MINK_results_dir}/mutation_report_abspaths.md 2>> {log}

        pandoc {output.MINK_results_dir}/mutation_report_abspaths.md \
            -V linkcolor:blue \
            -V geometry:a4paper \
            -V geometry:margin=2cm \
            -V mainfont="Helvetica Neue" \
            -V monofont="Helvetica Neue" \
            -V fontsize=10pt \
            --latex-engine=pdflatex \
            -o {output.MINK_report} 2>> {log}
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
        config["output_path"] + "/logs/1b_publish_updated_global_lineages.log"
    resources: mem_per_cpu=20000
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

rule publish_cog_gisaid_data_for_lineage_release_work:
    input:
        combined_fasta = rules.combine_cog_gisaid.output.fasta,
        combined_metadata = rules.combine_cog_gisaid.output.metadata,
    output:
        fasta = config["publish_path"] + "/lineage_release/cog_gisaid.fasta",
        metadata = config["publish_path"] + "/lineage_release/cog_gisaid.csv",
    log:
        config["output_path"] + "/logs/1b_publish_cog_gisaid_data_for_lineage_release_work.log"
    resources: mem_per_cpu=20000
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
        fasta = rules.uk_remove_duplicates_root_biosample_by_gaps.output.fasta,
        metadata = rules.uk_type_dels.output.metadata,
        unmasked_alignment = rules.uk_alignment.output.fasta,
        alignment = rules.uk_trim_alignment.output.fasta,
    output:
        fasta = config["export_path"] + "/public/cog_" + config["date"] + "_all.fasta",
        metadata = config["export_path"] + "/public/cog_" + config["date"] + "_metadata.csv",
        alignment = config["export_path"] + "/public/cog_" + config["date"] + "_alignment.fasta",
        unmasked_alignment = config["export_path"] + "/public/cog_" + config["date"] + "_unmasked_alignment.fasta",
    log:
        config["output_path"] + "/logs/1b_publish_public_cog_data.log"
    resources: mem_per_cpu=20000
    shell:
        """
        cp {input.alignment} {output.alignment} &>> {log}
        cp {input.unmasked_alignment} {output.unmasked_alignment} &>> {log}

        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name country adm1 \
                          sample_date epi_week lineage lineage_support \
                          d614g n439k p323l a222v y453f n501y t1001i p681h q27stop del_21765_6 del_1605_3  \
          --where-column epi_week=edin_epi_week country=adm0 \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --restrict &>> {log}
        """

rule summarize_publish:
    input:
        COG_meta_master = rules.publish_COG_master_metadata.output.metadata_master,
        COG_meta_report = rules.publish_COG_master_metadata.output.metadata_report,

        export_meta_variants = rules.make_metadata_dir_outputs.output.variants_metadata,
        export_meta_geog = rules.make_metadata_dir_outputs.output.geography_metadata,
        export_meta_public = rules.make_metadata_dir_outputs.output.public_metadata,
        export_meta_consort = rules.make_metadata_dir_outputs.output.consortium_metadata,

        COG_seq_all = rules.publish_unaligned_cog_sequences.output.fasta,
        COG_seq_all_aligned = rules.publish_full_aligned_cog_data.output.fasta,
        COG_meta_all_aligned = rules.publish_full_aligned_cog_data.output.metadata,

        COG_seq_all_aligned_filtered = rules.publish_filtered_aligned_cog_data.output.fasta,
        COG_meta_all_aligned_filtered = rules.publish_filtered_aligned_cog_data.output.metadata,

        updated_global_lineages = rules.publish_updated_global_lineages.output.metadata,

        public_COG_GISAID_seq_all = rules.publish_unaligned_cog_sequences.output.fasta,
        public_COG_meta = rules.publish_public_cog_data.output.metadata,

        lineage_report_fasta = rules.publish_cog_gisaid_data_for_lineage_release_work.output.fasta,
        lineage_report_metadata = rules.publish_cog_gisaid_data_for_lineage_release_work.output.metadata,

    params:
        date = config["date"],
        grapevine_webhook = config["grapevine_webhook"],
        export_path = config["export_path"],
        json_path = config["json_path"],
    log:
        config["output_path"] + "/logs/1b_summarize_publish.log"
    shell:
        """
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
        echo "> Data (there is an alignment too) for Pangolin lineage releases published to {input.lineage_report_metadata}\\n" >> {log}
        echo "> \\n" >> {log}
        
        mkdir -p {params.json_path}  
        echo '{{"text":"' > {params.json_path}/1b_data.json
        echo "*DATAPIPE: Publish data to {params.export_path} complete*\\n" >> {params.json_path}/1b_data.json
        cat {log} >> {params.json_path}/1b_data.json
        echo '"}}' >> {params.json_path}/1b_data.json
        echo 'webhook {params.grapevine_webhook}'
        curl -X POST -H "Content-type: application/json" -d @{params.json_path}/1b_data.json {params.grapevine_webhook}
        """

#rule postpublish_rsync_phylogenetics_data:
#    input:
#        publishdone = rules.summarize_publish.log
#    params:
#        date = config["date"],
#        parsed_date = config["date"].replace('-', ''),
#        export_path = config["export_path"],
#    log:
#        config["output_path"] + "/logs/1b_postpublish_rsync_phylogenetics_data.log"
#    shell:
#        """
#        rsync -r {params.export_path}/ /cephfs/covid/bham/results/datapipe/{params.parsed_date}/
#        ln -sfn /cephfs/covid/bham/results/datapipe/{params.parsed_date} /cephfs/covid/bham/results/datapipe/latest
#        """

