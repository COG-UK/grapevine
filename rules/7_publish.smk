


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
          --where-column lineage_support=UFbootstrap \
          --out-metadata {output.metadata_temp2} &>> {log}

        fastafunk add_columns \
          --in-metadata {output.metadata_temp2} \
          --in-data {input.new_global_lineages} \
          --index-column sequence_name \
          --join-on taxon  \
          --new-columns lineage lineage_support lineages_version \
          --where-column lineage_support=UFbootstrap \
          --out-metadata {output.metadata} &>> {log}
        """


rule publish_COG_master_metadata:
    input:
        fasta = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.trimmed.low_covg_filtered.fasta",
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
                          sequencing_org_code \
          --where-column epi_week=edin_epi_week country=adm0 \
                         sample_date=received_date sample_date=collection_date \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata_report} \
          --restrict &>> {log}
        """


rule gisaid_add_lineage_information_back_to_master_metadata:
    input:
        metadata = config["output_path"] + "/0/gisaid.all.csv",
        global_lineage_data = config["global_lineages"]
    output:
        metadata = config["output_path"] + "/7/gisaid.master.csv"
    log:
        config["output_path"] + "/logs/7_gisaid_add_lineage_information_back_to_master_metadata.log"
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
        config["output_path"] + "/logs/7_publish_gisaid_master_metadata.log"
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
        config["output_path"] + "/logs/7_publish_filtered_aligned_cog_data.log"
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
        gisaid_fasta = config["output_path"] + "/0/gisaid.RD.UH.filt.mapped.filt2.masked.filt3.fasta",
        gisaid_metadata = rules.gisaid_add_lineage_information_back_to_master_metadata.output.metadata
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
          --filter-column covv_accession_id \
                          sequence_name sample_date epi_week \
                          country adm1 adm2 submission_org_code \
                          is_surveillance is_community is_hcw \
                          is_travel_history travel_history lineage \
                          lineage_support uk_lineage microreact_lineage acc_lineage del_lineage acc_introduction del_introduction phylotype d614g \
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
                          lineage_support uk_lineage microreact_lineage acc_lineage del_lineage acc_introduction del_introduction phylotype d614g \
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
          --filter-column sequence_name sample_date epi_week \
                          country adm1 adm2 outer_postcode \
                          is_surveillance is_community is_hcw \
                          is_travel_history travel_history lineage \
                          lineage_support uk_lineage acc_lineage del_lineage acc_introduction del_introduction phylotype \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --restrict &>> {log}
        """


rule publish_uk_lineage_specific_fasta_and_metadata_files:
    input:
        tree = config["export_path"] + "/trees/uk_lineages/uk_lineage_UK{i}.tree",
        metadata = rules.publish_full_annotated_tree_and_metadata.output.metadata,
        fasta = rules.publish_full_annotated_tree_and_metadata.output.fasta
    params:
        temp_fasta = temp(config["output_path"] + "/7/uk_lineage_UK{i}.fasta")
    output:
        fasta = config["export_path"] + "/trees/uk_lineages/uk_lineage_UK{i}.fasta",
        metadata = config["export_path"] + "/trees/uk_lineages/uk_lineage_UK{i}.csv"
    log:
        config["output_path"] + "/logs/7_publish_uk_lineage_specific_fasta_and_metadata_files_uk{i}.log"
    shell:
        """
        fastafunk extract \
          --in-fasta {input.fasta} \
          --in-tree {input.tree} \
          --out-fasta {output.fasta} &>> {log}

          fastafunk fetch \
          --in-fasta {output.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name sample_date epi_week \
                          country adm1 adm2 outer_postcode \
                          is_surveillance is_community is_hcw \
                          is_travel_history travel_history lineage \
                          lineage_support uk_lineage acc_lineage del_lineage acc_introduction del_introduction phylotype \
          --out-fasta {params.temp_fasta} \
          --out-metadata {output.metadata} \
          --restrict &>> {log}
        """


LIN,X = glob_wildcards(config["output_path"] + "/5/{lineage}/trees/uk_lineage_UK{i}.tree")

rule summarize_publish_uk_lineage_specific_fasta_and_metadata_files:
    input:
        expand(config["export_path"] + "/trees/uk_lineages/uk_lineage_UK{i}.csv", i = X)
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
        fasta = config["export_path"] + "/lineage_release/cog_gisaid.fasta",
        metadata = config["export_path"] + "/lineage_release/cog_gisaid.csv",
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
        metadata = rules.uk_add_lineage_information_back_to_master_metadata.output.metadata
    output:
        public_tree = config["export_path"] + "/public/cog_global_" + config["date"] + "_tree.newick",
        fasta = config["export_path"] + "/public/cog_" + config["date"] + ".fasta",
        metadata = config["export_path"] + "/public/cog_" + config["date"] + "_metadata.csv"
    log:
        config["output_path"] + "/logs/7_publish_public_cog_data.log"
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
                          lineage_support uk_lineage primary_uk_lineage d614g \
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
                          lineage lineage_support uk_lineage primary_uk_lineage d614g \
          --where-column primary_uk_lineage=microreact_lineage \
          --out-fasta {output.fasta2} \
          --out-metadata {output.private_metadata} \
          --restrict &>> {log}
        """


rule publish_time_trees:
    input:
        config["output_path"] + "/logs/6_summarize_treetime.log"
    params:
        treedir = config["output_path"] + "/6/",
        outdir = config["export_path"] + "/trees/uk_lineages/"
    log:
        config["output_path"] + "/logs/7_publish_time_trees.log"
    shell:
        """
        for DIR in {params.treedir}*/trees/*timetree/
        do
            FILECOUNT=$(ls $DIR | wc -l)
            if [ $FILECOUNT -gt 0 ]
            then
                cp -r $DIR {params.outdir}
            fi
        done &>> {log}
        """


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
    conda: "/cephfs/covid/bham/climb-covid19-jacksonb/git/report/environment.yml"
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
            -V fontsize=10pt \
            --template={params.template} \
            --latex-engine=pdflatex \
            -o {output.outdir}UK_report.pdf &>> {log}

        mkdir -p {params.publish_dir} &>> {log}
        cp {output.outdir}UK_report.pdf {output.outdir}UK_report.md {params.publish_dir} &>> {log}
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
    conda: "/cephfs/covid/bham/climb-covid19-jacksonb/git/report/environment.yml"
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
            -V fontsize=10pt \
            --template={params.template} \
            --latex-engine=pdflatex \
            -o {output.outdir}{wildcards.adm1}.pdf &>> {log}

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
        sc = "[A-Z]{4}"
    output:
        outdir = directory(config["output_path"] + "/7/regional_reports/{sc}/"),
    log:
        config["output_path"] + "/logs/7_publish_report_{sc}.log",
    conda: "/cephfs/covid/bham/climb-covid19-jacksonb/git/report/environment.yml"
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
            -V fontsize=10pt \
            --template={params.template} \
            --latex-engine=pdflatex \
            -o {output.outdir}report_{wildcards.sc}.pdf &>> {log}

        mkdir -p {params.publish_dir}results/results_{wildcards.sc} &>> {log}
        cp {output.outdir}report_{wildcards.sc}.pdf {output.outdir}report_{wildcards.sc}.md {params.publish_dir} &>> {log}
        cp -r {output.outdir}/figures/ {params.publish_dir}results/results_{wildcards.sc}/ &>> {log}
        cp -r {output.outdir}summary_files/ {params.publish_dir}results/results_{wildcards.sc}/ &>> {log}
        """


ADM1 = ['England', 'Scotland', 'Wales', 'Northern_Ireland']
SC = ['LIVE', 'PHWC', 'CAMB', 'NORW', 'GLAS', 'EDIN', 'SHEF', 'EXET', 'NOTT', 'PORT', 'OXON', 'NORT', 'NIRE', 'LOND', 'SANG', 'BIRM', 'PHEC']
REPORTS = ['full'] + ADM1 + SC

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

        GISAID_meta_master = rules.publish_gisaid_master_metadata.output.metadata,
        COG_meta_master = rules.publish_COG_master_metadata.output.metadata_master,

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

        microreact_public_tree = rules.publish_microreact_specific_output.output.public_tree,
        microreact_private_tree = rules.publish_microreact_specific_output.output.private_tree,
        microreact_public_metadata = rules.publish_microreact_specific_output.output.public_metadata,
        microreact_private_metadata = rules.publish_microreact_specific_output.output.private_metadata,

        lineage_report_fasta = rules.publish_cog_gisaid_data_for_lineage_release_work.output.fasta,
        lineage_report_metadata = rules.publish_cog_gisaid_data_for_lineage_release_work.output.metadata,

        uk_lineage_fasta_csv_summary = rules.summarize_publish_uk_lineage_specific_fasta_and_metadata_files.log,

        log_uk_lineage_timetrees = rules.publish_time_trees.log
    params:
        webhook = config["webhook"],
        uk_trees_path = config["export_path"] + "/trees/uk_lineages/",
        reports_path = config["export_path"] + "/reports/",
    log:
        config["output_path"] + "/logs/7_summarize_publish.log"
    shell:
        """
        echo "> Reports published to {params.reports_path}\\n" >> {log}
        echo "> \\n" >> {log}
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
        echo "> UK lineage timetrees published in {params.uk_trees_path}\\n" >> {log}
        echo "> \\n" >> {log}
        echo "> Public tree published to {input.public_COG_GISAID_newick_tree}\\n" >> {log}
        echo "> Associated unaligned sequences published to {input.public_COG_GISAID_seq_all}\\n" >> {log}
        echo "> Matching metadata with public fields only published to {input.public_COG_meta}\\n" >> {log}
        echo "> \\n" >> {log}
        echo "> Public tree for microreact published to {input.microreact_public_tree}\\n" >> {log}
        echo "> Public metadata for microreact published to {input.microreact_public_metadata}\\n" >> {log}
        echo "> Private metadata for microreact published to {input.microreact_private_metadata}\\n" >> {log}
        echo "> \\n" >> {log}

        echo '{{"text":"' > 7_data.json
        echo "*Step 7: publish data complete*\\n" >> 7_data.json
        cat {log} >> 7_data.json
        echo '"}}' >> 7_data.json
        echo 'webhook {params.webhook}'
        curl -X POST -H "Content-type: application/json" -d @7_data.json {params.webhook}
        """










#
