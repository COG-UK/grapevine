import pandas as pd
from Bio import SeqIO


rule uk_strip_header_digits:
    input:
        fasta = config["latest_uk_fasta"],
    output:
        fasta = config["output_path"] + "/1/uk_latest.headerstripped.fasta",
    log:
        config["output_path"] + "/logs/1_uk_strip_header_digits.log"
    run:
        fasta_in = SeqIO.parse(str(input.fasta), "fasta")
        with open(str(output.fasta), 'w') as f:
            for record in fasta_in:
                ID = record.description.split("|")[0]
                f.write(">" + ID + "\n")
                f.write(str(record.seq) + "\n")


rule uk_add_sample_date:
    input:
        metadata = config["latest_uk_metadata"],
    output:
        metadata = config["output_path"] + "/1/uk_latest.add_sample_date.csv",
    log:
        config["output_path"] + "/logs/1_uk_add_sample_date.log"
    resources: mem_per_cpu=20000
    run:
        df = pd.read_csv(input.metadata, sep = "\t")

        sample_date = []

        for i,row in df.iterrows():

            if not pd.isnull(row['collection_date']) and row['collection_date'] != "None":
                sample_date.append(row['collection_date'])
            elif not pd.isnull(row['received_date']) and row['received_date'] != "None":
                sample_date.append(row['received_date'])
            else:
                sample_date.append("")

        df['sample_date'] = sample_date
        df.to_csv(output.metadata, index=False, sep = ",")


rule uk_add_pillar_2:
    input:
        metadata = rules.uk_add_sample_date.output.metadata,
    output:
        metadata = config["output_path"] + "/1/uk_latest.add_pillar_2.csv",
    log:
        config["output_path"] + "/logs/1_uk_add_pillar_2.log"
    resources: mem_per_cpu=20000
    run:
        df = pd.read_csv(input.metadata, sep = ",")

        pillar_2 = []

        for i,row in df.iterrows():

            if row['collection_pillar'] == 2 or row['central_sample_id'][0:4] in ["ALDP", "CAMC", "MILK", "QEUH"]:
                pillar_2.append(True)
            else:
                pillar_2.append(False)

        df['pillar_2'] = pillar_2
        df.to_csv(output.metadata, index=False, sep = ",")


rule uk_make_sequence_name:
    input:
        metadata = rules.uk_add_pillar_2.output.metadata,
    output:
        metadata = config["output_path"] + "/1/uk_latest.add_sample_date.add_sequence_name.csv",
    log:
        config["output_path"] + "/logs/1_uk_make_sequence_name.log"
    resources: mem_per_cpu=20000
    run:
        df = pd.read_csv(input.metadata)

        adm1_to_country = {"UK-SCT": "Scotland",
                           "UK-WLS": "Wales",
                           "UK-ENG": "England",
                           "UK-NIR": "Northern_Ireland"}

        sequence_name = []

        for i,row in df.iterrows():
            country = adm1_to_country[row['adm1']]
            id = row['central_sample_id']
            year = str(row['sample_date']).split("-")[0]
            name = country + "/" + id + "/" + year

            sequence_name.append(name)

        df['sequence_name'] = sequence_name
        df.to_csv(output.metadata, index=False)


rule uk_add_gisaid_accession:
    input:
        metadata = rules.uk_make_sequence_name.output.metadata,
        accessions_table = config["latest_uk_accessions"],
    output:
        metadata = config["output_path"] + "/1/uk_latest.add_sample_date.add_sequence_name.accessions.csv"
    log:
        config["output_path"] + "/logs/1_uk_add_gisaid_accession.log"
    resources: mem_per_cpu=20000
    run:
        logfile = open(str(log), "w")
        df_acc = pd.read_csv(input.accessions_table, sep='\t')
        accessions_dict = {}
        for i,row in df_acc.iterrows():
            central_sample_id = row["central_sample_id"]
            run_name = row["run_name"]
            gisaid_accession = row["gisaid.accession"]

            if central_sample_id in accessions_dict:
                if run_name in accessions_dict[central_sample_id]:
                    logfile.write(f'duplicate central_sample_id * run_name in accessions list: {central_sample_id} {run_name}\n')
                    continue
                accessions_dict[central_sample_id][run_name] = gisaid_accession
            else:
                accessions_dict[central_sample_id] = {run_name: gisaid_accession}

        df = pd.read_csv(input.metadata)
        covv_accession_id = []
        for i,row in df.iterrows():
            acc = ""
            if row["central_sample_id"] in accessions_dict:
                if row["run_name"] in accessions_dict[row["central_sample_id"]]:
                    acc = accessions_dict[row["central_sample_id"]][row["run_name"]]

            covv_accession_id.append(acc)

        df['covv_accession_id'] = covv_accession_id
        df.to_csv(output.metadata, index=False)
        logfile.close()


rule uk_annotate_to_remove_duplicates:
    input:
        fasta = rules.uk_strip_header_digits.output.fasta,
        metadata = rules.uk_add_gisaid_accession.output.metadata,
    output:
        metadata = config["output_path"] + "/1/uk_latest.add_sample_date.add_sequence_name.accessions.annotated.csv"
    log:
        config["output_path"] + "/logs/1_uk_annotate_to_remove_duplicates.log"
    resources: mem_per_cpu=20000
    shell:
        """
        fastafunk annotate \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --out-metadata {output.metadata} \
          --log-file {log} \
          --index-column fasta_header &> {log}
        """
          # --add-cov-id \


rule uk_remove_duplicates_COGID_by_gaps:
    input:
        fasta = rules.uk_strip_header_digits.output.fasta,
        metadata = rules.uk_annotate_to_remove_duplicates.output.metadata
    output:
        fasta = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated_cov_id.fasta",
        metadata = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated_cov_id.csv"
    log:
        config["output_path"] + "/logs/1_uk_filter_duplicates_bycovid.log"
    resources: mem_per_cpu=20000
    shell:
        """
        fastafunk subsample \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --group-column central_sample_id \
          --index-column fasta_header \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --sample-size 1 \
          --select-by-min-column gaps &> {log}
        """


rule uk_update_sample_dates:
    input:
        metadata = rules.uk_remove_duplicates_COGID_by_gaps.output.metadata,
        updated_dates = config["uk_updated_dates"],
    output:
        metadata = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated_cov_id.sample_date.updated_sample_date.csv",
    log:
        config["output_path"] + "/logs/1_uk_update_sample_dates.log"
    resources: mem_per_cpu=20000
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.updated_dates} \
          --index-column central_sample_id \
          --join-on central_sample_id \
          --new-columns sample_date \
          --out-metadata {output.metadata} &>> {log}
        """


rule uk_add_epi_week:
    input:
        metadata = rules.uk_update_sample_dates.output.metadata
    output:
        metadata = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated_cov_id.sample_date.updated_sample_date.epi_week.csv",
    log:
        config["output_path"] + "/logs/1_uk_add_epi_week.log"
    resources: mem_per_cpu=20000
    shell:
        """
        datafunk add_epi_week \
        --input-metadata {input.metadata} \
        --output-metadata {output.metadata} \
        --date-column sample_date \
        --epi-week-column-name edin_epi_week \
        --epi-day-column-name edin_epi_day &> {log}
        """


rule uk_annotate_to_remove_duplicates_by_biosample:
    input:
        metadata = rules.uk_add_epi_week.output.metadata,
    output:
        metadata = config["output_path"] + "/1/uk_latest.epi_week.annotated2.csv",
    log:
        config["output_path"] + "/logs/1_uk_annotate_to_remove_duplicates_by_biosample.log",
    resources: mem_per_cpu=20000
    run:
        df = pd.read_csv(input.metadata)

        edin_biosample = []
        edin_biosample_root = []

        for i,row in df.iterrows():
            if pd.isnull(row['biosample_source_id']):
                edin_biosample.append(row['central_sample_id'])
            else:
                edin_biosample.append(row['biosample_source_id'])

            if pd.isnull(row['root_biosample_source_id']):
                edin_biosample_root.append(row['central_sample_id'])
            else:
                edin_biosample_root.append(row['root_biosample_source_id'])

        df['edin_biosample'] = edin_biosample
        df['edin_biosample_root'] = edin_biosample_root
        df.to_csv(output.metadata, index=False)


rule uk_remove_duplicates_biosamplesourceid_by_date:
    input:
        fasta = rules.uk_remove_duplicates_COGID_by_gaps.output.fasta,
        metadata = rules.uk_annotate_to_remove_duplicates_by_biosample.output.metadata
    output:
        fasta = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated_cov_id_biosample_source_id.fasta",
        metadata = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated_cov_id_biosample_source_id.csv"
    log:
        config["output_path"] + "/logs/1_uk_filter_duplicates_by_biosample.log"
    resources: mem_per_cpu=20000
    shell:
        """
        fastafunk subsample \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --group-column edin_biosample \
          --index-column fasta_header \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --sample-size 1 \
          --select-by-min-column edin_epi_day &> {log}
        """


rule uk_remove_duplicates_root_biosample_by_gaps:
    input:
        fasta = rules.uk_remove_duplicates_biosamplesourceid_by_date.output.fasta,
        metadata = rules.uk_remove_duplicates_biosamplesourceid_by_date.output.metadata
    output:
        fasta = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated_rootbiosample.fasta",
        metadata = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated_rootbiosample.csv"
    log:
        config["output_path"] + "/logs/1_uk_filter_duplicates_root_biosample_by_gaps.log"
    resources: mem_per_cpu=20000
    shell:
        """
        fastafunk subsample \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --group-column edin_biosample_root \
          --index-column fasta_header \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --sample-size 1 \
          --select-by-min-column gaps &> {log}
        """


rule uk_unify_headers:
    input:
        fasta = rules.uk_remove_duplicates_root_biosample_by_gaps.output.fasta,
        metadata = rules.uk_remove_duplicates_root_biosample_by_gaps.output.metadata
    output:
        fasta = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated.unify_headers.fasta",
        metadata = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated.unify_headers.csv"
    log:
        config["output_path"] + "/logs/1_uk_unify_headers.log"
    resources: mem_per_cpu=20000
    run:
        df = pd.read_csv(input.metadata)
        header_dict = {}

        for i,row in df.iterrows():
            header_dict[row['fasta_header']] = row['sequence_name']

        fasta_in = SeqIO.parse(str(input.fasta), "fasta")
        with open(str(output.fasta), 'w') as f:
            for record in fasta_in:
                new_ID = header_dict[record.description]
                f.write(">" + new_ID + "\n")
                f.write(str(record.seq) + "\n")

        df.to_csv(output.metadata, index=False)


rule uk_sed_United_Kingdom_to_UK:
    input:
        metadata = rules.uk_unify_headers.output.metadata
    output:
        metadata = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated.unify_headers.uk_sedded.csv"
    log:
        config["output_path"] + "/logs/1_uk_sed_United_Kingdom_to_UK.log"
    resources: mem_per_cpu=20000
    shell:
        """
        sed 's/United Kingdom/UK/g' {input.metadata} > {output.metadata}
        """


rule uk_minimap2_to_reference:
    input:
        fasta = rules.uk_unify_headers.output.fasta,
        reference = config["reference_fasta"]
    output:
        sam = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.mapped.sam"
    log:
        config["output_path"] + "/logs/1_uk_minimap2_to_reference.log"
    threads: 16
    resources: mem_per_cpu=2000
    shell:
        """
        minimap2 -t {threads} -a -x asm5 {input.reference} {input.fasta} > {output.sam} 2> {log}
        """


rule uk_get_variants:
    input:
        sam = rules.uk_minimap2_to_reference.output.sam,
        reference = config["reference_fasta"],
        genbank_anno = config["reference_genbank_annotation"],
    output:
        variants = config["output_path"] + "/1/uk.variants.csv",
    log:
        config["output_path"] + "/logs/1_uk_get_variants.log"
    threads: 12
    shell:
        """
        /cephfs/covid/bham/raccoon-dog/programs/gofasta/gofasta sam variants -t {threads} \
            --samfile {input.sam} \
            --reference {input.reference} \
            --genbank {input.genbank_anno} \
            --outfile {output.variants} &>> {log}
        """


rule uk_remove_insertions_and_trim_and_pad:
    input:
        sam = rules.uk_minimap2_to_reference.output.sam,
        reference = config["reference_fasta"]
    params:
        trim_start = config["trim_start"],
        trim_end = config["trim_end"],
        insertions = config["output_path"] + "/1/uk_insertions.txt",
        deletions = config["output_path"] + "/1/uk_deletions.txt"
    output:
        fasta = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.alignment.trimmed.fasta",
        insertions = config["export_path"] + "/metadata/uk_insertions.txt",
        deletions = config["export_path"] + "/metadata/uk_deletions.txt",
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
          --log-inserts \
          --log-deletions &> {log}
        mv insertions.txt {params.insertions}
        mv deletions.txt {params.deletions}
        cp {params.insertions} {output.insertions}
        cp {params.deletions} {output.deletions}
        """


rule uk_mask_1:
    input:
        fasta = rules.uk_remove_insertions_and_trim_and_pad.output.fasta,
        mask = config["uk_mask_file"]
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
    resources: mem_per_cpu=20000
    shell:
        """
        datafunk filter_fasta_by_covg_and_length \
          -i {input.fasta} \
          -o {output.fasta} \
          --min-covg {params.min_covg} &> {log}
        """


rule uk_filter_omitted_sequences:
    input:
        fasta = rules.uk_filter_low_coverage_sequences.output.fasta,
        omissions = config["uk_omissions"]
    output:
        fasta = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.trimmed.low_covg_filtered.omissions_filtered.fasta"
    log:
        config["output_path"] + "/logs/1_uk_filter_omitted_sequences.log"
    resources: mem_per_cpu=20000
    shell:
        """
        datafunk remove_fasta \
          -i {input.fasta} \
          -f {input.omissions} \
          -o {output.fasta}  &> {log}
        """


rule uk_full_untrimmed_alignment:
    input:
        sam = rules.uk_minimap2_to_reference.output.sam,
        reference = config["reference_fasta"],
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
        mask = config["uk_mask_file"]
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


# get the same alignment as we use for tree building but with no mask
rule uk_get_unmasked_alignment:
    input:
        fasta_template = rules.uk_filter_omitted_sequences.output.fasta,
        fasta = rules.uk_full_untrimmed_alignment.output.fasta,
    output:
        fasta = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.alignment.full.unmasked.fasta",
    log:
        config["output_path"] + "/logs/1_uk_get_unmasked_alignment.log"
    run:
        fasta_template = SeqIO.index(str(input.fasta_template), "fasta")
        fasta_in = SeqIO.index(str(input.fasta), "fasta")

        with open(str(output.fasta), 'w') as fasta_out:
            for record in fasta_in:
                if record in fasta_template:
                    fasta_out.write('>' + fasta_in[record].id + '\n')
                    fasta_out.write(str(fasta_in[record].seq) + '\n')


rule UK_AA_finder:
    input:
        fasta = rules.uk_full_untrimmed_alignment.output.fasta,
        AAs = config["AAs"]
    output:
        found = config["output_path"] + "/1/cog.AA_finder.csv",
    log:
        config["output_path"] + "/logs/1_UK_AA_finder.log"
    shell:
        """
        datafunk AA_finder -i {input.fasta} --codons-file {input.AAs} --genotypes-table {output.found} &> {log}
        """


rule add_AA_finder_result_to_metadata:
    input:
        AAs = config["AAs"],
        metadata = rules.uk_sed_United_Kingdom_to_UK.output.metadata,
        new_data = rules.UK_AA_finder.output.found
    output:
        metadata = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.with_AA_finder.csv"
    log:
        config["output_path"] + "/logs/1_add_AA_finder_result_to_metadata.log"
    resources: mem_per_cpu=20000
    shell:
        """
        columns=$(head -n1 {input.new_data} | cut -d',' -f2- | tr ',' ' ')
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.new_data} \
          --index-column sequence_name \
          --join-on sequence_name \
          --new-columns $columns \
          --out-metadata {output.metadata} &>> {log}
        """


rule uk_del_finder:
    input:
        fasta = rules.uk_full_untrimmed_alignment.output.fasta,
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
        metadata = rules.add_AA_finder_result_to_metadata.output.metadata,
        new_data = rules.uk_del_finder.output.metadata
    output:
        metadata = config["output_path"] + "/1/uk_latest.unify_headers.epi_week.deduplicated.with_AA_finder.with_del_finder.csv"
    log:
        config["output_path"] + "/logs/1_add_del_finder_result_to_metadata.log"
    resources: mem_per_cpu=20000
    shell:
        """
        columns=$(head -n1 {input.new_data} | cut -d',' -f2- | tr ',' ' ')
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.new_data} \
          --index-column sequence_name \
          --join-on sequence_name \
          --new-columns $columns \
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
    resources: mem_per_cpu=20000
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.previous_metadata} \
          --index-column sequence_name \
          --join-on sequence_name \
          --new-columns uk_lineage edin_date_stamp \
          --out-metadata {output.metadata_temp} &>> {log}

        fastafunk add_columns \
          --in-metadata {output.metadata_temp} \
          --in-data {input.global_lineages} \
          --index-column sequence_name \
          --join-on taxon \
          --new-columns lineage lineage_support lineages_version \
          --where-column lineage_support=probability lineages_version=pangoLEARN_version \
          --out-metadata {output.metadata} &>> {log}
        """


rule uk_extract_lineageless:
    input:
        fasta = rules.uk_filter_omitted_sequences.output,
        metadata = rules.uk_add_previous_lineages_to_metadata.output.metadata,
    output:
        fasta = config["output_path"] + "/1/uk.new.pangolin_lineages.fasta",
    log:
        config["output_path"] + "/logs/1_extract_lineageless.log"
    resources: mem_per_cpu=20000
    run:
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


rule uk_add_dups_to_lineageless:
    input:
        master_fasta = rules.uk_filter_omitted_sequences.output.fasta,
        lineageless_fasta = rules.uk_extract_lineageless.output.fasta,
        metadata = rules.uk_add_previous_lineages_to_metadata.output.metadata,
        dup_cogid_log = rules.uk_remove_duplicates_COGID_by_gaps.log,
        # dup_biosample_log = rules.uk_remove_duplicates_biosamplesourceid_by_date.log,
        # dup_rootbio_log = rules.uk_remove_duplicates_root_biosample_by_gaps.log,
    output:
        fasta = config["output_path"] + "/1/uk.new.dedupes.pangolin_lineages.fasta",
    log:
        config["output_path"] + "/logs/1_uk_add_dups_to_lineageless.log"
    resources: mem_per_cpu=20000
    run:
        master_fasta_in = SeqIO.index(str(input.master_fasta), "fasta")

        df = pd.read_csv(input.metadata)

        lineageless_sequence_record = set()
        dup_record = set()

        lineageless_fasta_in = open(str(input.lineageless_fasta), "r")
        pangolin_fasta_out = open(str(output.fasta), "w")

        for record in SeqIO.parse(lineageless_fasta_in, "fasta"):
            pangolin_fasta_out.write(">" + record.id + "\n")
            pangolin_fasta_out.write(str(record.seq) + "\n")
            lineageless_sequence_record.add(record.id)

        with open(str(input.dup_cogid_log), "r") as cogid_log:
            for line in cogid_log:
                l = line.rstrip()
                if l.startswith("COGUK/"):
                    dup_record.add(l)

        # # can uncomment these if necessary, but it shouldn't be
        # with open(str(input.dup_biosample_log), "r") as biosample_log:
        #     for line in dup_biosample_log:
        #         l = line.rstrip()
        #         if l.startswith("COGUK/"):
        #             dup_record.add(l)
        #
        # with open(str(input.dup_rootbio_log), "r") as rootbio_log:
        #     for line in rootbio_log:
        #         l = line.rstrip()
        #         if l.startswith("COGUK/"):
        #             dup_record.add(l)

        for i,row in df.iterrows():
            if row["fasta_header"] in dup_record:
                if row["sequence_name"] in lineageless_sequence_record:
                    continue
                if row["sequence_name"] in master_fasta_in:
                    record = master_fasta_in[row["sequence_name"]]
                    pangolin_fasta_out.write(">" + record.id + "\n")
                    pangolin_fasta_out.write(str(record.seq) + "\n")

        lineageless_fasta_in.close()
        pangolin_fasta_out.close()


rule summarize_preprocess_uk:
    input:
        raw_fasta = config["latest_uk_fasta"],
        deduplicated_fasta_by_covid = rules.uk_remove_duplicates_COGID_by_gaps.output.fasta,
        deduplicated_fasta_by_biosampleid = rules.uk_remove_duplicates_biosamplesourceid_by_date.output.fasta,
        deduplicated_fasta_by_rootbiosample = rules.uk_remove_duplicates_root_biosample_by_gaps.output.fasta,
        unify_headers_fasta = rules.uk_unify_headers.output.fasta,
        removed_low_covg_fasta = rules.uk_filter_low_coverage_sequences.output.fasta,
        removed_omitted_fasta = rules.uk_filter_omitted_sequences.output.fasta,
        full_unmasked_alignment = rules.uk_full_untrimmed_alignment.output.fasta,
        full_metadata = rules.uk_add_previous_lineages_to_metadata.output.metadata,
        pangolin_fasta = rules.uk_add_dups_to_lineageless.output.fasta,
        variants = rules.uk_get_variants.output.variants,
    params:
        grapevine_webhook = config["grapevine_webhook"],
        json_path = config["json_path"],
        date=config["date"]
    log:
        config["output_path"] + "/logs/1_summarize_preprocess_uk.log"
    shell:
        """
        echo "> Number of sequences in raw UK fasta: $(cat {input.raw_fasta} | grep ">" | wc -l)\\n" &> {log}
        echo "> Number of sequences after deduplication by central_sample_id: $(cat {input.deduplicated_fasta_by_covid} | grep ">" | wc -l)\\n" &>> {log}
        echo "> Number of sequences after deduplication by bio_sample_id: $(cat {input.deduplicated_fasta_by_biosampleid} | grep ">" | wc -l)\\n" &>> {log}
        echo "> Number of sequences after deduplication by root_source_biosample_id: $(cat {input.deduplicated_fasta_by_rootbiosample} | grep ">" | wc -l)\\n" &>> {log}
        echo "> Number of sequences after unifying headers: $(cat {input.unify_headers_fasta} | grep ">" | wc -l)\\n" &>> {log}
        echo "> Number of sequences after mapping: $(cat {input.removed_low_covg_fasta} | grep ">" | wc -l)\\n" &>> {log}
        echo "> Number of sequences after removing those in omissions file: $(cat {input.removed_omitted_fasta} | grep ">" | wc -l)\\n" &>> {log}
        echo "> Number of new/deduped sequences passed to Pangolin for typing: $(cat {input.pangolin_fasta} | grep ">" | wc -l)\\n" &>> {log}
        echo ">\\n" >> {log}

        echo '{{"text":"' > {params.json_path}/1_data.json
        echo "*Step 1: {params.date} COG-UK preprocessing complete*\\n" >> {params.json_path}/1_data.json
        cat {log} >> {params.json_path}/1_data.json
        echo '"}}' >> {params.json_path}/1_data.json
        echo "webhook {params.grapevine_webhook}"
        curl -X POST -H "Content-type: application/json" -d @{params.json_path}/1_data.json {params.grapevine_webhook}
        """
