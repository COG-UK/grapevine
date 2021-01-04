rule gisaid_process_json:
    input:
        json = config["latest_gisaid_json"],
        omitted = config["previous_omitted_file"],
    output:
        fasta = config["output_path"] + "/0/gisaid.fasta",
        metadata = config["output_path"] + "/0/gisaid.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_process_json.log"
    resources: mem_per_cpu=20000
    shell:
        """
        datafunk process_gisaid_data \
          --input-json {input.json} \
          --input-metadata False \
          --exclude-file {input.omitted} \
          --output-fasta {output.fasta} \
          --output-metadata {output.metadata} \
          --exclude-undated &> {log}
        """


rule gisaid_unify_headers:
    input:
        fasta = rules.gisaid_process_json.output.fasta,
        metadata = rules.gisaid_process_json.output.metadata,
    output:
        fasta = config["output_path"] + "/0/gisaid.UH.fasta",
        metadata = config["output_path"] + "/0/gisaid.UH.csv",
    log:
        config["output_path"] + "/logs/0_gisaid_unify_headers.log"
    run:
        import pandas as pd
        from Bio import SeqIO

        fasta_in = SeqIO.index(str(input.fasta), "fasta")
        df = pd.read_csv(input.metadata, sep=',')

        sequence_name = []

        with open(str(output.fasta), 'w') as fasta_out:
            for i,row in df.iterrows():
                edin_header = row["edin_header"]
                new_header = edin_header.split("|")[0]
                sequence_name.append(new_header)

                try:
                    record = fasta_in[edin_header]
                    fasta_out.write(">" + new_header + "\n")
                    fasta_out.write(str(record.seq) + "\n")
                except:
                    continue

        df['sequence_name'] = sequence_name
        df.to_csv(output.metadata, index=False, sep = ",")


rule gisaid_add_previous_lineages:
    input:
        metadata = rules.gisaid_unify_headers.output.metadata,
        previous_lineages = config["previous_gisaid_lineages"],
    output:
        metadata = config["output_path"] + "/0/gisaid.lineages.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_add_previous_lineages.log"
    resources: mem_per_cpu=20000
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.previous_lineages} \
          --index-column sequence_name \
          --join-on sequence_name \
          --new-columns lineage lineage_support lineages_version \
          --out-metadata {output.metadata} &> {log}
        """


rule gisaid_remove_duplicates:
    input:
        fasta = rules.gisaid_unify_headers.output.fasta,
        metadata = rules.gisaid_add_previous_lineages.output.metadata
    output:
        fasta = config["output_path"] + "/0/gisaid.UH.RD.fasta",
        metadata = config["output_path"] + "/0/gisaid.UH.RD.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_filter_duplicates.log"
    resources: mem_per_cpu=20000
    shell:
        """
        fastafunk subsample \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --group-column sequence_name \
          --index-column sequence_name \
          --out-fasta {output.fasta} \
          --sample-size 1 \
          --out-metadata {output.metadata} \
          --select-by-min-column edin_epi_day &> {log}
        """


rule gisaid_counts_by_country:
    input:
        metadata = rules.gisaid_remove_duplicates.output.metadata
    output:
        counts = config["output_path"] + "/0/gisaid_counts_by_country.csv",
    log:
        config["output_path"] + "/logs/0_gisaid_counts_by_country.log"
    resources: mem_per_cpu=20000
    shell:
        """
        fastafunk count \
          --in-metadata {input.metadata} \
          --group-column edin_admin_0 \
          --log-file {output.counts} &> {log}
        """


rule gisaid_filter_1:
    input:
        fasta = rules.gisaid_remove_duplicates.output.fasta
    params:
        min_covg = config["min_covg"],
        min_length = config["min_length"]
    output:
        fasta = config["output_path"] + "/0/gisaid.RD.UH.filt1.fasta"
    log:
        config["output_path"] + "/logs/0_gisaid_filter_1.log"
    shell:
        """
        datafunk filter_fasta_by_covg_and_length \
          -i {input.fasta} \
          -o {output} \
          --min-length {params.min_length} &> {log}
        """


rule gisaid_minimap2_to_reference:
    input:
        fasta = rules.gisaid_filter_1.output.fasta,
        reference = config["reference_fasta"]
    output:
        sam = config["output_path"] + "/0/gisaid.mapped.sam"
    log:
        config["output_path"] + "/logs/0_gisaid_minimap2_to_reference.log"
    threads: 8
    shell:
        """
        minimap2 -t8 -a -x asm5 {input.reference} {input.fasta} -o {output} &> {log}
        """


rule gisaid_get_variants:
    input:
        sam = rules.gisaid_minimap2_to_reference.output.sam,
        reference = config["reference_fasta"],
        genbank_anno = config["reference_genbank_annotation"],
    output:
        variants = config["output_path"] + "/0/gisaid.variants.csv",
        global_variants = config["output_path"] + "/0/gisaid.global.variants.csv",
    log:
        config["output_path"] + "/logs/0_gisaid_get_variants.log"
    threads: 8
    shell:
        """
        /cephfs/covid/bham/climb-covid19-jacksonb/programs/gofasta/gofasta sam variants -t {threads} \
            --samfile {input.sam} \
            --reference {input.reference} \
            --genbank {input.genbank_anno} \
            --outfile {output.variants} &>> {log}

        head -n1 {output.variants} > {output.global_variants}
        tail -n+2 {output.variants} | grep -v -E "^England|^Northern_Ireland|^Wales|^Scotland" >> {output.global_variants}
        """


rule gisaid_remove_insertions_and_pad:
    input:
        sam = rules.gisaid_minimap2_to_reference.output.sam,
        reference = config["reference_fasta"]
    params:
        trim_start = config["trim_start"],
        trim_end = config["trim_end"],
        insertions = config["output_path"] + "/0/gisaid_insertions.txt",
        deletions = config["output_path"] + "/0/gisaid_deletions.txt",
    output:
        fasta = config["output_path"] + "/0/gisaid.RD.UH.filt1.mapped.fasta"
    threads: 8
    log:
        config["output_path"] + "/logs/0_gisaid_remove_insertions_and_pad.log"
    shell:
        """
        /cephfs/covid/bham/climb-covid19-jacksonb/programs/gofasta/gofasta sam toMultiAlign -t {threads} -s {input.sam} -o {output.fasta} --trim --trimstart 265 --trimend 29674 --pad &> {log}
        """


rule gisaid_filter_2:
    input:
        fasta = rules.gisaid_remove_insertions_and_pad.output.fasta
    output:
        fasta = config["output_path"] + "/0/gisaid.RD.UH.filt.mapped.filt2.fasta"
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


rule gisaid_mask:
    input:
        fasta = rules.gisaid_filter_2.output.fasta,
        mask = config["gisaid_mask_file"]
    output:
        fasta = config["output_path"] + "/0/gisaid.RD.UH.filt.mapped.filt2.masked.fasta",
    shell:
        """
        datafunk mask \
          --input-fasta {input.fasta} \
          --output-fasta {output.fasta} \
          --mask-file \"{input.mask}\"
        """


rule gisaid_distance_QC:
    input:
        fasta = rules.gisaid_mask.output.fasta,
        metadata = rules.gisaid_remove_duplicates.output.metadata
    log:
        config["output_path"] + "/logs/0_gisaid_distance_QC.log"
    output:
        table = config["output_path"] + "/0/QC_distances.tsv",
    resources: mem_per_cpu=20000
    shell:
        """
        datafunk distance_to_root \
          --input-fasta {input.fasta} \
          --input-metadata {input.metadata} &> {log}

        mv distances.tsv {output.table}
        """


rule gisaid_filter_on_distance_to_WH04:
    input:
        fasta = rules.gisaid_mask.output.fasta,
        table = rules.gisaid_distance_QC.output.table,
    output:
        fasta = config["output_path"] + "/0/gisaid.RD.UH.filt.mapped.filt2.masked.filt3.fasta",
    log:
        config["output_path"] + "/logs/0_gisaid_filter_on_distance_to_WH04.log"
    resources: mem_per_cpu=20000
    run:
        from Bio import SeqIO
        import pandas as pd

        fasta_in = SeqIO.index(str(input.fasta), "fasta")
        df = pd.read_csv(input.table, sep='\t')

        with open(str(output.fasta), 'w') as fasta_out, open(str(log), 'w') as log_out:
            for i,row in df.iterrows():
                sequence_name = row['sequence_name']
                distance = row['distance_stdevs']
                if distance < 4.0:
                    if sequence_name in fasta_in:
                        record = fasta_in[sequence_name]
                        fasta_out.write('>' + record.id + '\n')
                        fasta_out.write(str(record.seq) + '\n')
                else:
                    log_out.write(sequence_name + ' was filtered for having too high a distance to WH04 (' + str(distance) + ' epi-week std devs)\n')



rule gisaid_AA_finder:
    input:
        fasta = rules.gisaid_filter_on_distance_to_WH04.output.fasta,
        AAs = config["AAs"]
    output:
        found = config["output_path"] + "/0/gisaid.AA_finder.csv",
    log:
        config["output_path"] + "/logs/0_gisaid_AA_finder.log"
    shell:
        """
        datafunk AA_finder -i {input.fasta} --codons-file {input.AAs} --genotypes-table {output.found} &> {log}
        """


rule gisaid_add_AA_finder_result_to_metadata:
    input:
        snps = config["snps"],
        metadata = rules.gisaid_remove_duplicates.output.metadata,
        new_data = rules.gisaid_AA_finder.output.found
    output:
        metadata = config["output_path"] + "/0/gisaid.RD.UH.AAfinder.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_add_AA_finder_result_to_metadata.log"
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


rule gisaid_extract_lineageless:
    input:
        fasta = rules.gisaid_filter_on_distance_to_WH04.output,
        metadata = rules.gisaid_add_AA_finder_result_to_metadata.output.metadata,
    output:
        fasta = config["output_path"] + "/0/gisaid.new.pangolin_lineages.fasta",
    log:
        config["output_path"] + "/logs/0_extract_lineageless.log"
    resources: mem_per_cpu=20000
    run:
        from Bio import SeqIO
        import pandas as pd

        fasta_in = SeqIO.index(str(input.fasta), "fasta")
        df = pd.read_csv(input.metadata)

        sequence_record = []

        with open(str(output.fasta), 'w') as fasta_out:
            if 'lineage' in df.columns:
                for i,row in df.iterrows():
                    if pd.isnull(row['lineage']):
                        sequence_name = row['sequence_name']
                        if sequence_name in fasta_in:
                            if sequence_name not in sequence_record:
                                record = fasta_in[sequence_name]
                                fasta_out.write('>' + record.id + '\n')
                                fasta_out.write(str(record.seq) + '\n')
                                sequence_record.append(sequence_name)
            else:
                for i,row in df.iterrows():
                    sequence_name = row['sequence_name']
                    if sequence_name in fasta_in:
                        if sequence_name not in sequence_record:
                            record = fasta_in[sequence_name]
                            fasta_out.write('>' + record.id + '\n')
                            fasta_out.write(str(record.seq) + '\n')
                            sequence_record.append(sequence_name)


rule gisaid_normal_pangolin:
    input:
        fasta = rules.gisaid_extract_lineageless.output.fasta
    params:
        outdir = config["output_path"] + "/0/normal_pangolin",
        tmpdir = config["output_path"] + "/0/normal_pangolin/tmp"
    output:
        lineages = config["output_path"] + "/0/normal_pangolin/lineage_report.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_normal_pangolin.log"
    resources: mem_per_cpu=20000
    shell:
        """
        pangolin {input.fasta} --tempdir {params.tmpdir} --outdir {params.outdir} --verbose > {log} 2>&1
        """


rule gisaid_add_pangolin_lineages_to_metadata:
    input:
        metadata = rules.gisaid_add_AA_finder_result_to_metadata.output.metadata,
        normal_lineages = rules.gisaid_normal_pangolin.output.lineages
    output:
        metadata = config["output_path"] + "/0/gisaid.RD.UH.SNPfinder.lineages.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_add_pangolin_lineages_to_metadata.log"
    resources: mem_per_cpu=20000
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.normal_lineages} \
          --index-column sequence_name \
          --join-on taxon \
          --new-columns lineage lineage_support lineages_version \
          --where-column lineage_support=probability lineages_version=pangoLEARN_version \
          --out-metadata {output.metadata} &> {log}
        """


rule gisaid_del_finder:
    input:
        fasta = rules.gisaid_filter_on_distance_to_WH04.output.fasta,
        dels = config["dels"]
    output:
        metadata = config["output_path"] + "/0/gisaid.del_finder.csv",
    log:
        config["output_path"] + "/logs/0_gisaid_del_finder.log"
    shell:
        """
        datafunk del_finder \
            -i {input.fasta} \
            --deletions-file {input.dels} \
            --genotypes-table {output.metadata} &> {log}
        """


rule gisaid_add_del_finder_result_to_metadata:
    input:
        dels = config["dels"],
        metadata = rules.gisaid_add_pangolin_lineages_to_metadata.output.metadata,
        new_data = rules.gisaid_del_finder.output.metadata
    output:
        metadata = config["output_path"] + "/0/gisaid.RD.UH.SNPfinder.lineages.del_finder.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_add_del_finder_result_to_metadata.log"
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


rule gisaid_output_all_matched_metadata:
    input:
        fasta = rules.gisaid_filter_on_distance_to_WH04.output.fasta,
        metadata = rules.gisaid_add_del_finder_result_to_metadata.output.metadata
    output:
        fasta = config["output_path"] + "/0/gisaid.all.fasta",
        metadata = config["output_path"] + "/0/gisaid.all.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_output_all_matched_metadata.log"
    resources: mem_per_cpu=20000
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column \
                         sequence_name \
                         country \
                         edin_admin_1 \
                         edin_admin_2 \
                         edin_travel \
                         edin_date_stamp \
                         sample_date \
                         epi_week \
                         lineage \
                         lineages_version \
                         lineage_support \
                         d614g \
                         n439k \
                         p323l \
                         a222v \
                         y453f \
                         n501y \
                         t1001i \
                         p681h \
                         q27stop \
                         del_21765_6 \
                         del_1605_3 \
                         covv_accession_id \
                         covv_virus_name \
                         covv_location \
                         covv_add_host_info \
                         covv_assembly_method \
                         covv_gender \
                         covv_host \
                         covv_passage \
                         covv_patient_age \
                         covv_seq_technology \
                         covv_specimen \
                         covv_subm_date \
          --where-column sample_date=covv_collection_date \
                                 epi_week=edin_epi_week country=edin_admin_0 \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --log-file {log} \
          --restrict
        """


rule gisaid_exclude_uk_seqs:
    input:
        fasta = rules.gisaid_filter_on_distance_to_WH04.output.fasta,
    output:
        fasta = config["output_path"] + "/0/gisaid.global.fasta"
    run:
        from Bio import SeqIO

        out_handle = open(str(output.fasta), 'w')
        with open(str(input.fasta), 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                id = record.id
                seq = str(record.seq)
                if id.split('/')[0] in ['England', 'Wales', 'Scotland', 'Northern_Ireland', 'NorthernIreland']:
                    continue
                else:
                    out_handle.write('>' + id + '\n')
                    out_handle.write(seq + '\n')

        out_handle.close()


rule gisaid_output_global_matched_metadata:
    input:
        fasta = rules.gisaid_exclude_uk_seqs.output.fasta,
        metadata = rules.gisaid_add_del_finder_result_to_metadata.output.metadata
    output:
        fasta = temp(config["output_path"] + "/0/gisaid.global.temp.fasta"),
        metadata = config["output_path"] + "/0/gisaid.global.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_output_global_matched_metadata.log"
    resources: mem_per_cpu=20000
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column \
                         sequence_name \
                         country \
                         edin_admin_1 \
                         edin_admin_2 \
                         edin_travel \
                         edin_date_stamp \
                         sample_date \
                         epi_week \
                         lineage \
                         lineages_version \
                         lineage_support \
                         d614g \
                         n439k \
                         p323l \
                         a222v \
                         y453f \
                         n501y \
                         t1001i \
                         p681h \
                         q27stop \
                         del_21765_6 \
                         del_1605_3 \
                         covv_accession_id \
                         covv_virus_name \
                         covv_location \
                         covv_add_host_info \
                         covv_assembly_method \
                         covv_gender \
                         covv_host \
                         covv_passage \
                         covv_patient_age \
                         covv_seq_technology \
                         covv_specimen \
                         covv_subm_date \
          --where-column sample_date=covv_collection_date \
                                 epi_week=edin_epi_week country=edin_admin_0 \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --log-file {log} \
          --restrict
        """


rule gisaid_collapse:
    input:
        fasta = rules.gisaid_exclude_uk_seqs.output.fasta,
    output:
        fasta = config["output_path"] + "/0/gisaid.global.collapsed.fasta",
        tip_to_redudants = config["output_path"] + "/0/tip_to_redundants.csv",
        redundant_to_tips = config["output_path"] + "/0/redundant_to_tips.csv",
    log:
        config["output_path"] + "/logs/0_gisaid_collapse.log"
    resources: mem_per_cpu=8000
    threads: 8
    shell:
        """
        /cephfs/covid/bham/climb-covid19-jacksonb/programs/julia-1.5.0/bin/julia --threads {threads} /cephfs/covid/bham/climb-covid19-jacksonb/programs/julialign/src/collapse.jl \
        	-i {input.fasta} \
        	--retain Wuhan/WH04/2020 \
        	--retain Wuhan/WHU01/2020 \
        	--retain Italy/ABR-IZSGC-TE5166/2020 \
            --retain Germany/BY-MVP-V2010837/2020 \
            --retain Spain/VC-IBV-98006461/2020 \
	        -o {output.fasta} &> {log}

            mv tip_to_redundants.csv {output.tip_to_redudants} &>> {log}
            mv redundant_to_tips.csv {output.redundant_to_tips} &>> {log}
        """


rule gisaid_get_unique_redundants:
    input:
        fasta = rules.gisaid_exclude_uk_seqs.output.fasta,
        redundant_to_tips = rules.gisaid_collapse.output.redundant_to_tips,
    output:
        fasta = config["output_path"] + "/0/unique_redundants.fasta",
    log:
        config["output_path"] + "/logs/0_gisaid_get_unique_redundants.log"
    run:
        from Bio import SeqIO

        fasta_in = SeqIO.index(str(input.fasta), "fasta")

        firstline = True
        with open(str(input.redundant_to_tips), "r") as f_in, open(str(output.fasta), "w") as f_out:
            for line in f_in:
                if firstline:
                    firstline = False
                    continue

                l = line.rstrip().split(",")
                redundant = l[0]
                possible_tips = l[1].split("|")
                if len(possible_tips) > 1:
                    continue

                record = fasta_in[redundant]

                f_out.write(">" + record.id + "\n")
                f_out.write(str(record.seq) + "\n")


rule gisaid_get_collapsed_metadata:
    input:
        fasta = rules.gisaid_collapse.output.fasta,
        metadata = rules.gisaid_add_del_finder_result_to_metadata.output.metadata,
    output:
        fasta = temp(config["output_path"] + "/0/gisaid.global.collapsed.temp.fasta"),
        metadata = config["output_path"] + "/0/gisaid.global.collapsed.csv",
    log:
        config["output_path"] + "/logs/0_gisaid_get_collapsed_metadata.log"
    resources: mem_per_cpu=20000
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column \
                         sequence_name \
                         country \
                         edin_admin_1 \
                         edin_admin_2 \
                         edin_travel \
                         edin_date_stamp \
                         sample_date \
                         epi_week \
                         lineage \
                         lineages_version \
                         lineage_support \
                         d614g \
                         n439k \
                         p323l \
                         a222v \
                         y453f \
                         n501y \
                         t1001i \
                         p681h \
                         q27stop \
                         del_21765_6 \
                         del_1605_3 \
                         covv_accession_id \
                         covv_virus_name \
                         covv_location \
                         covv_add_host_info \
                         covv_assembly_method \
                         covv_gender \
                         covv_host \
                         covv_passage \
                         covv_patient_age \
                         covv_seq_technology \
                         covv_specimen \
                         covv_subm_date \
          --where-column sample_date=covv_collection_date \
                         epi_week=edin_epi_week country=edin_admin_0 \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --log-file {log} \
          --restrict
        """


rule gisaid_get_collapsed_expanded_metadata:
    input:
        collapsed_fasta = rules.gisaid_collapse.output.fasta,
        unique_fasta = rules.gisaid_get_unique_redundants.output.fasta,
        metadata = rules.gisaid_add_del_finder_result_to_metadata.output.metadata,
    output:
        fasta = config["output_path"] + "/0/gisaid.global.collapsed.unique_expanded.fasta",
        metadata = config["output_path"] + "/0/gisaid.global.collapsed.unique_expanded.csv",
    log:
        config["output_path"] + "/logs/0_gisaid_get_collapsed_expanded_metadata.log"
    resources: mem_per_cpu=20000
    shell:
        """
        cat {input.collapsed_fasta} {input.unique_fasta} > {output.fasta} 2> {log}

        fastafunk fetch \
          --in-fasta {output.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column \
                         sequence_name \
                         country \
                         edin_admin_1 \
                         edin_admin_2 \
                         edin_travel \
                         edin_date_stamp \
                         sample_date \
                         epi_week \
                         lineage \
                         lineages_version \
                         lineage_support \
                         d614g \
                         n439k \
                         p323l \
                         a222v \
                         y453f \
                         n501y \
                         t1001i \
                         p681h \
                         q27stop \
                         del_21765_6 \
                         del_1605_3 \
                         covv_accession_id \
                         covv_virus_name \
                         covv_location \
                         covv_add_host_info \
                         covv_assembly_method \
                         covv_gender \
                         covv_host \
                         covv_passage \
                         covv_patient_age \
                         covv_seq_technology \
                         covv_specimen \
                         covv_subm_date \
          --where-column sample_date=covv_collection_date \
                         epi_week=edin_epi_week country=edin_admin_0 \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --restrict 2>> {log}
        """


rule summarize_preprocess_gisaid:
    input:
        latest_fasta = rules.gisaid_process_json.output.fasta,
        deduplicated_fasta = rules.gisaid_remove_duplicates.output.fasta,
        unify_headers_fasta = rules.gisaid_unify_headers.output.fasta,
        removed_short_fasta = rules.gisaid_filter_1.output,
        removed_low_covg_fasta = rules.gisaid_filter_2.output.fasta,
        removed_distance_to_root_fasta = rules.gisaid_filter_on_distance_to_WH04.output.fasta,
        all_fasta = rules.gisaid_output_all_matched_metadata.output.fasta,
        all_metadata = rules.gisaid_output_all_matched_metadata.output.metadata,
        global_fasta = rules.gisaid_exclude_uk_seqs.output.fasta,
        global_metadata = rules.gisaid_output_global_matched_metadata.output.metadata,
        collapsed_fasta = rules.gisaid_collapse.output.fasta,
        collapsed_metadata = rules.gisaid_get_collapsed_metadata.output.metadata,
        collapsed_expanded_metadata = rules.gisaid_get_collapsed_expanded_metadata.output.metadata,
        counts = rules.gisaid_counts_by_country.output.counts,
        variants = rules.gisaid_get_variants.output.variants,
    params:
        publish_path = config["publish_path"] + "/GISAID/",
        published_counts = config["publish_path"] + "/GISAID/gisaid_counts_by_country.csv",
        published_all_fasta = config["publish_path"] + "/GISAID/gisaid.all.fasta",
        published_all_metadata = config["publish_path"] + "/GISAID/gisaid.all.csv",
        grapevine_webhook = config["grapevine_webhook"],
        date = config["date"],
    log:
        config["output_path"] + "/logs/0_summarize_preprocess_gisaid.log"
    shell:
        """
        mkdir -p {params.publish_path}

        cp {input.counts} {params.published_counts}
        cp {input.all_fasta} {params.published_all_fasta}
        cp {input.all_metadata} {params.published_all_metadata}

        echo "Number of sequences in latest GISAID download: $(cat {input.latest_fasta} | grep '>' | wc -l)\\n" >> {log}
        echo "Number of sequences after unifying headers: $(cat {input.unify_headers_fasta} | grep '>' | wc -l)\\n" >> {log}
        echo "Number of sequence after deduplicating: $(cat {input.deduplicated_fasta} | grep '>' | wc -l)\\n" >> {log}
        echo "Number of sequences after removing sequences <29000bps: $(cat {input.removed_short_fasta} | grep '>' | wc -l)\\n" >> {log}
        echo "Number of sequences after mapping and removing those with <95% coverage: $(cat {input.removed_low_covg_fasta} | grep '>' | wc -l)\\n" >> {log}
        echo "Number of sequences after removing those >4 epi-week stddevs to WH04: $(cat {input.removed_distance_to_root_fasta} | grep '>' | wc -l)\\n" >> {log}
        echo "Number of non-UK sequences: $(cat {input.global_fasta} | grep '>' | wc -l)\\n" >> {log}
        echo "Number of non-UK sequences after collapsing: $(cat {input.collapsed_fasta} | grep '>' | wc -l)\\n" >> {log}
        echo "> \\n" >> {log}
        echo "> Counts by country published to {params.published_counts}\\n" >> {log}
        echo "> Full alignment published to {params.published_all_fasta}\\n" >> {log}
        echo "> Matched metadata published to {params.published_all_metadata}\\n" >> {log}
        echo '{{"text":"' > 0_data.json
        echo "*Step 0: GISAID processing complete*\\n" >> 0_data.json
        cat {log} >> 0_data.json
        echo '"}}' >> 0_data.json
        echo 'webhook {params.grapevine_webhook}'

        curl -X POST -H "Content-type: application/json" -d @0_data.json {params.grapevine_webhook}
        """

rule alert_sam:
    input: rules.summarize_preprocess_gisaid.log,
    log: config["output_path"] + "/logs/0_alert_sam.log"
    params:
        date = config["date"],
    shell:
        """
        ln -sfn /cephfs/covid/bham/raccoon-dog/{params.date}_gisaid /cephfs/covid/bham/raccoon-dog/gisaid-latest 2> {log}

        ~/.conda/envs/ben-ipc/bin/python /cephfs/covid/software/sam/public/mqtt-message.py -t 'COGUK/infrastructure/housekeeping/gisaid/status' --attr status finished &>> {log}
        """
