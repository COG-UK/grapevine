
rule gisaid_process_json:
    input:
        json = config["latest_gisaid_json"],
        omitted = config["previous_omitted_file"]
        # metadata = config["previous_gisaid_metadata"]
    output:
        fasta = config["output_path"] + "/0/gisaid.fasta",
        metadata = config["output_path"] + "/0/gisaid.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_process_json.log"
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
        metadata = rules.gisaid_process_json.output.metadata
    output:
        fasta = config["output_path"] + "/0/gisaid.UH.fasta",
        metadata = config["output_path"] + "/0/gisaid.UH.csv"
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


rule gisaid_remove_duplicates:
    input:
        fasta = rules.gisaid_unify_headers.output.fasta,
        metadata = rules.gisaid_unify_headers.output.metadata
    output:
        fasta = config["output_path"] + "/0/gisaid.UH.RD.fasta",
        metadata = config["output_path"] + "/0/gisaid.UH.RD.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_filter_duplicates.log"
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
        fasta = config["output_path"] + "/0/gisaid.RD.UH.filt1.mapped.fasta"
    log:
        config["output_path"] + "/logs/0_gisaid_remove_insertions_and_pad.log"
    run:
        shell(
        """
        datafunk sam_2_fasta \
          -s {input.sam} \
          -r {input.reference} \
          -o {output} \
          -t [{params.trim_start}:{params.trim_end}] \
          --pad \
          --log-inserts &> {log} """)
        if os.path.exists("insertions.txt"):
            shell("""
            mv insertions.txt {params.insertions}
            """)



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



rule gisaid_snp_finder:
    input:
        fasta = rules.gisaid_filter_on_distance_to_WH04.output.fasta,
        snps = config["snps"]
    output:
        found = config["output_path"] + "/0/gisaid.snp_finder.csv",
    log:
        config["output_path"] + "/logs/0_gisaid_snp_finder.log"
    shell:
        """
        datafunk snp_finder -a {input.fasta} -o {output.found} --snp-csv {input.snps} &> {log}
        """


rule gisaid_add_snp_finder_result_to_metadata:
    input:
        snps = config["snps"],
        metadata = rules.gisaid_remove_duplicates.output.metadata,
        new_data = rules.gisaid_snp_finder.output.found
    output:
        metadata = config["output_path"] + "/0/gisaid.RD.UH.SNPfinder.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_add_snp_finder_result_to_metadata.log"
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


rule gisaid_extract_lineageless:
    input:
        fasta = rules.gisaid_filter_on_distance_to_WH04.output,
        metadata = rules.gisaid_add_snp_finder_result_to_metadata.output.metadata,
    output:
        fasta = config["output_path"] + "/0/gisaid.new.pangolin_lineages.fasta",
    log:
        config["output_path"] + "/logs/0_extract_lineageless.log"
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
    threads: 40
    shell:
        """
        pangolin {input.fasta} -p -t{threads} --tempdir {params.tmpdir} --outdir {params.outdir} --verbose > {log} 2>&1
        """


rule gisaid_add_pangolin_lineages_to_metadata:
    input:
        metadata = rules.gisaid_add_snp_finder_result_to_metadata.output.metadata,
        normal_lineages = rules.gisaid_normal_pangolin.output.lineages
    output:
        metadata = config["output_path"] + "/0/gisaid.RD.UH.SNPfinder.lineages.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_add_pangolin_lineages_to_metadata.log"
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
        metadata = config["output_path"] + "/0/gisaid.all.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_add_del_finder_result_to_metadata.log"
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


rule gisaid_output_lineage_table:
    input:
        fasta = rules.gisaid_filter_on_distance_to_WH04.output.fasta,
        metadata = rules.gisaid_add_del_finder_result_to_metadata.output.metadata
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
                          lineage uk_lineage \
          --where-column uk_omit=is_uk sample_date=covv_collection_date \
                                 epi_week=edin_epi_week country=edin_admin_0 \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --log-file {log} \
          --restrict
        """


rule gisaid_output_matched_fasta_and_metadata_table:
    input:
        fasta = rules.gisaid_filter_on_distance_to_WH04.output.fasta,
        metadata = rules.gisaid_add_del_finder_result_to_metadata.output.metadata
    output:
        published_fasta = config["publish_path"] + "/GISAID/gisaid.trimmed_alignment.fasta",
        published_metadata = config["publish_path"] + "/GISAID/gisaid.metadata.csv",
    log:
        config["output_path"] + "/logs/0_gisaid_output_matched_fasta_and_metadata_table.log"
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
          --out-fasta {output.published_fasta} \
          --out-metadata {output.published_metadata} \
          --log-file {log} \
          --restrict
        """


rule summarize_preprocess_gisaid:
    input:
        latest_fasta = rules.gisaid_process_json.output.fasta,
        deduplicated_fasta = rules.gisaid_remove_duplicates.output.fasta,
        unify_headers_fasta = rules.gisaid_unify_headers.output.fasta,
        removed_short_fasta = rules.gisaid_filter_1.output,
        removed_low_covg_fasta = rules.gisaid_filter_2.output.fasta,
        removed_distance_to_root_fasta = rules.gisaid_filter_on_distance_to_WH04.output.fasta,
        full_metadata = rules.gisaid_add_pangolin_lineages_to_metadata.output.metadata,
        matched_fasta = rules.gisaid_output_lineage_table.output.fasta,
        matched_lineage_table = rules.gisaid_output_lineage_table.output.metadata,
        counts = rules.gisaid_counts_by_country.output.counts,
        published_fasta = rules.gisaid_output_matched_fasta_and_metadata_table.output.published_fasta,
        published_metadata = rules.gisaid_output_matched_fasta_and_metadata_table.output.published_metadata,
    params:
        publish_path = config["publish_path"] + "/GISAID/",
        published_counts = config["publish_path"] + "/GISAID/gisaid_counts_by_country.csv",
        #grapevine_webhook = config["grapevine_webhook"]
    log:
        config["output_path"] + "/logs/0_summarize_preprocess_gisaid.log"
    shell:
        """
        mkdir -p {params.publish_path}

        cp {input.counts} {params.published_counts}

        echo "Number of sequences in latest GISAID download: $(cat {input.latest_fasta} | grep '>' | wc -l)\\n" >> {log}
        echo "Number of sequences after unifying headers: $(cat {input.unify_headers_fasta} | grep '>' | wc -l)\\n" >> {log}
        echo "Number of sequence after deduplicating: $(cat {input.deduplicated_fasta} | grep '>' | wc -l)\\n" >> {log}
        echo "Number of sequences after removing sequences <29000bps: $(cat {input.removed_short_fasta} | grep '>' | wc -l)\\n" >> {log}
        echo "Number of sequences after mapping and removing those with <95% coverage: $(cat {input.removed_low_covg_fasta} | grep '>' | wc -l)\\n" >> {log}
        echo "Number of sequences after removing those >4 epi-week stddevs to WH04: $(cat {input.removed_distance_to_root_fasta} | grep '>' | wc -l)\\n" >> {log}
        echo "> \\n" >> {log}
        echo "> Counts by country published to {params.published_counts}\\n" >> {log}
        echo "> Trimmed alignment published to {input.published_fasta}\\n" >> {log}
        echo "> Matched metadata published to {input.published_metadata}\\n" >> {log}
        echo '{{"text":"' > 0_data.json
        echo "*Step 0: GISAID preprocessing complete*\\n" >> 0_data.json
        cat {log} >> 0_data.json
        echo '"}}' >> 0_data.json
        # rm 0_data.json
        """

rule reference_tree_mid:
    input:
        fasta = rules.gisaid_output_matched_fasta_and_metadata_table.output.published_fasta
    output:
        tree = config["output_path"] + "/0.5/gisaid.trimmed_alignment.1.multi.fasttree"
    shell:
        """
        fasttree \
        -nt -gamma \
        -nosupport \
        -sprlength 500 \
        -nni 0 -spr 5 \
        -refresh 0.8 \
        -topm 1.5 \
        -close 0.75 -noml \
        {input.fasta:q} > {output.tree:q}
        """

rule reference_tree:
    input:
        tree = rules.reference_tree_mid.output.tree,
        fasta = rules.gisaid_output_matched_fasta_and_metadata_table.output.published_fasta
    output:
        tree = config["output_path"] + "/0.5/gisaid.trimmed_alignment.multi.fasttree"
    shell:
        """
        fasttree \
        -nt -gamma -sprlength 200 -spr 5 \
        -intree {input.tree} \
        {input.fasta:q} > {output.tree:q}
        """

rule build_sequence_bootstraps:
    input:
        fasta = rules.gisaid_output_matched_fasta_and_metadata_table.output.published_fasta
    params:
        bootprefix = config["output_path"] + "/0.5/gisaid.trimmed_alignment.boot{bs}"
    output:
        bootstrap = config["output_path"] + "/0.5/gisaid.trimmed_alignment.boot{bs}0.fa"
    shell:
        """
        goalign build seqboot -i {input.fasta:q} -t 1 -n 1 -S -o {params.bootprefix}
        """

rule run_bootstraps:
    input:
        bootstrap = config["output_path"] + "/0.5/gisaid.trimmed_alignment.boot{bs}0.fa",
    output:
        tree = config["output_path"] + "/0.5/gisaid.trimmed_alignment.boot{bs}.unrooted.tree"
    shell:
        """
        fasttree -nosupport -nt -fastest {input.bootstrap} > {output.tree:q}
        """

rule gather_bootstraps:
    input:
        expand(config["output_path"] + "/0.5/gisaid.trimmed_alignment.boot{bs}.unrooted.tree", bs = range(100))
    output:
        trees = config["output_path"] + "/0.5/gisaid.trimmed_alignment.ft_replicates.multi.tree"
    shell:
        """
        cat {input} > {output.trees:q}
        """

rule compute_tbe:
    input:
        trees = rules.gather_bootstraps.output.trees,
        tree = rules.reference_tree.output.tree
    threads: workflow.cores
    output:
        tree = config["output_path"] + "/0.5/gisaid.trimmed_alignment.TBE.unrooted.tree"
    shell:
        """
        gotree compute support tbe \
        -i {input.tree:q} \
        -b {input.trees:q} \
        -t {threads} \
        -o {output.tree}
        """

rule reroot:
    input:
        tree = rules.compute_tbe.output.tree
    params:
        outgroup = config["outgroup"]
    output:
        tree = config["output_path"] + "/0.5/gisaid.trimmed_alignment.TBE.tree"
    shell:
        """
        clusterfunk root -i {input.tree} --in-format newick --out-format newick -o {output.tree} --outgroup {params.outgroup}
        """

