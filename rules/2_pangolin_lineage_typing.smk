import pandas as pd
import sys

rule uk_extract_new:
    input:
        previous_stage = config["output_path"] + "/logs/1_summarize_preprocess_uk.log",
        fasta = rules.uk_filter_low_coverage_sequences.output,
        metadata = rules.add_snp_finder_result_to_metadata.output.metadata,
        previous_metadata = config["previous_uk_metadata"]
    output:
        fasta = config["output_path"] + "/2/uk.new.fasta",
        metadata = config["output_path"] + "/2/uk.new.csv",
    log:
        config["output_path"] + "/logs/2_extract_new.log"
    shell:
        """
        fastafunk new \
          --in-fasta {input.fasta} \
          --in-metadata {input.previous_metadata} {input.metadata} \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --date-column 'edin_date_stamp' \
          --log {log}
        """

# rule update_pangolin_lineages:
#     log:
#         config["output_path"] + "/logs/2_update_pangolin_lineages.log"
#     shell:
#         """
#         #pip install --upgrade git+https://github.com/hCoV-2019/spangolin.git
#         """

rule uk_pangolin:
    input:
        fasta = rules.uk_extract_new.output.fasta,
        pangolin_updated = rules.update_pangolin_lineages.log
    params:
        outdir = config["output_path"] + "/2/pangolin",
        tmpdir = config["output_path"] + "/2/pangolin/tmp"
    output:
        lineages = protected(config["output_path"] + "/2/pangolin/lineage_report.csv")
    log:
        config["output_path"] + "/logs/2_uk_pangolin.log"
    threads: 16
    shell:
        """
        pangolin {input.fasta} \
        --threads {threads} \
        --outdir {params.outdir} \
        --tempdir {params.tmpdir}  >> {log} 2>&1
        """

        # pangolin --lineages-version >> {log}
        # pangolin --version >> {log}

rule uk_add_previous_uk_lineages_to_metadata:
    input:
        previous_metadata = config["previous_uk_metadata"],
        metadata = rules.add_snp_finder_result_to_metadata.output.metadata,
    output:
        metadata = config["output_path"] + "/2/uk.with_previous_lineages.csv",
    log:
        config["output_path"] + "/logs/2_uk_add_previous_uk_lineages_to_metadata.log"
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.previous_metadata} \
          --index-column sequence_name \
          --join-on sequence_name \
          --new-columns uk_lineage special_lineage lineage lineage_support edin_date_stamp \
          --out-metadata {output.metadata} &>> {log}
        """

rule uk_add_pangolin_lineages_to_metadata:
    input:
        metadata = rules.uk_add_previous_uk_lineages_to_metadata.output.metadata,
        new_metadata = rules.uk_extract_new.output.metadata,
        lineages = rules.uk_pangolin.output.lineages
    output:
        metadata = config["output_path"] + "/2/uk.with_new_lineages.csv",
        tmp_metadata = temp(config["output_path"] + "/2/uk.with_new_lineages.csv.tmp")
    log:
        config["output_path"] + "/logs/2_uk_add_pangolin_lineages_to_metadata.log"
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.new_metadata} \
          --index-column sequence_name \
          --join-on sequence_name \
          --new-columns edin_date_stamp \
          --out-metadata {output.tmp_metadata} &> {log}

        fastafunk add_columns \
          --in-metadata {output.tmp_metadata} \
          --in-data {input.lineages} \
          --index-column sequence_name \
          --join-on taxon \
          --new-columns special_lineage lineage_support \
          --where-column lineage_support=UFbootstrap special_lineage=lineage\
          --out-metadata {output.metadata} &>> {log}
        """

rule uk_update_metadata_lineages:
    input:
        metadata = rules.uk_add_pangolin_lineages_to_metadata.output.metadata
    output:
        metadata = config["output_path"] + "/2/uk.with_new_lineages.special.csv",
    log:
        config["output_path"] + "/logs/2_uk_update_metadata_lineages.log"
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


"""
Let's output a full-width metadata table here, but restrict it (in length) to
sequences in the fasta file that we're going to build trees on

NB do we need to check that every row has 'special lineage' at this point?
"""
rule uk_output_lineage_table:
    input:
        fasta = rules.uk_filter_low_coverage_sequences.output.fasta,
        metadata = rules.uk_update_metadata_lineages.output.metadata
    output:
        fasta = config["output_path"] + "/2/uk.matched.fasta",
        metadata = config["output_path"] + "/2/uk.matched.special_linages.csv"
    log:
        config["output_path"] + "/logs/2_uk_output_lineage_table.log"
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name  special_lineage \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --log-file {log} \
          --restrict
        """


rule summarize_pangolin_lineage_typing:
    input:
        fasta = rules.uk_output_cog.output.fasta,
        metadata = rules.uk_output_cog.output.metadata,
        # public_fasta = rules.uk_output_cog_public.output.fasta,
        # public_metadata = rules.uk_output_cog_public.output.metadata,
        # full_metadata = rules.uk_update_metadata_lineages.output.metadata
    params:
        webhook = config["webhook"],
        outdir = config["publish_path"] + "/COG",
        prefix = config["publish_path"] + "/COG/cog",
        export_dir1 = config["export_path"] + "/public",
        export_prefix1 = config["export_path"] + "/public/cog_" + config["date"],
        export_dir2 = config["export_path"] + "/alignments",
        export_prefix2 = config["export_path"] + "/alignments/cog_" + config["date"]
    log:
        config["output_path"] + "/logs/2_summarize_pangolin_lineage_typing.log"
    shell:
        """
        echo '{{"text":"' > 2_data.json
        echo "*Step 2: COG-UK pangolin typing complete*\\n" >> 2_data.json
        cat {log} >> 2_data.json
        echo '"}}' >> 2_data.json
        echo "webhook {params.webhook}"
        curl -X POST -H "Content-type: application/json" -d @2_data.json {params.webhook}
        #rm 2_data.json
        """

        # mkdir -p {params.outdir}
        # mkdir -p {params.export_dir1}
        # mkdir -p {params.export_dir2}
        #
        # cp {input.full_metadata} {params.prefix}_metadata.full.csv
        # echo "> Full COG metadata published to _{params.prefix}_metadata.csv_\\n" >> {log}
        # echo ">\\n" >> {log}

        # cp {input.fasta} {params.prefix}_alignment.trimmed.matched.fasta
        # cp {input.metadata} {params.prefix}_metadata.matched.csv
        # cp {input.fasta} {params.export_prefix2}_alignment.fasta
        # cp {input.metadata} {params.export_prefix2}_metadata.csv
        # echo "> Matched COG trimmed alignment and restricted metadata published to _{params.prefix}_alignment.trimmed.matched.fasta_ and _{params.prefix}_metadata.matched.csv_\\n" >> {log}
        # echo "> and to _{params.export_prefix2}_alignment.fasta_ and _{params.export_prefix2}_metadata.csv_\\n" >> {log}
        # echo "> Number of sequences in matched COG files: $(cat {input.fasta} | grep ">" | wc -l)\\n" &>> {log}
        # echo ">\\n" >> {log}

        # cp {input.public_fasta} {params.prefix}_sequences.public.fasta
        # cp {input.public_metadata} {params.prefix}_metadata.public.csv
        # cp {input.public_fasta} {params.export_prefix1}_sequences.fasta
        # cp {input.public_metadata} {params.export_prefix1}_metadata.csv
        # echo "> Public unaligned COG fasta and restricted metadata published to _{params.prefix}_sequences.public.fasta_ and _{params.prefix}_metadata.public.csv_\\n" >> {log}
        # echo "> and to _{params.export_prefix1}_sequences.fasta_ and _{params.export_prefix1}_metadata.csv_\\n" >> {log}
        # echo "> Number of sequences in public COG files: $(cat {input.public_fasta} | grep ">" | wc -l)\\n" &>> {log}
        # echo ">\\n" >> {log}
