import datetime

date = datetime.date.today()

rule update_pangolin:
    output:
        temp("updated_pangolin")
    shell:
        """
        pip install --upgrade git+https://github.com/hCoV-2019/pangolin.git
        touch updated_pangolin
        """

rule uk_extract_new:
    input:
        fasta = rules.uk_filter_low_coverage_sequences.output,
        metadata = rules.uk_remove_duplicates.output.metadata,
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

rule uk_pangolin:
    input:
        fasta = rules.uk_extract_new.output.fasta,
        update = rules.update_pangolin.output
    params:
        outdir = config["output_path"] + "/2/pangolin"
    output:
        lineages = config["output_path"] + "/2/pangolin/lineage_report.csv"
    log:
        config["output_path"] + "/logs/2_uk_pangolin.log"
    threads: 32
    shell:
        """
        pangolin {input.fasta} \
        --threads {threads} \
        --outdir {params.outdir} > {log} 2>&1
        """

rule uk_add_pangolin_lineages_to_metadata:
    input:
        metadata = rules.uk_extract_new.output.metadata,
        lineages = rules.uk_pangolin.output.lineages
    output:
        metadata = config["output_path"] + "/2/uk_with_lineages.csv"
    log:
        config["output_path"] + "/logs/2_uk_add_pangolin_lineages_to_metadata.log"
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.lineages} \
          --index-column sequence_name \
          --join-on taxon \
          --new-columns lineage UFbootstrap \
          --out-metadata {output.metadata} &> {log}
        """

rule uk_combine_previous_and_new:
    input:
        previous_fasta = config["previous_uk_fasta"],
        previous_metadata = config["previous_uk_metadata"],
        new_fasta = rules.uk_extract_new.output.fasta,
        new_metadata = rules.uk_add_pangolin_lineages_to_metadata.output.metadata
    output:
        fasta = config["output_path"] + "/2/uk.combined.fasta",
        metadata = config["output_path"] + "/2/uk.combined.csv"
    log:
        config["output_path"] + "/logs/2_uk_combine_previous_and_new.log"
    shell:
        """
        fastafunk merge \
          --in-fasta {input.previous_fasta} {input.new_fasta} \
          --in-metadata {input.previous_metadata} {input.new_metadata} \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --index-column sequence_name \
          --log-file {log}
        """

rule uk_output_cog:
    input:
        fasta = rules.uk_combine_previous_and_new.output.fasta,
        metadata = rules.uk_combine_previous_and_new.output.metadata
    params:
        outdir = config["publish_path"] + "/COG",
        prefix = config["publish_path"] + "/COG/cog_%s" %s
    output:
        fasta = config["output_path"] + "/2/uk.combined.regularized.fasta",
        metadata = config["output_path"] + "/2/uk.combined.regularized.csv"
    log:
        config["output_path"] + "/logs/2_uk_output_cog.log"
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name collection_date epi_week \
                          country adm1 adm2 outer_postcode \
                          is_surveillance is_community is_hcw \
                          is_travel_history travel_history lineage \
                          lineage_support uk_lineage \
          --where-column epi_week=edin_epi_week country=adm0 lineage_support=ufbootstrap \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --log-file {log} \
          --restrict

        mkdir -p {params.outdir}
        cp {output.fasta} {params.prefix}_alignment.fasta
        cp {output.metadata} {params.prefix}_metadata.fasta
        """