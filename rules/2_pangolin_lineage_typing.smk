rule update_pangolin:
    output:
        temp("updated_pangolin")
    shell:
        """
        pip install --upgrade git+https://github.com/hCoV-2019/pangolin.git
        touch updated_pangolin
        """

rule uk_pangolin:
    input:
        fasta = rules.uk_filter_low_coverage_sequences.output,
        update = rules.update_pangolin.output
    params:
        outdir = config["output_path"] + "/2/pangolin"
    output:
        config["output_path"] + "/2/pangolin/lineage_report.csv"
    log:
        config["output_path"] + "/logs/2_uk_pangolin.log"
    threads: 8
    shell:
        """
        pangolin {input.fasta} \
        --threads {threads} \
        --outdir {params.outdir}
        """

# rule uk_add_pangolin_lineages_to_metadata:
#