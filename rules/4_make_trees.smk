import datetime

date = datetime.date.today()

rule split_based_on_lineages:
    input:
        fasta = rules.combine_gisaid_and_cog.output.fasta,
        metadata = rules.combine_gisaid_and_cog.output.metadata,
        lineage = config["lineage_splits"]
    params:
        prefix = config["output_path"] + "/4/lineage_"
    output:
        temp(config["output_path"] + "/4/split_done")
    log:
        config["output_path"] + "/logs/4_split_based_on_lineages.log"
    shell:
        """
        lineages=$(cat {input.lineage} | cut -f1 --delim "," | tr '\n' '  ')
        fastafunk split \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --index-field lineage \
          --lineage $lineages \
          --out-folder {params.prefix} &> {log}

        touch {output}
        """

rule run_subroutine_on_lineage:
    input:
        split_done = rules.split_based_on_lineages.output,
        metadata = rules.combine_gisaid_and_cog.output.metadata,
        lineage = config["lineage_splits"]
    params:
        path_to_script = workflow.current_basedir,
        output_path = config["output_path"],
        prefix = config["output_path"] + "/4/lineage_"

    log:
        config["output_path"] + "/logs/4_run_subroutine_on_lineage.log"
    shell:
        """
        while IFS=, read -r lineage lineage_specific_outgroup
        do
          snakemake --nolock \
            --snakefile {params.path_to_script}/rules/4_subroutine/process_lineage.smk \
            --cores 8 \
            --configfile {params.path_to_script}/rules/4_subroutine/ \
            --config \
            output_path={params.output_path} \
            lineage_fasta={params.prefix}$lineage.fasta \
            lineage=$lineage \
            lineage_specific_outgroup=$lineage_specific_outgroup \
            metadata={input.metadata}
        done &> {log}
        """