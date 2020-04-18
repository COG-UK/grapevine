import datetime

date = datetime.date.today()

rule split_based_on_lineages:
    input:
        fasta = rules.combine_gisaid_and_cog.output.fasta,
        metadata = rules.combine_gisaid_and_cog.output.metadata,
        lineage = config["lineage_splits"]
    params:
        path_to_script = workflow.current_basedir,
        output_path = config["output_path"],

    output:
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
          --lineage $lineages &> {log}


        snakemake --nolock --snakefile {params.path_to_script}/rules/4_subroutine/process_lineage.smk \
        --cores 8
        --configfile {params.path_to_script}/rules/4_subroutine/ \
        --config \
        output_path={params.output_path} \
        lineage_fasta={params.min_length} \
        lineage={params.max_length} \
        lineage_specific_outgroup={params.kraken_fasta} \
        metadata={input.metadata} &>> {log}
        """