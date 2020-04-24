import datetime

date = datetime.date.today()

rule split_based_on_lineages:
    input:
        fasta = rules.combine_gisaid_and_cog.output.fasta,
        metadata = rules.combine_gisaid_and_cog.output.metadata,
        lineage = config["lineage_splits"]
    params:
        prefix = config["output_path"] + "/4/lineage_",
        webhook = config["webhook"]
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

        echo '{{"text":"' > 4_data.json
        echo "*Step 4: Ready for tree building*\n" >> 4_data.json
        num_lineages=$(cat {input.lineage} | wc -l)
        num_lineages=$((num_lineages+1))
        tail -n$num_lineages {log} >> 4_data.json
        echo '"}}' >> 4_data.json
        curl -X POST -H "Content-type: application/json" -d @4_data.json {params.webhook}
        rm 4_data.json

        touch {output}
        """

rule run_subroutine_on_lineages:
    input:
        split_done = rules.split_based_on_lineages.output,
        metadata = rules.combine_gisaid_and_cog.output.metadata,
        lineage = config["lineage_splits"]
    params:
        path_to_script = workflow.current_basedir,
        output_path = config["output_path"],
        publish_path = config["publish_path"],
        prefix = config["output_path"] + "/4/lineage_"
    output:
        config["output_path"] + "/4/trees_done"
    log:
        config["output_path"] + "/logs/4_run_subroutine_on_lineages.log"
    threads: 16
    shell:
        """
        lineages=$(cat {input.lineage} | cut -f1 -d"," | tr '\n' '  ')
        outgroups=$(cat {input.lineage} | cut -f2 -d"," | tr '\n' '  ')
        snakemake --nolock \
          --snakefile {params.path_to_script}/4_subroutine/process_lineages.smk \
          --cores {threads} \
          --configfile {params.path_to_script}/4_subroutine/config.yaml \
          --config \
          output_path={params.output_path} \
          publish_path={params.publish_path} \
          lineages="$lineages" \
          lineage_specific_outgroups="$outgroups" \
          metadata={input.metadata} &> {log}

        touch {output}
        """

rule summarize_make_trees:
    input:
        lineage = config["lineage_splits"],
        trees_done = rules.run_subroutine_on_lineages.output,
    params:
        webhook = config["webhook"],
        outdir = config["publish_path"] + "/COG_GISAID",
    log:
        config["output_path"] + "/logs/4_summarize_make_trees.log"
    shell:
        """
        echo "> Trees have been published in {params.outdir}\n" >> {log}

        echo '{{"text":"' > 4_data.json
        echo "*Step 4: Construct and annotate trees completed*\n" >> 4_data.json
        cat {log} >> 4_data.json
        echo '"}}' >> 4_data.json
        curl -X POST -H "Content-type: application/json" -d @4_data.json {params.webhook}
        rm 4_data.json
        """