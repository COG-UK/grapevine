rule split_based_on_lineages:
    input:
        previous_stage = config["output_path"] + "/logs/3_summarize_combine_gisaid_and_cog.log",
        fasta = config["output_path"] + "/3/cog_gisaid.fasta",
        metadata = config["output_path"] + "/3/cog_gisaid.lineages.csv",
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
        fastafunk split \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --index-field special_lineage \
          --lineage-csv {input.lineage} \
          --out-folder {params.prefix} &> {log}

        echo {params.webhook}

        echo '{{"text":"' > 4a_data.json
        echo "*Step 4: Ready for tree building*\\n" >> 4a_data.json
        num_lineages=$(cat {input.lineage} | wc -l)
        range={{$num_lineages..1}}
        for i in $(eval echo ${{range}})
        do
          line=$(tail -n$i {log} | head -n1)
          echo ">$line\\n" >> 4a_data.json
        done
        echo '"}}' >> 4a_data.json
        echo "webhook {params.webhook}"
        curl -X POST -H "Content-type: application/json" -d @4a_data.json {params.webhook}
        #rm 4a_data.json

        touch {output}
        """

rule run_4_subroutine_on_lineages:
    input:
        split_done = rules.split_based_on_lineages.output,
        metadata = config["output_path"] + "/3/cog_gisaid.csv",
        lineage = config["lineage_splits"]
    params:
        path_to_script = workflow.current_basedir,
        output_path = config["output_path"],
        publish_path = config["publish_path"],
        export_path = config["export_path"],
        date = config["date"],
        guide_tree = config["guide_tree"],
        prefix = config["output_path"] + "/4/lineage_"
    output:
        traits = config["output_path"] + "/4/all_traits.csv",
        public_tree = config["output_path"] + "/4/cog_gisaid_full.tree.public.newick"
        private_tree = config["output_path"] + "/4/cog_gisaid_full.tree.private.nexus"
    log:
        config["output_path"] + "/logs/4_run_4_subroutine_on_lineages.log"
    threads: 16
    shell:
        """
        lineages=$(cat {input.lineage} | cut -f1 -d"," | tr '\\n' '  ')
        outgroups=$(cat {input.lineage} | cut -f2 -d"," | tr '\\n' '  ')
        snakemake --nolock \
          --snakefile {params.path_to_script}/4_subroutine/4_process_lineages.smk \
          --cores {threads} \
          --configfile {params.path_to_script}/4_subroutine/config.yaml \
          --config \
          output_path={params.output_path} \
          publish_path={params.publish_path} \
          lineages="$lineages" \
          lineage_specific_outgroups="$outgroups" \
          guide_tree="{params.guide_tree}" \
          metadata={input.metadata} &> {log}
        """

rule summarize_make_trees:
    input:
        traits = rules.run_4_subroutine_on_lineages.output.traits,
        public_tree = rules.run_4_subroutine_on_lineages.output.public_tree
        private_tree = rules.run_4_subroutine_on_lineages.output.private_tree
    output:
        published_tree = config["publish_path"] + "/COG_GISAID/cog_gisaid_full_tree.nexus",
        exported_tree1 = config["export_path"] + "/trees/cog_global_" + config["date"] + "_tree.nexus"
        exported_tree2 = config["export_path"] + "/public/cog_global_" + config["date"] + "_tree.newick",
    params:
        webhook = config["webhook"],
        outdir = config["publish_path"] + "/COG_GISAID",
    log:
        config["output_path"] + "/logs/4_summarize_make_trees.log"
    shell:
        """
        echo "> Lineage trees have been published in _{params.outdir}_\\n" >> {log}
        echo ">\\n" >> {log}

        cp {input.private_tree} {output.published_tree}
        cp {input.private_tree} {output.exported_tree1}
        cp {input.public_tree} {output.exported_tree2}
        echo "> Full annotated tree has been published in _{output.published_tree}_\\n" >> {log}
        echo "> and _{output.exported_tree1}_\\n" >> {log}
        echo "> Full unannotated tree has been published in _{output.exported_tree2}_\\n" >> {log}

        echo '{{"text":"' > 4b_data.json
        echo "*Step 4: Construct and annotate lineage trees completed*\\n" >> 4_data.json
        cat {log} >> 4b_data.json
        echo '"}}' >> 4b_data.json
        echo "webhook {params.webhook}"
        curl -X POST -H "Content-type: application/json" -d @4b_data.json {params.webhook}
        #rm 4b_data.json
        """
