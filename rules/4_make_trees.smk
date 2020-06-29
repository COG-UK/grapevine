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
          --index-field lineage \
          --lineage-csv {input.lineage} \
          --out-folder {params.prefix} &> {log}

        # echo {params.webhook}

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
        # echo "webhook {params.webhook}"

        touch {output}
        curl -X POST -H "Content-type: application/json" -d @4a_data.json {params.webhook}
        """


rule run_4_subroutine_on_lineages:
    input:
        split_done = rules.split_based_on_lineages.log,
        metadata = config["output_path"] + "/3/cog_gisaid.lineages.csv",
        lineage = config["lineage_splits"],
    params:
        path_to_script = workflow.current_basedir,
        output_path = config["output_path"],
        publish_path = config["publish_path"],
        export_path = config["export_path"],
        date = config["date"],
        guide_tree = config["guide_tree"],
        prefix = config["output_path"] + "/4/lineage_"
    output:
        grafted_tree = config["output_path"] + '/4/cog_gisaid_grafted.tree'
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
          metadata="{input.metadata}" &> {log}
        """

rule sort:
    input:
        grafted_tree = rules.run_4_subroutine_on_lineages.output.grafted_tree,
    output:
        sorted_tree = config["output_path"] + '/4/cog_gisaid_full.tree.public.newick',
    log:
        config["output_path"] + "/logs/4_sort.log",
    shell:
        """
        gotree rotate sort -i {input.grafted_tree} -o {output.sorted_tree} &> {log}
        """

rule step_4_annotate_tree:
    input:
        tree = rules.sort.output.sorted_tree,
        metadata = config["output_path"] + "/3/cog_gisaid.lineages.csv",
    params:
        collapse=0.000005,
    output:
        tree = config["output_path"] + "/4/cog_gisaid_grafted.annotated.tree"
    log:
        config["output_path"] + "/logs/4_annotate_tree.log"
    shell:
        """
        clusterfunk annotate_tips \
          --in-metadata {input.metadata} \
          --trait-columns country lineage uk_lineage \
          --index-column sequence_name \
          --boolean-for-trait country='UK' country='UK' country='UK' country='UK' \
          --boolean-trait-names country_uk country_uk_acctran country_uk_deltran \
          --in-format newick \
          --out-format nexus \
          --collapse_to_polytomies {params.collapse} \
          --input {input.tree} \
          --output {output.tree} &> {log}
        """

rule acctran_ancestral_reconstruction:
    input:
        tree = rules.step_4_annotate_tree.output.tree
    output:
        tree = config["output_path"] + "/4/cog_gisaid_grafted.annotated.acc.tree"
    log:
        config["output_path"] + "/logs/4_acctran_ancestral_reconstruction.log"
    shell:
        """
        clusterfunk ancestral_reconstruction \
        --traits country_uk_acctran \
        --acctran \
        --ancestral-state False \
        --input {input.tree} \
        --output {output.tree} &> {log}
        """

rule deltran_ancestral_reconstruction:
    input:
        tree = rules.acctran_ancestral_reconstruction.output.tree
    output:
        tree = config["output_path"] + "/4/cog_gisaid_grafted.annotated.acc.del.tree"
    log:
        config["output_path"] + "/logs/4_deltran_ancestral_reconstruction.log"
    shell:
        """
        clusterfunk ancestral_reconstruction \
        --traits country_uk_deltran \
        --deltran \
        --ancestral-state False \
        --input {input.tree} \
        --output {output.tree} &> {log}
        """

rule label_acctran_introductions:
    input:
        tree = rules.deltran_ancestral_reconstruction.output.tree
    output:
        tree = config["output_path"] + "/4/cog_gisaid_grafted.annotated.acc.del.acc_labelled.tree"
    log:
        config["output_path"] + "/logs/4_label_acctran_introductions.log"
    shell:
        """
        clusterfunk label_transitions \
          --trait country_uk_acctran \
          --to True \
          --transition-name acc_introduction \
          --transition-prefix acc_trans_ \
          --include_root \
          --stubborn \
          --input {input.tree} \
          --output {output.tree} &> {log}
        """

rule label_deltran_introductions:
    input:
        tree = rules.label_acctran_introductions.output.tree
    output:
        tree = config["output_path"] + "/4/cog_gisaid_grafted.annotated.acc.del.acc_labelled.del_labelled.tree"
    log:
        config["output_path"] + "/logs/4_label_deltran_introductions.log"
    shell:
        """
        clusterfunk label_transitions \
          --trait country_uk_deltran \
          --to True \
          --transition-name del_introduction \
          --transition-prefix del_trans_ \
          --stubborn \
          --input {input.tree} \
          --output {output.tree} &> {log}
        """

rule merge_sibling_acc_introduction:
    input:
        tree = rules.label_deltran_introductions.output.tree
    output:
        tree = config["output_path"] + "/4/cog_gisaid_grafted.annotated.acc.del.acc_labelled.del_labelled.acc_merged.tree"
    log:
        config["output_path"] + "/logs/4_merge_acctran_introductions.log"
    shell:
        """
        clusterfunk merge_transitions \
          --trait-to-merge acc_introduction \
          --merged-trait-name acc_lineage \
          --merge-siblings \
          --input {input.tree} \
          --output {output.tree} &> {log}
        """

rule merge_sibling_del_introduction:
    input:
        tree = rules.merge_sibling_acc_introduction.output.tree
    params:
        outdir = config["publish_path"] + "/COG_GISAID",
    output:
        tree = config["output_path"] + "/4/cog_gisaid_grafted.annotated.acc.del.acc_labelled.del_labelled.acc_merged.del_merged.tree"
    log:
        config["output_path"] + "/logs/4_merge_deltran_introductions.log"
    shell:
        """
        clusterfunk merge_transitions \
          --trait-to-merge del_introduction \
          --merged-trait-name del_lineage \
          --merge-siblings \
          --input {input.tree} \
          --output {output.tree} &> {log}

        mkdir -p {params.outdir}
        """


rule output_annotations:
    input:
        tree = rules.merge_sibling_del_introduction.output.tree,
    output:
        traits = config["output_path"] + "/4/all_traits.csv"
    log:
        config["output_path"] + "/logs/4_output_annotations.log"
    shell:
        """
        clusterfunk extract_tip_annotations \
          --traits country uk_lineage acc_introduction acc_lineage del_introduction del_lineage \
          --input {input.tree} \
          --output {output.traits} &> {log}
        """


rule summarize_make_trees:
    input:
        traits = rules.output_annotations.output.traits,
        public_tree = rules.sort.output.sorted_tree,
    params:
        webhook = config["webhook"],
    log:
        config["output_path"] + "/logs/4_summarize_make_trees.log"
    shell:
        """
        echo '{{"text":"' > 4b_data.json
        echo "*Step 4: Construct and annotate lineage trees completed*\\n" >> 4b_data.json
        echo '"}}' >> 4b_data.json
        echo "webhook {params.webhook}"
        curl -X POST -H "Content-type: application/json" -d @4b_data.json {params.webhook}
        """
