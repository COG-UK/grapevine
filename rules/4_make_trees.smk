import os

##### Configuration from original subsnake #####
config["publish_path"] = os.path.abspath(config["publish_path"])
lineage_to_outgroup_map = {}
LINEAGES = []
OUTGROUPS = []
with open(config["lineage_splits"]) as lineage_fh:
    # Read the lineage file into LINEAGES and OUTGROUPS
    for i, line in enumerate(lineage_fh):
        lineage, outgroup = line.strip().split(',')

        if i == 0:
            continue
        LINEAGES.append(lineage)
        OUTGROUPS.append(outgroup)
        lineage_to_outgroup_map[lineage] = outgroup

rule split_based_on_lineages:
    input:
        previous_stage = config["output_path"] + "/logs/3_summarize_combine_gisaid_and_cog.log",
        fasta = config["output_path"] + "/3/cog_gisaid.fasta",
        metadata = config["output_path"] + "/3/cog_gisaid.lineages.csv",
        lineage = config["lineage_splits"],
        aliases = config["lineage_aliases"]
    params:
        prefix = config["output_path"] + "/4/lineage_",
        grapevine_webhook = config["grapevine_webhook"],
        json_path = config["json_path"],
    output:
        done=temp(config["output_path"] + "/4/split_done"),
        lineage_fastas = [config["output_path"] + "/4/lineage_%s.fasta" % l for l in LINEAGES],
    resources: mem_per_cpu=20000
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
          --aliases {input.aliases} \
          --out-folder {params.prefix} &> {log}

        echo '{{"text":"' > {params.json_path}/4a_data.json
        echo "*Step 4: Ready for tree building*\\n" >> {params.json_path}/4a_data.json
        num_lineages=$(cat {input.lineage} | wc -l)
        range={{$num_lineages..1}}
        for i in $(eval echo ${{range}})
        do
          line=$(tail -n$i {log} | head -n1)
          echo ">$line\\n" >> {params.json_path}/4a_data.json
        done
        echo '"}}' >> {params.json_path}/4a_data.json
        echo "webhook {params.grapevine_webhook}"

        touch {output.done}
        curl -X POST -H "Content-type: application/json" -d @{params.json_path}/4a_data.json {params.grapevine_webhook}
        """


rule fasttree:
    input:
        split_done = rules.split_based_on_lineages.log,
        lineage_fasta = config["output_path"] + "/4/lineage_{lineage}.fasta"
    params:
        lineage = "{lineage}",
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.unrooted.tree",
    log:
        config["output_path"] + "/logs/4_fasttree_{lineage}.log"
    resources: mem_per_cpu=20000 # 10gb * 4 threads * 4 lineages = 160gb < 192gb
    threads: 3
    shell:
        """
        echo "{input.lineage_fasta} {params.lineage}"

        export OMP_NUM_THREADS={threads}

        FastTreeMP -nosupport -nt {input.lineage_fasta} > {output.tree} 2> {log}
        """


rule root_tree:
    input:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.unrooted.tree",
    params:
        lineage = "{lineage}",
        outgroup = lambda wildcards: lineage_to_outgroup_map[wildcards.lineage]
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.tree",
    log:
        config["output_path"] + "/logs/4_root_tree_{lineage}.log",
    resources: mem_per_cpu=20000
    shell:
        """
        echo "{params.lineage} {params.outgroup}"

        clusterfunk root \
          --in-format newick \
          -i {input.tree} \
          --out-format newick \
          -o {output.tree} \
          --outgroup {params.outgroup} &>> {log}
        """


rule graft:
    input:
        # not sure how to pass this as a space separated list below. Also assuming the order here matches lineages
        scions = expand(config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.tree", lineage=sorted(LINEAGES)),
        guide_tree = config["guide_tree"]
    params:
        lineages = sorted(LINEAGES),
    output:
        grafted_tree = config["output_path"] + '/4/cog_gisaid_grafted.tree',
    log:
        config["output_path"] + "/logs/4_graft.log",
    resources: mem_per_cpu=20000
    run:
        if len(input.scions) > 1:
            shell("""
                clusterfunk graft \
                  --scions {input.scions} \
                  --scion-annotation-name scion_lineage \
                  --annotate-scions {params.lineages} \
                  --input {input.guide_tree} \
                  --in-format newick \
                  --out-format newick \
                  --output {output.grafted_tree} &> {log}
                """)

        else:
            shell("""
                clusterfunk reformat \
                  -i {input.scions} \
                  --in-format newick \
                  --out-format newick \
                  -o {output.grafted_tree} &> {log}
            """)


rule insert_cog_seqs:
    input:
        grafted_tree = rules.graft.output.grafted_tree,
        metadata = rules.cog_hash_seqs.output.metadata,
    output:
        grafted_tree_expanded_cog = config["output_path"] + '/4/cog_gisaid_grafted.cog_expanded.tree',
    log:
        config["output_path"] + "/logs/4_replace_cog_seqs.log",
    resources: mem_per_cpu=20000
    shell:
        """
        sed -i.bak "s/'//g" {input.grafted_tree} 2> {log}

        /cephfs/covid/bham/climb-covid19-jacksonb/programs/jclusterfunk_v0.0.4pre/jclusterfunk insert \
            -i {input.grafted_tree} \
            --metadata {input.metadata} \
            --unique-only \
            --format newick \
            -o {output.grafted_tree_expanded_cog} 2>> {log}
        """


rule sort_collapse:
    input:
        expanded_tree = rules.insert_cog_seqs.output.grafted_tree_expanded_cog,
    output:
        sorted_tree = config["output_path"] + '/4/cog_gisaid_grafted.sorted.tree',
        sorted_collapsed_tree = config["output_path"] + '/4/cog_gisaid_full.tree.public.newick',
    params:
        collapse=0.000005,
    log:
        config["output_path"] + "/logs/4_sort_collapse.log",
    resources: mem_per_cpu=20000
    shell:
        """
        gotree rotate sort -i {input.expanded_tree} -o {output.sorted_tree} &> {log}
        gotree collapse length --length {params.collapse} -i {output.sorted_tree} -o {output.sorted_collapsed_tree} &>> {log}
        """


rule step_4_annotate_tree:
    input:
        tree = rules.sort_collapse.output.sorted_collapsed_tree,
        metadata = config["output_path"] + "/3/cog_gisaid.lineages.expanded.csv",
    output:
        tree = config["output_path"] + "/4/cog_gisaid_grafted.annotated.tree"
    log:
        config["output_path"] + "/logs/4_annotate_tree.log"
    resources: mem_per_cpu=20000
    shell:
        """
        clusterfunk annotate_tips \
          --in-metadata {input.metadata} \
          --trait-columns country lineage uk_lineage \
          --index-column sequence_name \
          --boolean-for-trait country='UK' country='UK' \
          --boolean-trait-names country_uk country_uk_deltran \
          --in-format newick \
          --out-format nexus \
          --input {input.tree} \
          --output {output.tree} &> {log}
        """


# rule acctran_ancestral_reconstruction:
#     input:
#         tree = rules.step_4_annotate_tree.output.tree
#     output:
#         tree = config["output_path"] + "/4/cog_gisaid_grafted.annotated.acc.tree"
#     log:
#         config["output_path"] + "/logs/4_acctran_ancestral_reconstruction.log"
#     shell:
#         """
#         clusterfunk ancestral_reconstruction \
#         --traits country_uk_acctran \
#         --acctran \
#         --ancestral-state False \
#         --input {input.tree} \
#         --output {output.tree} &> {log}
#         """

rule deltran_ancestral_reconstruction:
    input:
        tree = rules.step_4_annotate_tree.output.tree
    output:
        tree = config["output_path"] + "/4/cog_gisaid_grafted.annotated.del.tree"
    log:
        config["output_path"] + "/logs/4_deltran_ancestral_reconstruction.log"
    resources: mem_per_cpu=20000
    shell:
        """
        clusterfunk ancestral_reconstruction \
        --traits country_uk_deltran \
        --deltran \
        --ancestral-state False \
        --input {input.tree} \
        --output {output.tree} &> {log}
        """

# rule label_acctran_introductions:
#     input:
#         tree = rules.deltran_ancestral_reconstruction.output.tree
#     output:
#         tree = config["output_path"] + "/4/cog_gisaid_grafted.annotated.acc.del.acc_labelled.tree"
#     log:
#         config["output_path"] + "/logs/4_label_acctran_introductions.log"
#     shell:
#         """
#         clusterfunk label_transitions \
#           --trait country_uk_acctran \
#           --to True \
#           --transition-name acc_introduction \
#           --transition-prefix acc_trans_ \
#           --include_root \
#           --stubborn \
#           --input {input.tree} \
#           --output {output.tree} &> {log}
#         """

rule label_deltran_introductions:
    input:
        tree = rules.deltran_ancestral_reconstruction.output.tree
    output:
        tree = config["output_path"] + "/4/cog_gisaid_grafted.annotated.del.del_labelled.tree"
    log:
        config["output_path"] + "/logs/4_label_deltran_introductions.log"
    resources: mem_per_cpu=20000
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

# rule merge_sibling_acc_introduction:
#     input:
#         tree = rules.label_deltran_introductions.output.tree
#     output:
#         tree = config["output_path"] + "/4/cog_gisaid_grafted.annotated.acc.del.acc_labelled.del_labelled.acc_merged.tree"
#     log:
#         config["output_path"] + "/logs/4_merge_acctran_introductions.log"
#     shell:
#         """
#         clusterfunk merge_transitions \
#           --trait-to-merge acc_introduction \
#           --merged-trait-name acc_lineage \
#           --merge-siblings \
#           --input {input.tree} \
#           --output {output.tree} &> {log}
#         """

rule merge_sibling_del_introduction:
    input:
        tree = rules.label_deltran_introductions.output.tree
    params:
        outdir = config["publish_path"] + "/COG_GISAID",
    output:
        tree = config["output_path"] + "/4/cog_gisaid_grafted.annotated.del.del_labelled.del_merged.tree"
    log:
        config["output_path"] + "/logs/4_merge_deltran_introductions.log"
    resources: mem_per_cpu=20000
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
    resources: mem_per_cpu=20000
    shell:
        """
        clusterfunk extract_tip_annotations \
          --traits country uk_lineage del_introduction del_lineage \
          --input {input.tree} \
          --output {output.traits} &> {log}
        """


rule summarize_make_trees:
    input:
        traits = rules.output_annotations.output.traits,
        public_tree = rules.sort_collapse.output.sorted_collapsed_tree,
        grafted_tree = rules.graft.output.grafted_tree,
        expanded_tree = rules.insert_cog_seqs.output.grafted_tree_expanded_cog,
    params:
        grapevine_webhook = config["grapevine_webhook"],
    log:
        config["output_path"] + "/logs/4_summarize_make_trees.log"
    shell:
        """
        echo "> Total number of sequences in tree before expanding: $(gotree stats tips -i {input.grafted_tree} | tail -n+2 | wc -l)\\n" > {log}
        echo "> Number of UK sequences in tree before expanding: $(gotree stats tips -i {input.grafted_tree} | tail -n+2 | grep -E "England|Wales|Scotland|Northern_Ireland" | wc -l)\\n" >> {log}
        echo "> Total number of sequences in tree after expanding: $(gotree stats tips -i {input.expanded_tree} | tail -n+2 | wc -l)\\n" >> {log}
        echo "> Number of UK sequences in tree after expanding: $(gotree stats tips -i {input.expanded_tree} | tail -n+2 | grep -E "England|Wales|Scotland|Northern_Ireland" | wc -l)\\n" >> {log}

        echo '{{"text":"' > 4b_data.json
        echo "*Step 4: Construct and annotate lineage trees completed*\\n" >> 4b_data.json
        cat {log} >> 4b_data.json
        echo '"}}' >> 4b_data.json
        echo "webhook {params.grapevine_webhook}"
        curl -X POST -H "Content-type: application/json" -d @4b_data.json {params.grapevine_webhook}
        """
