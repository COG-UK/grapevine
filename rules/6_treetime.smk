rule treetime:
    input:
        tree=config["output_path"] + "/5/trees/uk_lineage_UK{i}.tree",
        metadata = config["output_path"] + "/5/cog_gisaid.lineages.with_all_traits.with_phylotype_traits.csv",
        fasta = config["output_path"] + "/3/cog_gisaid.fasta",
    params:
        i="{i}"
    output:
        temp_fasta = temp(config["output_path"] + "/6/trees/uk_lineage_UK{i}.temp.fasta"),
        fasta = config["output_path"] + "/6/trees/uk_lineage_UK{i}.fasta",
        metadata = config["output_path"] + "/6/trees/uk_lineage_UK{i}.timetree.csv",
        tree = config["output_path"] + "/6/trees/uk_lineage_UK{i}.nexus",
        directory = directory(config["output_path"] + "/6/trees/uk_lineage_UK{i}_timetree/")
    log:
        config["output_path"] + "/logs/6_timetree_run_uk{i}.log"
    resources: mem_per_cpu=48000
    shell:
        """
        sed "s/'//g" {input.tree} > {output.tree}

        fastafunk extract \
          --in-fasta {input.fasta} \
          --in-tree {input.tree} \
          --out-fasta {output.fasta} &>> {log}

        fastafunk fetch \
          --in-fasta {output.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name sample_date \
          --out-fasta {output.temp_fasta} \
          --out-metadata {output.metadata} \
          --restrict &>> {log}


        set +e

        treetime \
          --aln {output.fasta} \
          --clock-rate 0.00100 \
          --clock-filter 0 \
          --keep-root \
          --max-iter 10 \
          --tree {output.tree} \
          --dates {output.metadata} \
          --name-column sequence_name \
          --date-column sample_date \
          --outdir {output.directory} &>> {log}

        exitcode=$?
        if [ $exitcode -eq 0 ]
        then
            exit 0
        else
            echo "timetree failed for this sample" &>> {log}
            exit 0
        fi
        """


def aggregate_input_treetime_logs(wildcards):
    checkpoint_output_directory = checkpoints.get_uk_lineage_samples.get(**wildcards).output[0]
    print(checkpoints.get_uk_lineage_samples.get(**wildcards).output[0])
    required_files = expand( "%s/logs/6_timetree_run_uk{i}.log" %(config["output_path"]),
                            i=glob_wildcards(os.path.join(checkpoint_output_directory, "UK{i}.samples.txt")).i)
    return (sorted(required_files))

rule summarize_treetime:
    input:
        logs=aggregate_input_treetime_logs
    params:
        grapevine_webhook = config["grapevine_webhook"],
        json_path = config["json_path"],
        date = config["date"]
    log:
        config["output_path"] + "/logs/6_summarize_treetime.log"
    shell:
        """
        echo '{{"text":"' > {params.json_path}/6_data.json
        echo "*Step 6: treetime for {params.date} UK lineages complete*\\n" >> {params.json_path}/6_data.json
        echo '"}}' >> {params.json_path}/6_data.json
        echo 'webhook {params.grapevine_webhook}'
        curl -X POST -H "Content-type: application/json" -d @{params.json_path}/6_data.json {params.grapevine_webhook}
        """










#
