
rule treetime:
    input:
        tree = config["output_path"] + "/5/{lineage}/trees/uk_lineage_UK{i}.tree",
        metadata = config["output_path"] + "/5/cog_gisaid.lineages.with_all_traits.with_phylotype_traits.csv",
        fasta = config["output_path"] + "/3/cog_gisaid.fasta",
    params:
        lineage="{lineage}",
        i="{i}"
    output:
        temp_fasta = temp(config["output_path"] + "/6/{lineage}/trees/uk_lineage_UK{i}.temp.fasta"),
        fasta = config["output_path"] + "/6/{lineage}/trees/uk_lineage_UK{i}.fasta",
        metadata = config["output_path"] + "/6/{lineage}/trees/uk_lineage_UK{i}.timetree.csv",
        tree = config["output_path"] + "/6/{lineage}/trees/uk_lineage_UK{i}.nexus",
        directory = directory(config["output_path"] + "/6/{lineage}/trees/uk_lineage_UK{i}_timetree/")
    log:
        config["output_path"] + "/logs/6_timetree_run_lineage{lineage}_uk{i}.log"
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

        treetime \
          --aln {output.fasta} \
          --tree {output.tree} \
          --dates {output.metadata} \
          --name-column sequence_name \
          --date-column sample_date \
          --outdir {output.directory} &>> {log}
        """
