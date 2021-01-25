rule cog_hash_seqs:
    input:
        fasta = rules.uk_output_lineage_table.output.fasta,
        lineage_splits = config["lineage_splits"]
    output:
        fasta = config["output_path"] + "/3/uk.hashed.fasta",
        metadata = config["output_path"] + "/3/uk.hashmap.csv",
    log:
        config["output_path"] + "/logs/3_cog_hash_seqs.log",
    resources: mem_per_cpu=20000
    run:
        from Bio import SeqIO

        outgroups = []
        with open(input.lineage_splits, "r") as outgroup_handle:
            line = outgroup_handle.readline()
            while line:
                try:
                    outgroup = line.strip().split(",")[-1]
                    outgroups.append(outgroup)
                except:
                    continue
                line = outgroup_handle.readline()

        input_fasta = SeqIO.index(str(input.fasta), "fasta")

        hash_dict = {}

        with open(str(input.fasta), "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                id = record.id
                seq = str(record.seq)

                if seq in hash_dict:
                    hash_dict[seq] = hash_dict[seq] + [id]
                else:
                    hash_dict[seq] = [id]

        with open(str(output.fasta), "w") as fasta, open(str(output.metadata), "w") as metadata:
            metadata.write("tip,redundant\n")

            for key, value in hash_dict.items():
                if len(value) == 1:

                    r = input_fasta[value[0]]

                    fasta.write(">" + r.id + "\n")
                    fasta.write(str(r.seq) + "\n")

                elif len(value) > 1:
                    r = None

                    for id in value:
                        if id in outgroups:
                            r = input_fasta[id]
                            fasta.write(">" + r.id + "\n")
                            fasta.write(str(r.seq) + "\n")
                            value.remove(id)

                    if not r:
                        r = input_fasta[value[0]]
                        fasta.write(">" + r.id + "\n")
                        fasta.write(str(r.seq) + "\n")
                        value.remove(value[0])

                    metadata.write(r.id + ",")
                    metadata.write("|".join(value) + "\n")


rule uk_output_hashed_lineage_table:
    input:
        fasta = rules.cog_hash_seqs.output.fasta,
        metadata = rules.uk_output_lineage_table.output.metadata,
    output:
        fasta = temp(config["output_path"] + "/3/uk.hashed.temp.fasta"),
        metadata = config["output_path"] + "/3/uk.hashed.lineages.csv"
    log:
        config["output_path"] + "/logs/3_uk_output_hashed_lineage_table.log"
    resources: mem_per_cpu=20000
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name country adm1 adm2 \
                          sample_date epi_week \
                          lineage uk_lineage \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --log-file {log} \
          --restrict
        """


rule gisaid_output_lineage_table:
    input:
        fasta =  config["GISAID_background_fasta"],
        metadata = config["GISAID_background_metadata"],
    output:
        fasta = config["output_path"] + "/3/gisaid.matched.fasta",
        metadata = config["output_path"] + "/3/gisaid.matched.lineages.csv",
    log:
        config["output_path"] + "/logs/3_gisaid_output_lineage_table.log"
    resources: mem_per_cpu=20000
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column sequence_name \
          --filter-column sequence_name country adm1 adm2 \
                          sample_date epi_week \
                          lineage uk_lineage \
          --where-column adm1=edin_admin_1 adm2=edin_admin_2 \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --log-file {log} \
          --restrict
        """


rule combine_gisaid_and_cog:
    input:
        previous_stage = config["output_path"] + "/logs/2_summarize_pangolin_lineage_typing.log",
        gisaid_fasta = rules.gisaid_output_lineage_table.output.fasta,
        gisaid_metadata = rules.gisaid_output_lineage_table.output.metadata,
        uk_fasta = rules.cog_hash_seqs.output.fasta,
        uk_metadata = rules.uk_output_hashed_lineage_table.output.metadata
    output:
        fasta = config["output_path"] + "/3/cog_gisaid.fasta",
        metadata = config["output_path"] + "/3/cog_gisaid.lineages.csv"
    log:
        config["output_path"] + "/logs/3_combine_gisaid_and_cog.log"
    resources: mem_per_cpu=20000
    shell:
        """
        fastafunk merge \
          --in-fasta {input.gisaid_fasta} {input.uk_fasta} \
          --in-metadata {input.gisaid_metadata} {input.uk_metadata} \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --index-column sequence_name \
          --log-file {log}
        """


rule combine_gisaid_and_cog_expanded:
    input:
        previous_stage = config["output_path"] + "/logs/2_summarize_pangolin_lineage_typing.log",
        gisaid_fasta = rules.gisaid_output_lineage_table.output.fasta,
        gisaid_metadata = rules.gisaid_output_lineage_table.output.metadata,
        uk_fasta = rules.uk_output_lineage_table.output.fasta,
        uk_metadata = rules.uk_output_lineage_table.output.metadata
    output:
        fasta = config["output_path"] + "/3/cog_gisaid.expanded.fasta",
        metadata = config["output_path"] + "/3/cog_gisaid.lineages.expanded.csv"
    log:
        config["output_path"] + "/logs/3_combine_gisaid_and_cog.log"
    resources: mem_per_cpu=20000
    shell:
        """
        fastafunk merge \
          --in-fasta {input.gisaid_fasta} {input.uk_fasta} \
          --in-metadata {input.gisaid_metadata} {input.uk_metadata} \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --index-column sequence_name \
          --log-file {log}
        """


rule summarize_combine_gisaid_and_cog:
    input:
        cog_hashed_fasta = rules.cog_hash_seqs.output.fasta,
        fasta = rules.combine_gisaid_and_cog.output.fasta,
        metadata = rules.combine_gisaid_and_cog.output.metadata,
        full_fasta = rules.combine_gisaid_and_cog_expanded.output.fasta,
        full_metadata = rules.combine_gisaid_and_cog_expanded.output.metadata,
    params:
        grapevine_webhook = config["grapevine_webhook"],
        json_path = config["json_path"],
    log:
        config["output_path"] + "/logs/3_summarize_combine_gisaid_and_cog.log"
    shell:
        """
        echo "> Number of sequences in total COG and GISAID matched files for later steps: $(cat {input.full_fasta} | grep ">" | wc -l)\\n" &>> {log}
        echo "> Number of sequences in reduced COG fasta: $(cat {input.cog_hashed_fasta} | grep ">" | wc -l)\\n" &>> {log}
        echo "> Number of sequences in collapsed COG and GISAID matched files for tree building: $(cat {input.fasta} | grep ">" | wc -l)\\n" &>> {log}
        echo "> \\n" &>> {log}
        echo '{{"text":"' > {params.json_path}/3_data.json
        echo "*Step 3: Combine COG-UK and GISAID data complete*\\n" >> {params.json_path}/3_data.json
        cat {log} >> {params.json_path}/3_data.json
        echo '"}}' >> {params.json_path}/3_data.json
        echo "webhook {params.grapevine_webhook}"
        curl -X POST -H "Content-type: application/json" -d @{params.json_path}/3_data.json {params.grapevine_webhook}
        """
