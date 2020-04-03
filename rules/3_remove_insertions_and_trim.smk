rule remove_insertions_and_trim:
    input:
        sam = rules.minimap2_to_reference.output
        reference = config["reference_fasta"]
    params:
        trim_start = config["trim_start"]
        trim_end = config["trim_end"]
    output:
        config["output_path"] + "/filtered_fixed_consensus_sequences.fasta"
    shell:
        """
        datafunk sam_2_fasta \
          -s {input.sam} \
          -r {input.reference} \
          -o {output} \
          -t [params.trim_start, params.trim_end] \
          --prefix_ref
        """