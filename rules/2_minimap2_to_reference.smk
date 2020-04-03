rule minimap2_to_reference:
    input:
        consensus = rules.filter_low_coverage_sequences.output
        reference = config["reference_fasta"]
    output:
        config["output_path"] + "/filtered_consensus_sequences.mapped.sam"
    shell:
        """
        minimap2 -a -x asm5 {input.reference} {input.consensus} > {output}
        """