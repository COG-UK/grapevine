rule filter_low_coverage_sequences:
    input:
        consensus = config["input_fasta"]
        threshold = config["low_coverage_threshold"]
    output:
        config["output_path"] + "/filtered_consensus_sequences.fasta"
    shell:
        """
        datafunk filter_low_coverage -i {input.consensus} -o {output} -t {input.threshold}
        """