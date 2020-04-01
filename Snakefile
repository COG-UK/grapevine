configfile: workflow.current_basedir + "/config.yaml"

##### Configuration #####

if config.get("output_path"):
    config["output_path"] = config["output_path"].rstrip("/")
else:
    config["output_path"] = "analysis"

##### Target rules #####

rule all:
    input:
        config["output_path"] + "/filtered_consensus_sequences.fasta"

##### Modules #####

include: "rules/1_filter_low_coverage_sequences.smk"
include: "rules/2_minimap2_to_reference.smk"