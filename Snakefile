configfile: workflow.current_basedir + "/config.yaml"

##### Configuration #####

if config.get("output_path"):
    config["output_path"] = config["output_path"].rstrip("/")
else:
    config["output_path"] = "analysis"

##### Target rules #####

rule all:
    input:
        config["output_path"] + "/filtered_fixed_gisaid.fasta"

##### Modules #####
include: "rules/1_preprocess_gisaid.smk"
