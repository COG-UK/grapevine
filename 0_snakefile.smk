configfile: workflow.current_basedir + "/config.yaml"

##### Configuration #####

if config.get("output_path"):
    config["output_path"] = config["output_path"].rstrip('/')
config["output_path"] = os.path.abspath(config["output_path"])

if config.get("publish_path"):
    config["publish_path"] = config["publish_path"].rstrip('/')
config["publish_path"] = os.path.abspath(config["publish_path"])


##### Target rules #####

rule all:
    input:
        config["output_path"] + "/logs/0_summarize_preprocess_gisaid.log",
        config["output_path"] + "/0/gisaid.matched.fasta",
        config["output_path"] + "/0/gisaid.matched.lineages.csv"

##### Modules #####
include: "rules/0_preprocess_gisaid.smk"
