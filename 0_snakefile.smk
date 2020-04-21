configfile: workflow.current_basedir + "/config.yaml"

import datetime

date = datetime.date.today()

##### Configuration #####

if config.get("output_path"):
    config["output_path"] = config["output_path"].rstrip("/") + "/analysis"
else:
    config["output_path"] = "analysis"

if config.get("publish_path"):
    config["publish_path"] = config["publish_path"].rstrip("/") + "/publish"
else:
    config["publish_path"] = "publish"

##### Target rules #####

rule all:
    input:
        fasta = config["output_path"] + "/0/gisaid.regularized.fasta",
        metadata = config["output_path"] + "/0/gisaid.regularized.csv",
        summary = config["output_path"] + "/logs/0_summary_preprocess_gisaid.log",
        counts = config["output_path"] + "/0/gisaid_counts_by_country.csv",


##### Modules #####
include: "rules/0_preprocess_gisaid.smk"
