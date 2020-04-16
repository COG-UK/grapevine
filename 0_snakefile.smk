configfile: workflow.current_basedir + "/config.yaml"

import datetime

date = datetime.date.today()

##### Configuration #####

if config.get("output_path"):
    config["output_path"] = config["output_path"].rstrip("/") + "/"
config["output_path"] += "analysis"

##### Target rules #####

rule all:
    input:
        fasta = config["output_path"] + "/0/gisaid_%s.fasta" %date,
        metadata = config["output_path"] + "/0/gisaid_%s.csv" %date,
        counts = config["output_path"] + "/0/gisaid_counts_by_country.csv",
        lineages = config["output_path"] + "/0/pangolin/lineage_report.csv"


##### Modules #####
include: "rules/0_preprocess_gisaid.smk"
