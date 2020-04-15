configfile: workflow.current_basedir + "/config.yaml"

import datetime

date = datetime.date.today()

##### Configuration #####

if config.get("output_path"):
    config["output_path"] = config["output_path"].rstrip("/")
config["output_path"] += "/%s/grapevine" %date 

##### Target rules #####

rule all:
    input:
        fasta = config["output_path"] + "/0/gisaid_%s.fasta" %date,
        metadata = config["output_path"] + "/0/gisaid_%s.csv" %date,
        lineages = config["output_path"] + "/0/pangolin/lineage_report.csv"


##### Modules #####
include: "rules/0_preprocess_gisaid.smk"
