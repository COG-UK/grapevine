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
        fasta_cog = config["publish_path"] + "/GISAID/gisaid.matched.fasta",
        metadata_cog = config["publish_path"] + "/GISAID/gisaid.matched.csv",
        fasta_gisaid = config["publish_path"] + "/GISAID/gisaid.full.fasta",
        metadata_gisaid = config["publish_path"] + "/GISAID/gisaid.full.csv",
        fasta_gisaid_new = config["publish_path"] + "/GISAID/gisaid.new.fasta",
        counts = config["output_path"] + "/0/gisaid_counts_by_country.csv",
        QC_table = config["output_path"] + "/0/QC_distances.tsv",
        QC_plot = config["output_path"] + "/0/QC_distances.png",
        summary = config["output_path"] + "/logs/0_summary_preprocess_gisaid.log"

##### Modules #####
include: "rules/0_preprocess_gisaid.smk"
