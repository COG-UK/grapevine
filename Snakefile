configfile: workflow.current_basedir + "/config.yaml"

import datetime

date = datetime.date.today()

##### Configuration #####

if config.get("output_path"):
    config["output_path"] = config["output_path"].rstrip("/")
else:
    config["output_path"] = "analysis"

##### Target rules #####

rule all:
    input:
        config["output_path"] + "/gisaid_%s_filtered.fasta" %date,
        config["output_path"] + "/uk_%s_filtered.fasta" %date

##### Modules #####
include: "rules/1_preprocess_gisaid.smk"
include: "rules/2_preprocess_uk.smk"
