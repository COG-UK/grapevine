configfile: workflow.current_basedir + "/config.yaml"

import datetime

date = datetime.date.today()

##### Configuration #####

if config.get("output_path"):
    config["output_path"] = config["output_path"].rstrip("/") + "/analysis"
else:
    config["output_path"] = "analysis"

##### Target rules #####

rule all:
    input:
        config["output_path"] + "/logs/1_summary_preprocess_uk.log"


##### Modules #####
include: "rules/1_preprocess_uk.smk"
