configfile: workflow.current_basedir + "/config.yaml"

import datetime

date = datetime.date.today()

##### Configuration #####

if config.get("output_path"):
    config["output_path"] = config["output_path"].rstrip("/")
else:
    config["output_path"] = "analysis"
config["output_path"] += "/%s" %date

##### Target rules #####

rule all:
    input:
        rules.gisaid_filter_low_coverage_sequences.output,
        rules.gisaid_process_json.output.metadata,
        rules.uk_filter_low_coverage_sequences.output,
        #rules.uk_fix_headers.output.metadata

##### Modules #####
include: "rules/0_preprocess_gisaid.smk"
include: "rules/1_preprocess_uk.smk"
