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
        config["output_path"] + "/3/cog_gisaid_%s.fasta" %date,
        config["output_path"] + "/3/cog_gisaid_%s.csv" %date,
        config["output_path"] + "/logs/4_run_subroutine_on_lineage.log"

##### Modules #####
include: "rules/0_preprocess_gisaid.smk"
include: "rules/1_preprocess_uk.smk"
include: "rules/2_pangolin_lineage_typing.smk"
include: "rules/3_combine_gisaid_and_uk.smk"
include: "rules/4_make_trees.smk"

