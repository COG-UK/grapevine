# configfile: workflow.current_basedir + "/config.yaml"

import os

##### Configuration #####

if config.get("output_path"):
    config["output_path"] = config["output_path"].rstrip('/')
config["output_path"] = os.path.abspath(config["output_path"])

if config.get("publish_path"):
    config["publish_path"] = config["publish_path"].rstrip('/')
config["publish_path"] = os.path.abspath(config["publish_path"])

if config.get("export_path"):
    config["export_path"] = config["export_path"].rstrip('/')
config["export_path"] = os.path.abspath(config["export_path"])

if not config.get("cwd"):
    config["cwd"] = os.getcwd()

if not config.get("date"):
    config["date"] = os.path.basename(config["cwd"])[:10]

##### Target rules #####

rule all:
    input:
        config["output_path"] + "/logs/7_summarize_postpublish.log",


##### Modules #####
include: "rules/1_preprocess_uk.smk"
include: "rules/2_pangolin_lineage_typing.smk"
include: "rules/3_combine_gisaid_and_uk.smk"
include: "rules/4_make_trees.smk"
include: "rules/5_define_uk_lineages_and_cut_out_trees.smk"
include: "rules/6_treetime.smk"
include: "rules/7_publish.smk"
