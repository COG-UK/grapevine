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


##### Target rules #####

rule all:
    input:
        config["output_path"] + "/logs/1a_summarize_preprocess_uk.log",
        config["output_path"] + "/logs/1b_summarize_publish.log"


##### Modules #####
include: "rules/1a_preprocess_uk.smk"
include: "rules/1b_publish.smk"
