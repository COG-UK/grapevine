import os

configfile: workflow.current_basedir + "/config.yaml"

##### Configuration #####

if config.get("output_path"):
    config["output_path"] = config["output_path"].rstrip('/')
config["output_path"] = os.path.abspath(config["output_path"])

if config.get("publish_path"):
    config["publish_path"] = config["publish_path"].rstrip('/')
config["publish_path"] = os.path.abspath(config["publish_path"])

if not config.get("cwd"):
    config["cwd"] = os.getcwd()

if not config.get("date"):
    config["date"] = os.path.basename(config["cwd"])[:10]

##### Target rules #####

rule all:
    input:
        config["output_path"] + "/logs/0_alert_sam.log",

##### Modules #####
include: "rules/0_preprocess_gisaid.smk"
