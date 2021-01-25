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
        config["output_path"] + "/logs/1_summarize_preprocess_uk.log",


##### Modules #####
include: "/cephfs/covid/bham/climb-covid19-jacksonb/git/grapevine/rules/1_preprocess_uk.smk"
