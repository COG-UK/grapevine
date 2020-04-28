configfile: workflow.current_basedir + "/config.yaml"

import os

##### Configuration #####

if config.get("output_path"):
    config["output_path"] = config["output_path"].rstrip("/") + "/analysis"
else:
    config["output_path"] = "analysis"

if config.get("publish_path"):
    config["publish_path"] = config["publish_path"].rstrip("/") + "/publish"
else:
    config["publish_path"] = "publish"
config["publish_path"] = os.path.abspath(config["publish_path"])

##### Target rules #####

rule all:
    input:
        config["output_path"] + "/logs/4_summarize_make_trees.log",
        config["output_path"] + "/logs/3_summarize_combine_gisaid_and_cog.log",
        config["output_path"] + "/logs/2_summarize_pangolin_lineage_typing.log",
        config["output_path"] + "/logs/1_summarize_preprocess_uk.log",
        config["output_path"] + "/snakejunk"

rule clean_up:
    input:
        config["output_path"] + "/logs/4_summarize_make_trees.log",
    output:
        config["output_path"] + "/snakejunk"
    shell:
        """
        mkdir -p {output}
        mv slurm-*.out *_data.json {output}/
        for file in pre trace default.profraw
        do
          if [ -f "$file" ]
          then
            rm $file
          fi
        done
        """

##### Modules #####
include: "rules/1_preprocess_uk.smk"
include: "rules/2_pangolin_lineage_typing.smk"
include: "rules/3_combine_gisaid_and_uk.smk"
include: "rules/4_make_trees.smk"