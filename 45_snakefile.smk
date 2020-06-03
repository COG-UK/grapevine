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

if not config.get("date"):
    cwd = os.getcwd()
    config["date"] = os.path.basename(cwd)[:10]

##### Target rules #####

rule all:
    input:
        config["output_path"] + "/logs/5_define_uk_lineages_and_cut_out_trees.log",
        config["output_path"] + "/logs/4_summarize_make_trees.log",
        # config["output_path"] + "/snakejunk/45"
#
# rule clean_up:
#     input:
#         config["output_path"] + "/logs/5_summarize_generate_report_and_cut_out_trees.log",
#     output:
#         config["output_path"] + "/snakejunk/45"
#     shell:
#         """
#         mkdir -p {output}
#         mv slurm-*.out *_data.json {output}/
#         for file in pre trace default.profraw
#         do
#           if [ -f "$file" ]
#           then
#             rm $file
#           fi
#         done
#         """

##### Modules #####
include: "rules/4_make_trees.smk"
include: "rules/5_define_uk_lineages_and_cut_out_trees.smk"
