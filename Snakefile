configfile: workflow.current_basedir + "/config.yaml"

##### Configuration #####

if config.get("output_path"):
    config["output_path"] = config["output_path"].rstrip('/')

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
       config["output_path"] + "/logs/5_summarize_generate_report_and_cut_out_trees.log",
        config["output_path"] + "/logs/4_summarize_make_trees.log",
        config["output_path"] + "/logs/3_summarize_combine_gisaid_and_cog.log",
        config["output_path"] + "/logs/2_summarize_pangolin_lineage_typing.log",
        config["output_path"] + "/logs/1_summarize_preprocess_uk.log",
        config["output_path"] + "/logs/0_gisaid_summarize_preprocess.log",
        config["output_path"] + "/snakejunk/all"

rule clean_up:
    input:
        config["output_path"] + "/logs/5_summarize_generate_report_and_cut_out_trees.log",
    output:
        config["output_path"] + "/snakejunk/all"
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
include: "rules/0_preprocess_gisaid.smk"
include: "rules/1_preprocess_uk.smk"
include: "rules/2_pangolin_lineage_typing.smk"
include: "rules/3_combine_gisaid_and_uk.smk"
include: "rules/4_make_trees.smk"
include: "rules/5_generate_report_and_cut_out_trees.smk"
