configfile: workflow.current_basedir + "/config.yaml"

##### Configuration #####

if config.get("output_path"):
    config["output_path"] = config["output_path"].rstrip('/')

if config.get("publish_path"):
    config["publish_path"] = config["publish_path"].rstrip('/')
config["publish_path"] = os.path.abspath(config["publish_path"])


##### Target rules #####

rule all:
    input:
        config["output_path"] + "/logs/0_summarize_preprocess_gisaid.log",
        config["output_path"] + "/0/gisaid_counts_by_country.csv",
        config["output_path"] + "/0/QC_distances.tsv",
        config["output_path"] + "/0/QC_distances.png",
        config["output_path"] + "/snakejunk/0"

rule clean_up:
    input:
        config["output_path"] + "/logs/0_summarize_preprocess_gisaid.log",
    output:
        config["output_path"] + "/snakejunk/0"
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
