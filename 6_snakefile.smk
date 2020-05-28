import os
import pandas as pd

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

LINEAGES = []
df = pd.read_csv(config["lineage_splits"])
for i,row in df.iterrows():
    LINEAGES.append(row['lineage'])

UK = []
for l in LINEAGES:
    UK.append(glob_wildcards(config["output_path"] + "/5/" + l + "/trees/uk_lineage_UK{i}.tree").i)

LINEAGES_REP = []
for i,x in enumerate(UK):
    LINEAGES_REP = LINEAGES_REP + [LINEAGES[i]] * len(x)

UK = [item for sublist in UK for item in sublist]

rule all:
    input:
        expand(config["output_path"] + "/logs/6_timetree_run_lineage{lineage}_uk{i}.log", zip, lineage = LINEAGES_REP, i = UK)


##### Modules #####
include: "rules/6_treetime.smk"
