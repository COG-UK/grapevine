configfile: workflow.current_basedir + "/config.yaml"

import os
import pandas as pd

##### Configuration #####

config["publish_path"] = os.path.abspath(config["publish_path"])

LINEAGES = config["lineages"].split()[1:]
OUTGROUPS = config["lineage_specific_outgroups"].split()[1:]

lineage_to_outgroup_map = {}
for i,lin in enumerate(LINEAGES):
    lineage_to_outgroup_map[lin] = OUTGROUPS[i]

print("lineages", LINEAGES)
print("outgroups", OUTGROUPS)

##### Target rules #####

rule all:
    input:
        grafted_tree=config["output_path"] + '/4/cog_gisaid_grafted.tree'


rule fasttree:
    input:
        lineage_fasta = config["output_path"] + "/4/lineage_{lineage}.fasta"
    params:
        lineage = "{lineage}",
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.unrooted.tree",
    log:
        config["output_path"] + "/logs/4_fasttree_{lineage}.log"
    resources: mem_per_cpu=20000
    threads: 3
    shell:
        """
        echo "{input.lineage_fasta} {params.lineage}"

        export OMP_NUM_THREADS={threads}

        FastTreeMP -nosupport -nt {input.lineage_fasta} > {output.tree} 2> {log}
        """


rule root_tree:
    input:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.unrooted.tree",
    params:
        lineage = "{lineage}",
        outgroup = lambda wildcards: lineage_to_outgroup_map[wildcards.lineage]
    output:
        tree = config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.tree",
    log:
        config["output_path"] + "/logs/4_root_tree_{lineage}.log",
    shell:
        """
        echo "{params.lineage} {params.outgroup}"

        clusterfunk root \
          --in-format newick \
          -i {input.tree} \
          --out-format newick \
          -o {output.tree} \
          --outgroup {params.outgroup} &>> {log}
        """


rule graft:
    input:
        # not sure how to pass this as a space separated list below. Also assuming the order here matches lineages
        scions = expand(config["output_path"] + "/4/{lineage}/cog_gisaid_{lineage}.tree", lineage=sorted(LINEAGES)),
        guide_tree = config["guide_tree"]
    params:
        lineages = sorted(LINEAGES),
    output:
        grafted_tree = config["output_path"] + '/4/cog_gisaid_grafted.tree',
    log:
        config["output_path"] + "/logs/4_graft.log",
    run:
        if len(input.scions) > 1:
            shell("""
                clusterfunk graft \
                  --scions {input.scions} \
                  --scion-annotation-name scion_lineage \
                  --annotate-scions {params.lineages} \
                  --input {input.guide_tree} \
                  --in-format newick \
                  --out-format newick \
                  --output {output.grafted_tree} &> {log}
                """)

        else:
            shell("""
                clusterfunk reformat \
                  -i {input.scions} \
                  --in-format newick \
                  --out-format newick \
                  -o {output.grafted_tree} &> {log}
            """)
