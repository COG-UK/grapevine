// 1) get samples from metadata
// 2) cut out trees
// 3) extract fasta & metadata

/*
 * Define the default parameters
 */

bigalignment = "/cephfs/covid/bham/results/phylogenetics/20201026/alignments/cog_2020-10-26_alignment.fasta"
bigtree = "/cephfs/covid/bham/results/phylogenetics/20201026/trees/cog_global_2020-10-26_tree.newick"
bigmetadata = "/cephfs/covid/bham/results/phylogenetics/20201026/trees/cog_global_2020-10-26_metadata.csv"

/*
 * Define the output directory
 */

results_path = "$PWD/uk_lineages"

/*
 * Define the processes
 */

process get_samples {
  cpus 1

  output:
    file "UK*.samples" into sample_path_channel

  script:
  """
  grep -oh -E ",UK[0-9]+," $bigmetadata | cut -d, -f2 | sort | uniq | while read LINE
  do
    I=`echo \${LINE} | cut -d"K" -f2`
    if [ \$(grep -c ",\${LINE}," $bigmetadata) -ge 3 ]
    then
      grep ",\${LINE}," $bigmetadata | cut -d"," -f1 > UK\${I}.samples.temp
      sed -i.bak "s/^/'/g" UK\${I}.samples.temp
      sed "s/\$/'/g" UK\${I}.samples.temp > UK\${I}.samples
    fi
  done
  """
}

process cut_trees {
  publishDir "$results_path", mode: "copy"
  cpus 1

  input:
   each sample_file from sample_path_channel

  output:
    file "${sample_file.baseName}.newick" into tree_path_channel

  script:
  """
  gotree prune -r \
      -i $bigtree \
      --tipfile ${sample_file} \
      -o ${sample_file.baseName}.newick
  """
}

process extract {
  publishDir "$results_path", mode: "copy"

  input:
    each tree from tree_path_channel

  output:
  tuple \
    file("${tree.baseName}.fasta"), \
    file("${tree.baseName}.csv")

  script:
  """
  fastafunk extract \
    --in-fasta $bigalignment \
    --in-tree ${tree} \
    --out-fasta ${tree.baseName}.fasta

    fastafunk fetch \
    --in-fasta ${tree.baseName}.fasta \
    --in-metadata $bigmetadata \
    --index-column sequence_name \
    --filter-column sequence_name secondary_identifier sample_date epi_week \
                    country adm1 \
                    is_surveillance is_community is_hcw \
                    is_travel_history travel_history lineage \
                    lineage_support uk_lineage acc_lineage del_lineage acc_introduction del_introduction phylotype \
    --out-fasta ${tree.baseName}.fasta.junk \
    --out-metadata ${tree.baseName}.csv \
    --restrict
  """
}
