// 1) cut out trees
// 2) make treetime inputs
// 3) run treetime

/*
 * Define the default parameters
 */

params.alignment = "$PWD/analysis/3/cog_gisaid.fasta"
params.tree = "$PWD/analysis/5/cog_gisaid_full.tree.nexus"
params.metadata = "$PWD/analysis/5/cog_gisaid.lineages.with_all_traits.with_phylotype_traits.csv"

/*
 * Define the output directory
 */

results_path = "$PWD/deltrans_intro_treetime_results"

/*
 * Define the processes
 */

process 'cut_out_trees' {
  cpus 40

  output:
    path 'del_introduction*' into cut_trees

  script:
  """
  clusterfunk prune \
    --extract \
    --trait del_introduction \
    --input ${params.tree} \
    --threads 40 \
    --output .
  """
}

process 'get_treetime_input' {

  input:
   each cut_tree from cut_trees

  output:
    tuple \
      val(cut_tree.baseName), \
      file("${cut_tree.baseName}.noquotes"), \
      file("${cut_tree.baseName}.fasta"), \
      file("${cut_tree.baseName}.csv") into treetime_inputs

  script:
  """
  sed "s/'//g" $cut_tree > ${cut_tree.baseName}.noquotes

  fastafunk extract \
    --in-fasta ${params.alignment} \
    --in-tree ${cut_tree.baseName}.noquotes \
    --out-fasta ${cut_tree.baseName}.fasta

  fastafunk fetch \
    --in-fasta ${cut_tree.baseName}.fasta \
    --in-metadata ${params.metadata} \
    --index-column sequence_name \
    --filter-column sequence_name sample_date \
    --out-fasta ${cut_tree.baseName}.fasta.junk \
    --out-metadata ${cut_tree.baseName}.csv \
    --restrict
  """
}

process treetime {
  publishDir "$results_path/$datasetID"
  memory '12 GB'

  input:
    tuple val(datasetID), file('sedded_tree.nexus'), file('alignment.fasta'), file('metadata.csv') from treetime_inputs

  output:
    tuple file('ancestral_sequences.fasta'), \
    file('dates.tsv'), \
    file('divergence_tree.nexus'), \
    file('molecular_clock.txt'), \
    file('root_to_tip_regression.pdf'), \
    file('sequence_evolution_model.txt'), \
    file('timetree.nexus'), \
    file('timetree.pdf') optional true

  script:
  """
  set +e

  treetime \
    --aln alignment.fasta \
    --clock-rate 0.00100 \
    --clock-filter 0 \
    --tree sedded_tree.nexus \
    --keep-root \
    --max-iter 10 \
    --dates metadata.csv \
    --name-column sequence_name \
    --date-column sample_date \
    --outdir .

  exitcode=\$?
  if [ \$exitcode -eq 0 ]
  then
      exit 0
  else
      echo "timetree failed for this sample"
      exit 0
  fi
  """
}





















//
