# Grapevine

This document provides an overview of the complete Grapevine pipeline, beginning with unaligned sequence data (from GISAID and CoG-UK) and ending up in annotated phylogenetic trees, lineage designations and other published output.

## Gisaid processing

### Relevant snakemake files

`rules/0_preprocess_gisaid.smk`

### Input

`export.json`
> *A dump of GISAID data in .json format*

### Steps

1. parse GISAID dump (`export.json`):
> Extract:
>
>	a) Fasta file of sequences
>
>	b) Associated metadata
>
> Filter sequences:
>
> a) That are present in a predefined file of known problematic sequences
>
> b) `covv_host.lower() != 'human'`
>
> c) With a mal-formed (not `YYYY-MM-DD`) or impossible (earlier than `2019-11-30` or later than today) date in `covv_collection_date`
>

2. Deduplicate samples, set a Grapevine-style header (format: `location/sequence_dentifier/year`) and filter sequences with fewer than 29000 sites

3. Map to `Wuhan-Hu-1` using Minimap:
	`minimap2 -a -x asm5 reference.fasta gisaid.fasta -o output.sam`

4. Remove insertions relative to the reference and pad sites outside of the CDS with `N`
5. Filter sequences with more than 5% `N`s (excluding the padding in step 4)
6. Mask predefined problematic sites

	These are defined by column and regex match to sequence ID:
	```
	11083,?,\w
	```
	So (1-based) column `11083` is masked in all sequences
	See this [Virological post](https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473) for more information

7. Calculate distance to `WH04` - the sequence closest to the 'root' of the outbreak
8. Filter on excessive distance to `WH04` (anything more than 4 standard deviations from each epi-week's mean distance is removed)
9. Genotype each sequence for predefined SNPs of interest
10. Assign a pangolin lineage to sequences that lack one
11. Output complete dataset with associated metadata, and non-UK dataset and associated metadata, which is carried forward to be combined with the CoG-UK data

### Outputs

`config["output_path"] + /0/gisaid.matched.fasta`
> *Non-UK sequences passing Grapevine QC*
>

`config["output_path"] + /0/gisaid.matched.lineages.csv`  
> *Metadata with Grapevine-specific information matching the fasta file the above*
>
> Columns: `sequence_name,country,adm1,adm2,sample_date,epi_week,lineage,uk_lineage`
>

`config["publish_path"] + "GISAID/gisaid.trimmed_alignment.fasta"`
> *All sequences passing QC, mapped with regions outside of the CDS padded with `N`*
>

`config["publish_path"] + "GISAID/gisaid.metadata.csv`
> *Metadata matching the fasta file above*
> Columns: `sequence_name,country,edin_admin_1,edin_admin_2,edin_travel,edin_date_stamp,sample_date,epi_week,special_lineage,lineage,lineages_version,lineage_support,d614g,covv_accession_id,covv_virus_name,covv_location,covv_add_host_info,covv_assembly_method,covv_gender,covv_host,covv_passage,covv_patient_age,covv_seq_technology,covv_specimen,covv_subm_date`

`config["publish_path"] + "GISAID/gisaid_counts_by_country.csv`
> *Sequence counts by country*
>


## CoG-UK processing

### Relevant snakemake files

`1_preprocess_uk.smk`

`2_pangolin_lineage_typing.smk`

`3_combine_gisaid_and_uk.smk`

### Input

`/cephfs/covid/bham/artifacts/published/latest/elan.<DATE>.consensus.fasta`
> *All unaligned CoG-UK sequences*
>

`/cephfs/covid/bham/artifacts/published/latest/majora.<DATE>.metadata.tsv`
> *CoG-UK metadata associated with the fasta file above*
>

### Steps

1. Set a Grapevine-style header (format: `location/sequence_dentifier/year`)

2. Add annotations to metadata:
> epi-week from `collection_date`. If not `collection_date`, then epi-week from `received_date`
>
> Proportion of missing data for each sequence
>
> `covid` from sequence name

3. Deduplicate samples:
> By `covid` - keep the sample with the lowest proportion of missing data
>
> By `biosample_source_id` - keep the sample with the earliest collection date

4. Map to `Wuhan-Hu-1` using Minimap:
	`minimap2 -a -x asm5 reference.fasta gisaid.fasta -o output.sam`

5. Remove insertions relative to the reference and pad sites outside of the CDS with `N`
6. Filter sequences with more than 5% `N`s (excluding the padding in step 5)
7. Mask predefined problematic sites

	These are defined by column and regex match to sequence ID:
	```
	11083,?,\w
	```
	So (1-based) column `11083` is masked in all sequences
	See this [Virological post](https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473) for more information
8. Genotype each sequence for predefined SNPs of interest
9. Add previous Pangolin and UK lineage designations to metadata
10. Assign a pangolin lineage to sequences that lack one
11. Combine CoG-UK dataset with non-UK GISAID data

### Output

` config["output_path"] +/2/uk.matched.fasta`
> *CoG-UK sequences passing Grapevine QC*
>

`config["output_path"] + /2/uk.matched.lineages.csv`  
> *Metadata with Grapevine-specific information matching the fasta file the above*
>
> Columns: `sequence_name,country,adm1,adm2,sample_date,epi_week,lineage,uk_lineage`

`config["output_path"] + /3/cog_gisaid.fasta`
> *GISAID + CoG-UK sequences passing Grapevine QC*
>

`config["output_path"] + /3/cog_gisaid.lineages.csv`  
> *Metadata with Grapevine-specific information matching the fasta file the above*
>
> Columns: `sequence_name,country,adm1,adm2,sample_date,epi_week,lineage,uk_lineage`

## Tree building and ancestral reconstruction

### Relevant snakemake files

`4_make_trees.smk`

`4_subroutine/4_process_lineages.smk`

### Input

`config['output_path'] + /3/cog_gisaid.fasta`
> *GISAID + CoG-UK sequences passing Grapevine QC*
>

`config['output_path'] + /3/cog_gisaid.lineages.csv`
> *CoG-UK metadata associated with the fasta file above*
>

`config["lineage_splits"]`
> *Definitions of (Pangolin) lineage splits with an associated outgroup for each sub-lineage. Minimally, this can be Lineage A and "Wuhan/WH04/2020". Eg:*
>
> ```
> lineage,outgroup
> A,Wuhan/WH04/2020
> ```
> *Defining more splits will result in more (smaller) subtrees being built. Eg:*
>
> ```
> lineage,outgroup
> A,Wuhan/WH04/2020
> B,Wuhan/WHU01/2020
> B.1,Italy/TE5166/2020
> B.1.1,Germany/BAV-V2010837/2020
> B.2,Norway/1379/2020
> ```
> *Subtrees are grafted together into one complete tree after they are built*

### Steps

1. Split the complete dataset into the subsets defined in
`config["lineage_splits"]`

	Then, for each subset:

	1. Run FastTree:
		`FastTreeMP -nosupport -nt <subset.fasta> > <subset.fasta>`
	2. Root the tree with its designationed outgroup

2. Graft together the subtrees into one complete tree
3. Annotate the tree with country and previous UK lineage designation
4. **Designate lineages** - See ([_CoG-UK Lineage Designation_](https://github.com/COG-UK/grapevine/blob/master/docs/lineages.md)) for more detail
5. Output a full, unannotated public tree and lineage designation information to pass to the next step, which is (re-)defining UK lineages

### Output

`config["output_path"] + "/4/all_traits.csv"`
> *step-specific metadata containing information about (previous) UK lineage designations, with new acctrans and deltrans lineage designations*
>

`["output_path"] + "/4/cog_gisaid_full.tree.public.newick"`
> *phylogenetic tree containing all GISAID and CoG-UK sequences*
>

##  UK lineage designation and operations

### Relevant snakemake files

`5_define_uk_lineages_and_cut_out_trees.smk`

`6_treetime.smk`

### Input

`config["output_path"] + "/4/all_traits.csv`
> *step-specific metadata containing information about (previous) UK lineage designations, with new deltrans lineage designations*
>


### Steps

1. Merge and create new UK lineage definitions and apply them to all CoG-UK samples:

	>
	> See ([_CoG-UK Lineage Designation_](https://github.com/COG-UK/grapevine/blob/master/docs/lineages.md)) for more detail
	>

2. Annotate the full tree with UK lineage and deltrans information

3. Cut out UK lineage-specific trees from the full tree
4. Phylotyping

	>
	> See ([_CoG-UK Lineage Designation_](https://github.com/COG-UK/grapevine/blob/master/docs/lineages.md)) for more detail
	>
	
5. Run `TreeTime` on UK lineages:

	```
	treetime \
  	  --aln UK_lineage.fasta \
  	  --clock-rate 0.00100 \
  	  --clock-filter 0 \
  	  --tree UK_lineage.tree \
  	  --dates UK_lineage.csv \
  	  --name-column sequence_name \
  	  --date-column sample_date \
  	  --outdir UK_lineage/    
	```

### Output

`config["output_path"] + "/5/cog_gisaid_full.tree.nexus"`
> *phylogenetic tree annotated with location and lineage information*

`config["output_path"] + "/5/updated_traits.csv"`
> *metadata containing information about new UK lineage designations and current deltrans lineage designations*

`config["output_path"] + "/5/trees/"`
and
`config["output_path"] + "/5/phylotyped_trees/"`
> *Directories of trees which correspond to UK lineages, as defined in step 1. above, with and without phylotype annotation, respectively*

`config["output_path"] + "/5/UK_phylotypes.csv"`
> *metadata containing information about new UK phylotype designations*

`config["output_path"] + "/6/trees/"`
> *Directory containing TimeTree analyses in subdirectories*


## Publishing

Consortium-wide outputs are published to `config["export_path"]` (and `rsync`ed to `/cephfs/covid/bham/artifacts/published/[DATE]/phylogenetics/` on `climb`).

See the [Grapevine README](https://github.com/COG-UK/grapevine/blob/master/README.md) for details of these outputs.
