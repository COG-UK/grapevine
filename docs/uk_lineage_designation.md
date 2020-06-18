## UK lineage designation and operations

### Rules

`5_define_uk_lineages_and_cut_out_trees.smk`

`5_subroutine/5_process_lineages.smk`

`6_treetime.smk`

### Input

`config["output_path"] + "/4/all_traits.csv`
> *step-specific metadata containing information about (previous) UK lineage designations, with new acctrans and deltrans lineage designations*
> 


### Steps

1. Merge and create new UK lineage definitions and apply them to all CoG-UK samples:

	>
	> TO DO - more detail
	>

2. Annotate the full tree with UK lineage and acctrans/deltrans information

3. Cut out UK lineage-specific trees from the full tree
4. Phylotyping - TO DO - more detail
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

`...`
> *...*
> 
