## Tree building and ancestral reconstruction

### Rules

`4_make_trees.smk`

`4_subroutine/4_process_lineages.smk`

### Input

`config['output_path'] + /3/cog_gisaid.fasta`
> *All unaligned CoG-UK sequences*
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
3. Annotate the tree with country and previous uk lineage designation
4. **Designate Lineages**
	
	
	
5. Output a full, unannotated public tree and lineage designation information to pass to the next step, defining UK lineages

### Output

`config["output_path"] + "/4/all_traits.csv"`
> *step-specific metadata containing information about (previous) UK lineage designations, with new acctrans and deltrans lineage designations*
> 

`["output_path"] + "/4/cog_gisaid_full.tree.public.newick"`
> *phylogenetic tree containing all GISAID and CoG-UK sequences*
> 
 