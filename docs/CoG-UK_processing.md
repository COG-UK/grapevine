## CoG-UK processing

### Rules

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
