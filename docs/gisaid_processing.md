## Gisaid processing

### Rules

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
