# grapevine
snakemake pipelines to run all the funky funcs

## Install
```
git clone https://github.com/cov-ert/grapevine.git
conda install -f environment.yml
conda activate grapevine
```

## The pipeline

1) Import CoG-UK consensus fasta file

2) Clean and filter

3) Align to reference

4) Lineage type using PANGOLIN

5) Merge with GISAID alignment

6) Split into lineages of approximately equal size - A, B, B.1

7) Build IQ-Tree maximum likelihood trees for each

```
iqtree -m GTR+G -nt AUTO -czb -bb 1000 -s <sequences.fasta>
```
> uses `iq-brownie` to maintain taxon names unscathed

8) Clip out UK clusters using `cluster-funk`

9) Build report

10) trotters up

## Outputs

### `phylogenetics/alignment/`
`cog_2020-05-08_all.fasta`
> *all CoG-UK sequences passing Majora's QC*

`cog_2020-05-08_all_alignment.fasta`  
> *all CoG-UK sequences passing Majora's QC aligned to reference*

> `cog_2020-05-08_all_metadata.csv`
> *metadata file for above*
> Columns: `adm0,adm1,adm2,biosample_source_id,central_sample_id,collection_date,flowcell_id,flowcell_type,instrument_make,instrument_model,layout_insert_length,l
ayout_read_length,library_layout_config,library_name,library_primers,library_protocol,library_selection,library_seq_kit,library_seq_protocol,library_s
ource,library_strategy,meta.artic.primers,meta.artic.protocol,metadata,published_as,received_date,root_sample_id,run_group,run_name,sample_type_collec
ted,sample_type_received,sampling_strategy,secondary_accession,secondary_identifier,sequencing_org,sequencing_org_code,sequencing_submission_date,sequ
encing_uuid,source_age,source_sex,start_time,submission_org,submission_org_code,submission_user,swab_site,header,sequence_name,length,missing,gaps,cov
_id,subsample_omit,edin_epi_week,d614g`

`cog_2020-05-08_alignment.fasta`
> *phylogenetic subset (padded to only CDS)*

> `cog_2020-05-08_metadata.csv`
> *metadata file for above*
> Columns: `sequence_name,sample_date,epi_week,country,adm1,adm2,outer_postcode,is_surveillance,is_community,is_hcw,is_travel_history,travel_history,lineage,linea
ge_support,uk_lineage,acc_lineage,del_lineage,phylotype`

### `phylogenetics/trees/`

`cog_global_2020-05-08_tree.nexus`
> tree of CoG and Global genomes

`cog_global_2020-05-08_metadata.csv`
> *metadata file for above*
> Columns: `sequence_name,sample_date,epi_week,country,adm1,adm2,outer_postcode,is_surveillance,is_community,is_hcw,is_travel_history,travel_history,lineage,linea
ge_support,uk_lineage,acc_lineage,del_lineage,phylotype`


### `phylogenetics/public/`

`cog_2020-05-08_sequences.fasta`
> *all CoG-UK sequences passing Majora's QC*

`cog_2020-05-08_alignment.fasta`
> *all CoG-UK sequences as alignment to reference*

`cog_2020-05-08_metadata.csv`
> *metadata file for above*
> Columns: `sequence_name,country,adm1,sample_date,epi_week,lineage,lineage_support`

`cog_global_2020-05-08_tree.newick`
> *tree of UK phylogenetic subset and GISAID*
