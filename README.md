# grapevine
Phylogenetic pipeline for the COG-UK project

## Install
```
git clone https://github.com/COG-UK/grapevine.git
cd grapevine
conda env create -f environment.yml
conda activate grapevine
```

## The pipeline

See [the documentation](https://github.com/COG-UK/grapevine/tree/master/docs) for more detailed methods

## The pipeline - in brief

0) Process latest GISAID data 

1) Import CoG-UK consensus fasta file

2) Clean and filter

3) Align to reference

4) Lineage type using PANGOLIN

5) Merge with GISAID alignment

6) (Optionally) split into lineages based on Pangolin typing to relieve some of the burden on tree building

7) Build a Fasttree maximum likelihood tree (for each lineage if so defined):

```
FastTreeMP -nosupport -nt <sequences.fasta> > <output.tree>
```

8) If split, graft together the subtrees into a complete tree

9) Define and extract UK clusters using `cluster-funk`

10) Build reports


## Outputs

### `phylogenetics/alignment/`
`cog_2020-06-12_all.fasta`
> *all CoG-UK sequences passing Majora's QC*

`cog_2020-06-12_all_alignment.fasta`  
> *all CoG-UK sequences passing Majora's QC aligned to reference*

`cog_2020-06-12_all_metadata.csv`
> *metadata file for above*

> Columns: `adm0,adm1,adm2,adm2_private,biosample_source_id,central_sample_id,collected_by,collection_date,end_time,flowcell_id,flowcell_type,instrument_make,instrument_model,layout_insert_length,layout_read_length,library_adaptor_barcode,library_layout_config,library_name,library_primers,library_protocol,library_selection,library_seq_kit,library_seq_protocol,library_source,library_strategy,meta.artic.primers,meta.artic.protocol,meta.epi.cluster,meta.investigation.cluster,meta.investigation.name,meta.investigation.site,metric.ct.1.ct_value,metric.ct.1.test_kit,metric.ct.1.test_platform,metric.ct.1.test_target,metric.ct.2.ct_value,metric.ct.2.test_kit,metric.ct.2.test_platform,metric.ct.2.test_target,metric.ct.max_ct,metric.ct.min_ct,metric.ct.num_tests,published_as,received_date,root_sample_id,run_group,run_name,sample_type_collected,sample_type_received,secondary_accession,secondary_identifier,sequencing_org,sequencing_org_code,sequencing_submission_date,sequencing_uuid,source_age,source_sex,start_time,submission_org,submission_org_code,submission_user,swab_site,header,sequence_name,length,missing,gaps,cov_id,subsample_omit,edin_epi_week,d614g`

`cog_2020-06-12_alignment.fasta`
> *phylogenetic subset (padded with Ns to only CDS)*

`cog_2020-06-12_metadata.csv`
> *metadata file for above*
> 
> Columns: `sequence_name,sample_date,epi_week,country,adm1,adm2,outer_postcode,is_surveillance,is_community,is_hcw,is_travel_history,travel_history,lineage,linea
ge_support,uk_lineage,acc_lineage,del_lineage,phylotype`

### `phylogenetics/trees/`

`cog_global_2020-06-12_tree.newick`
> unannotated tree of CoG and Global genomes

`cog_global_2020-06-12_tree.nexus`
> annotated tree of CoG and Global genomes

`cog_global_2020-06-12_metadata.csv`
> *metadata file for above*
> 
> Columns: `sequence_name,sample_date,epi_week,country,adm1,adm2,outer_postcode,is_surveillance,is_community,is_hcw,is_travel_history,travel_history,lineage,lineage_support,uk_lineage,acc_lineage,del_lineage,acc_introduction,del_introduction,phylotype`

`uk_lineages/`
> containing, for each non-singleton UK lineage:
> 
> `uk_lineage_UK1003.csv`
> `uk_lineage_UK1003.fasta`
> `uk_lineage_UK1003_timetree/`
> 
> an alignment, a metadata csv, and a directory with the output of `TreeTime` (where applicable)

### `phylogenetics/public/`

`cog_2020-06-12.fasta`
> *all CoG-UK sequences passing Majora's QC*

`cog_2020-06-12_alignment.fasta`
> *phylogenetic subset (padded with Ns to only CDS)*

`cog_2020-06-12_metadata.csv`
> *metadata file for above*
> 
> Columns: `sequence_name,country,adm1,sample_date,epi_week,lineage,lineage_support`

`cog_global_2020-06-12_tree.newick`
> *tree of UK phylogenetic subset and GISAID*

### `phylogenetics/reports/`

`UK_report.pdf`
> *Phylogenetic report for the UK*
> 
> `figures/` and `summary_files/` contain figures and data files associated with this report

`adm1_reports/`
> *Reports for the four countries of the UK: Wales, Scotland, Northern Ireland and England (in subdirectories)*

`regional_reports/`
>*Reports by sequencing centre*
>
> Associated figures and summary files at in subdirectories under `results/`

### `phylogenetics/microreact/`

`cog_global_2020-06-12_tree_public.newick`
> *tree of UK phylogenetic subset and GISAID, with all CoG samples given an anonymous ID*

`cog_global_2020-06-12_metadata_public.csv`
> *metadata file for above with matching anonymized sample names. `adm2` values represented by <5 samples are blanked*
> 
> Columns:
> `sequence_name,sample_date,epi_week,country,adm1,adm2,submission_org_code,lineage,lineage_support,uk_lineage,primary_uk_lineage,d614g`

`cog_global_2020-06-12_tree_private.newick`
> *tree of UK phylogenetic subset and GISAID*

`cog_global_2020-06-12_metadata_private.csv`
> *metadata file for above*
>
> Columns:
> `sequence_name,sample_date,epi_week,country,adm1,adm2,submission_org_code,is_hcw,travel_history,lineage,lineage_support,uk_lineage,primary_uk_lineage,d614g`

### `phylogenetics/civet/`

`cog_global_2020-06-12_tree.nexus`
> *annotated tree of UK phylogenetic subset and GISAID*

`cog_2020-06-12_alignment_all.fasta`
> *all CoG-UK sequences passing Majora's QC aligned to reference*

`cog_2020-06-12_metadata_all.csv`  
> *metadata file for above*
> 
> Columns:
> `central_sample_id,biosample_source_id,sequence_name,sample_date,epi_week,country,adm1,adm2,outer_postcode,is_surveillance,is_community,is_hcw,is_travel_history,travel_history,lineage,lineage_support,uk_lineage,acc_lineage,del_lineage,phylotype`

`cog_global_2020-06-12_alignment.fasta`
> *phylogenetic subset and GISAID aligned to reference*

`cog_global_2020-06-12_metadata.csv`
> *metadata file for above*
> 
> Columns:
> `central_sample_id,biosample_source_id,sequence_name,sample_date,epi_week,country,adm1,adm2,outer_postcode,is_surveillance,is_community,is_hcw,is_travel_history,travel_history,lineage,lineage_support,uk_lineage,acc_lineage,del_lineage,phylotype`

`cog_2020-06-12_alignment.fasta`
> *phylogenetic subset aligned to reference*

`cog_2020-06-12_metadata.csv`
> *metadata file for above*
> 
> Columns:
> `central_sample_id,biosample_source_id,sequence_name,sample_date,epi_week,country,adm1,adm2,outer_postcode,is_surveillance,is_community,is_hcw,is_travel_history,travel_history,lineage,lineage_support,uk_lineage,acc_lineage,del_lineage,phylotype`

## Acknowledgements

Grapevine uses the following tools:

```
bioconda
conda
biopython
minimap2
python
snakemake
pandoc
fasttree
gotree
TreeTime (https://github.com/neherlab/treetime)
Pangolin (https://github.com/hCoV-2019/pangolin)
```
