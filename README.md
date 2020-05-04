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
