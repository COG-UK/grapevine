## How to run `deltrans_intro_treetime.nf` on `climb`

If you don't have one already, you need a Grapevine conda environment. This could go somewhere in your home directory (`cd ~`):
```
git clone https://github.com/COG-UK/grapevine.git
cd grapevine
conda env create -f environment.yml
conda activate grapevine
```

Navigate to the basal directory of a Grapevine run, eg.

`cd /cephfs/covid/bham/raccoon-dog/2020-06-26`

**Open a new [Screen](https://linuxize.com/post/how-to-use-linux-screen/) session**, with a sensible name

`screen -S timetree_nextflow`

**Activate the grapevine environment**

`conda activate grapevine`

Then if you **run this**:

```
nextflow -C /cephfs/covid/nextflow/conf/nextflow.config \
  run $PATH/$TO/detrans_intro_treetime.nf \
  -profile slurm
```

(for example, for Ben, this looks like):

```
climb> nextflow -C /cephfs/covid/nextflow/conf/nextflow.config \
  run /cephfs/covid/bham/climb-covid19-jacksonb/git/grapevine/utilities/detrans_intro_treetime.nf \
  -profile slurm
```

**The pipeline should begin**. At the moment the required input files are hardcoded into `detrans_intro_treetime.nf`, and are relative to `$PWD`:

```
params.alignment = "$PWD/analysis/3/cog_gisaid.fasta"
params.tree = "$PWD/analysis/5/cog_gisaid_full.tree.nexus"
params.metadata = "$PWD/analysis/5/cog_gisaid.lineages.with_all_traits.with_phylotype_traits.csv"
```

The pipeline takes a while to run, so you should detach the current Screen session so that it carries on happily in the background, even if you log out of Climb or are disconnected. To do this, you need to issue the following commands using your keyboard:

`ctrl` + `a`, then `d`

If this is the only Screen session you have running, you can reattach the session at any time by typing:

`screen -r`

Otherwise you can reattach Screen sessions by their ID:

```
climb> screen -ls
There are screens on:
	94676.some_other_session	(Detached)
	12360.nextflow_running	(Detached)
2 Sockets in /var/run/screen/S-climb-covid19-jacksonb.


climb> screen -r 12360
(grapevine) climb> nextflow -C /cephfs/covid/nextflow/conf/nextflow.config \
>   run /cephfs/covid/bham/climb-covid19-jacksonb/git/grapevine/utilities/detrans_intro_treetime.nf \
>   -profile slurm
N E X T F L O W  ~  version 20.01.0
Launching `/cephfs/covid/bham/climb-covid19-jacksonb/git/grapevine/utilities/detrans_intro_treetime.nf` [cheesy_stonebraker] - revision: 1bb4a07438
executor >  slurm (1)
[b1/fa90ae] process > cut_out_trees      [  0%] 0 of 1
[-        ] process > get_treetime_input -
[-        ] process > treetime           -
...
```






