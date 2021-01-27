## How to run `cut_uk_lineages.nf` on `climb`

If you don't have one already, you need a Grapevine conda environment. This could go somewhere in your home directory (`cd ~`):

```
git clone https://github.com/COG-UK/grapevine.git
cd grapevine
conda env create -f environment.yml
conda activate grapevine
```

Edit `cut_uk_lineages.nf` on lines `9-11` so that it points to the phylopipe output for the date of your choosing.

Optionally edit line `17` to define where the output will be written (default is a subdirectory called `uk_lineages` in the directory that you run the script from).

**Open a new [Screen](https://linuxize.com/post/how-to-use-linux-screen/) session**, with a sensible name

`screen -S uk_lineages_nextflow`

**Activate the grapevine environment**

`conda activate grapevine`

Then if you **run this**:

```
nextflow -C /cephfs/covid/nextflow/conf/nextflow.config \
  run $PATH/$TO/cut_uk_lineages.nf \
  -profile slurm
```

(for example, for Ben, this looks like):

```
climb> nextflow -C /cephfs/covid/nextflow/conf/nextflow.config \
  run /cephfs/covid/bham/climb-covid19-jacksonb/git/grapevine/utilities/cut_uk_lineages.nf \
  -profile slurm
```

**The workflow should begin**. 

The pipeline takes a few minutes to run, so you should detach the current Screen session so that it carries on happily in the background, even if you log out of Climb or are disconnected. To do this, you need to issue the following commands using your keyboard:

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
>   run /cephfs/covid/bham/climb-covid19-jacksonb/git/grapevine/utilities/cut_uk_lineages.nf \
>   -profile slurm
N E X T F L O W  ~  version 20.01.0
Launching `/cephfs/covid/bham/climb-covid19-jacksonb/git/grapevine/utilities/cut_uk_lineages.nf` [cheesy_stonebraker] - revision: 1bb4a07438
executor >  slurm (1)
[b1/fa90ae] process > get_samples      [  0%] 0 of 1
[-        ] process > cut_trees -
[-        ] process > extract           -
...
```






