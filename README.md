[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
# PISPINO (PIpits SPIN-Off tools)
> A bioinformatics toolkit for processing NGS data


A bunch of tools which were originally part of [PIPITS](https://github.com/hsgweon/pipits) (fungal ITS pipeline), but now found a new nest to be used more generically and independently managed. For Linux and Mac OS only.


## Installation

It is recommended that you use [Conda](https://conda.io/) environment so that you install tools and dependencies in this "sandbox" environment without messing with your system default. Don't worry, it's easy - just follow the instruction.

> create a Conda environment named "pispino" with python 3.6 as the default

```shell
$ conda create -n pispino python=3.6 conda
```

> get into your environment you just created

```shell
$ source activate pispino
```

> now, you can install pispino and its dependencies in the environment. "-c biobuilds" and "-c bioconda" is there to specify the locations of pispino and and its dependencies

```shell
$ conda install pispino -c biobuilds -c bioconda
```

> once finishes, leave the environment

```shell
$ source deactivate
```


<br>


## Tool 1: SEQPREP
For preparing (Quality filter, reindex, join, merge etc.) raw data from Illumina sequencing platform for further processing by PIPITS, QIIME etc.


### How to run

Illumina reads are generally provided as demultiplexed FASTQ files where BASESPACE (Illumina software) splits the reads into separate files, one for each barcode.


> get into the Conda environment (above)

```shell
$ source activate pispino
```

> First create a list file to specify (1) sample names, (2) file names of forward reads and (3) file names of reverse reads. This can be done with a script called ***pispino_createreadpairslist*** which will generate a tab-delimited text file for all read-pairs from the directory containing your fresh raw sequences from sequencer.

```shell
$ pispino_createreadpairslist -i rawdata -o readpairslist.txt
```


> Then, process the data. Simple!

```shell
$ pispino_seqprep -i rawdata -o pispino_seqprep_output -t 30 -l readpairslist.txt
```


### (Optional) Test with a test dataset

Download a test dataset from [here]() and save it in a directory of your choice (yourdirectory), uncompress it then run:

```shell
$ cd yourdirectory
$ wget __
$ unzip __
$ pispino_createreadpairslist -i rawdata -o readpairslist.txt
$ pispino_seqprep -i rawdata -o pispino_seqprep_output -t 30 -l readpairslist.txt
```
