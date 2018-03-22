[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
# PISPINO (PIpits SPIN-Off tools)
### A bioinformatics toolkit for processing NGS data
###### Tools which were originally part of [PIPITS](https://github.com/hsgweon/pipits) (a fungal ITS pipeline), but moved here for more generic use as well as to be managed separately. For Linux and Mac OS only.


## Installation

### Prerequisite: set up conda channels (only for the first time)
> add the [Bioconda](https://bioconda.github.io/index.html) channel as well as the other channels bioconda depends on. It is important to add them in this order (this needs to be done once)

```shell
$ conda config --add channels defaults
$ conda config --add channels conda-forge
$ conda config --add channels bioconda
```

### Create a Conda environment with PISPINO
It is recommended that you use [Conda](https://conda.io/) environment so that you install tools and dependencies in this "sandbox" environment without messing with your system. Don't worry, it's easy - just type the following command.

> create a Conda environment (here named "myNGSenv" but you can choose any name)

```shell
$ conda create -n myNGSenv pispino
```



<br>


## Tool 1: SEQPREP
###For preparing (quality filter, reindex, join, merge etc.) raw data from Illumina sequencing platform for further processing by PIPITS, QIIME etc.


**Prerequisite**

All you need is a directory with your raw FASTQ sequences (can be compressed with .gz or .bzip2 or uncompressed).


**Usage**

Illumina reads are generally provided as demultiplexed FASTQ files where BASESPACE (Illumina software) splits the reads into separate files, one for each barcode.


> first of all, get into your environment you just created

```shell
$ source activate myNGSenv
```

> create a list file to specify (1) sample names, (2) file names of forward reads and (3) file names of reverse reads. This can be done with ***pispino_createreadpairslist*** which will generate a tab-delimited text file for all read-pairs from the directory containing your fresh raw sequences from sequencer. "rawdata" is the directory with your FASTQ sequences

```shell
$ pispino_createreadpairslist -i rawdata -o readpairslist.txt
```

> inspect "readpairslist.txt" to see everything looks right, and once happy, process the data with the following: (see more options by "pispino_seqprep -h")

```shell
$ pispino_seqprep -i rawdata -o pispino_seqprep_output -l readpairslist.txt
```

> to leave the environment

```shell
$ source deactivate
```
