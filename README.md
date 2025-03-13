# Sort and Assemble
[![Open Source Starter Files](https://github.com/nrminor/sort-and-assemble/actions/workflows/open-source-starter.yaml/badge.svg)](https://github.com/nrminor/sort-and-assemble/actions/workflows/open-source-starter.yaml) [![Docker CI](https://github.com/nrminor/sort-and-assemble/actions/workflows/docker-image.yaml/badge.svg)](https://github.com/sort-and-assemble/ALPINE/actions/workflows/docker-image.yaml) [![Go](https://github.com/nrminor/sort-and-assemble/actions/workflows/go.yml/badge.svg)](https://github.com/nrminor/sort-and-assemble/actions/workflows/go.yml)

## Overview
This pipeline takes Oxford Nanopore reads from antibody chain amplicons, sorts out the reads from each amplicon, and assembles them into the most likely haplotypes using somatic variants.

In short, the pipeline does the following:
1. Merges compressed Nanopore FASTQs based on sample barcodes with `seqkit scat`.
2. Discovers adapter sequences with `bbmerge` for usage in `bbduk` in step 4.
3. Splits each sample's FASTQ into one FASTQ for each primer of interest, as well as an unprimed FASTQ. To do this, it uses `seqkit grep` to search for the forward or reverse primer sequence in each read.
4. Uses `bbduk` to trim to a minimum read length, minimum quality, and to remove adapter and primer sequences.
5. Uses `seqkit stats` to output statistics for each sample's amplicons.
6. Uses `csvtk` to visualize some of those statistics.
7. Runs the `amplicon_sorter` asssembler on the reads for each primer of interest to provide a FASTA of antibody heavy and light haplotype consensus sequences.
8. Searches IgBLAST for each of those consensus sequences to find the closest matches. (*In development!*)

## Quick Start

If you already have [Docker](https://www.docker.com/get-started/) and [NextFlow](https://www.nextflow.io/) installed on your system, simply run the following command in the directory of your choice:

```
nextflow run fruggles11/sort-and-assemble -latest --fastq_dir /path/to/fastq_files.fastq.gz
```

This command automatically pulls the workflow from GitHub and runs it. If you do not have Docker and NextFlow installed, or want to tweak any of the default configurations in the workflow, proceed to the following sections.

## Detailed Setup Instructions

To run this workflow on a MacOS or Linux machine, make sure you have `git` installed, and then simply `git clone` it into your working directory of choice:

```
git clone https://github.com/nrminor/sort-and-assemble.git .
```

You will also need to install the Docker engine if you haven't already. The workflow pulls all the software it needs automatically from Docker Hub, which means you will never need to permanently install that software on your system. To install Docker, simply visit the Docker installation page at [https://docs.docker.com/get-docker/](https://docs.docker.com/get-docker/).

### Nextflow Installation

This workflow uses the [NextFlow](https://www.nextflow.io/) workflow manager. We recommend you install NextFlow to your system in one of the two following ways:

#### 1) Installation with Conda

1. Install the miniconda python distribution, if you haven't already: [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)
2. Install the `mamba` package installation tool in the command line:
   `conda install -y -c conda-forge mamba`
3. Install Nextflow to your base environment:
   `mamba install -c bioconda nextflow `

#### 2) Installation with curl

1. Run the following line in a directory where you'd like to install NextFlow, and run the following line of code:
   `curl -fsSL https://get.nextflow.io | bash`
2. Add this directory to your $PATH. If on MacOS, a helpful guide can be viewed [here](https://www.architectryan.com/2012/10/02/add-to-the-path-on-mac-os-x-mountain-lion/).

To double check that the installation was successful, type `nextflow -v` into the terminal. If it returns something like `nextflow version 21.04.0.5552`, you are set and ready to proceed.

To run the workflow, simply change into the workflow directory and run the following in your terminal:

```
nextflow run . \
--fastq_dir "~/Documents/fastq_pass" \
--primer_table "resources/primers.csv"
```

This workflow generates a large number of working directories. Unless you are debugging, it may be worth running the workflow together with a `nextflow clean` command, like so:

```
nextflow run . \
--fastq_dir "~/Documents/fastq_pass"
--primer_table "resources/primers.csv" \
&& nextflow clean -f 
```

If the workflow runs partway, but a computer outage or other issue interrupts its progress, no need to start over! Instead, run:

```
nextflow run . \
--fastq_dir "~/Documents/fastq_pass"
--primer_table "resources/primers.csv" \
-resume
```

Finally, the workflow's default configurations tell NextFlow to plot the workflow and record run statistics. However, the plot the workflow, note that NextFlow requires the package GraphViz, which is easiest to install via the intructions on [GraphViz's website](https://graphviz.org/download/).
