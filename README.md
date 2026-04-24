## DISCLAIMER
Repository moved to GitLab: https://gitlab.com/rc-reding/OxBreaker

# OxBreaker
Nextflow implementation of the pipeline described in [OxBreaker: species-agnostic
pipeline for the analysis of outbreaks using nanopore
sequencing](https://www.biorxiv.org/content/10.64898/2026.03.18.709804). The
pipeline takes a set of reads, sequenced using Oxford Nanopore Technologies
(ONT), and determine whether they are closely related.

The aim of this pipeline is to enable healthcare professionoals and
Infection Prevention and Control (IPC) teams idenfity incipient outbreaks
in their hospitals, and manage the available resources in case of no suspected
outbreak.

The implementation found here works in computers running GNU/Linux or macOS with
[Nextflow](https://www.nextflow.io/), and its associated dependencies installed,
and has the following syntax and options:

```
-entry oxbreaker
```
Nextflow program to be run.

```
-profile hpc|standard
```
Sets whether the pipeline uses Slurm workload manager, common
in research high performance computing (HPC) environments, or not. `standard` is
the preferred option if OxBreaker is run on a local workstation.

```
--input STR
```
Full path containing the location of the `.fastq.gz` files to analyse.

```
--output STR
```
Full path where the results from the pipeline will be saved.

```
--reference STR
```
Full path of the file containing the reference genome, in _genbank_ format. This
parameter is REQUIRED if the [Krake2 Microbial Database
(30GB)](https://lomanlab.github.io/mockcommunity/mc_databases.html) is not used.

```
--db STR
```
Full path of the directory containing the [Krake2 Microbial Database
(30GB)](https://lomanlab.github.io/mockcommunity/mc_databases.html). This
parameter is REQUIRED if no reference genome is provided.

```
--min_mq INT
```
Minimum mapping quality below which OxBreaker will discard the reads, in
PHRED format ranging 0-60 (defaults to `55`).

```
--min_read_number INT
```
Minimum supporting reads for OxBreaker to consider a variant with respect to
the reference genome (defaults to `10`).

```
--min_freq FLOAT
```
Minimum variant frequency for OxBreaker to consider a variant with respect to
the reference genome, ranging from 0.0 to 1.0 (defaults to `0.9`).

```
-phylogeny
```
Produce a phylogeny as part of the outputs of OxBreaker.

```
--debug
```
Root phylogeny on the reference given, useful to check if the results are
meaninful.

## Experimental
OxBreaker uses [bcftools](https://samtools.github.io/bcftools/bcftools.html) as
backend for variant-calling. But it can also use [Clair3](https://github.com/HKU-BAL/Clair3).
However, given the [current state of Clair3 for bacterial
genomes](), these settings remain experimental as further studies are run:

```
--clair3
```
Enables the use of Clair3 to find variants with respect the reference genome
provided.

```
--model STR
```
Full path of the directory containing the [model]() used to train Clair3 based
on chemistry and other ONT characteristics. This parameter is REQUIRED if
`-clair3` is used.

```
--chemistry R9|R10
```
Specify chemistry used for the ONT sequencing run (defaults to `R10`).


## Sample script
The code below is an example of how to run locallyl OxBreaker, with default
settings, adding a phylogeny rooted in `this_reference.gb` as part of the output
files:

```
#!/usr/bin/bash
READS_DIR=/directory/containing/data
OUT_DIR=/directory/project/output

    # Run OxBreaker
    nextflow run main.nf -resume \
    -entry outbreaker \
    -profile standard \
    --input $READS_DIR \
    --output $OUT_DIR \
    --reference /directory/references/this_reference.gb \
    --phylogeny --debug
```

Carlos Reding &copy; 2026. University of Oxford.
