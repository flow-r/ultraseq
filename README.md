
# Ultraseq

![](https://travis-ci.org/flow-r/ultraseq.svg?branch=master)

## Pipelines and the R package are still under active development

The idea is to create a framework addressing some of the following scenarios

**NOTE:** These have been tested limitedly, please file a bug report in-case of errors.



#### Installation

```
# install flowr
Rscript -e 'install.packages("flowr")'
Rscript -e 'library(flowr);flowr::setup()'

git clone https://github.com/flow-r/ultraseq.git

flowr devtools::install pkg=ultraseq/ultraseq
```



**scenario 1**

```
# create a table of fastq files to use
flowr create_sample_sheet
# create a set of commands to run
flowr run ultraseq
# execute commands on the cluster
flowr to_flow
```

**scenario 2**

start by using a single bam

```
flowr run -h
flowr run ultraseq
flowr run ultraseq tbam=$tbam nbam=$nbam oprefix=$oprefix
```

This repository has the folling structure. 
Folder `pipelines`, contains several pipelines used in cancer genome analysis.
The `ultraseq` folder contains a R package, with unit tests and vinettes and a stable source code
enabling the pipelines. Using a R package streamlines, tesinting, building and documentation.


```
├── pipelines
│   ├── bam_preprocess.R
│   ├── fastq_bam_bwa.R
│   ├── fastq_haplotyper.R
│   ├── fastq_mutect.R
│   ├── fastq_mutect_mem.R
│   ├── fastq_mutect_old.R
│   ├── fastq_preprocess.R
│   └── ultraseq.R
└── ultraseq
    ├── R
    │   ├── annotate.R
    │   ├── bwa.R
    │   ├── chipseq.R
    │   ├── depthofcoverage.R
    │   ├── fastqc.R
    │   ├── fetch_genomes.R
    │   ├── freebayes.R
    │   ├── generic.R
    │   ├── haplotyper.R
    │   ├── modules.R
    │   ├── mut2maf.R
    │   ├── mutect.R
    │   ├── parse_vcf.R
    │   ├── picard.R
    │   ├── preprocess.R
    │   ├── samblaster.R
    │   ├── samtools.R
    │   ├── scripture.R
    │   ├── sheets_fastqs.R
    │   ├── split_fq.R
    │   ├── star.R
    │   └── zzz.R
    ├── inst
    │   ├── conf
    │   └── scripts
    ├── man
    └── vignettes
```

## installation

```
# install dependencies
Rscript -e 'install.packages("flowr", repos = "http://cran.rstudio.com")'
# install ultraseq package, from this repository
devtools::install_github("flow-r/ultraseq", subdir = "ultraseq")

# get somatic variant calling workflows
git clone https://github.com/sahilseth/ultraseq.git
```



<!---- 
notes @sethsa

tree ultraseq/ -P *R 

cd ultraseq
chmod u+x ultraseq
bin/ultraseq
flowr should be available in all nodes
flowr ultraseq::merge_sheets


--->
