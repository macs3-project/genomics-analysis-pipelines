# Generic Pipeline for RNA-seq analysis using STAR/RSEM/BioC/DESeq2

Currently two conditions are supported. Including QC metrics, differential
genes calling, gene set enrichment analysis and KEGG analysis.

## Visualization of the workflow

![DAG of the pipeline](./dag.png)

## Step1: creating conda environment *da*
Using environment.yml file, which includes the packages needed, to create a conda environment, named *da*.
And install some other required packages, seperately. e.g. DO.db

```
conda env create -f environment.yml # you may use slurm to submit this command
conda install -c bioconda bioconductor-do.db 
```

## Step2: activate environment *da* and run snakemake workflow

```
conda activate da
snakemake -np
snakemake --cores 16  # you may use slurm to submit this command

```
