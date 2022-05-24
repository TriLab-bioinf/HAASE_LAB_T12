# Snakemake pipeline for processing bulk RNAseq sequencing data

## To load conda environment and install required tools
```
mamba env create -f=environment.yml -n rnaseq
conda activate rnaseq
```

## Run entire pipeline
```
snakemake --profile slurm
```

## To run specific rule
## First try a dry run
```
snakemake --profile slurm -R --until $MY_RULE -n
```
## Then do the actual run
```
snakemake --profile slurm -R --until $MY_RULE
```
## Use screen to check workflow progress when connected remotely in a terminal called workflow
```
screen -r workflow
```

