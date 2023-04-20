#!/bin/bash
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH -t 0-18:00:00
#SBATCH -n 1
#SBATCH -p general

module load mamba/latest
source activate wes_gatk_v2

nextflow run main.nf -profile cluster -with-report
