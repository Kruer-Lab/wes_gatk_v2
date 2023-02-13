#!/bin/bash
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH -t 0-48:00:00
#SBATCH -n 2

module load singularity/3.8.0
source activate wes_gatk_v2

nextflow run main.nf -profile cluster -with-report
