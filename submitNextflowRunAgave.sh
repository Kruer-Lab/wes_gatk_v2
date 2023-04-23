#!/bin/bash
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH -t 0-18:00:00
#SBATCH -n 1

module load singularity/3.8.0
source activate wes_gatk_v2

nextflow run main.nf -profile clusterAgave -with-report
