#!/bin/bash
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH -t 0-6:00:00
#SBATCH -n 1
#SBATCH -p general


source activate wes_gatk_v2

nextflow run main.nf -profile cluster -with-report
