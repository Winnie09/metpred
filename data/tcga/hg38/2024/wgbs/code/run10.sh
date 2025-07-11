#!/bin/bash -l
#SBATCH -A hji7_bigmem
#SBATCH --partition=bigmem
#SBATCH --time=47:58:59
#SBATCH --nodes=1
#SBATCH --ntasks=20

ml r/4.0.2
Rscript 10_CpG_location_dnasequence.R $1

