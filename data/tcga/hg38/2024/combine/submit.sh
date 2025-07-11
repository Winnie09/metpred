#!/bin/bash -l
#SBATCH -A hji7_bigmem
#SBATCH --partition=bigmem
#SBATCH --time=47:58:59
#SBATCH --nodes=1
#SBATCH --ntasks=22
#SBATCH -o %x.o%j

ml r/4.0.2
echo $1
Rscript $1
