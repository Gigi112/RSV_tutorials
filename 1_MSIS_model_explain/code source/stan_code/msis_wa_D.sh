#!/bin/bash

#SBATCH --time=10-
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zhe.zheng@yale.edu
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=30000


module load R/4.1.0-foss-2020b
Rscript ageMSIS_deterministric_annual.R