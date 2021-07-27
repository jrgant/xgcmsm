#!/bin/bash
#SBATCH -J Sim1-CleanLHS
#SBATCH--time=3:00:00
#SBATCH -p batch
#SBATCH--mem=150GB
#SBATCH -n 208
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jrgant@brown.edu

module load R/4.0.3

cd ~/data/jgantenb/xgcmsm/
Rscript ./inst/cal/sim1_01_clean.r --vanilla
