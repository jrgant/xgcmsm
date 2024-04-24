#!/bin/bash
#SBATCH -J sim1_maxrim10
#SBATCH --time=2:00:00
#SBATCH -p batch
#SBATCH --mem=3GB
#SBATCH -n 1
#SBATCH --array=1-10000
#SBATCH -o /users/jgantenb/data/xgcmsm-01-calibration/simulator_run/sim1_maxrim10/.logs/LHS_sim1_maxrim10_ARRAY-%A_JOB-%J_SIMNO-%4a.log
#SBATCH --export=ALL,NSIMS=1,NSTEPS=3120,ARRIVE_RATE_ADD_PER20K=1.285,SIMDIR=~/scratch/sim1_maxrim10,SLURM_ARRAY_MAP_OFFSET=$OFFSET
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jrgant@brown.edu
cd /users/jgantenb/data/xgcmsm-01-calibration
/users/jgantenb/build-sources/R-4.3.3/bin/Rscript /users/jgantenb/data/xgcmsm-01-calibration//users/jgantenb/data/xgcmsm-01-calibration/simulator_run/sim1_maxrim10_02_run_simulator.r --vanilla
