sbatch -J SENS_03.04_PrEP_3.0_SCREENTYPE_STI_BASE -o SENS_03.04_ARRAY-%A_JOB-%J_SIMNO-%4a.log --export=ALL,SIMDIR=~/scratch/SENS_03.04_PrEP_3.0_SCREEN_STI_BASE,STI_SCREEN_TYPE=base,STI_SCREEN_KISS_EXPOSURE=FALSE,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,PREP_SCALE_ALTPARAM=3 02_main.sh
