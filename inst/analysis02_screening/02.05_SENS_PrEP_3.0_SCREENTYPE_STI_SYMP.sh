sbatch -J SENS_02.05_PrEP_3.0_SCREENTYPE_STI_SYMP -o SENS_02.05_ARRAY-%A_JOB-%J_SIMNO-%4a.log -t 5:00:00  --export=ALL,SIMDIR=~/scratch/SENS_02.05_PrEP_3.0_SCREEN_STI_SYMP,EPI_RUN_TYPE=scenario,STI_SCREEN_TYPE=symptomatic,STI_SCREEN_KISS_EXPOSURE=FALSE,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,PREP_SCALE_ALTPARAM=3 02.01_sti_symp.sh
