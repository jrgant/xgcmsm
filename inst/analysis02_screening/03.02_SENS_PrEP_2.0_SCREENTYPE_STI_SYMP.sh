sbatch -J SENS_03.02_PrEP_2.0_SCREENTYPE_STI_SYMP -o SENS_03.02_ARRAY-%A_JOB-%J_SIMNO-%4a.log -t 5:00:00  --export=ALL,SIMDIR=~/scratch/SENS_03.02_PrEP_2.0_SCREEN_STI_SYMP,STI_SCREEN_TYPE=symptomatic,STI_SCREEN_KISS_EXPOSURE=FALSE,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,PREP_SCALE_ALTPARAM=2 02.02_sti_symp.sh
