sbatch -J SENS_03.03_PrEP_2.5_SCREENTYPE_STI_SYMP -o SENS_03.03_ARRAY-%A_JOB-%J_SIMNO-%4a.log --export=ALL,SIMDIR=~/scratch/SENS_03.03_PrEP_2.5_SCREEN_STI_SYMP,STI_SCREEN_TYPE=symptomatic,STI_SCREEN_KISS_EXPOSURE=FALSE,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,PREP_SCALE_ALTPARAM=2.5 02.02_sti_symp.sh
