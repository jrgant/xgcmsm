sbatch -J SENS_02.07_GCSYMPT_PHAR_0.056_RECT_0.435_URETH_0.802_SCREENTYPE_STI_CDC2 -o SENS_02.07_ARRAY-%A_JOB-%J_SIMNO-%4a.log -t 5:00:00  --export=ALL,SIMDIR=~/scratch/SENS_02.07_GCSYMPT_PHAR_0.056_RECT_0.435_URETH_0.802_SCREEN_STI_CDC2,EPI_RUN_TYPE=scenario,STI_SCREEN_TYPE=cdc,STI_SCREEN_KISS_EXPOSURE=TRUE,STARTTIME=3121,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,PHAR_GC_SYMPT_PROB_ALTPARAM=0.0557674131297335,RECT_GC_SYMPT_PROB_ALTPARAM=0.435270233639652,URETH_GC_SYMPT_PROB_ALTPARAM=0.802103788628546 02.01_sti_cdc2.sh
