sbatch -J SENS_02.09_GCSYMPT_PHAR_0.123_RECT_0.369_URETH_0.759_SCREENTYPE_STI_SYMP -o SENS_02.09_ARRAY-%A_JOB-%J_SIMNO-%4a.log -t 5:00:00  --export=ALL,SIMDIR=~/scratch/SENS_02.09_GCSYMPT_PHAR_0.123_RECT_0.369_URETH_0.759_SCREEN_STI_SYMP,EPI_RUN_TYPE=scenario,STI_SCREEN_TYPE=symptomatic,STI_SCREEN_KISS_EXPOSURE=FALSE,STARTTIME=3121,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,PHAR_GC_SYMPT_PROB_ALTPARAM=0.123042341732643,RECT_GC_SYMPT_PROB_ALTPARAM=0.368565309872326,URETH_GC_SYMPT_PROB_ALTPARAM=0.759238593831365 02.01_sti_symp.sh
