sbatch -J SENS_02.10_GCSYMPT_PHAR_0.05_RECT_0.283_URETH_0.57_SCREENTYPE_STI_UNIV -o SENS_02.10_ARRAY-%A_JOB-%J_SIMNO-%4a.log -t 5:00:00  --export=ALL,SIMDIR=~/scratch/SENS_02.10_GCSYMPT_PHAR_0.05_RECT_0.283_URETH_0.57_SCREEN_STI_UNIV,EPI_RUN_TYPE=scenario,STI_SCREEN_TYPE=universal,STI_SCREEN_KISS_EXPOSURE=FALSE,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,PHAR_GC_SYMPT_PROB_ALTPARAM=0.04998510913675,RECT_GC_SYMPT_PROB_ALTPARAM=0.282855302880887,URETH_GC_SYMPT_PROB_ALTPARAM=0.569794693487193 02.01_sti_univ.sh
