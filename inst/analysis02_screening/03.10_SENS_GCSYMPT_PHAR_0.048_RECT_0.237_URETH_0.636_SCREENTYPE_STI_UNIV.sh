sbatch -J SENS_03.10_GCSYMPT_PHAR_0.048_RECT_0.237_URETH_0.636_SCREENTYPE_STI_UNIV -o SENS_03.10_ARRAY-%A_JOB-%J_SIMNO-%4a.log --export=ALL,SIMDIR=~/scratch/SENS_03.10_GCSYMPT_PHAR_0.048_RECT_0.237_URETH_0.636_SCREEN_STI_UNIV,STI_SCREEN_TYPE=universal,STI_SCREEN_KISS_EXPOSURE=FALSE,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,PHAR_GC_SYMPT_PROB_ALTPARAM=0.0484890314460044,RECT_GC_SYMPT_PROB_ALTPARAM=0.237088431859587,URETH_GC_SYMPT_PROB_ALTPARAM=0.635944305395122 02.05_sti_univ.sh
