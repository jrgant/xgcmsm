sbatch -J SENS_03.13_GCSYMPT_PHAR_0.125_RECT_0.34_URETH_0.713_SCREENTYPE_STI_UNIV -o SENS_03.13_ARRAY-%A_JOB-%J_SIMNO-%4a.log --export=ALL,SIMDIR=~/scratch/SENS_03.13_GCSYMPT_PHAR_0.125_RECT_0.34_URETH_0.713_SCREEN_STI_UNIV,STI_SCREEN_TYPE=universal,STI_SCREEN_KISS_EXPOSURE=FALSE,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,PHAR_GC_SYMPT_PROB_ALTPARAM=0.125032290623446,RECT_GC_SYMPT_PROB_ALTPARAM=0.339845853011114,URETH_GC_SYMPT_PROB_ALTPARAM=0.713356923331563 02_main.sh
