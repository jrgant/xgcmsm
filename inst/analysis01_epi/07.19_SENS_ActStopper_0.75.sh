sbatch -J SENS_07.19_ActStopper_0.75 -o SENS_07.19_ARRAY-%A_JOB-%J_SIMNO-%4a.log --export=ALL,SIMDIR=~/scratch/SENS_07.19_ActStopper_0.75,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,ACT_STOPPER_PROB_ALTPARAM=0.75 02_main.sh
