sbatch -J SENS_07.09_KM_0_KC_0_KO_0_RM_5_RC_5_RO_0.5 -o SENS_07.09_ARRAY-%A_JOB-%J_SIMNO-%4a.log --export=ALL,SIMDIR=~/scratch/SENS_07.09_KM_0_KC_0_KO_0_RM_5_RC_5_RO_0.5,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,KISS_RATE_MAIN_ALTPARAM=0,KISS_RATE_CASUAL_ALTPARAM=0,KISS_PROB_ONETIME_ALTPARAM=0,RIM_RATE_MAIN_ALTPARAM=5,RIM_RATE_CASUAL_ALTPARAM=5,RIM_PROB_ONETIME_ALTPARAM=0.5 02_main.sh
