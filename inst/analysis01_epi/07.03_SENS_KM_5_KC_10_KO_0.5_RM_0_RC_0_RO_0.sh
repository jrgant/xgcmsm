sbatch -J SENS_07.03_KM_5_KC_10_KO_0.5_RM_0_RC_0_RO_0 -o SENS_07.03_ARRAY-%A_JOB-%J_SIMNO-%4a.log --export=ALL,SIMDIR=~/scratch/SENS_07.03_KM_5_KC_10_KO_0.5_RM_0_RC_0_RO_0,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,KISS_RATE_MAIN_ALTPARAM=5,KISS_RATE_CASUAL_ALTPARAM=10,KISS_PROB_ONETIME_ALTPARAM=0.5,RIM_RATE_MAIN_ALTPARAM=0,RIM_RATE_CASUAL_ALTPARAM=0,RIM_PROB_ONETIME_ALTPARAM=0 02_main.sh
