sbatch -J SENS_KM_0_KC_0_KO_0_RM_5_RC_10_RO_1 --export=ALL,SIMDIR=~/scratch/SENS_KM_0_KC_0_KO_0_RM_5_RC_10_RO_1,NSIMS=1,NSTEPS=3380,KISS_RATE_MAIN_ALTPARAM=0,KISS_RATE_CASUAL_ALTPARAM=0,KISS_PROB_ONETIME_ALTPARAM=0,RIM_RATE_MAIN_ALTPARAM=5,RIM_RATE_CASUAL_ALTPARAM=10,RIM_PROB_ONETIME_ALTPARAM=1 02_main.sh
