sbatch -J SENS_KM_5_KC_10_KO_1_RM_0_RC_0_RO_0 --export=ALL,SIMDIR=~/scratch/SENS_KM_5_KC_10_KO_1_RM_0_RC_0_RO_0,NSIMS=1,NSTEPS=3380,KISS_RATE_MAIN_ALTPARAM=5,KISS_RATE_CASUAL_ALTPARAM=10,KISS_PROB_ONETIME_ALTPARAM=1,RIM_RATE_MAIN_ALTPARAM=0,RIM_RATE_CASUAL_ALTPARAM=0,RIM_PROB_ONETIME_ALTPARAM=0 02_main.sh
