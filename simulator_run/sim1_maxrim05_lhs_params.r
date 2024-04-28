# GENERATE PRIOR RANGES --------------------------------------------------------
# Scenario: sim1_maxrim05

RIM_RATE_MAIN_MAX <- 5
RIM_RATE_CASL_MAX <- 5

ROOTDIR <- "/users/jgantenb/data/xgc/xgcmsm-01-calibration"
SIMRUN_DIR <- "simulator_run"
SCENARIO_NAME <- "sim1_maxrim05"

set.seed(5882300)
source(here::here(SIMRUN_DIR, "00_template_lhs_params.r"))
