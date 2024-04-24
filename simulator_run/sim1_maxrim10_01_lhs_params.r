# GENERATE PRIOR RANGES --------------------------------------------------------
# Scenario: sim1_maxrim10

RIM_RATE_MAIN_MAX <- 10
RIM_RATE_CASL_MAX <- 10

ROOTDIR <- "/users/jgantenb/data/xgcmsm-01-calibration"
SIMRUN_DIR <- "simulator_run"
SCENARIO_NAME <- "sim1_maxrim10"

set.seed(1971)
source(here::here(SIMRUN_DIR, "00_template_lhs_params.r"))
