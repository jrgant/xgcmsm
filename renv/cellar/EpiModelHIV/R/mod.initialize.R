
# MSM -----------------------------------------------------------------

#' @title Initialization Module
#'
#' @description This function initializes the master \code{dat} object on which
#'              data are stored, simulates the initial state of the network, and
#'              simulates disease status and other attributes.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netest}}.
#' @param param An \code{EpiModel} object of class \code{\link{param_msm}}.
#' @param init An \code{EpiModel} object of class \code{\link{init_msm}}.
#' @param control An \code{EpiModel} object of class \code{\link{control_msm}}.
#' @param s Simulation number, used for restarting dependent simulations.
#'
#' @return
#' This function returns the updated \code{dat} object with the initialized values
#' for demographics and disease-related variables.
#'
#' @export
#' @keywords module msm
#'
initialize_msm <- function(x, param, init, control, s) {

  ## Master Data List Setup ##
  dat <- list()
  dat$param <- param
  dat$init <- init
  dat$control <- control

  ## Network Setup ##
  # Initial network simulations
  dat$nw <- list()
  for (i in 1:3) {
    dat$nw[[i]] <- simulate(x[[i]]$fit, basis = x[[i]]$fit$newnetwork)
  }
  nw <- dat$nw

  # Pull Network parameters
  dat$nwparam <- list()
  for (i in 1:3) {
    dat$nwparam[i] <- list(x[[i]][-which(names(x[[i]]) == "fit")])
  }

  # Convert to tergmLite method
  dat <- init_tergmLite(dat)

  ## Nodal Attributes Setup ##
  dat$attr <- param$netstats$attr

  num <- network.size(nw[[1]])
  dat$attr$active <- rep(1, num)
  dat$attr$arrival.time <- rep(1, num)
  dat$attr$uid <- 1:num

  # Circumcision
  rates <- param$circ.prob[dat$attr$race]
  dat$attr$circ <- rbinom(length(rates), 1, rates)

  # Insertativity Quotients

  ## Anal insertativity
  ins.quot <- rep(NA, num)
  role.class <- dat$attr$role.class
  ins.quot[role.class == 0]  <- 1
  ins.quot[role.class == 1]  <- 0
  ins.quot[role.class == 2]  <- runif(sum(role.class == 2))
  dat$attr$ins.quot <- ins.quot

  ## Oral insertativity
  ins.quot.oral <- rep(NA, num)
  ins.quot.oral <- runif(length(ins.quot.oral))
  dat$attr$ins.quot.oral <- ins.quot.oral

  ## Rimming insertativity
  ins.quot.rim <- rep(NA, num)
  ins.quot.rim <- runif(length(ins.quot.rim))
  dat$attr$ins.quot.rim <- ins.quot.rim

  # HIV-related attributes
  dat <- init_status_msm(dat)

  # STI Status
  dat <- init_sti_msm(dat)

  # PrEP-related attributes
  dat$attr$prepClass <- rep(NA, num)
  dat$attr$prepElig <- rep(NA, num)
  dat$attr$prepStat <- rep(0, num)
  dat$attr$prepStartTime <- rep(NA, num)
  dat$attr$prepLastRisk <- rep(NA, num)
  dat$attr$prepLastStiScreen <- rep(NA, num)

  ## Other Setup ##
  dat$stats <- list()
  dat$stats$nwstats <- list()
  dat$temp <- list()
  dat$epi <- list()

  # Prevalence Tracking
  dat$temp$max.uid <- num
  dat <- prevalence_msm(dat, at = 1)

  # Setup Partner List
  plist <- cbind(dat$el[[1]], ptype = 1)
  plist <- rbind(plist, cbind(dat$el[[2]], ptype = 2))
  plist <- cbind(plist, start = 1, stop = NA)
  colnames(plist)[1:2] <- c("p1", "p2")
  dat$temp$plist <- plist

  # Clinical history
  if (dat$control$save.clin.hist == TRUE) {
    # assigned in mod.hivvl.R (not exported to namespace)
    dat <- save_clin_hist(dat, at = 1)
  }

  # Network statistics
  if (dat$control$save.nwstats == TRUE) {
    # assigned in mod.simnet.R (not exported to namespace)
    dat <- calc_nwstats(dat, at = 1)
  }

  # dat$param$netstats <- NULL
  class(dat) <- "dat"
  return(dat)
}


#' @title Initialize the HIV status of persons in the network
#'
#' @description Sets the initial individual-level disease status of persons
#'              in the network, as well as disease-related attributes for
#'              infected persons.
#'
#' @param dat Data object created in initialization module.
#'
#' @export
#' @keywords initiation utility msm
#'
init_status_msm <- function(dat) {

  num <- sum(dat$attr$active == 1)

  # Sub in diag.status from model for status
  status <- dat$attr$diag.status

  # Late (AIDS-stage) tester type
  rates <- dat$param$hiv.test.late.prob[dat$attr$race]
  dat$attr$late.tester <- rbinom(length(rates), 1, rates)

  # Treatment trajectory
  tt.traj <- rep(NA, num)
  races <- sort(unique(dat$attr$race))
  for (i in races) {
    ids.race <- which(dat$attr$race == i)

    tt.traj[ids.race] <-
      sample(
        x = 1:3,
        size = length(ids.race),
        replace = TRUE,
        prob = c(
          dat$param$tt.part.supp[i],
          dat$param$tt.full.supp[i],
          dat$param$tt.dur.supp[i]
        )
      )
  }
  dat$attr$tt.traj <- tt.traj

  ## Infection-related attributes
  dat$attr$status <- status
  idsInf <- which(status == 1)

  age <- dat$attr$age
  min.ages <- min(dat$param$netstats$demog$ages)
  time.sex.active <-
    pmax(1, round(52 * age[idsInf] - 52 * min.ages, 0))

  min.hiv.time <- round(
    dat$param$vl.acute.rise.int + dat$param$vl.acute.fall.int
  )

  max.hiv.time <- dat$param$vl.aids.onset.int

  time.infected <- round(
    pmax(
      min.hiv.time,
      pmin(
        time.sex.active,
        sample(min.hiv.time:max.hiv.time, length(idsInf), TRUE)
      )
    )
  )

  dat$attr$inf.time <- rep(NA, num)
  dat$attr$inf.time[idsInf] <- -time.infected

  dat$attr$stage <- rep(NA, num)
  dat$attr$stage.time <- rep(NA, num)
  dat$attr$aids.time <- rep(NA, num)
  dat$attr$stage[idsInf] <- 3
  dat$attr$stage.time[idsInf] <- time.infected - min.hiv.time

  dat$attr$diag.stage <- rep(NA, num)
  dat$attr$diag.stage[idsInf] <- dat$attr$stage[idsInf]

  dat$attr$vl <- rep(NA, num)
  dat$attr$vl[idsInf] <- dat$param$vl.set.point
  dat$attr$vl.last.usupp <- rep(NA, num)
  dat$attr$vl.last.supp <- rep(NA, num)

  dat$attr$diag.time <- rep(NA, num)
  dat$attr$diag.time[idsInf] <-
    dat$attr$inf.time[idsInf] + round(mean(1 / dat$param$hiv.test.rate))
  dat$attr$last.neg.test <- rep(NA, num)

  dat$attr$tx.status <- rep(NA, num)
  dat$attr$tx.status[idsInf] <- 0
  dat$attr$cuml.time.on.tx <- rep(NA, num)
  dat$attr$cuml.time.on.tx[idsInf] <- 0
  dat$attr$cuml.time.off.tx <- rep(NA, num)
  dat$attr$cuml.time.off.tx[idsInf] <- time.infected
  dat$attr$tx.period.first <- rep(NA, num)
  dat$attr$tx.period.last <- rep(NA, num)
  dat$attr$tx.init.time <- rep(NA, num)

  dat$attr$count.trans <- rep(0, num)
  dat$attr$num.neg.tests <- rep(0, length(status))

  return(dat)
}



#' @title Initialize the STI status of persons in the network
#'
#' @description Sets the initial individual-level disease status of persons
#'              in the network, as well as disease-related attributes for
#'              infected persons.
#'
#' @param dat Data object created in initialization module.
#'
#' @export
#' @keywords initiation utility msm
#'
init_sti_msm <- function(dat) {

  role.class <- dat$attr$role.class
  num <- length(role.class)

  idsUreth <- which(role.class %in% c(0, 2))
  idsRect <- which(role.class %in% c(1, 2))
  idsPhar <- which(role.class %in% 0:2)

  uGC <- rGC <- pGC <- rep(0, num)

  # Initialize GC infection at all anatomic sites

  ## urethral
  idsUGC <- sample(
    idsUreth,
    size = round(dat$init$prev.ugc * num),
    FALSE
  )

  uGC[idsUGC] <- 1

  ## rectal
  idsRGC <- sample(
    setdiff(idsRect, idsUGC),
    size = round(dat$init$prev.rgc * num),
    FALSE
  )

  rGC[idsRGC] <- 1

  ## pharyngeal
  idsPGC <- sample(
    setdiff(idsRect, idsUGC),
    size = round(dat$init$prev.pgc * num),
    FALSE
  )

  pGC[idsPGC] <- 1

  dat$attr$rGC <- rGC
  dat$attr$uGC <- uGC
  dat$attr$pGC <- pGC

  # Set GC gonorrhea symptom status
  dat$attr$rGC.sympt <- dat$attr$uGC.sympt <- dat$attr$pGC.sympt <- rep(NA, num)

  dat$attr$rGC.sympt[rGC == 1] <-
    rbinom(sum(rGC == 1), 1, dat$param$rgc.sympt.prob)

  dat$attr$uGC.sympt[uGC == 1] <-
    rbinom(sum(uGC == 1), 1, dat$param$ugc.sympt.prob)

  dat$attr$pGC.sympt[pGC == 1] <-
    rbinom(sum(pGC == 1), 1, dat$param$pgc.sympt.prob)

  # Set GC infection time
  dat$attr$rGC.infTime <-
    dat$attr$uGC.infTime <-
      dat$attr$pGC.infTime <-
        rep(NA, length(dat$attr$active))

  dat$attr$rGC.infTime[rGC == 1] <- 1
  dat$attr$uGC.infTime[uGC == 1] <- 1
  dat$attr$pGC.infTime[pGC == 1] <- 1

  dat$attr$rGC.timesInf <- rep(0, num)
  dat$attr$rGC.timesInf[rGC == 1] <- 1

  dat$attr$uGC.timesInf <- rep(0, num)
  dat$attr$uGC.timesInf[uGC == 1] <- 1

  dat$attr$pGC.timesInf <- rep(0, num)
  dat$attr$pGC.timesInf[pGC == 1] <- 1

  ## keep track of most recent STI test
  dat$attr$rGC.tx <-
    dat$attr$uGC.tx <-
      dat$attr$pGC.tx <-
        dat$attr$anyGC.tx <- rep(NA, num)

  dat$attr$rGC.tx.prep <-
    dat$attr$uGC.tx.prep <-
      dat$attr$pGC.tx.prep <- rep(NA, num)

  dat$attr$last.rGC.test <-
    dat$attr$last.uGC.test <-
      dat$attr$last.pGC.test <-
        rep(NA, num)

  ## keep track of most recent anatomic site exposure
  dat$attr$last.rectal.exp <-
    dat$attr$last.ureth.exp <-
      dat$attr$last.phar.exp <-
        rep(NA, num)

  return(dat)

}


#' @title Re-Initialization Module
#'
#' @description This function reinitializes an epidemic model to restart at a
#'              specified time step given an input \code{netsim} object.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netsim}}.
#' @inheritParams initialize_msm
#'
#' @details
#' Currently, the necessary components that must be on \code{x} for a simulation
#' to be restarted must be: param, control, nwparam, epi, attr, temp, el, p.
#' TODO: describe this more.
#'
#' @return
#' This function resets the data elements on the \code{dat} master data object
#' in the needed ways for the time loop to function.
#'
#' @export
#' @keywords module msm
#'
reinit_msm <- function(x, param, init, control, s) {

  need.for.reinit <- c(
    "param", "control", "nwparam",
    "epi", "attr", "temp", "el", "p"
  )

  if (!all(need.for.reinit %in% names(x))) {
    stop("x must contain the following elements for restarting: ",
         "param, control, nwparam, epi, attr, temp, el, p",
         call. = FALSE)
  }

  if (length(x$el) == 1) {
    s <- 1
  }

  dat <- list()

  dat$param <- param
  dat$param$modes <- 1
  dat$control <- control
  dat$nwparam <- x$nwparam

  dat$epi <- sapply(x$epi, function(var) var[s])
  names(dat$epi) <- names(x$epi)

  dat$el <- x$el[[s]]
  dat$p <- x$p[[s]]

  dat$attr <- x$attr[[s]]

  if (!is.null(x$stats)) {
    dat$stats <- list()
    if (!is.null(x$stats$nwstats)) {
      dat$stats$nwstats <- x$stats$nwstats[[s]]
    }
  }

  dat$temp <- x$temp[[s]]

  class(dat) <- "dat"
  return(dat)
}
