
#' @title Prevalence Calculations within Time Steps
#'
#' @description This module calculates demographic, transmission, and clinical
#'              statistics at each time step within the simulation.
#'
#' @inheritParams aging_msm
#'
#' @details
#' Summary statistic calculations are of two broad forms: prevalence and
#' incidence. This function establishes the summary statistic vectors for both
#' prevalence and incidence at time 1, and then calculates the prevalence
#' statistics for times 2 onward. Incidence statistics (e.g., number of new
#' infections or deaths) are calculated within the modules as they depend on
#' vectors that are not stored external to the module.
#'
#' @return
#' This function returns the \code{dat} object with an updated summary of current
#' attributes stored in \code{dat$epi}.
#'
#' @keywords module msm
#'
#' @export
#'
prevalence_msm <- function(dat, at) {

  active <- dat$attr$active
  status <- dat$attr$status
  diag.status <- dat$attr$diag.status
  diag.stage <- dat$attr$diag.stage
  diag.time <- dat$attr$diag.time
  aids.time <- dat$attr$aids.time
  inf.time <- dat$attr$inf.time
  race <- dat$attr$race
  age <- dat$attr$age
  tx.init.time <- dat$attr$tx.init.time
  vl <- dat$attr$vl
  vl.last.usupp <- dat$attr$vl.last.usupp
  last.neg.test <- dat$attr$last.neg.test
  stage <- dat$attr$stage

  prepElig <- dat$attr$prepElig
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass

  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  pGC <- dat$attr$pGC

  # Pop Size / Demog

  ## Total N
  dat$epi$num[at] <- sum(active == 1, na.rm = TRUE)
  ##  Race/ethnicity (n)
  dat$epi$num.B[at] <- sum(race == 1, na.rm = TRUE)
  dat$epi$num.H[at] <- sum(race == 2, na.rm = TRUE)
  dat$epi$num.O[at] <- sum(race == 3, na.rm = TRUE)
  dat$epi$num.W[at] <- sum(race == 4, na.rm = TRUE)
  dat$epi$age.mean[at] <- mean(age, na.rm = TRUE)
  dat$epi$s.num[at] <- sum(status == 0, na.rm = TRUE)

  # HIV infection (N)
  dat$epi$i.num[at] <- sum(status == 1, na.rm = TRUE)
  dat$epi$i.num.B[at] <- sum(status == 1 & race == 1, na.rm = TRUE)
  dat$epi$i.num.H[at] <- sum(status == 1 & race == 2, na.rm = TRUE)
  dat$epi$i.num.O[at] <- sum(status == 1 & race == 3, na.rm = TRUE)
  dat$epi$i.num.W[at] <- sum(status == 1 & race == 4, na.rm = TRUE)
  dat$epi$i.num.dx[at] <- sum(diag.status == 1, na.rm = TRUE)

  # Calculate prevalence and incidence

  dat$epi$i.prev[at] <- dat$epi$i.num[at] / dat$epi$num[at]

  dat$epi$i.prev.B[at] <-
    sum(race == 1 & status == 1, na.rm = TRUE) / sum(race == 1, na.rm = TRUE)

  dat$epi$i.prev.H[at] <-
    sum(race == 2 & status == 1, na.rm = TRUE) / sum(race == 2, na.rm = TRUE)

  dat$epi$i.prev.O[at] <-
    sum(race == 3 & status == 1, na.rm = TRUE) / sum(race == 3, na.rm = TRUE)

  dat$epi$i.prev.W[at] <-
    sum(race == 3 & status == 1, na.rm = TRUE) / sum(race == 4, na.rm = TRUE)

  dat$epi$i.prev.dx[at] <- sum(diag.status == 1, na.rm = TRUE) / dat$epi$num[at]
  dat$epi$i.prev.dx.B[at] <-
    sum(race == 1 & diag.status == 1, na.rm = TRUE) / sum(race == 1, na.rm = TRUE)
  dat$epi$i.prev.dx.H[at] <-
    sum(race == 2 & diag.status == 1, na.rm = TRUE) / sum(race == 2, na.rm = TRUE)
  dat$epi$i.prev.dx.O[at] <-
    sum(race == 3 & diag.status == 1, na.rm = TRUE) / sum(race == 3, na.rm = TRUE)
  dat$epi$i.prev.dx.W[at] <-
    sum(race == 4 & diag.status == 1, na.rm = TRUE) / sum(race == 4, na.rm = TRUE)

  dat$epi$ir100[at] <- (dat$epi$incid[at] / sum(status == 0, dat$epi$incid[at], na.rm = TRUE)) * 5200
  dat$epi$ir100.B[at] <-
    (dat$epi$incid.B[at] /
       sum(status == 0 & race == 1, dat$epi$incid.B[at], na.rm = TRUE)) * 5200
  dat$epi$ir100.H[at] <-
    (dat$epi$incid.H[at] /
       sum(status == 0 & race == 2, dat$epi$incid.H[at], na.rm = TRUE)) * 5200
  dat$epi$ir100.O[at] <-
    (dat$epi$incid.O[at] /
       sum(status == 0 & race == 3, dat$epi$incid.O[at], na.rm = TRUE)) * 5200
  dat$epi$ir100.W[at] <-
    (dat$epi$incid.W[at] /
       sum(status == 0 & race == 4, dat$epi$incid.W[at], na.rm = TRUE)) * 5200


  # Care continuum stats (primary)

  ## HIV diagnosis
  dat$epi$cc.dx[at] <-
    sum(diag.status == 1 & inf.time >= 2, na.rm = TRUE) /
      sum(status == 1 & inf.time >= 2, na.rm = TRUE)

  dat$epi$cc.dx.B[at] <-
    sum(diag.status == 1 & inf.time >= 2 & race == 1, na.rm = TRUE) /
      sum(status == 1 & inf.time >= 2 & race == 1, na.rm = TRUE)

  dat$epi$cc.dx.H[at] <-
    sum(diag.status == 1 & inf.time >= 2 & race == 2, na.rm = TRUE) /
      sum(status == 1 & inf.time >= 2 & race == 2, na.rm = TRUE)

  dat$epi$cc.dx.O[at] <-
    sum(diag.status == 1 & inf.time >= 2 & race == 3, na.rm = TRUE) /
      sum(status == 1 & inf.time >= 2 & race == 3, na.rm = TRUE)

  dat$epi$cc.dx.W[at] <-
    sum(diag.status == 1 & inf.time >= 2 & race == 4, na.rm = TRUE) /
      sum(status == 1 & inf.time >= 2 & race == 4, na.rm = TRUE)

  ## AIDS diagnosed
  dat$epi$cc.dx.aids[at] <-
    sum(diag.status == 1 & stage == 4 & inf.time >= 2 &
          aids.time - diag.time <= 52, na.rm = TRUE) /
    sum(diag.status == 1 & inf.time >= 2, na.rm = TRUE)

  dat$epi$cc.dx.aids.B[at] <-
    sum(diag.status == 1 & stage == 4 & inf.time >= 2 &
          aids.time - diag.time <= 52 & race == 1, na.rm = TRUE) /
    sum(diag.status == 1 & inf.time >= 2 & race == 1, na.rm = TRUE)

  dat$epi$cc.dx.aids.H[at] <-
    sum(diag.status == 1 & stage == 4 &  inf.time >= 2 &
          aids.time - diag.time <= 52 & race == 2, na.rm = TRUE) /
    sum(diag.status == 1 & inf.time >= 2 & race == 2, na.rm = TRUE)

  dat$epi$cc.dx.aids.O[at] <-
    sum(diag.status == 1 & stage == 4 &  inf.time >= 2 &
          aids.time - diag.time <= 52 & race == 3, na.rm = TRUE) /
    sum(diag.status == 1 & inf.time >= 2 & race == 2, na.rm = TRUE)

  dat$epi$cc.dx.aids.W[at] <-
    sum(diag.status == 1 & stage == 4 & inf.time >= 2 &
          aids.time - diag.time <= 52 & race == 3, na.rm = TRUE) /
    sum(diag.status == 1 & inf.time >= 2 & race == 3, na.rm = TRUE)

  ## Linked to care
  dat$epi$cc.linked1m[at] <-
    sum(tx.init.time - diag.time <= 4 & diag.time >= 2, na.rm = TRUE) /
    sum(dat$attr$diag.status == 1 & diag.time >= 2, na.rm = TRUE)

  dat$epi$cc.linked1m.B[at] <-
    sum(tx.init.time - diag.time <= 4 & diag.time >= 2 & race == 1, na.rm = TRUE) /
    sum(diag.status == 1 & diag.time >= 2 & race == 1, na.rm = TRUE)

  dat$epi$cc.linked1m.H[at] <-
    sum(tx.init.time - diag.time <= 4 & diag.time >= 2 & race == 2, na.rm = TRUE) /
    sum(diag.status == 1 & diag.time >= 2 & race == 2, na.rm = TRUE)

  dat$epi$cc.linked1m.O[at] <-
    sum(tx.init.time - diag.time <= 4 & diag.time >= 2 & race == 3, na.rm = TRUE) /
    sum(diag.status == 1 & diag.time >= 2 & race == 3, na.rm = TRUE)

  dat$epi$cc.linked1m.W[at] <-
    sum(tx.init.time - diag.time <= 4 & diag.time >= 2 & race == 4, na.rm = TRUE) /
    sum(diag.status == 1 & diag.time >= 2 & race == 4, na.rm = TRUE)

  dat$epi$cc.linked1m.int[at] <-
    sum(tx.init.time - diag.time <= 4 & diag.time >= 3380, na.rm = TRUE) /
    sum(dat$attr$diag.status == 1 & diag.time >= 3380, na.rm = TRUE)

  dat$epi$cc.linked1m.int.B[at] <-
    sum(tx.init.time - diag.time <= 4 & diag.time >= 3380 & race == 1, na.rm = TRUE) /
    sum(diag.status == 1 & diag.time >= 3380 & race == 1, na.rm = TRUE)

  dat$epi$cc.linked1m.int.H[at] <-
    sum(tx.init.time - diag.time <= 4 & diag.time >= 3380 & race == 2, na.rm = TRUE) /
    sum(diag.status == 1 & diag.time >= 3380 & race == 2, na.rm = TRUE)

  dat$epi$cc.linked1m.int.H[at] <-
    sum(tx.init.time - diag.time <= 4 & diag.time >= 3380 & race == 3, na.rm = TRUE) /
    sum(diag.status == 1 & diag.time >= 3380 & race == 3, na.rm = TRUE)

  dat$epi$cc.linked1m.int.O[at] <-
    sum(tx.init.time - diag.time <= 4 & diag.time >= 3380 & race == 3, na.rm = TRUE) /
    sum(diag.status == 1 & diag.time >= 3380 & race == 3, na.rm = TRUE)

  dat$epi$cc.linked1m.int.W[at] <-
    sum(tx.init.time - diag.time <= 4 & diag.time >= 3380 & race == 4, na.rm = TRUE) /
    sum(diag.status == 1 & diag.time >= 3380 & race == 4, na.rm = TRUE)

  ## Viral Suppression

  dat$epi$cc.vsupp[at] <-
    sum(vl <= log10(200) & diag.status == 1 & diag.time >= 2, na.rm = TRUE) /
    sum(diag.status == 1 & diag.time >= 2, na.rm = TRUE)

  dat$epi$cc.vsupp.B[at] <-
    sum(vl <= log10(200) & diag.status == 1 &
          diag.time >= 2 & race == 1, na.rm = TRUE) /
    sum(diag.status == 1 & diag.time >= 2 & race == 1, na.rm = TRUE)

  dat$epi$cc.vsupp.H[at] <-
    sum(vl <= log10(200) & diag.status == 1 &
          diag.time >= 2 & race == 2, na.rm = TRUE) /
    sum(diag.status == 1 & diag.time >= 2 & race == 2, na.rm = TRUE)

  dat$epi$cc.vsupp.O[at] <-
    sum(vl <= log10(200) & diag.status == 1 &
          diag.time >= 2 & race == 3, na.rm = TRUE) /
    sum(diag.status == 1 & diag.time >= 2 & race == 3, na.rm = TRUE)

  dat$epi$cc.vsupp.W[at] <-
    sum(vl <= log10(200) & diag.status == 1 &
          diag.time >= 2 & race == 4, na.rm = TRUE) /
    sum(diag.status == 1 & diag.time >= 2 & race == 4, na.rm = TRUE)

  dat$epi$cc.vsupp.all[at] <-
    sum(vl <= log10(200) & status == 1 & inf.time >= 2, na.rm = TRUE) /
    sum(status == 1 & inf.time >= 2, na.rm = TRUE)

  dat$epi$cc.vsupp.all.B[at] <-
    sum(vl <= log10(200) & status == 1 & inf.time >= 2 & race == 1, na.rm = TRUE) /
    sum(status == 1 & inf.time >= 2 & race == 1, na.rm = TRUE)

  dat$epi$cc.vsupp.all.H[at] <-
    sum(vl <= log10(200) & status == 1 & inf.time >= 2 & race == 2, na.rm = TRUE) /
    sum(status == 1 & inf.time >= 2 & race == 2, na.rm = TRUE)

  dat$epi$cc.vsupp.all.O[at] <-
    sum(vl <= log10(200) & status == 1 & inf.time >= 2 & race == 3, na.rm = TRUE) /
    sum(status == 1 & inf.time >= 2 & race == 3, na.rm = TRUE)

  dat$epi$cc.vsupp.all.W[at] <-
    sum(vl <= log10(200) & status == 1 & inf.time >= 2 & race == 4, na.rm = TRUE) /
    sum(status == 1 & inf.time >= 2 & race == 4, na.rm = TRUE)

  ## Viral suppression duration
  dat$epi$cc.vsupp.dur1y[at] <-
    1 - (sum((at - vl.last.usupp) <= 52 & diag.time >= 2 &
               diag.status == 1, na.rm = TRUE) /
           sum(diag.status == 1 & diag.time >= 2, na.rm = TRUE))
  dat$epi$cc.vsupp.dur1y.B[at] <-
    1 - (sum((at - vl.last.usupp) <= 52 & diag.time >= 2 &
               diag.status == 1 & race == 1, na.rm = TRUE) /
           sum(diag.status == 1 & diag.time >= 2 & race == 1, na.rm = TRUE))
  dat$epi$cc.vsupp.dur1y.H[at] <-
    1 - (sum((at - vl.last.usupp) <= 52 & diag.time >= 2 &
               diag.status == 1 & race == 2, na.rm = TRUE) /
           sum(diag.status == 1 & diag.time >= 2 & race == 2, na.rm = TRUE))

  dat$epi$cc.vsupp.dur1y.O[at] <-
    1 - (sum((at - vl.last.usupp) <= 52 & diag.time >= 2 &
               diag.status == 1 & race == 3, na.rm = TRUE) /
           sum(diag.status == 1 & diag.time >= 2 & race == 3, na.rm = TRUE))

  dat$epi$cc.vsupp.dur1y.W[at] <-
    1 - (sum((at - vl.last.usupp) <= 52 & diag.time >= 2 &
               diag.status == 1 & race == 4, na.rm = TRUE) /
           sum(diag.status == 1 & diag.time >= 2 & race == 4, na.rm = TRUE))

  dat$epi$cc.HIV.mr[at] <- (dat$epi$dep.HIV[at]/dat$epi$i.num[at]) * 52

  # Care continuum stats (secondary)
  dat$epi$cc.test.int[at] <-
    mean(at - last.neg.test[diag.status == 0], na.rm = TRUE)

  dat$epi$cc.test.int.B[at] <-
    mean(at - last.neg.test[diag.status == 0 & race == 1], na.rm = TRUE)

  dat$epi$cc.test.int.H[at] <-
    mean(at - last.neg.test[diag.status == 0 & race == 2], na.rm = TRUE)

  dat$epi$cc.test.int.O[at] <-
    mean(at - last.neg.test[diag.status == 0 & race == 3], na.rm = TRUE)

  dat$epi$cc.test.int.W[at] <-
    mean(at - last.neg.test[diag.status == 0 & race == 4], na.rm = TRUE)

  dat$epi$cc.dx.delay[at] <-
    mean(diag.time[diag.time >= 2] - inf.time[diag.time >= 2], na.rm = TRUE)

  dat$epi$cc.dx.delay.B[at] <- mean(
    diag.time[diag.time >= 2 & race == 1] -
      inf.time[diag.time >= 2 & race == 1],
    na.rm = TRUE
  )

  dat$epi$cc.dx.delay.H[at] <- mean(
    diag.time[diag.time >= 2 & race == 2] -
      inf.time[diag.time >= 2 & race == 2],
    na.rm = TRUE
  )

  dat$epi$cc.dx.delay.O[at] <- mean(
    diag.time[diag.time >= 2 & race == 3] -
      inf.time[diag.time >= 2 & race == 3],
    na.rm = TRUE
  )

  dat$epi$cc.dx.delay.W[at] <- mean(
    diag.time[diag.time >= 2 & race == 4] -
      inf.time[diag.time >= 2 & race == 4],
    na.rm = TRUE
  )

  dat$epi$cc.dx.delay.int[at] <- mean(
    diag.time[diag.time >= 3380] - inf.time[diag.time >= 3380],
    na.rm = TRUE
  )

  dat$epi$cc.dx.delay.int.B[at] <- mean(
    diag.time[diag.time >= 3380 & race == 1] -
      inf.time[diag.time >= 3380 & race == 1],
    na.rm = TRUE
  )

  dat$epi$cc.dx.delay.int.H[at] <- mean(
    diag.time[diag.time >= 3380 & race == 2] -
      inf.time[diag.time >= 3380 & race == 2],
    na.rm = TRUE
  )

  dat$epi$cc.dx.delay.int.O[at] <- mean(
    diag.time[diag.time >= 3380 & race == 3] -
      inf.time[diag.time >= 3380 & race == 3],
    na.rm = TRUE
  )

  dat$epi$cc.dx.delay.int.W[at] <- mean(
    diag.time[diag.time >= 3380 & race == 4] -
    inf.time[diag.time >= 3380 & race == 4],
    na.rm = TRUE
  )

  # dat$epi$cc.tx.any1y[at] <- sum((at - dat$attr$tx.period.last <= 52), na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1, na.rm = TRUE)
  # dat$epi$cc.tx.any1y.B[at] <- sum((at - dat$attr$tx.period.last <= 52) & race == 1, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1 & race == 1, na.rm = TRUE)
  # dat$epi$cc.tx.any1y.H[at] <- sum((at - dat$attr$tx.period.last <= 52) & race == 2, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1 & race == 2, na.rm = TRUE)
  # dat$epi$cc.tx.any1y.W[at] <- sum((at - dat$attr$tx.period.last <= 52) & race == 3, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1 & race == 3, na.rm = TRUE)

  # dat$epi$cc.dx.delay[at] <- mean(dat$attr$diag.time - dat$attr$inf.time, na.rm = TRUE)
  # dat$epi$cc.testpy[at] <- 1-sum((at - dat$attr$last.neg.test) > 52 & status == 0,
  #     is.na(dat$attr$last.neg.test) & status == 0, na.rm = TRUE) /
  #   sum(status == 0)
  # dat$epi$cc.linked[at] <- sum(dat$attr$cuml.time.on.tx > 0, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1, na.rm = TRUE)
  # dat$epi$cc.tx[at] <- sum(dat$attr$tx.status == 1, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1, na.rm = TRUE)
  # dat$epi$cc.tx.ret3m[at] <- sum((at - dat$attr$tx.period.last) <= 52 &
  #       (dat$attr$tx.period.last - dat$attr$tx.period.first) > 13, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1, na.rm = TRUE)
  # dat$epi$cc.vsupp.tt1[at] <- sum(dat$attr$vl <= log10(200) & dat$attr$tt.traj == 1, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1 & dat$attr$tt.traj == 1, na.rm = TRUE)
  # dat$epi$cc.vsupp.tt2[at] <- sum(dat$attr$vl <= log10(200) & dat$attr$tt.traj == 2, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1 & dat$attr$tt.traj == 2, na.rm = TRUE)
  # dat$epi$cc.vsupp.tt3[at] <- sum(dat$attr$vl <= log10(200) & dat$attr$tt.traj == 3, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1 & dat$attr$tt.traj == 3, na.rm = TRUE)


  # HIV screening outcomes
  dat$epi$mean.neg.tests[at] <-
    mean(dat$attr$num.neg.tests[diag.status == 0], na.rm = TRUE)

  dat$epi$mean.neg.tests.B[at] <-
    mean(dat$attr$num.neg.tests[diag.status == 0 & race == 1], na.rm = TRUE)

  dat$epi$mean.neg.tests.H[at] <-
    mean(dat$attr$num.neg.tests[diag.status == 0 & race == 2], na.rm = TRUE)

  dat$epi$mean.neg.tests.W[at] <-
    mean(dat$attr$num.neg.tests[diag.status == 0 & race == 3], na.rm = TRUE)

  dat$epi$test.past.year[at] <-
    sum(at - dat$attr$last.neg.test <= 52 & diag.status == 0, na.rm = TRUE) /
      sum(diag.status == 0, na.rm = TRUE)

  dat$epi$test.past.year.B[at] <-
    sum(at - dat$attr$last.neg.test <= 52 & diag.status == 0 & race == 1,
        na.rm = TRUE) / sum(diag.status == 0 & race == 1, na.rm = TRUE)

  dat$epi$test.past.year.H[at] <-
    sum(at - dat$attr$last.neg.test <= 52 & diag.status == 0 & race == 2,
        na.rm = TRUE) / sum(diag.status == 0 & race == 2, na.rm = TRUE)

  dat$epi$test.past.year.O[at] <-
    sum(at - dat$attr$last.neg.test <= 52 & diag.status == 0 & race == 3,
        na.rm = TRUE) / sum(diag.status == 0 & race == 3, na.rm = TRUE)

  dat$epi$test.past.year.W[at] <-
    sum(at - dat$attr$last.neg.test <= 52 & diag.status == 0 & race == 4,
        na.rm = TRUE) / sum(diag.status == 0 & race == 4, na.rm = TRUE)

  # HIV stage
  dat$epi$hstage.acute[at] <-
    sum(stage %in% 1:2 & diag.time >= 2, na.rm = TRUE) /
    sum(status == 1 & diag.time >= 2, na.rm = TRUE)

  dat$epi$hstage.chronic[at] <-
    sum(stage == 3 & diag.time >= 2, na.rm = TRUE) /
    sum(status == 1 & diag.time >= 2, na.rm = TRUE)

  dat$epi$hstage.aids[at] <-
    sum(stage == 4 & diag.time >= 2, na.rm = TRUE) /
    sum(status == 1 & diag.time >= 2, na.rm = TRUE)

  dat$epi$prepElig[at] <- sum(prepElig == 1, na.rm = TRUE)
  dat$epi$prepElig.B[at] <- sum(prepElig == 1 & race == 1, na.rm = TRUE)
  dat$epi$prepElig.H[at] <- sum(prepElig == 1 & race == 2, na.rm = TRUE)
  dat$epi$prepElig.O[at] <- sum(prepElig == 1 & race == 3, na.rm = TRUE)
  dat$epi$prepElig.W[at] <- sum(prepElig == 1 & race == 4, na.rm = TRUE)

  dat$epi$prepCurr[at] <- sum(prepStat == 1, na.rm = TRUE)
  dat$epi$prepCurr.B[at] <- sum(prepStat == 1 & race == 1, na.rm = TRUE)
  dat$epi$prepCurr.H[at] <- sum(prepStat == 1 & race == 2, na.rm = TRUE)
  dat$epi$prepCurr.O[at] <- sum(prepStat == 1 & race == 3, na.rm = TRUE)
  dat$epi$prepCurr.W[at] <- sum(prepStat == 1 & race == 4, na.rm = TRUE)

  dat$epi$prepCurr.hadr[at] <- sum(prepStat == 1 & prepClass == 3, na.rm = TRUE)

  # STIs
  dat$epi$i.num.gc[at] <- sum((rGC == 1 | uGC == 1 | pGC == 1), na.rm = TRUE)
  dat$epi$i.num.rgc[at] <- sum(rGC == 1, na.rm = TRUE)
  dat$epi$i.num.ugc[at] <- sum(uGC == 1, na.rm = TRUE)
  dat$epi$i.num.pgc[at] <- sum(pGC == 1, na.rm = TRUE)

  dat$epi$prev.gc[at] <- dat$epi$i.num[at] / dat$epi$num[at]
  dat$epi$prev.rgc[at] <- dat$epi$i.num.rgc[at] / dat$epi$num[at]
  dat$epi$prev.ugc[at] <- dat$epi$i.num.ugc[at] / dat$epi$num[at]
  dat$epi$prev.pgc[at] <- dat$epi$i.num.pgc[at] / dat$epi$num[at]

  rgc.anatsite.wks.atrisk <- sum(rGC == 0, dat$epi$incid.rgc[at], na.rm = TRUE)
  ugc.anatsite.wks.atrisk <- sum(uGC == 0, dat$epi$incid.ugc[at], na.rm = TRUE)
  pgc.anatsite.wks.atrisk <- sum(pGC == 0, dat$epi$incid.pgc[at], na.rm = TRUE)

  dat$epi$ir100.rgc[at] <-
    dat$epi$incid.rgc[at] / rgc.anatsite.wks.atrisk * 5200

  dat$epi$ir100.ugc[at] <-
    dat$epi$incid.ugc[at] / ugc.anatsite.wks.atrisk * 5200

  dat$epi$ir100.pgc[at] <-
    dat$epi$incid.pgc[at] / pgc.anatsite.wks.atrisk * 5200

  dat$epi$ir100.gc[at] <- sum(
    dat$epi$ir100.rgc[at],
    dat$epi$ir100.ugc[at],
    dat$epi$ir100.pgc[at]
  )

  dat$epi$ir100.anatsite.yrs.gc[at] <-
    dat$epi$incid.gc[at] /
    sum(
      rgc.anatsite.wks.atrisk,
      ugc.anatsite.wks.atrisk,
      pgc.anatsite.wks.atrisk
    ) * 5200

  return(dat)
}


whichVlSupp <- function(attr, param) {
  which(attr$status == 1 &
        attr$vlLevel <= log10(50) &
        (attr$age - attr$ageInf) * (365 / param$time.unit) >
        (param$vl.acute.topeak + param$vl.acute.toset))
}
