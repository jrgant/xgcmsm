
#' @title STI Recovery Module
#'
#' @description Stochastically simulates GC/CT recovery.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
stirecov_msm <- function(dat, at) {

  # Parameters ----------------------------------------------------------

  rgc.ntx.int <- dat$param$rgc.ntx.int
  ugc.ntx.int <- dat$param$ugc.ntx.int
  pgc.ntx.int <- dat$param$pgc.ntx.int
  gc.tx.int <- dat$param$gc.tx.int


  # GC Recovery ----------------------------------------------------------------

  # Untreated (asymptomatic and symptomatic)
  idsRGC_ntx <- which(
    dat$attr$rGC == 1 &
    dat$attr$rGC.infTime < at &
    (is.na(dat$attr$rGC.tx) | dat$attr$rGC.tx == 0) &
    (is.na(dat$attr$rGC.tx.prep) | dat$attr$rGC.tx.prep == 0)
  )

  idsUGC_ntx <- which(
    dat$attr$uGC == 1 &
    dat$attr$uGC.infTime < at &
    (is.na(dat$attr$uGC.tx) | dat$attr$uGC.tx == 0) &
    (is.na(dat$attr$uGC.tx.prep) | dat$attr$uGC.tx.prep == 0)
  )

  idsPGC_ntx <- which(
    dat$attr$pGC == 1 &
    dat$attr$pGC.infTime < at &
    (is.na(dat$attr$pGC.tx) | dat$attr$pGC.tx == 0) &
    (is.na(dat$attr$pGC.tx.prep) | dat$attr$pGC.tx.prep == 0)
  )

  if (dat$control$gcUntreatedRecovDist == "geom") {

    # translate median duration into a weekly probability of
    # GC resolution
    geom_pr <- function(durat) 1 - 2 ^ (-1 / durat)

    recovRGC_ntx <- idsRGC_ntx[
      which(rbinom(length(idsRGC_ntx), 1, geom_pr(rgc.ntx.int)) == 1)]

    recovUGC_ntx <- idsUGC_ntx[
      which(rbinom(length(idsUGC_ntx), 1, geom_pr(ugc.ntx.int)) == 1)]

    recovPGC_ntx <- idsPGC_ntx[
      which(rbinom(length(idsPGC_ntx), 1, geom_pr(pgc.ntx.int)) == 1)]

  } else if (dat$control$gcUntreatedRecovDist == "pois") {

    # rectal recovery
    pr.rgc.recov <- ppois(
      q = at - dat$attr$rGC.infTime[idsRGC_ntx],
      lambda = rgc.ntx.int,
      lower.tail = TRUE
    )

    recovRGC_ntx <- idsRGC_ntx[
      which(rbinom(length(idsRGC_ntx), 1, pr.rgc.recov) == 1)]

    # urethral recovery
    pr.ugc.recov <- ppois(
      q = at - dat$attr$uGC.infTime[idsUGC_ntx],
      lambda = ugc.ntx.int,
      lower.tail = TRUE
    )

    recovUGC_ntx <- idsUGC_ntx[
      which(rbinom(length(idsUGC_ntx), 1, pr.ugc.recov) == 1)]

    # pharyngeal recovery
    pr.pgc.recov <- ppois(
      q = at - dat$attr$pGC.infTime[idsPGC_ntx],
      lambda = pgc.ntx.int,
      lower.tail = TRUE
    )

    recovPGC_ntx <- idsPGC_ntx[
      which(rbinom(length(idsPGC_ntx), 1, pr.pgc.recov) == 1)]

  }

  # Treated (asymptomatic and symptomatic)
  idsRGC_tx <- which(
    dat$attr$rGC == 1 &
    dat$attr$rGC.infTime < at &
    (dat$attr$rGC.tx == 1 | dat$attr$rGC.tx.prep == 1)
  )

  idsUGC_tx <- which(
    dat$attr$uGC == 1 &
    dat$attr$uGC.infTime < at &
    (dat$attr$uGC.tx == 1 | dat$attr$uGC.tx.prep == 1)
  )

  idsPGC_tx <- which(
    dat$attr$pGC == 1 &
    dat$attr$pGC.infTime < at &
    (dat$attr$pGC.tx == 1 | dat$attr$pGC.tx.prep == 1)
  )

  # recovRGC_tx <- idsRGC_tx[which(rbinom(length(idsRGC_tx), 1,
  #                                       1/gc.tx.int) == 1)]
  # recovUGC_tx <- idsUGC_tx[which(rbinom(length(idsUGC_tx), 1,
  #                                       1/gc.tx.int) == 1)]

  recovRGC_tx_draw <- which(
    rbinom(
      length(idsRGC_tx), 1,
      dat$param$rgc.tx.recov.pr[at - dat$attr$rGC.infTime[idsRGC_tx]]
    ) == 1
  )

  recovRGC_tx <- idsRGC_tx[recovRGC_tx_draw]

  recovUGC_tx_draw <- which(
    rbinom(
      length(idsUGC_tx), 1,
      dat$param$ugc.tx.recov.pr[at - dat$attr$uGC.infTime[idsUGC_tx]]
    ) == 1
  )

  recovUGC_tx <- idsUGC_tx[recovUGC_tx_draw]


  recovPGC_tx_draw <- which(
    rbinom(
      length(idsPGC_tx), 1,
      dat$param$pgc.tx.recov.pr[at - dat$attr$pGC.infTime[idsPGC_tx]]
    ) == 1
  )

  recovPGC_tx <- idsPGC_tx[recovPGC_tx_draw]

  ## recovRGC_tx <-
  ##   idsRGC_tx[at - dat$attr$rGC.infTime[idsRGC_tx] >= gc.tx.int]

  ## recovUGC_tx <-
  ##   idsUGC_tx[at - dat$attr$uGC.infTime[idsUGC_tx] >= gc.tx.int]

  ## recovPGC_tx <-
  ##   idsPGC_tx[at - dat$attr$pGC.infTime[idsPGC_tx] >= gc.tx.int]

  recovRGC <- c(recovRGC_ntx, recovRGC_tx)
  recovUGC <- c(recovUGC_ntx, recovUGC_tx)
  recovPGC <- c(recovPGC_ntx, recovPGC_tx)

  dat$attr$rGC[recovRGC] <- 0
  dat$attr$rGC.sympt[recovRGC] <- NA
  dat$attr$rGC.infTime[recovRGC] <- NA
  dat$attr$rGC.tx[recovRGC] <- NA
  dat$attr$rGC.tx.prep[recovRGC] <- NA

  dat$attr$uGC[recovUGC] <- 0
  dat$attr$uGC.sympt[recovUGC] <- NA
  dat$attr$uGC.infTime[recovUGC] <- NA
  dat$attr$uGC.tx[recovUGC] <- NA
  dat$attr$uGC.tx.prep[recovUGC] <- NA

  dat$attr$pGC[recovPGC] <- 0
  dat$attr$pGC.sympt[recovPGC] <- NA
  dat$attr$pGC.infTime[recovPGC] <- NA
  dat$attr$pGC.tx[recovPGC] <- NA
  dat$attr$pGC.tx.prep[recovPGC] <- NA

  dat$attr$anyGC.tx[c(recovRGC, recovUGC, recovPGC)] <- NA

  # Output ------------------------------------------------------------------

  dat$epi$recov.rgc[at] <- length(recovRGC)
  dat$epi$recov.pgc[at] <- length(recovPGC)
  dat$epi$recov.ugc[at] <- length(recovUGC)

  return(dat)
}
