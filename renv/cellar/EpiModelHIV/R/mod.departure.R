
#' @title Depature Module
#'
#' @description Module function for simulting both general and disease-related
#'              departures, including deaths, among population members.
#'
#' @inheritParams aging_msm
#'
#' @details
#' Deaths are divided into two categories: general deaths, for which demographic
#' data on age-specific mortality rates applies; and disease-related diseases,
#' for which the rate of death is a function of progression to end-stage AIDS.
#'
#' @return
#' This function returns the updated \code{dat} object accounting for deaths.
#' The deaths are deactivated from the main and casual networks, as those are in
#' \code{networkDynamic} class objects; dead nodes are not deleted from the
#' instant network until the \code{\link{simnet_msm}} module for bookkeeping
#' purposes.
#'
#' @keywords module msm
#' @export
#'
departure_msm <- function(dat, at) {

  ## General departures
  active <- dat$attr$active
  age <- floor(dat$attr$age)
  race <- dat$attr$race
  status <- dat$attr$status
  stage <- dat$attr$stage
  tx.status <- dat$attr$tx.status

  aids.mr <- dat$param$aids.mr
  asmr <- dat$param$netstats$demog$asmr

  idsElig <- which(active == 1)
  rates <- rep(NA, length(idsElig))
  races <- sort(unique(race))

  for (i in seq_along(races)) {
    ids.race <- which(race == races[i])
    ages <- age[ids.race]
    rates[ids.race] <- asmr[match(ages, asmr$age), i + 1, with = FALSE][[1]]
  }

  idsDep <- idsElig[rbinom(length(rates), 1, rates) == 1]

  ## HIV-related deaths
  idsEligAIDS <- which(stage == 4)
  idsDepAIDS <- idsEligAIDS[rbinom(length(idsEligAIDS), 1, aids.mr) == 1]

  idsDepAll <- unique(c(idsDep, idsDepAIDS))
  depHIV <- intersect(idsDepAll, which(status == 1))
  depHIV.old <- intersect(depHIV, which(age >= 65))

  if (length(idsDepAll) > 0) {
    dat$attr$active[idsDepAll] <- 0
    for (i in 1:3) {
      dat$el[[i]] <- tergmLite::delete_vertices(dat$el[[i]], idsDepAll)
    }


    dat$attr <- deleteAttr(dat$attr, idsDepAll)

    if (unique(sapply(dat$attr, length)) != attributes(dat$el[[1]])$n) {
      stop("mismatch between el and attr length in departures mod")
    }
  }

  # Update clinical history
  if (dat$control$save.clin.hist == TRUE & length(idsDepAll) > 0) {
    m <- dat$temp$clin.hist
    for (i in 1:length(m)) {
      m[[i]] <- m[[i]][-idsDepAll, ]
    }
    dat$temp$clin.hist <- m
  }

  ## Summary Output
  ## FIXME should this be using idsDep where not in idsDepAIDS? Check and update
  ##       techdocs if needed.
  dat$epi$dep.gen[at] <- length(idsDep)
  dat$epi$dep.AIDS[at] <- length(idsDepAIDS)
  dat$epi$dep.HIV[at] <- length(depHIV)
  dat$epi$dep.HIV.old[at] <- length(depHIV.old)

  return(dat)
}
