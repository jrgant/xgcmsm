
#' @title Arrivals Module
#'
#' @description Module function for arrivals into the sexually active
#'              population.
#'
#' @inheritParams aging_msm
#'
#' @details
#' New population members are added based on expected numbers of entries,
#' stochastically determined with draws from Poisson distributions. For each new
#' entry, a set of attributes is added for that node, and the nodes are added onto
#' the network objects. Only attributes that are a part of the network model
#' formulae are updated as vertex attributes on the network objects.
#'
#' @return
#' This function updates the \code{attr} list with new attributes for each new
#' population member, and the \code{nw} objects with new vertices.
#'
#' @keywords module msm
#' @export
#'
arrival_msm <- function(dat, at) {

  ## Variables

  # Parameters
  a.rate <- dat$param$a.rate

  ## Process
  num <- dat$epi$num[1]

  nNew <- rpois(1, a.rate * num)

  ## Update Attr
  if (nNew > 0) {
    dat <- setNewAttr_msm(dat, at, nNew)
  }

  # Update Networks
  if (nNew > 0) {
    for (i in 1:3) {
      dat$el[[i]] <- tergmLite::add_vertices(dat$el[[i]], nNew)
    }
  }

  ## Output
  dat$epi$nNew[at] <- nNew
  return(dat)
}


setNewAttr_msm <- function(dat, at, nNew) {

  # Set all attributes NA by default
  dat$attr <- lapply(dat$attr, {
    function(x)
      c(x, rep(NA, nNew))
  })

  newIds <- which(is.na(dat$attr$active))

  # Demographic
  dat$attr$active[newIds] <- rep(1, nNew)
  dat$attr$uid[newIds] <- dat$temp$max.uid + (1:nNew)
  dat$temp$max.uid <- dat$temp$max.uid + nNew

  dat$attr$arrival.time[newIds] <- rep(at, nNew)

  race.dist <- prop.table(table(dat$param$netstats$attr$race))
  race <- sample(sort(unique(dat$attr$race)), nNew, TRUE, race.dist)
  dat$attr$race[newIds] <- race

  dat$attr$age.wk[newIds] <- rep(dat$param$arrival.age * 52, nNew)
  dat$attr$age[newIds] <- rep(dat$param$arrival.age, nNew)
  dat$attr$age.grp[newIds] <- rep(1, nNew)

  # Assign HIV status and related
  dat$attr$status[newIds] <- rep(0, nNew)
  dat$attr$diag.status[newIds] <- rep(0, nNew)

  # Assign GC status
  dat$attr$rGC[newIds] <- dat$attr$GC.timesInf[newIds] <- 0
  dat$attr$uGC[newIds] <- dat$attr$GC.timesInf[newIds] <- 0
  dat$attr$pGC[newIds] <- dat$attr$GC.timesInf[newIds] <- 0

  dat$attr$count.trans[newIds] <- 0

  rates <- dat$param$hiv.test.late.prob[race]
  dat$attr$late.tester[newIds] <- rbinom(length(rates), 1, rates)

  races <- sort(unique(dat$attr$race[newIds]))
  tt.traj <- rep(NA, nNew)
  for (i in races) {
    ids.race <- which(dat$attr$race[newIds] == i)

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

  dat$attr$tt.traj[newIds] <- tt.traj

  # Circumcision
  circ <- rep(NA, nNew)
  for (i in races) {
    ids.race <- which(dat$attr$race[newIds] == i)
    circ[ids.race] <- rbinom(length(ids.race), 1, dat$param$circ.prob[i])
  }
  dat$attr$circ[newIds] <- circ

  # Role
  ns <- dat$param$netstats$attr
  role.class <- rep(NA, nNew)
  for (i in races) {
    ids.race <- which(dat$attr$race[newIds] == i)
    rc.probs <- prop.table(table(ns$role.class[ns$race == i]))
    role.class[ids.race] <- sample(0:2, length(ids.race), TRUE, rc.probs)
  }
  dat$attr$role.class[newIds] <- role.class

  # Insertativity quotients

  ## Anal insertativity
  ins.quot <- rep(NA, nNew)
  ins.quot[dat$attr$role.class[newIds] == 0]  <- 1
  ins.quot[dat$attr$role.class[newIds] == 1]  <- 0
  ins.quot[dat$attr$role.class[newIds] == 2] <- runif(
    sum(dat$attr$role.class[newIds] == 2)
  )
  dat$attr$ins.quot[newIds] <- ins.quot

  ## Oral insertativity
  ins.quot.oral <- rep(NA, nNew)
  ins.quot.oral <- runif(length(ins.quot.oral))
  dat$attr$ins.quot.oral[newIds] <- ins.quot.oral

  ## Rimming insertativity
  ins.quot.rim <- rep(NA, nNew)
  ins.quot.rim <- runif(length(ins.quot.rim))
  dat$attr$ins.quot.rim[newIds] <- ins.quot.rim

  # Degree
  dat$attr$deg.main[newIds] <- 0
  dat$attr$deg.casl[newIds] <- 0

  # One-off risk group
  dat$attr$risk.grp[newIds] <- sample(1:5, nNew, TRUE)

  # PrEP
  dat$attr$prepStat[newIds] <- 0

  # HIV screening
  dat$attr$num.neg.tests[newIds] <- 0

  # Update clinical history
  if (dat$control$save.clin.hist == TRUE & length(newIds) > 0) {
    m <- dat$temp$clin.hist
    for (i in 1:length(m)) {
      new.m <- array(dim = c(length(newIds), dat$control$nsteps))
      m[[i]] <- rbind(m[[i]], new.m)
    }
    dat$temp$clin.hist <- m
  }

  ## Check attributes written as expected
  # cbind(sapply(dat$attr, function(x) is.na(tail(x, 1))))
  return(dat)
}
