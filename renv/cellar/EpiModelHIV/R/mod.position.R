
#' @title Position Module
#'
#' @description Module function for establishing sexual role or position in each
#'              act on the discordant edgelist.
#'
#' @inheritParams aging_msm
#'
#' @details
#' The sexual role within each act is determined by each nodes "role identity"
#' as exclusively receptive, exclusively insertive, or versatile. This function
#' determines whether the infected or the susceptible partner is the insertive
#' partner for that act. For the first two role identity types, that is
#' deterministic based on identity. For versatile-versatile pairs, this is
#' determined stochastically for each act.
#'
#' @return
#' This function returns the updated discordant edgelist with a \code{ins}
#' attribute for values of whether the infected node is insertive or the
#' susceptible node is insertive for that act.
#'
#' @keywords module msm
#'
#' @export
#'
position_msm <- function(dat, at) {

  al <- dat$temp$al
  ol <- dat$temp$ol

  if (nrow(al) == 0) return(dat)

  # Attributes
  role.class <- dat$attr$role.class
  ins.quot <- dat$attr$ins.quot
  ins.quot.oral <- dat$attr$ins.quot.oral


  # Parameters

  ## Process
  p1.role.class <- role.class[al[, "p1"]]
  p2.role.class <- role.class[al[, "p2"]]

  ins <- rep(NA, length(p1.role.class))
  ins[which(p1.role.class == 0)] <- 1
  ins[which(p1.role.class == 1)] <- 0
  ins[which(p2.role.class == 0)] <- 0
  ins[which(p2.role.class == 1)] <- 1

  # Versatile MSM
  vv <- which(p1.role.class == 2 & p2.role.class == 2)

  p1.ins.prob <-
    ins.quot[al[, 1][vv]] /
    (ins.quot[al[, 1][vv]] + ins.quot[al[, 2][vv]])

  p1.ins <- rbinom(length(vv), 1, p1.ins.prob)

  ins[vv[p1.ins == 1]] <- 1
  ins[vv[p1.ins == 0]] <- 0

  dat$temp$al <- cbind(al, ins)

  ## Oral sex position
  if (nrow(ol) == 0) return(dat)

  p1.ins.oral.prob <-
    ins.quot.oral[ol[, 1]] /
    (ins.quot.oral[ol[, 1]] + ins.quot.oral[ol[, 2]])

  ins.oral <- rbinom(length(p1.ins.oral.prob), 1, p1.ins.oral.prob)

  dat$temp$ol <- cbind(ol, ins.oral)

  ## Rimming position

  ri <- dat$temp$ri
  ins.quot.rim <- dat$attr$ins.quot.rim

  p1.ins.rim.prob <- ins.quot.rim[ri[, 1]] /
    (ins.quot.rim[ri[, 1]] + ins.quot.rim[ri[, 2]])

  ins.rim <- rbinom(length(p1.ins.rim.prob), 1, p1.ins.rim.prob)

  dat$temp$ri <- cbind(ri, ins.rim)


  # Record time of anatomic site exposure.
  # Exposure simply defined as sexual event involving anatomic site.
  ids_exposedRectum <- union(
    dat$temp$al[, "p1"][which(dat$temp$al[, "ins"] == 0)],
    dat$temp$al[, "p2"][which(dat$temp$al[, "ins"] == 1)]
  )

  ids_exposedUrethra <- union(
    c(
      dat$temp$al[, "p1"][which(dat$temp$al[, "ins"] == 1)],
      dat$temp$al[, "p2"][which(dat$temp$al[, "ins"] == 0)]
    ),
    c(
      dat$temp$ol[, "p1"][which(dat$temp$ol[, "ins.oral"] == 1)],
      dat$temp$ol[, "p2"][which(dat$temp$ol[, "ins.oral"] == 0)]
    )
  )

  ids_exposedPharynx <- union(
    dat$temp$ol[, "p1"][which(dat$temp$ol[, "ins.oral"] == 0)],
    dat$temp$ol[, "p2"][which(dat$temp$ol[, "ins.oral"] == 1)]
  )

  ## Add exposures due to kissing a rimming if applicable

  ids_exposedRectum <- union(
    ids_exposedRectum,
    c(
      dat$temp$ri[, "p1"][which(dat$temp$ri[, "ins.rim"] == 0)],
      dat$temp$ri[, "p2"][which(dat$temp$ri[, "ins.rim"] == 1)]
    )
  )

  ids_exposedPharynx <- union(
    ids_exposedPharynx,
    c(
      dat$temp$ri[, "p1"][which(dat$temp$ri[, "ins.rim"] == 1)],
      dat$temp$ri[, "p2"][which(dat$temp$ri[, "ins.rim"] == 0)]
    )
  )

  ids_exposedPharynx <- union(
    ids_exposedPharynx,
    c(
      dat$temp$kiss[, "p1"],
      dat$temp$kiss[, "p2"]
    )
  )

  dat$attr$last.rectal.exp[ids_exposedRectum] <- at
  dat$attr$last.ureth.exp[ids_exposedUrethra] <- at
  dat$attr$last.phar.exp[ids_exposedPharynx] <- at

  ## Output
  return(dat)
}
