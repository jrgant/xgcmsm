
#' @title Condom Use Module
#'
#' @description Module function stochastically simulates potential condom use
#'              for each act on the discordant edgelist.
#'
#' @inheritParams aging_msm
#'
#' @details
#' For each act on the discordant edgelist, condom use is stochastically simulated
#' based on the partnership type and racial combination of the dyad. Other
#' modifiers for the probability of condom use in that pair are diagnosis of
#' disease, and full or partial HIV viral suppression
#' given HIV anti-retroviral therapy.
#'
#' @return
#' Updates the discordant edgelist with a \code{uai} variable indicating whether
#' condoms were used in that act.
#'
#' @import data.table
#' @keywords module msm
#' @export
#'
condoms_msm <- function(dat, at) {

  # Attributes
  race <- dat$attr$race
  age <- dat$attr$age
  diag.status <- dat$attr$diag.status
  prepStat <- dat$attr$prepStat

  # Condom Use Models (anal sex acts only)
  mc_qts3 <- dat$param$epistats$mc_qts$q3
  mc_qts5 <- dat$param$epistats$mc_qts$q5

  cond.mc.mod <- dat$param$epistats$cond.mc
  attr(cond.mc.mod$terms, ".Environment") <- environment()

  cond.oo.mod <- dat$param$epistats$cond.oo
  attr(cond.oo.mod$terms, ".Environment") <- environment()

  cond.scale <- dat$param$cond.scale

  # Temp edgelist
  el <- dat$temp$el


  ## Main/casual partnerships ##
  mc.parts <- which(el[, "ptype"] != 3)
  el.mc <- el[mc.parts, ]


  race.combo <- data.table(
    r1 = race[el.mc[, 1]],
    r2 = race[el.mc[, 2]]
  )

  race.combo <- race.combo[,
    rc := paste(sort(c(r1, r2)), collapse = ""),
    by = seq_len(NROW(race.combo))
    ][, rc]

  diag.status.i <- diag.status.j <- rep(NA, nrow(el.mc))
  diag.status.i <- diag.status[el.mc[, 1]]
  diag.status.j <- diag.status[el.mc[, 2]]

  hiv.concord <- rep(NA, nrow(el.mc))
  hiv.concord[which(diag.status.i == 0 & diag.status.j == 0)] <- 1
  hiv.concord[which(diag.status.i != diag.status.j)] <- 2
  hiv.concord[which(diag.status.i == 1 & diag.status.j == 1)] <- 3

  # Prediction dataset for main/casual condom use
  pred_df <- data.frame(
    ptype = el.mc[, "ptype"],
    durat_wks = el.mc[, "durations"],
    hiv.concord = hiv.concord,
    race.combo = race.combo,
    age.i = age[el.mc[, 1]],
    age.j = age[el.mc[, 2]],
    any.prep = as.numeric((prepStat[el.mc[, 1]] + prepStat[el.mc[, 2]]) > 0)
  )

  cond.prob <- unname(
    predict(cond.mc.mod, newdata = pred_df, type = "response")
  )

  el.mc <- cbind(el.mc, cond.prob)

  ## Clean up
  race.combo <- hiv.concord <- pred_df <- diag.status.i <- diag.status.j <- NULL

  ## One-off partnerships ##
  oo.parts <- which(el[, "ptype"] == 3)
  el.oo <- el[oo.parts, ]

  race.combo <- data.table(
    r1 = race[el.oo[, 1]],
    r2 = race[el.oo[, 2]]
  )

  race.combo <- race.combo[,
    rc := paste(sort(c(r1, r2)), collapse = ""),
    by = seq_len(NROW(race.combo))
    ][, rc]

  diag.status.i <- diag.status.j <- rep(NA, nrow(el.oo))
  diag.status.i <- diag.status[el.oo[, 1]]
  diag.status.j <- diag.status[el.oo[, 2]]

  hiv.concord <- rep(NA, nrow(el.oo))
  hiv.concord[which(diag.status.i == 0 & diag.status.j == 0)] <- 1
  hiv.concord[which(diag.status.i != diag.status.j)] <- 2
  hiv.concord[which(diag.status.i == 1 & diag.status.j == 1)] <- 3

  # Prediction dataset for one-time contact condom use
  pred_df <- data.table(
    age.i = age[el.oo[, 1]],
    age.j = age[el.oo[, 2]],
    hiv.concord = hiv.concord,
    abs_sqrt_agediff = abs(sqrt(age[el.oo[, 1]]) - sqrt(age[el.oo[, 2]])),
    any.prep = as.numeric((prepStat[el.oo[, 1]] + prepStat[el.oo[, 2]]) > 0)
  )[, abs_sqrt_agediff := abs(sqrt(age.i) - sqrt(age.j))]

  cond.prob.oo <- unname(
    predict(cond.oo.mod, newdata = pred_df, type = "response")
  )

  el.oo <- cbind(el.oo, cond.prob.oo)

  ## Bind el together
  el <- rbind(el.mc, el.oo)

  # anal sex acts
  ai.vec <- el[, "ai"]
  pid <- rep(seq_len(length(ai.vec)), ai.vec)
  p1 <- rep(el[, "p1"], ai.vec)
  p2 <- rep(el[, "p2"], ai.vec)
  ptype <- rep(el[, "ptype"], ai.vec)
  cond.prob <- rep(el[, "cond.prob"], ai.vec)

  cond.prob <- cond.prob * cond.scale

  ## Draw for whether anal sex act is unprotected.
  uai <- rbinom(length(cond.prob), 1, 1 - cond.prob)

  # Anal act list construction
  al <- cbind(p1, p2, ptype, uai, pid)
  dat$temp$al <- al

  # Oral act list construction
  pid <- p1 <- p2 <- ptype <- NULL
  oi.vec <- el[, "oi"]
  pid <- rep(seq_len(length(oi.vec)), oi.vec)
  p1 <- rep(el[, "p1"], oi.vec)
  p2 <- rep(el[, "p2"], oi.vec)
  ptype <- rep(el[, "ptype"], oi.vec)

  ol <- cbind(p1, p2, ptype, pid)
  dat$temp$ol <- ol

  # Kissing list construction

  pid <- p1 <- p2 <- ptype <- NULL
  kiss.vec <- el[, "kiss"]
  pid <- rep(seq_len(length(kiss.vec)), kiss.vec)
  p1 <- rep(el[, "p1"], kiss.vec)
  p2 <- rep(el[, "p2"], kiss.vec)
  ptype <- rep(el[, "ptype"], kiss.vec)

  kiss <- cbind(p1, p2, ptype, pid)
  dat$temp$kiss <- kiss


  # Rimming list construction
  pid <- p1 <- p2 <- ptype <- NULL

  rim.vec <- el[, "ri"]

  pid <- rep(seq_len(length(rim.vec)), rim.vec)
  p1 <- rep(el[, "p1"], rim.vec)
  p2 <- rep(el[, "p2"], rim.vec)
  ptype <- rep(el[, "ptype"], rim.vec)

  ri <- cbind(p1, p2, ptype, pid)
  dat$temp$ri <- ri

  return(dat)
}
