
#' @title Sexual Acts Module
#'
#' @description Module function for setting the number of sexual acts on the
#'              discordant edgelist.
#'
#' @inheritParams aging_msm
#'
#' @details
#' The number of acts at each time step is specified as a function of the race of
#' both members in a pair and the expected values within black-black, black-white,
#' and white-white combinations. For one-off partnerships, this is deterministically
#' set at 1, whereas for main and causal partnerships it is a stochastic draw
#' from a negative binomial distribution. The number of total acts may further be modified
#' by the level of HIV viral suppression in an infected person.
#'
#' @importFrom rms rcs predictrms
#' @import data.table
#' @keywords module msm
#' @export
#'
acts_msm <- function(dat, at) {

  ## Attributes
  status <- dat$attr$status
  diag.status <- dat$attr$diag.status
  race <- dat$attr$race
  age <- dat$attr$age
  stage <- dat$attr$stage
  vl <- dat$attr$vl
  uid <- dat$attr$uid
  prepStat <- dat$attr$prepStat

  plist <- dat$temp$plist

  # Parameter setup ------------------------------------------------------------
 
  ## Parameters from Epistats object
  ## NOTE:
  ##   We assign the function's environment to the fit$terms environment so
  ##   that predictions relying on spline terms will work.
  mc_qts3 <- dat$param$epistats$mc_qts$q3
  mc_qts5 <- dat$param$epistats$mc_qts$q5

  ai.acts.mc <- dat$param$epistats$ai.acts.mc
  attr(ai.acts.mc$terms, ".Environment") <- environment()
  ai.acts.mc.theta <- dat$param$epistats$ai.acts.mc.theta

  oi.acts.mc <- dat$param$epistats$oi.acts.mc
  attr(oi.acts.mc$terms, ".Environment") <- environment()
  oi.acts.mc.theta <- dat$param$epistats$oi.acts.mc.theta

  otp_qts5 <- dat$param$epistats$otp_qts$q5
  otp_qts3 <- dat$param$epistats$otp_qts$q3

  ai.acts.oo <- dat$param$epistats$ai.acts.oo
  attr(ai.acts.oo$terms, ".Environment") <- environment()

  oi.acts.oo <- dat$param$epistats$oi.acts.oo
  attr(oi.acts.oo$terms, ".Environment") <- environment()

  ## Other parameters
  acts.aids.vl <- dat$param$acts.aids.vl
  ai.acts.scale <- dat$param$ai.acts.scale
  oi.acts.scale <- dat$param$oi.acts.scale

  ## Construct edgelist
  el <- rbind(dat$el[[1]], dat$el[[2]], dat$el[[3]])
  ptype  <- rep(1:3, times = c(nrow(dat$el[[1]]),
                               nrow(dat$el[[2]]),
                               nrow(dat$el[[3]])))

  st1 <- status[el[, 1]]
  st2 <- status[el[, 2]]

  el <- cbind(el, st1, st2, ptype)
  colnames(el) <- c("p1", "p2", "st1", "st2", "ptype")


  # Main/casual partnerships ---------------------------------------------------

  ## Subset to main/casual
  el.mc <- el[el[, "ptype"] != 3, ]

  ## Align partner characteristics
  age.i <- age.j <- rep(NA, nrow(el.mc))
  age.i <- age[el.mc[, 1]]
  age.j <- age[el.mc[, 2]]

  abs_sqrt_agediff <- rep(NA, nrow(el.mc))
  abs_sqrt_agediff <- abs(sqrt(age.i) - sqrt(age.j))

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

  ## Current partnership durations
  pid_plist <- plist[, 1] * 1e7 + plist[, 2]
  pid_el <- uid[el.mc[, 1]] * 1e7 + uid[el.mc[, 2]]
  matches <- match(pid_el, pid_plist)

  # NOTE This section necessary only if statistical models used splines
  #      to model the effect of partnership duration.
  if (at == 2) {
    md <- dat$param$netstats$netmain$durat_wks
    cd <- dat$param$netstats$netcasl$durat_wks
    me <- dat$param$netstats$netmain$edges
    ce <- dat$param$netstats$netcasl$edges
    mw <- me / sum(me, ce)
    cw <- ce / sum(me, ce)

    init_pstart <- weighted.mean(c(md, cd), w = c(mw, cw))
    plist[, "start"] <- rnorm(length(plist[, "start"]), -init_pstart, 100)
  }

  durations <- (at - plist[, "start"])[matches]

  ## Simulate anal acts in main/casual partnerships.
  pred_aimc <- data.table(
    ptype = el.mc[, "ptype"],
    durat_wks = durations,
    race.combo = race.combo,
    age.i = age.i,
    age.j = age.j,
    abs_sqrt_agediff, abs_sqrt_agediff,
    hiv.concord = hiv.concord
  )

  ai.acts <- predict(
    ai.acts.mc,
    newdata = pred_aimc,
    type = "response"
  )

  ai.acts.sim <- MASS::rnegbin(ai.acts, theta = ai.acts.mc.theta)
  ai.rates <- ai.acts.sim / 52
  ai <- round(ai.rates * ai.acts.scale)

  ## Simulate oral acts in main/casual partnerships.

  any.prep <- rep(0, nrow(el.mc))
  any.prep[which(prepStat[el.mc[, 1]] == 1 | prepStat[el.mc[, 2]] == 1)] <- 1

  pred_oimc <- data.table(
    any.prep = any.prep,
    race.combo = race.combo,
    ptype = el.mc[, "ptype"],
    age.i = age.i,
    age.j = age.j,
    abs_sqrt_agediff = abs_sqrt_agediff,
    durat_wks = durations
  )

  oi.acts <- predict(
    oi.acts.mc,
    newdata = pred_oimc,
    type = "response"
  )

  oi.acts.sim <- MASS::rnegbin(oi.acts, theta = oi.acts.mc.theta)
  oi.rates <- oi.acts.sim / 52
  oi <- round(oi.rates * oi.acts.scale)

  el.mc <- cbind(el.mc, durations, ai, oi)

  ## Clean up.
  age.i <- age.j <- abs_sqrt_agediff <- NULL
  race.i <- race.j <- NULL
  diag.status.i <- diag.status.j <- NULL
  race.combo <- hiv.concord <- any.prep <- NULL
  pred_aimc <- pred_oimc <- NULL
  ai.acts <- ai.acts.sim <- ai.rates <- ai <- NULL
  oi.acts <- oi.acts.sim <- oi.rates <- oi <- NULL


  ## Simulate rimming acts (governed by flag)
  rrmain <- dat$param$rim.rate.main
  rrcasl <- dat$param$rim.rate.casl

  rim.acts <-
    I(el.mc[, "ptype"] == 1) * rrmain +
    I(el.mc[, "ptype"] == 2) * rrcasl

  ri <- rpois(nrow(el.mc), rim.acts)

  el.mc <- cbind(el.mc, ri)


  ## Simulate kissing acts (governed by flag)

  kiss.rates <- ifelse(
    el.mc[, "ptype"] == 1,
    dat$param$kiss.rate.main,
    dat$param$kiss.rate.casl
  )

  kiss <- rep(NA, nrow(el.mc))
  kiss  <- rpois(nrow(el.mc), kiss.rates)

  el.mc <- cbind(el.mc, kiss)


  # One-time sexual contacts ---------------------------------------------------

  el.oo <- el[el[, "ptype"] == 3, ]

  age.i <- age.j <- rep(NA, nrow(el.oo))
  age.i <- age[el.oo[, 1]]
  age.j <- age[el.oo[, 2]]

  race.i <- race.j <- rep(NA, nrow(el.oo))
  race.i <- race[el.oo[, 1]]
  race.j <- race[el.oo[, 2]]

  # Clear existing race.combo object first.
  race.combo <- data.table(r1 = race.i, r2 = race.j)
  race.combo[,
    rc := paste(sort(c(r1, r2)), collapse = ""),
    by = seq_len(NROW(race.combo))
    ]
  race.combo <- race.combo[, rc]

  any.prep <- rep(0, nrow(el.oo))
  any.prep[which(prepStat[el.oo[, 1]] == 1 | prepStat[el.oo[, 2]] == 1)] <- 1

  diag.status.i <- diag.status.j <- rep(NA, nrow(el.oo))
  diag.status.i <- diag.status[el.oo[, 1]]
  diag.status.j <- diag.status[el.oo[, 2]]

  hiv.concord <- rep(NA, nrow(el.oo))
  hiv.concord[which(diag.status.i == 0 & diag.status.j == 0)] <- 1
  hiv.concord[which(diag.status.i != diag.status.j)] <- 2
  hiv.concord[which(diag.status.i == 1 & diag.status.j == 1)] <- 3

  ## Simulate anal acts.
  pred_oo <- data.table(
    any.prep = any.prep,
    hiv.concord = hiv.concord,
    race.combo = race.combo,
    age.i = age.i,
    age.j = age.j
  )

  ai.acts <- predict(
    ai.acts.oo,
    newdata = pred_oo,
    type = "response"
  )

  ai.acts.sim <- rbinom(length(ai.acts), 1, prob = ai.acts)
  ai <- ai.acts.sim

  ## Simulate oral acts.
  ## Add variable for the oral act model.
  pred_oo_oi <- copy(pred_oo)
  pred_oo_oi[, abs_sqrt_agediff := abs(sqrt(age.i) - sqrt(age.j))]

  oi.acts <- predict(oi.acts.oo, newdata = pred_oo_oi, type = "response")
  oi.acts.sim <- rbinom(length(oi.acts), 1, prob = oi.acts)
  oi <- oi.acts.sim

  durations <- rep(1, nrow(el.oo))
  el.oo <- cbind(el.oo, durations, ai, oi)


  ## Simulate rimming acts.
  ri <- NULL
  ri <- rep(NA, nrow(el.oo))
  ri <- rbinom(nrow(el.oo), 1, dat$param$rim.prob.oo)
  el.oo <- cbind(el.oo, ri)


  ## Simulate kissing during one-time contacts.
  kiss <- NULL
  kiss <- rep(NA, nrow(el.oo))
  kiss  <- rbinom(nrow(el.oo), 1, dat$param$kiss.prob.oo)
  el.oo <- cbind(el.oo, kiss)


  # Bind el back together-------------------------------------------------------

  el <- rbind(el.mc, el.oo)


  # Post-processing for AIDS ---------------------------------------------------

  # For AIDS cases with VL above acts.aids.vl, reduce their their anal acts to 0
  p1HIV <- which(el[, "st1"] == 1)
  p1AIDS <- stage[el[p1HIV, "p1"]] == 4 & vl[el[p1HIV, "p1"]] >= acts.aids.vl
  el[p1HIV[p1AIDS == TRUE], "ai"] <- 0

  p2HIV <- which(el[, "st2"] == 1)
  p2AIDS <- stage[el[p2HIV, "p2"]] == 4 & vl[el[p2HIV, "p2"]] >= acts.aids.vl
  el[p2HIV[p2AIDS == TRUE], "ai"] <- 0

  # Flip order of discordant edges
  disc <- abs(el[, "st1"] - el[, "st2"]) == 1
  disc.st2pos <- which(disc == TRUE & el[, "st2"] == 1)
  el[disc.st2pos, 1:4] <- el[disc.st2pos, c(2, 1, 4, 3)]

  # Remove inactive edges from el
  el <- el[-which(
              el[, "ai"] == 0 &
              el[, "oi"] == 0 &
              el[, "ri"] == 0 &
              el[, "kiss"] == 0
            ), ]

  # Save edgelist
  dat$temp$el <- el

  return(dat)
}
