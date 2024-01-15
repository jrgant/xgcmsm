#' @title STI Transmission Module
#'
#' @description Stochastically simulates GC transmission given the current
#'              state of the edgelist. This version of the STI transmission
#'              function randomizes the transmission pathway responsible
#'              for a transmission event when multiple such events occur
#'              at a single anatomic site. Randomization addresses
#'              the depletion of susceptible anatomic sites that would
#'              occur if transmissions were assigned in the same order
#'              of transmission pathways every time, artificially depressing
#'              incidence due to transmission pathways that are processed
#'              later in the sequence.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @import data.table
#' @import magrittr
#' @export
#'
stitrans_msm_rand <- function(dat, at) {

  # Parameters ----------------------------------------------------------

  # Acquisition probabilities given contact with infected man
  u2rgc.tprob <- dat$param$u2rgc.tprob
  r2ugc.tprob <- dat$param$r2ugc.tprob
  u2pgc.tprob <- dat$param$u2pgc.tprob
  p2ugc.tprob <- dat$param$p2ugc.tprob

  # Probability of symptoms given infection
  rgc.sympt.prob <- dat$param$rgc.sympt.prob
  ugc.sympt.prob <- dat$param$ugc.sympt.prob
  pgc.sympt.prob <- dat$param$pgc.sympt.prob

  # Relative risk of infection given condom use during act
  sti.cond.eff <- dat$param$sti.cond.eff
  sti.cond.fail <- dat$param$sti.cond.fail


  # Attributes ----------------------------------------------------------

  # Demographics
  race <- dat$attr$race
  age.grp <- dat$attr$age.grp

  # Current infection state
  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  pGC <- dat$attr$pGC

  # n Times infected
  rGC.timesInf <- dat$attr$rGC.timesInf
  uGC.timesInf <- dat$attr$uGC.timesInf
  pGC.timesInf <- dat$attr$pGC.timesInf

  # Infection time
  rGC.infTime <- dat$attr$rGC.infTime
  uGC.infTime <- dat$attr$uGC.infTime
  pGC.infTime <- dat$attr$pGC.infTime

  # Infection symptoms (non-varying)
  rGC.sympt <- dat$attr$rGC.sympt
  uGC.sympt <- dat$attr$uGC.sympt
  pGC.sympt <- dat$attr$pGC.sympt


  # Pull act lists
  al <- dat$temp$al
  ol <- dat$temp$ol

  ## ins variable coding
  # ins = 0 : p2 is insertive
  # ins = 1 : p1 is insertive

  ##############################################################################
  ## ANAL SEX ##
  ##############################################################################

  # Urethral-to-Rectal GC ------------------------------------------------------

  # P1 infects P2
  p1Inf_u2rgc <- which(
    uGC[al[, "p1"]] == 1 & uGC.infTime[al[, "p1"]] < at &
    rGC[al[, "p2"]] == 0 & al[, "ins"] == 1
  )

  # P2 infects P1
  p2Inf_u2rgc <- which(
    uGC[al[, "p2"]] == 1 & uGC.infTime[al[, "p2"]] < at &
    rGC[al[, "p1"]] == 0 & al[, "ins"] == 0
  )

  allActs_u2rgc <- c(p1Inf_u2rgc, p2Inf_u2rgc)

  # UAI modifier
  uai_u2rgc <- al[allActs_u2rgc, "uai"]
  tprob_u2rgc <- rep(u2rgc.tprob, length(allActs_u2rgc))

  # Transform to log odds
  tlo_u2rgc <- log(tprob_u2rgc / (1 - tprob_u2rgc))

  # Modify log odds by race-specific condom effectiveness
  races <- c(race[al[p1Inf_u2rgc, "p1"]], race[al[p2Inf_u2rgc, "p2"]])
  condom.rr <- rep(NA, length(races))
  for (i in sort(unique(races))) {
    ids.race <- which(races == i)
    condom.rr[ids.race] <- 1 - (sti.cond.eff - sti.cond.fail[i])
  }

  tlo_u2rgc[uai_u2rgc == 0] <-
    tlo_u2rgc[uai_u2rgc == 0] +
    ifelse(length(uai_u2rgc == 0), log(condom.rr[uai_u2rgc == 0]), 0)

  # Back-transform to probability
  tprob_u2rgc <- plogis(tlo_u2rgc)

  # Stochastic transmission
  trans_u2rgc <- rbinom(length(allActs_u2rgc), 1, tprob_u2rgc)

  # Determine the infected partner
  idsInf_u2rgc <- NULL
  if (sum(trans_u2rgc, na.rm = TRUE) > 0) {
    transAL_u2rgc <- al[allActs_u2rgc[trans_u2rgc == 1], , drop = FALSE]
    idsInf_u2rgc <- c(
      intersect(al[p1Inf_u2rgc, "p2"], transAL_u2rgc[, "p2"]),
      intersect(al[p2Inf_u2rgc, "p1"], transAL_u2rgc[, "p1"])
    )
    stopifnot(all(rGC[idsInf_u2rgc] == 0))
  }

  # Rectal-to-Urethral GC ----------------------------------------------------

  # P1 infects P2
  p1Inf_r2ugc <- which(
    rGC[al[, "p1"]] == 1 & rGC.infTime[al[, "p1"]] < at &
    uGC[al[, "p2"]] == 0 & al[, "ins"] == 0
  )

  # P2 infects P1
  p2Inf_r2ugc <- which(
    rGC[al[, "p2"]] == 1 & rGC.infTime[al[, "p2"]] < at &
    uGC[al[, "p1"]] == 0 & al[, "ins"] == 1
  )

  allActs_r2ugc <- c(p1Inf_r2ugc, p2Inf_r2ugc)

  # UAI modifier
  uai_r2ugc <- al[allActs_r2ugc, "uai"]
  tprob_r2ugc <- rep(r2ugc.tprob, length(allActs_r2ugc))

  # Transform to log odds
  tlo_r2ugc <- log(tprob_r2ugc / (1 - tprob_r2ugc))

  # Modify log odds by race-specific condom effectiveness
  races <- c(race[al[p1Inf_r2ugc, "p2"]], race[al[p2Inf_r2ugc, "p1"]])
  condom.rr <- rep(NA, length(races))
  for (i in sort(unique(races))) {
    ids.race <- which(races == i)
    condom.rr[ids.race] <- 1 - (sti.cond.eff - sti.cond.fail[i])
  }

  tlo_r2ugc[uai_r2ugc == 0] <-
    tlo_r2ugc[uai_r2ugc == 0] +
    ifelse(length(uai_r2ugc == 0), log(condom.rr[uai_r2ugc == 0]), 0)

  # Back-transform to probability
  tprob_r2ugc <- plogis(tlo_r2ugc)

  # Stochastic transmission
  trans_r2ugc <- rbinom(length(allActs_r2ugc), 1, tprob_r2ugc)

  # Determine the newly infected partner
  idsInf_r2ugc <- NULL
  if (sum(trans_r2ugc, na.rm = TRUE) > 0) {
    transAL_r2ugc <- al[allActs_r2ugc[trans_r2ugc == 1], , drop = FALSE]
    idsInf_r2ugc <- c(
      intersect(al[p1Inf_r2ugc, "p2"], transAL_r2ugc[, "p2"]),
      intersect(al[p2Inf_r2ugc, "p1"], transAL_r2ugc[, "p1"])
    )
    stopifnot(all(uGC[idsInf_r2ugc] == 0))
  }


  ##############################################################################
  ## ORAL SEX ##
  ##############################################################################

  # Urethral-to-Pharyngeal GC --------------------------------------------------

  # P1 infects P2
  p1Inf_u2pgc <- which(
    uGC[ol[, "p1"]] == 1 & uGC.infTime[ol[, "p1"]] < at &
    pGC[ol[, "p2"]] == 0 & ol[, "ins.oral"] == 1
  )

  # P2 infects P1
  p2Inf_u2pgc <- which(
    uGC[ol[, "p2"]] == 1 & uGC.infTime[ol[, "p2"]] < at &
    pGC[ol[, "p1"]] == 0 & ol[, "ins.oral"] == 0
  )

  allActs_u2pgc <- c(p1Inf_u2pgc, p2Inf_u2pgc)

  # REVIEW: This step may not be necessary and may simply use up memory.
  #         Check for each transmission pathway.
  # Pathway-specific transmission probability
  tprob_u2pgc <- rep(u2pgc.tprob, length(allActs_u2pgc))

  # Stochastic transmission
  trans_u2pgc <- rbinom(length(allActs_u2pgc), 1, tprob_u2pgc)

  # Determine the newly infected partner
  idsInf_u2pgc <- NULL
  if (sum(trans_u2pgc, na.rm = TRUE) > 0) {
    transOL_u2pgc <- ol[allActs_u2pgc[trans_u2pgc == 1], , drop = FALSE]
    idsInf_u2pgc <- c(
      intersect(ol[p1Inf_u2pgc, "p2"], transOL_u2pgc[, "p2"]),
      intersect(ol[p2Inf_u2pgc, "p1"], transOL_u2pgc[, "p1"])
    )
    stopifnot(all(pGC[idsInf_u2pgc] == 0))
  }

  # Pharyngeal-to-Urethral -----------------------------------------------------

  # P1 infects P2
  p1Inf_p2ugc <- which(
    pGC[ol[, "p1"]] == 1 & pGC.infTime[ol[, "p1"]] < at &
    uGC[ol[, "p2"]] == 0 & ol[, "ins.oral"] == 0
  )

  # P2 infects P1
  p2Inf_p2ugc <- which(
    pGC[ol[, "p2"]] == 1 & pGC.infTime[ol[, "p2"]] < at &
    uGC[ol[, "p1"]] == 0 & ol[, "ins.oral"] == 1
  )

  allActs_p2ugc <- c(p1Inf_p2ugc, p2Inf_p2ugc)
  tprob_p2ugc <- rep(p2ugc.tprob, length(allActs_p2ugc))

  # Stochastic transmission with pathway-specific probability
  trans_p2ugc <- rbinom(length(allActs_p2ugc), 1, dat$param$p2ugc.tprob)

  # Determine the newly infected partner
  idsInf_p2ugc <- NULL
  if (sum(trans_p2ugc, na.rm = TRUE) > 0) {
    transOL_p2ugc <- ol[allActs_p2ugc[trans_p2ugc == 1], , drop = FALSE]
    idsInf_p2ugc <- c(
      intersect(ol[p1Inf_p2ugc, "p2"], transOL_p2ugc[, "p2"]),
      intersect(ol[p2Inf_p2ugc, "p1"], transOL_p2ugc[, "p1"])
    )

    stopifnot(all(uGC[idsInf_p2ugc] == 0))
  }


  ##############################################################################
  ## RIMMING ##
  ##############################################################################

  ri <- dat$temp$ri

  # Pharyngeal-to-Rectal GC --------------------------------------------------

  # P1 infects P2
  p1Inf_p2rgc <- which(
    pGC[ri[, "p1"]] == 1 & pGC.infTime[ri[, "p1"]] < at &
    rGC[ri[, "p2"]] == 0 & ri[, "ins.rim"] == 1
  )
  # P2 infects P1
  p2Inf_p2rgc <- which(
    pGC[ri[, "p2"]] == 1 & pGC.infTime[ri[, "p2"]] < at &
    rGC[ri[, "p1"]] == 0 & ri[, "ins.rim"] == 0
  )

  allActs_p2rgc <- c(p1Inf_p2rgc, p2Inf_p2rgc)

  # REVIEW: This step may not be necessary and may simply use up memory.
  #         Check for each transmission pathway.
  # Pathway-specific transmission probability
  tprob_p2rgc <- rep(dat$param$p2rgc.tprob, length(allActs_p2rgc))

  # Stochastic transmission
  trans_p2rgc <- rbinom(length(allActs_p2rgc), 1, tprob_p2rgc)

  # Determine the newly infected partner
  idsInf_p2rgc <- NULL
  if (sum(trans_p2rgc, na.rm = TRUE) > 0) {
    transRI_p2rgc <- ri[allActs_p2rgc[trans_p2rgc == 1], , drop = FALSE]
    idsInf_p2rgc <- c(
      intersect(ri[p1Inf_p2rgc, "p2"], transRI_p2rgc[, "p2"]),
      intersect(ri[p2Inf_p2rgc, "p1"], transRI_p2rgc[, "p1"])
    )
    stopifnot(all(rGC[idsInf_p2rgc] == 0))
  }

  # Rectal-to-Pharyngeal GC --------------------------------------------------

  # P1 infects P2
  p1Inf_r2pgc <- which(
    rGC[ri[, "p1"]] == 1 & rGC.infTime[ri[, "p1"]] < at &
    pGC[ri[, "p2"]] == 0 & ri[, "ins.rim"] == 0
  )

  # P2 infects P1
  p2Inf_r2pgc <- which(
    rGC[ri[, "p2"]] == 1 & rGC.infTime[ri[, "p2"]] < at &
    pGC[ri[, "p1"]] == 0 & ri[, "ins.rim"] == 1
  )

  allActs_r2pgc <- c(p1Inf_r2pgc, p2Inf_r2pgc)

  # Pathway-specific transmission probability
  tprob_r2pgc <- rep(dat$param$r2pgc.tprob, length(allActs_r2pgc))

  # Stochastic transmission with pathway-specific probability
  trans_r2pgc <- rbinom(length(allActs_r2pgc), 1, tprob_p2rgc)

  # Determine the newly infected partner
  idsInf_r2pgc <- NULL
  if (sum(trans_r2pgc, na.rm = TRUE) > 0) {
    transRI_r2pgc <- ri[allActs_r2pgc[trans_r2pgc == 1], , drop = FALSE]
    idsInf_r2pgc <- c(
      intersect(ri[p1Inf_r2pgc, "p2"], transRI_r2pgc[, "p2"]),
      intersect(ri[p2Inf_r2pgc, "p1"], transRI_r2pgc[, "p1"])
    )

    stopifnot(all(pGC[idsInf_r2pgc] == 0))
  }


  ##############################################################################
  ## KISSING ##
  ##############################################################################

  # Pharyngeal-to-Pharyngeal GC ------------------------------------------------

  kiss <- dat$temp$kiss

  # P1 infects P2
  p1Inf_p2pgc <- which(
    pGC[kiss[, "p1"]] == 1 &
    pGC.infTime[kiss[, "p1"]] < at &
    pGC[kiss[, "p2"]] == 0
  )

  # P2 infects P1
  p2Inf_p2pgc <- which(
    pGC[kiss[, "p2"]] == 1 &
    pGC.infTime[kiss[, "p2"]] < at &
    pGC[kiss[, "p1"]] == 0
  )

  allActs_p2pgc <- c(p1Inf_p2pgc, p2Inf_p2pgc)
  tprob_p2pgc <- rep(dat$param$p2pgc.tprob, length(allActs_p2pgc))
  trans_p2pgc <- rbinom(length(allActs_p2pgc), 1, dat$param$p2pgc.tprob)

  # Determine the newly infected partner
  idsInf_p2pgc <- NULL
  if (sum(trans_p2pgc, na.rm = TRUE) > 0) {
    transKL_p2pgc <- kiss[allActs_p2pgc[trans_p2pgc == 1], , drop = FALSE]
    idsInf_p2pgc <- c(
      intersect(kiss[p1Inf_p2pgc, "p2"], transKL_p2pgc[, "p2"]),
      intersect(kiss[p2Inf_p2pgc, "p1"], transKL_p2pgc[, "p1"])
    )
    }


  ##############################################################################
  ## RANDOMIZE ACT ORDER ##
  ##############################################################################

  idsInf_objects <- ls(pattern = "idsInf")
  idsInf_list <- lapply(
    setNames(idsInf_objects, idsInf_objects), function(x) {
      data.table(ids = get(x))
    }) %>%
    rbindlist(idcol = "pathway")

  ## Find IDs for those with multiple transmission events to an anatomic site.
  if (nrow(idsInf_list) > 0) {

    inf_dest_multiEvents <- idsInf_list[,
      dest := substring(pathway, 10)
      ][, .N, .(ids, dest)
      ][N > 1, .(ids, N)]

    ## If any IDs were stochastically assigned more than one transmission event
    ## to a single anatomic site, randomly select one of the events.
    if (any(inf_dest_multiEvents[, N] > 1)) {
      idsInf_list <- inf_dest_multiEvents[idsInf_list, on = "ids"]

      multiEvents_select <- split(
        idsInf_list[N > 1],
        by = c("ids", "dest")
      ) %>%
        lapply(., function(x) x[sample(1:.N)][1]) %>%
        rbindlist()

      ## Drop rows with multi-event IDs and then re-add only the selected
      ## transmission event.
      idsInf_list <- idsInf_list[!ids %in% multiEvents_select$ids]
      idsInf_list <- rbind(idsInf_list, multiEvents_select)
      idsInf_list <- idsInf_list[, .(ids, pathway)]
    }


  ##############################################################################
  ## UPDATE ATTRIBUTES ##
  ##############################################################################

    idsInf_destR <- idsInf_list[pathway %like% "rgc", ids]
    rGC[idsInf_destR] <- 1
    rGC.infTime[idsInf_destR] <- at
    rGC.sympt[idsInf_destR] <- rbinom(length(idsInf_destR), 1, rgc.sympt.prob)
    rGC.timesInf[idsInf_destR] <- rGC.timesInf[idsInf_destR] + 1

    idsInf_destU <- idsInf_list[pathway %like% "ugc", ids]
    uGC[idsInf_destU] <- 1
    uGC.infTime[idsInf_destU] <- at
    uGC.sympt[idsInf_destU] <- rbinom(length(idsInf_destU), 1, ugc.sympt.prob)
    uGC.timesInf[idsInf_destU] <- uGC.timesInf[idsInf_destU] + 1

    idsInf_destP <- idsInf_list[pathway %like% "pgc", ids]
    pGC[idsInf_destP] <- 1
    pGC.infTime[idsInf_destP] <- at
    pGC.sympt[idsInf_destP] <- rbinom(length(idsInf_destP), 1, pgc.sympt.prob)
    pGC.timesInf[idsInf_destP] <- pGC.timesInf[idsInf_destP] + 1


  ##############################################################################
  ## OUTPUT ##
  ##############################################################################

    # Attributes
    dat$attr$rGC <- rGC
    dat$attr$uGC <- uGC
    dat$attr$pGC <- pGC

    dat$attr$rGC.infTime <- rGC.infTime
    dat$attr$uGC.infTime <- uGC.infTime
    dat$attr$pGC.infTime <- pGC.infTime

    dat$attr$rGC.sympt <- rGC.sympt
    dat$attr$uGC.sympt <- uGC.sympt
    dat$attr$pGC.sympt <- pGC.sympt

    dat$attr$rGC.timesInf <- rGC.timesInf
    dat$attr$uGC.timesInf <- uGC.timesInf
    dat$attr$pGC.timesInf <- pGC.timesInf

    # Tally incidence
    incData <- idsInf_list[, ":="(
      pathway = substring(pathway, 8),
      race = race[ids],
      age.grp = age.grp[ids]
    )]

    cat("\ndata.table: incData ----------------------------------------------\n")
    print(incData)

    ## Incidence by transmission pathway, race/ethniciy, and age group.
    incid.byDemog <- incData[, .(incid = .N), keyby = .(race, age.grp, pathway)]
    incid.byDemog <- incid.byDemog[, anatsite := substring(pathway, 3)]
    setkeyv(incid.byDemog, c("race", "age.grp"))

    cat("\ndata.table: incid.byDemog -----------------------------------------\n")
    print(incid.byDemog)
    # label vectors for incidence assignments below
    anatsites <- c("rgc", "ugc", "pgc")
    races <- c("B", "H", "O", "W")

    # store overall GC incidence
    dat$epi$incid.gc[at] <- incid.byDemog[, sum(incid)]

    # store incidence by anatomic site
    lapply(anatsites, function(x) {
      dat$epi[[paste0("incid.", x)]][at] <<-
        incid.byDemog[pathway %like% paste0(x, "$"), sum(incid)]
    })

    # NOTE
    # In the section below, the joins and and ifelse() statements ensure
    # the vector lengths for incidence by sub-stratum are the same across
    # simulations. Prevents errors due to timesteps in which no one in a
    # given small stratum was infected.

    # Incidence by anatomic site and race
    lu <- as.data.table(expand.grid(
      anatsite = anatsites,
      race = match(races, races)
    ))

    cat("\nANATSITExRACE LOOKUP ------------ \n")
    print(lu)

    incid.ar <- incid.byDemog[, .(incid = sum(incid)), .(anatsite, race)]
    cat("\nPREJOIN.INCID.AR --------------- \n")
    print(incid.ar)

    incid.ar <- incid.ar[lu, on = c("anatsite", "race")]
    cat("\nPOSTJOIN.INCID.AR -------------- \n")
    print(incid.ar)

    lapply(seq_len(nrow(incid.ar)), function(x, racelabs = races) {
      rslug <- racelabs[incid.ar[x, race]]
      aslug <- incid.ar[x, anatsite]

      incid.curr <- incid.ar[x, incid]
      dat$epi[[paste0("incid.", rslug, ".", aslug)]][at] <<-
        ifelse(!is.na(incid.curr), incid.curr, 0)
    })

    # Incidence by anatomic site and age group
    lu <- as.data.table(expand.grid(
      anatsite = anatsites,
      age.grp = 1:5
    ))

    cat("\nANATSITExAGE.GRP LOOKUP ------------ \n")
    print(lu)

    incid.aa <- incid.byDemog[, .(incid = sum(incid)), .(anatsite, age.grp)]
    cat("\nPREJOIN.INCID.AA --------------- \n")
    print(incid.aa)

    incid.aa <- incid.aa[lu, on = c("anatsite", "age.grp")]
    cat("\nPOSTJOIN.INCID.AA --------------- \n")
    print(incid.aa)

    lapply(seq_len(nrow(incid.aa)), function(x) {
      gslug <- incid.aa[x, age.grp]
      aslug <- incid.aa[x, anatsite]

      incid.curr <- incid.aa[x, incid]
      dat$epi[[paste0("incid.age.", gslug, ".", aslug)]][at] <<-
        ifelse(!is.na(incid.curr), incid.curr, 0)
    })

    # Incidence by anatomic site, race, and age group
    lu <- as.data.table(expand.grid(
      anatsite = anatsites,
      race = match(races, races),
      age.grp = 1:5
    ))
    cat("\nANATSITExRACExAGR.GRP LOOKUP ------------ \n")
    print(lu)

    incid.ara <- incid.byDemog[, .(
      incid = sum(incid)
    ), .(anatsite, race, age.grp)]

    cat("\nPREJOIN.INCID.ARA --------------- \n")
    print(incid.ara)

    incid.ara <- incid.ara[lu, on = c("anatsite", "race", "age.grp")]
    cat("\nPREJOIN.INCID.ARA --------------- \n")
    print(incid.ara)

    lapply(seq_len(nrow(incid.ara)), function(x, racelabs = races) {
      rslug <- racelabs[incid.ara[x, race]]
      gslug <- incid.ara[x, age.grp]
      aslug <- incid.ara[x, anatsite]

      incid.curr <- incid.ara[x, incid]
      dat$epi[[paste0("incid.", rslug, ".age", gslug, ".", aslug)]][at] <<-
        ifelse(!is.na(incid.curr), incid.curr, 0)
    })

    # incidence by transmission pathway
    incid.tp <- incid.byDemog[, .(incid = sum(incid)), pathway]
    lapply(seq_len(nrow(incid.tp)), function(x) {
      incid.curr <- incid.tp[x, incid]
      dat$epi[[paste0("incid.", incid.tp[x, pathway])]][at] <<-
        ifelse(!is.na(incid.curr), incid.curr, 0)
    })

    # Check all infected have all STI attributes
    stopifnot(
      all(!is.na(rGC.infTime[rGC == 1])),
      all(!is.na(rGC.sympt[rGC == 1])),
      all(!is.na(uGC.infTime[uGC == 1])),
      all(!is.na(uGC.sympt[uGC == 1])),
      all(!is.na(pGC.infTime[pGC == 1])),
      all(!is.na(pGC.sympt[pGC == 1]))
    )
  }

  return(dat)

}
