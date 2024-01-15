#' @title STI Treatment Module
#'
#' @description Stochastically simulates GC/CT diagnosis and treatment.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @import pscl
#' @export
#'
stitx_msm <- function(dat, at) {

  # Parameters
  gc.sympt.prob.tx <- dat$param$gc.sympt.prob.tx
  gc.asympt.prob.tx <- dat$param$gc.asympt.prob.tx

  prep.sti.screen.int <- dat$param$prep.sti.screen.int
  prep.sti.prob.tx <- dat$param$prep.sti.prob.tx

  stitest.mod <- dat$param$epistats$stitest

  # Attributes
  race <- dat$attr$race
  age <- dat$attr$age

  ## Identify agents with symptoms by anatomic site

  ### Rectal: infected, symptomatic, no treatment
  idsRGC_sympt <- which(
    dat$attr$rGC == 1 &
    dat$attr$rGC.infTime < at &
    dat$attr$rGC.sympt == 1 &
    is.na(dat$attr$rGC.tx)
  )

  ### Urethral: infected, symptomatic, not on treatment
  idsUGC_sympt <- which(
    dat$attr$uGC == 1 &
    dat$attr$uGC.infTime < at &
    dat$attr$uGC.sympt == 1 &
    is.na(dat$attr$uGC.tx)
  )

  ### Pharyngeal: infected, symptomatic, not on treatment
  idsPGC_sympt <- which(
    dat$attr$pGC == 1 &
    dat$attr$pGC.infTime < at &
    dat$attr$pGC.sympt == 1 &
    is.na(dat$attr$pGC.tx)
  )

  ### GC symptoms at any site
  idsGC_sympt <- union(c(idsRGC_sympt, idsUGC_sympt), idsPGC_sympt)

  ## Identify infected agents asymptomatic at all sites

  ### Rectal: infected, asymptomatic, no treatment
  idsRGC_asympt <- which(
    dat$attr$rGC == 1 &
    dat$attr$rGC.infTime < at &
    dat$attr$rGC.sympt == 0 &
    is.na(dat$attr$rGC.tx) &
    dat$attr$prepStat == 0
  )

  ### Urethral: infected, asymptomatic, no treatment
  idsUGC_asympt <- which(
    dat$attr$uGC == 1 &
    dat$attr$uGC.infTime < at &
    dat$attr$uGC.sympt == 0 &
    is.na(dat$attr$uGC.tx) &
    dat$attr$prepStat == 0
  )

  ### Pharyngeal: infected, asymptomatic, no treatment
  idsPGC_asympt <- which(
    dat$attr$pGC == 1 &
    dat$attr$pGC.infTime < at &
    dat$attr$pGC.sympt == 0 &
    is.na(dat$attr$pGC.tx) &
    dat$attr$prepStat == 0
  )

  stopifnot(
    dat$control$stiScreeningProtocol %in%
    c("base", "symptomatic", "universal", "cdc")
  )

  sti_clinic_scenarios <- c("base", "symptomatic", "universal")

  if (dat$control$stiScreeningProtocol %in% sti_clinic_scenarios) {

    idsGC_asympt <- union(c(idsRGC_asympt, idsUGC_asympt), idsPGC_asympt)

    # Simulate STI test seeking
    pred_df_sti_seek <- data.table(
      race.cat = race,
      age = age
    )

    # predict yearly test rate and divide by 52 to get weekly rate
    prob_sti_test <- predict(
      stitest.mod,
      newdata = pred_df_sti_seek,
      type = "response"
    ) / 52

    # Scale probability of testing in week among those with GC symptoms
    # at any anatomic site.
    prob_sti_test[idsGC_sympt] <-
      prob_sti_test[idsGC_sympt] * dat$param$gc.sympt.seek.test.scale

    # use weekly rate as the probability of testing within that week
    seekTest <- rbinom(length(prob_sti_test), 1, prob_sti_test)

    # Set seek test to 0 for those on treatment
    idsTx <- which(
      !is.na(dat$attr$rGC.tx) |
      !is.na(dat$attr$uGC.tx) |
      !is.na(dat$attr$pGC.tx)
    )

    seekTest[idsTx] <- 0
    ids_seekTest <- which(seekTest == 1)

    gc.sympt.prob.test <- dat$param$gc.sympt.prob.test
    gc.asympt.prob.test <- dat$param$gc.asympt.prob.test

    # Determine who is seeking an STI test
    ## Symptomatic, seeking test
    idsRGC_sympt_seekingTest <- intersect(idsRGC_sympt, ids_seekTest)
    idsUGC_sympt_seekingTest <- intersect(idsUGC_sympt, ids_seekTest)
    idsPGC_sympt_seekingTest <- intersect(idsPGC_sympt, ids_seekTest)

    ## Asymptomatic, seeking test
    idsRGC_asympt_seekingTest <- setdiff(ids_seekTest, idsRGC_sympt_seekingTest)
    idsUGC_asympt_seekingTest <- setdiff(ids_seekTest, idsUGC_sympt_seekingTest)
    idsPGC_asympt_seekingTest <- setdiff(ids_seekTest, idsPGC_sympt_seekingTest)

    # Simulate test receipt within the clinic

    if (dat$control$stiScreeningProtocol == "base") {

      ## Symptomatic anatomic sites
      rgc.sympt.test <-
        rbinom(length(idsRGC_sympt_seekingTest), 1, gc.sympt.prob.test[1])

      ugc.sympt.test <-
        rbinom(length(idsUGC_sympt_seekingTest), 1, gc.sympt.prob.test[2])

      pgc.sympt.test <-
        rbinom(length(idsPGC_sympt_seekingTest), 1, gc.sympt.prob.test[3])

      ## Asymptomatic anatomic sites
      rgc.asympt.test <-
        rbinom(length(idsRGC_asympt_seekingTest), 1, gc.asympt.prob.test[1])

      ugc.asympt.test <-
        rbinom(length(idsUGC_asympt_seekingTest), 1, gc.asympt.prob.test[2])

      pgc.asympt.test <-
        rbinom(length(idsPGC_asympt_seekingTest), 1, gc.asympt.prob.test[3])

      # Assign test receipt
      idsRGC_getTest <- union(
        idsRGC_sympt_seekingTest[rgc.sympt.test == 1],
        idsRGC_asympt_seekingTest[rgc.asympt.test == 1]
      )

      idsUGC_getTest <- union(
        idsUGC_sympt_seekingTest[ugc.sympt.test == 1],
        idsUGC_asympt_seekingTest[ugc.asympt.test == 1]
      )

      idsPGC_getTest <- union(
        idsPGC_sympt_seekingTest[pgc.sympt.test == 1],
        idsPGC_asympt_seekingTest[pgc.asympt.test == 1]
      )
    }

    ## Anyone seeking a test receives a test at all symptomatic sites.
    ## - No asymptomatic testing
    if (dat$control$stiScreeningProtocol == "symptomatic") {

      idsRGC_getTest <- idsRGC_sympt_seekingTest
      idsUGC_getTest <- idsUGC_sympt_seekingTest
      idsPGC_getTest <- idsPGC_sympt_seekingTest

    }

    ## Anyone seeking a test receives a test at all anatomic sites.
    if (dat$control$stiScreeningProtocol == "universal") {

      idsRGC_getTest <-
        union(idsRGC_sympt_seekingTest, idsRGC_asympt_seekingTest)

      idsUGC_getTest <-
        union(idsUGC_sympt_seekingTest, idsUGC_asympt_seekingTest)

      idsPGC_getTest <-
        union(idsPGC_sympt_seekingTest, idsPGC_asympt_seekingTest)

    }

  }

  ## ===========================================================================
  ## CDC Guidelines
  ## ===========================================================================

  if (dat$control$stiScreeningProtocol == "cdc") {

    ## General annual screening for all MSM

    cdc.sti.int <- dat$param$cdc.sti.int

    idsRGC_getTest_int <- which(
      (at - dat$attr$last.rGC.test >= cdc.sti.int &
       at - dat$attr$last.rectal.exp <= cdc.sti.int) |
      is.na(dat$attr$last.rGC.test)
    )

    idsUGC_getTest_int <- which(
      (at - dat$attr$last.uGC.test >= cdc.sti.int &
       at - dat$attr$last.ureth.exp <= cdc.sti.int) |
      is.na(dat$attr$last.uGC.test)
    )

    idsPGC_getTest_int <- which(
      (at - dat$attr$last.pGC.test >= cdc.sti.int &
       at - dat$attr$last.phar.exp <= cdc.sti.int) |
      is.na(dat$attr$last.pGC.test)
    )

    ## More frequent screening for high-risk MSM
    ## - Defined as those with multiple partners, OR
    ## - Those with partners who have multiple partners
    cdc.sti.hr.int <- dat$param$cdc.sti.hr.int

    ids_multiPart <- which(
      dat$attr$deg.main + dat$attr$deg.casl > 1
    )

    # those with at least one main/casual partners
    ids_mcPart <- which(dat$attr$deg.main + dat$attr$deg.casl > 0)
    el.mc <- rbind(dat$el[[1]], dat$el[[1]])

    part_multiBool <- sapply(ids_mcPart, function(x) {

      pids <- union(el.mc[, 2][el.mc[, 1] == x], el.mc[, 1][el.mc[, 2] == x])

      part_pnums <- sapply(pids, function(x) {
        pnums <- dat$attr$deg.main[x] + dat$attr$deg.casl[x]
        pnums
      })

      # check if any of the agent's partners has multiple partners
      part_multiBool <- any(part_pnums > 1)
      part_multiBool
    })

    ids_highRiskElig <- ids_mcPart[part_multiBool]

    idsRGC_getTest_hr <- intersect(
      ids_highRiskElig,
      which(
        at - dat$attr$last.rGC.test >= cdc.sti.hr.int |
        is.na(dat$attr$last.rGC.test)
      )
    )

    idsUGC_getTest_hr <- intersect(
      ids_highRiskElig,
      which(
        at - dat$attr$last.uGC.test >= cdc.sti.hr.int |
        is.na(dat$attr$last.uGC.test)
      )
    )

    idsPGC_getTest_hr <- intersect(
      ids_highRiskElig,
      which(
        at - dat$attr$last.pGC.test >= cdc.sti.hr.int |
        is.na(dat$attr$last.pGC.test)
      )
    )

    idsRGC_getTest <- union(idsRGC_getTest_int, idsRGC_getTest_hr)
    idsUGC_getTest <- union(idsUGC_getTest_int, idsUGC_getTest_hr)
    idsPGC_getTest <- union(idsPGC_getTest_int, idsPGC_getTest_hr)

  }


  ## ===========================================================================
  ## Record STI testing
  ## ===========================================================================

  dat$epi$num.rgc.test[at] <- length(idsRGC_getTest)
  dat$epi$num.ugc.test[at] <- length(idsUGC_getTest)
  dat$epi$num.pgc.test[at] <- length(idsPGC_getTest)

  dat$attr$last.rGC.test[idsRGC_getTest] <- at
  dat$attr$last.uGC.test[idsUGC_getTest] <- at
  dat$attr$last.pGC.test[idsPGC_getTest] <- at


  ## ===========================================================================
  ## DEBUG BLOCK
  ## ===========================================================================

  if (dat$control$debug_stitx) {
    cat("\nPROB_STI_TEST INITIAL MISSING")
    print(sum(is.na(prob_sti_test)))


    cat("\nidsGC_sympt\n")
    print(idsGC_sympt)

    cat("\nProb. testing among symptomatic\n", sep = "")
    print(prob_sti_test[idsGC_sympt])

    cat("\nPROB_STI_TEST AFTER SCALE MISSING")
    print(sum(is.na(prob_sti_test)))

    cat("\nTesting status RGC\n")
    print(table(dat$attr$rGC.tx, exclude = NULL))

    cat("\nIDs of those already on treatment\n")
    print(head(idsTx))

    cat("\nSeek testing (yes/no)\n")
    print(head(seekTest))

    cat("\nIDs of those seeking test\n")
    print(head(ids_seekTest))

    cat("\nIDs of those getting a rectal test\n")
    print(head(idsRGC_getTest))

    cat("\nIDs of those getting a urethral test\n")
    print(head(idsUGC_getTest))

    cat("\nIDs of those getting a pharyngeal test\n")
    print(head(idsPGC_getTest))
  }


  ## ===========================================================================
  ## Treatment
  ## ===========================================================================
  ## All Treated GC ##

  # IDs of men put on treatment
  txRGC <- intersect(idsRGC_getTest, which(dat$attr$rGC == 1))
  txUGC <- intersect(idsUGC_getTest, which(dat$attr$uGC == 1))
  txPGC <- intersect(idsPGC_getTest, which(dat$attr$pGC == 1))

  ## IDs of men eligible for treatment
  idsRGC_txElig <- union(idsRGC_sympt, idsRGC_asympt)
  idsUGC_txElig <- union(idsUGC_sympt, idsUGC_asympt)
  idsPGC_txElig <- union(idsPGC_sympt, idsPGC_asympt)

  ## Interval-based treatment for MSM on PrEP ##
  idsSTI_screen <- which(
    dat$attr$prepStartTime == at |
    (at - dat$attr$prepLastStiScreen >= prep.sti.screen.int)
  )

  dat$attr$prepLastStiScreen[idsSTI_screen] <- at

  idsRGC_prep_tx <- intersect(
    idsSTI_screen,
    which(
      dat$attr$rGC == 1 &
      dat$attr$rGC.infTime < at &
      is.na(dat$attr$rGC.tx.prep)
    )
  )

  idsUGC_prep_tx <- intersect(
    idsSTI_screen,
    which(
      dat$attr$uGC == 1 &
      dat$attr$uGC.infTime < at &
      is.na(dat$attr$uGC.tx.prep)
    )
  )

  idsPGC_prep_tx <- intersect(
    idsSTI_screen,
    which(
      dat$attr$pGC == 1 &
      dat$attr$pGC.infTime < at &
      is.na(dat$attr$pGC.tx.prep)
    )
  )

  txRGC_prep <- idsRGC_prep_tx[
    which(rbinom(length(idsRGC_prep_tx), 1, prep.sti.prob.tx) == 1)
  ]

  txUGC_prep <- idsUGC_prep_tx[
    which(rbinom(length(idsUGC_prep_tx), 1, prep.sti.prob.tx) == 1)
  ]

  txPGC_prep <- idsPGC_prep_tx[
    which(rbinom(length(idsPGC_prep_tx), 1, prep.sti.prob.tx) == 1)
  ]

  ## Update Attributes ##

  ## Not on PrEP
  dat$attr$rGC.tx[idsRGC_txElig] <- 0
  dat$attr$uGC.tx[idsUGC_txElig] <- 0
  dat$attr$pGC.tx[idsPGC_txElig] <- 0

  dat$attr$rGC.tx[txRGC] <- 1
  dat$attr$uGC.tx[txUGC] <- 1
  dat$attr$pGC.tx[txPGC] <- 1

  ## On PrEP
  dat$attr$rGC.tx.prep[txRGC_prep] <- 1
  dat$attr$uGC.tx.prep[txUGC_prep] <- 1
  dat$attr$pGC.tx.prep[txPGC_prep] <- 1

  dat$attr$rGC.tx.prep[idsRGC_prep_tx] <- 0
  dat$attr$uGC.tx.prep[idsUGC_prep_tx] <- 0
  dat$attr$pGC.tx.prep[idsPGC_prep_tx] <- 0


  # REVIEW I'm not sure why we need treatment variables for every anatomic site.
  #        If the person is treated, we assume treatment affects all sites.
  ## Add tx at other anatomical site ##
  dat$attr$rGC.tx[which((dat$attr$uGC.tx == 1 | dat$attr$uGC.tx.prep == 1) &
                        dat$attr$rGC == 1)] <- 1

  dat$attr$uGC.tx[which((dat$attr$rGC.tx == 1 | dat$attr$rGC.tx.prep == 1) &
                        dat$attr$uGC == 1)] <- 1

  dat$attr$pGC.tx[which((dat$attr$pGC.tx == 1 | dat$attr$pGC.tx.prep == 1) &
                        dat$attr$pGC == 1)] <- 1

  ## Summary treatment indicator
  idsGC_tx <- which(
    dat$attr$rGC.tx == 1 |
    dat$attr$uGC.tx == 1 |
    dat$attr$pGC.tx == 1 |
    dat$attr$rGC.tx.prep == 1 |
    dat$attr$uGC.tx.prep == 1 |
    dat$attr$pGC.tx.prep == 1
  )

  dat$attr$anyGC.tx[idsGC_tx] <- 1

  return(dat)
}
